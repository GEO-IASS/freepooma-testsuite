// -*- C++ -*-
//
// Copyright (C) 1998, 1999, 2000, 2002  Los Alamos National Laboratory,
// Copyright (C) 1998, 1999, 2000, 2002  CodeSourcery, LLC
//
// This file is part of FreePOOMA.
//
// FreePOOMA is free software; you can redistribute it and/or modify it
// under the terms of the Expat license.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Expat
// license for more details.
//
// You should have received a copy of the Expat license along with
// FreePOOMA; see the file LICENSE.
//

//-----------------------------------------------------------------------------
// array_test23.cpp more tests to verify correctnes of stencil objects
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Engine/BrickEngine.h"
#include "Engine/Stencil.h"
#include "Array/Array.h"
#include "Pooma/FunctorResult.h"
#include "Pooma/Indices.h"

#include <iostream>
#include <cmath>

template<class Array>
bool isSmall(const Array &a)
{
  double epsilon = 0.000000001;
  return (sum(a*a) < epsilon);
}

class AsymDoof
{
public:
  AsymDoof() { }
  AsymDoof(const AsymDoof &) { }

  template <class A>
  inline
  typename A::Element_t
  operator()(const A& x, int i, int j) const
  {
    return ( (1.0/15.0) *
             ( x(i+1,j+1) + 2*x(i+1,j  ) + 3*x(i+1,j-1) +
               3*x(i  ,j+1) + x(i  ,j  ) + 2*x(i  ,j-1) +
               4*x(i-1,j+1) + 3*x(i-1,j  ) + 5*x(i-1,j-1) ) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 1; }  

private:
};

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": More stencil tests.."
	       << std::endl;
  tester.out() << "------------------------------------------------"
	       << std::endl;


  Stencil<AsymDoof> doof;

  Interval<1> inew(10);
  Interval<2> d2(inew,inew);
  Interval<2> inset = doof.insetDomain(d2);

  Array<2,double,Brick> init(d2),version1(d2),version2(d2);

  init = iota(d2).comp(0) + sin(iota(d2).comp(1) * 0.4);
  version1 = 0.0;

  Interval<1> d1(4);

  Array<1,double,Brick> g(d1), h(d1);
  Array<1,Loc<2>,Brick> ind(d1);

  ind(0) = Loc<2>(3,4);
  ind(1) = Loc<2>(7,4);
  ind(2) = Loc<2>(4,4);
  ind(3) = Loc<2>(5,6);

  version1(inset) = doof(init);

  g = version1(inset)(ind);
  h = doof(init)(ind);

  tester.out() << version1 << std::endl;
  tester.out() << g << std::endl;
  tester.out() << h << std::endl;
  tester.check("indirection of stencil", isSmall(g - h));

  tester.out() << "------------------------------------------------"
	       << std::endl;

  int retval = tester.results("array_test23");

  Pooma::finalize();

  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test23.cpp,v $   $Author: richi $
// $Revision: 1.7 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
