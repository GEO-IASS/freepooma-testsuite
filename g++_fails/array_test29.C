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
// array_test29.cpp verify correctnes of stencil objects with expressions
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Arrays.h"
#include "Utilities/Tester.h"
#include <iostream>

class EvaluateExpr
{
public:
  template <class A>
  typename A::Element_t operator()(const A& x, int i) const
  {
	  return x.read(i);
  }

  int lowerExtent(int) const { return 0; }
  int upperExtent(int) const { return 0; }  
};

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> domain(8);
  UniformGridLayout<1> layout(domain, Loc<1>(2), GuardLayers<1>(1), DistributedTag());
  Array<1,int,MultiPatch<UniformTag, Remote<Brick> > > a(layout), b(layout), c(layout);

  a = 0;
  b = 1;
  c = 2;
  a(domain) = Stencil<EvaluateExpr>()(b(domain-1)+c(domain+1), domain);
  tester.check("a = b(I-1) + c(I+1)", all(a(domain) == 3));
  tester.out() << a(domain) << std::endl;

  a = 0;
  b = 2;
  c = 3;
  a(domain) = b(domain) + Stencil<EvaluateExpr>()(b(domain)+c(domain+1), domain);
  tester.check("a = b + b(I-1) + c", all(a(domain) == 7));
  tester.out() << a(domain) << std::endl;

  int retval = tester.results("array_test29");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test29.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
