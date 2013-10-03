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
// Miscellaneous tests of combinations of engines.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Engine/BrickEngine.h"
#include "Engine/ConstantFunctionEngine.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"

template<class Array>
bool isSmall(const Array &a)
{
  double epsilon = 0.000001;
  int i;
  int first = a.domain()[0].first();
  int last  = a.domain()[0].last();
  double sum = 0.0;
  for (i = first; i <= last; ++i)
  {
    sum += a.read(i)*a.read(i);
  }
  return (sum < epsilon);
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int from = 1;
  int to = 5;

  Interval<1> dom(from,to);

  ModelElement<Vector<3,double> > v(Vector<3,double>(1,2,3));

  Array<1,Vector<3,double>,ConstantFunction> a(dom,v);
  Array<1,Vector<3,double>,Brick> b(dom);
  Array<1,double,Brick> c(dom);

  b = 2*a;
  c = 2*a.comp(2);

  tester.out() 
    << b << '\n'
    << b.comp(2) << '\n'
    << c << std::endl;

  tester.check(isSmall(b.comp(2)-c));

  int ret = tester.results("array_test16");
  Pooma::finalize();
  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test16.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
