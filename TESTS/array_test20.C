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
// Array test 20: constant-function engine and index-function engine.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Engine/ConstantFunctionEngine.h"
#include "Engine/IndexFunctionEngine.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"

#include <iostream>

template<class EngineTag>
void test1(const Array<7, Vector<2, int>, EngineTag> &ca, 
  Pooma::Tester &t, int val)
{
  Interval<1> i1(5), i2(5,9);
  Range<1> r1(1,9,2), r2(0,8,4);
  
  int ans = sum(ca(i1, 3, r1, 4, r2, i2, 8).comp(1));
  t.check(ans == val);
  
  t.out() << ans << std::endl;
  t.out() << ca(i1, 3, r1, 4, r2, i2, 8).comp(1).domain() << std::endl;
}

template<class EngineTag>
void test2(const Array<1, Vector<2, int>, EngineTag> &ca, 
  Pooma::Tester &t, int val)
{
  Array<1, int> i(3);
  
  Pooma::blockAndEvaluate();
  i(0) = 3; i(1) = 11; i(2) = 16;
  
  int ans = sum(ca(i).comp(1));
  t.check(ans == val);
  
  t.out() << ans << std::endl;
  t.out() << ca(i).comp(1).domain() << std::endl;
}

struct VectorFunctor
{
  VectorFunctor() { }
  VectorFunctor(const VectorFunctor &) { }
  VectorFunctor &operator=(const VectorFunctor &) { return *this; }

  Vector<2, int> operator()(int i1) const
    {
      return Vector<2, int>(i1, 2 * i1);
    }
};

struct ConstantFunctor
{
  ConstantFunctor() { }
  ConstantFunctor(const ConstantFunctor &) { }
  ConstantFunctor &operator=(const ConstantFunctor &) { return *this; }

  Vector<2, int> operator()(int) const
    {
      return Vector<2, int>(1, 2);
    }
  Vector<2, int> operator()(int, int, int, int, int, int, int) const
    {
      return Vector<2, int>(1, 2);
    }
};

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  
  Array<7, Vector<2, int>, ConstantFunction> 
    cf(10, 10, 10, 10, 10, 10, 10);
  cf.engine().setConstant(Vector<2, int>(1, 2));
  
  test1(cf, tester, 750);

  Array<7, Vector<2, int>, IndexFunction<ConstantFunctor> > 
    ifa(10, 10, 10, 10, 10, 10, 10);
  ifa.engine().setFunctor(ConstantFunctor());
  
  test1(ifa, tester, 750);
  
  Array<1, Vector<2, int>, ConstantFunction> cf2(20);
  cf2.engine().setConstant(Vector<2, int>(1, 2));
  
  test2(cf2, tester, 6);
  
  Array<1, Vector<2, int>, IndexFunction<ConstantFunctor> > ifa2(20);
  ifa2.engine().setFunctor(ConstantFunctor());
  
  test2(ifa2, tester, 6);
  
  Array<1, Vector<2, int>, IndexFunction<VectorFunctor> > ifa3(20);
  ifa3.engine().setFunctor(VectorFunctor());
  
  test2(ifa3, tester, 60);
  
  int ret = tester.results("array_test20");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test20.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
