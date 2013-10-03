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
// evaluatorTest3 - a simple patch function
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Evaluator/PatchFunction.h"
#include "Utilities/Tester.h"
#include <iostream>

// This really simple example doesn't do anything a stencil couldn't do,
// but you could imagine writing something more complicated.

struct MyFunction
{
  template<class ArrayPatch>
  void apply(const ArrayPatch &a) const
  {
    int i;
    int from = a.domain()[0].first();
    int to = a.domain()[0].last();
    for (i = from; i <= to; ++i)
    {
      if (a(i)>5.0)
      {
	a(i) /= 4.0;
      }
    }
  }
};

// version that takes two arguments.  You are writing to the
// first argument and reading from the second.

struct MyFunction2
{
  MyFunction2(double v1,double v2) : v1_m(v1),v2_m(v2) { }

  template<class ArrayPatch1, class ArrayPatch2>
  void apply(const ArrayPatch1 &a1, const ArrayPatch2 &a2) const
  {
    int i;
    int from = a1.domain()[0].first();
    int to = a1.domain()[0].last();
    for (i = from; i <= to; ++i)
    {
      if (a2(i)>v1_m)
      {
	a1(i) /= v2_m;
      }
    }
  }
  double v1_m,v2_m;
};

struct MyFunction3
{
  MyFunction3(double v1,double v2) : v1_m(v1),v2_m(v2) { }

  template<class ArrayPatch1, class ArrayPatch2, class ArrayPatch3>
  void apply(const ArrayPatch1 &a1, const ArrayPatch2 &a2,
	     const ArrayPatch3 &a3) const
  {
    int i;
    int from = a1.domain()[0].first();
    int to = a1.domain()[0].last();
    for (i = from; i <= to; ++i)
    {
      if (a2(i)>v1_m)
      {
	a1(i) /= v2_m;
      }
      else
      {
	a1(i) += a3(i);
      }
    }
  }
  double v1_m,v2_m;
};

struct TestFunction
{
  TestFunction(Inform &inform) : inform_m(&inform) { }
  TestFunction(const TestFunction &tf) : inform_m(tf.inform_m) { }

  template<class Patch>
  void apply(const Patch &a)
  {
    (*inform_m) << "test:" << a << std::endl;
  }

  template<class Patch>
  void apply(const Patch &a, int node)
  {
    (*inform_m) << node << ":" << a << std::endl;
  }

  template<class P1, class P2>
  void apply(const P1 &a, const P2 &b, int node)
  {
    (*inform_m) << "a:" << node << ":" << a << std::endl;
    (*inform_m) << "b:" << ":" << b << std::endl;
  }

  template<class P1, class P2, class P3>
  void apply(const P1 &a, const P2 &b, const P3 &c, int node)
  {
    (*inform_m) << "a:" << node << ":" << a << std::endl;
    (*inform_m) << "b:" << ":" << b << std::endl;
    (*inform_m) << "c:" << ":" << c << std::endl;
  }

  mutable Inform *inform_m;
};


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Patch function test." << std::endl;
  tester.out() << "------------------------------------------------" << std::endl;

  int size = 120;

  Interval<1> domain(size);
  UniformGridPartition<1> partition(Loc<1>(10));
  UniformGridLayout<1> layout(domain, partition, ReplicatedTag());

  Array<1,double,MultiPatch<UniformTag,Brick> > a(layout),b(layout);

  int i;
  for (i = 0; i < size; ++i )
  {
    a(i) = i;
  }
  b = where(a>5.0,a/4.0,a);

  PatchFunction<MyFunction,PatchTag1>()(a);

  tester.out() << a << std::endl;
  tester.out() << b << std::endl;

  tester.check(sum((a-b)*(a-b))<0.001);

  UniformGridPartition<1> partition2(Loc<1>(12));
  UniformGridLayout<1> layout2(domain, partition2, ReplicatedTag());

  Array<1,double,MultiPatch<UniformTag,Brick> > a2(layout2),b2(layout2);

  for (i = 0; i < size; ++i )
  {
    a(i) = i;
    a2(i) = 3+i;
    b2(i) = 3+i;
  }
 
  PatchFunction<MyFunction2,PatchTag2>(5.0,4.0)(a2,a);
  b2 = where(a>5.0,b2/4.0,b2);

  tester.check(sum((a2-b2)*(a2-b2))<0.001);

  for (i = 0; i < size; ++i )
  {
    a(i) = i;
    a2(i) = 3+i;
    b2(i) = 3+i;
  }

  PatchFunction<MyFunction3,PatchTag3>(5.0,4.0)(a2,a,a);
  b2 = where(a>5.0,b2/4.0,b2+a);

  tester.out() << a << std::endl;
  tester.out() << a2 << std::endl;
  tester.out() << b2 << std::endl;

  tester.check(sum((a2-b2)*(a2-b2))<0.001);

  TestFunction tf(tester.out());
  PatchFunction<TestFunction, PatchParticle1<true> > test(tf);
  test(b2);
  PatchFunction<TestFunction, PatchParticle2<false, false> > test2(tf);
  test2(b2,a2);
  PatchFunction<TestFunction, PatchParticle3<false, false, false> > test3(tf);
  test3(b2,a2,a2);

  PatchFunction<TestFunction, PatchReadTag1> test4(tf);
  test4(b2*2+a);

  PatchFunction<TestFunction, PatchParticle1<false> > test5(tf);
  test5(b2*2+a2);

  tester.out() << "------------------------------------------------" << std::endl;
 
  int retval = tester.results("evaluatorTest3");

  Pooma::finalize();
  
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: evaluatorTest3.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
