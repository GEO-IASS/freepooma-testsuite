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
// evaluatorTest2 - a simple patch function using ScalarCode
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Pooma/Fields.h" // for PerformUpdateTag() only!
#include "Evaluator/ScalarCode.h"
#include "Utilities/Tester.h"
#include <iostream>

// This really simple example doesn't do anything a stencil couldn't do,
// but you could imagine writing something more complicated.

struct MyFunction
{
  MyFunction() {}

  template<class A>
  void operator()(const A &a, const Loc<1> &i) const
  {
    if (a(i)>5.0)
      {
	a(i) /= 4.0;
      }
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(1);
    i.dimensions(1);
    i.lowerExtent(0) = 0;
    i.upperExtent(0) = 0;
    i.write(0, true);
    i.useGuards(0, false);
  }
};

struct MyFunction2
{
  MyFunction2() {}

  template<class Array1, class Array2>
  void operator()(const Array1 &a, const Array2 &b, const Loc<2> &i) const
  {
    Loc<2> dx(1, 0), dy(0, 1);
    a(i) = 0.25 * (b(i-dx) + b(i+dx) + b(i-dy) + b(i+dy));
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(2);
    i.dimensions(2);
    i.lowerExtent(0) = 1;
    i.upperExtent(0) = 1;
    i.lowerExtent(1) = 1;
    i.upperExtent(1) = 1;
    i.write(0, true);
    i.write(1, false);
    i.useGuards(0, false);
    i.useGuards(1, true);
  }
};

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  {
    int size = 120;

    Interval<1> domain(size);
    UniformGridPartition<1> partition(Loc<1>(10));
    UniformGridLayout<1> layout(domain, partition, ReplicatedTag());

    Array<1,double,MultiPatch<UniformTag, Brick> > a(layout),b(layout);

    int i;
    for (i = 0; i < size; ++i )
      {
	a(i) = i;
      }
    b = where(a>5.0,a/4.0,a);

    ScalarCode<MyFunction>()(a);

    tester.out() << a << std::endl;
    tester.out() << b << std::endl;

    tester.check(sum((a-b)*(a-b))<0.001);
  }

  {
    Interval<2> domain(9, 9);
    UniformGridLayout<2> layout(domain, Loc<2>(3, 3),
				GuardLayers<2>(1), ReplicatedTag());

    Array<2, double, MultiPatch<UniformTag, Brick> > a(layout), b(layout), c(layout);

    a(a.domain()) = 1.0;
    b(b.domain()) = iota(b.domain()).comp(0) + iota(b.domain()).comp(1);
    c(domain) = 0.25 * (b(domain - Loc<2>(1, 0)) + b(domain + Loc<2>(1, 0))
			+ b(domain - Loc<2>(0, 1)) + b(domain + Loc<2>(0, 1)));
    ScalarCode<MyFunction2>()(a, b);

    tester.out() << a << std::endl;
    tester.out() << c << std::endl;

    tester.check("MultiPatch setup", all(a(domain) == c(domain)));
  }

  {
    Interval<2> domain(9, 9);
    UniformGridLayout<2> layout(domain, Loc<2>(3, 3),
				GuardLayers<2>(1), DistributedTag());

    Array<2, double, MultiPatch<UniformTag, Remote<Brick> > > a(layout), b(layout), c(layout);

    a(a.domain()) = 1.0;
    b(b.domain()) = iota(b.domain()).comp(0) + iota(b.domain()).comp(1);
    c(domain) = 0.25 * (b(domain - Loc<2>(1, 0)) + b(domain + Loc<2>(1, 0))
			+ b(domain - Loc<2>(0, 1)) + b(domain + Loc<2>(0, 1)));
    ScalarCode<MyFunction2>()(a, b);

    tester.out() << a << std::endl;
    tester.out() << c << std::endl;

    tester.check("Remote MultiPatch setup", all(a(domain) == c(domain)));
  }

  int retval = tester.results("evaluatorTest2 (ScalarCode)");

  Pooma::finalize();
  
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: evaluatorTest2.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
