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
// test of RNGSimple
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Utilities/Tester.h"
#include "Functions/RNGSimple.h"
#include "Functions/RNGComponent.h"

#include <iostream>

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  typedef Array<2, RNGSimple, MultiPatch<GridTag, Brick> > Array_t;
  typedef Array<2, double, MultiPatch<GridTag, Brick> > ArrayR_t;

  Interval<1> ii(20);
  Interval<2> dom(ii, ii);

  Loc<2> blocks(4,4);
  GridPartition<2> partition(blocks);
  GridLayout<2> layout(dom, partition,ReplicatedTag() );

  Array_t rng0(layout), rng1(layout);

  RNGValue value;
  RNGSeed seed;
  RNGAdvance advance;

  tester.out() << "----------------------------------" << std::endl;
  tester.out() << "some random numbers (all the same)" << std::endl;
  tester.out() << "----------------------------------" << std::endl;

  tester.out() << rng0.comp(value) << std::endl;
  tester.out() << rng1.comp(value) << std::endl;
  rng0.comp(advance) = 1;
  rng1.comp(advance) = 1;
  tester.out() << rng0.comp(value) << std::endl;
  tester.out() << rng1.comp(value) << std::endl;

  ArrayR_t a(layout);

  a = rng0.comp(value) - rng1.comp(value);

  tester.check("same values", sum(a) == 0.0);

  Vector<2, int> strides(1, dom[0].length());

  rng0.comp(seed) = dot(strides, iota(dom));
  rng1.comp(seed) = dot(strides, iota(dom));

  tester.out() << "-----------------------------------" << std::endl;
  tester.out() << "some random numbers (different now)" << std::endl;
  tester.out() << "-----------------------------------" << std::endl;

  tester.out() << rng0.comp(value) << std::endl;
  tester.out() << rng1.comp(value) << std::endl;

  a = rng0.comp(value) - rng1.comp(value);

  tester.check("same values", sum(a) == 0.0);

  rng0.comp(advance) = 20;
  rng1.comp(advance) = 10;

  tester.out() << "------------------------------------------" << std::endl;
  tester.out() << "some random numbers (completely different)" << std::endl;
  tester.out() << "------------------------------------------" << std::endl;

  tester.out() << rng0.comp(value) << std::endl;
  tester.out() << rng1.comp(value) << std::endl;

  a = rng0.comp(value) - rng1.comp(value);

  tester.check("different values", sum(a) != 0.0);

  tester.out() << "------------------------------------------" << std::endl;
  tester.out() << "finally the seeds:" << std::endl;
  tester.out() << "------------------------------------------" << std::endl;

  tester.out() << rng0.comp(seed) << std::endl;
  tester.out() << rng1.comp(seed) << std::endl;

  int ret = tester.results("rngArray");

  Pooma::finalize();

  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: rngArray.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:51 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
