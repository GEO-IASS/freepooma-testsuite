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
// array_test18.cpp multipatch expression tests.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/UniformGridLayout.h"
#include "Layout/GuardLayers.h"
#include "Engine/BrickEngine.h"
#include "Engine/RemoteEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "Tiny/Vector.h"
#include "Array/Array.h"
#include "Array/tests/ExpressionTest.h"

#include <iostream>

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] 
	       << ": Tests of expressions with multipatch." 
	       << std::endl;
  tester.out() << "------------------------------------------------" 
	       << std::endl;

  int size = 30;

  int from = 0;
  int to = size-1;
  int fromInterior = 0;
  int toInterior = size-1;

  int loc1 = 4;
  int loc2 = 14;
  int loc3 = 22;

  Interval<1> dom(from,to);
  Interval<1> I(fromInterior,toInterior);

  Interval<1> domain(size);

  UniformGridPartition<1> partition(Loc<1>(10),GuardLayers<1>(1));
  UniformGridLayout<1> layout(domain,partition,ReplicatedTag());

  Array<1,double,MultiPatch<UniformTag,Brick> >
    a1(layout), a2(layout), a3(layout), a4(layout), initial(layout);

  initial = 0.0;

  Pooma::blockAndEvaluate();

  initial(loc1) = 2.0;
  initial(loc2) = 3.0;
  initial(loc3) = 4.0;

  test1(tester, 1, a1, a2, a3, a4, initial, I);
  test2(tester, 2, a1, a2, a3, a4, initial, I);
  test3(tester, 3, a1, a2, a3, a4, initial, I);
  test4(tester, 4, a1, a2, a3, a4, initial, I);

  Array<1,Vector<2,double>,MultiPatch<UniformTag,Brick> >
    av1(layout),av2(layout),av3(layout),av4(layout),
    initialv(layout);

  initialv = Vector<2,double>(0.0,0.0);

  Pooma::blockAndEvaluate();

  initialv(loc1) = Vector<2,double>(2.0,3.0);
  initialv(loc2) = Vector<2,double>(3.0,-1.0);
  initialv(loc3) = Vector<2,double>(4.0,-5.0);

  test5(tester, 5, av1, av2, av3, av4, initialv, I);
  test1(tester, 6, av1, av2, av3, av4, initialv, I);
  test4(tester, 7, av1, av2, av3, av4, initialv, I);

  UniformGridLayout<1> layoutr(domain, partition, DistributedTag());

  Array<1, double, MultiPatch<UniformTag, Remote<Brick> > >
    ar1(layoutr), ar2(layoutr), ar3(layoutr), ar4(layoutr);

  test1(tester, 8,  ar1, ar2, ar3, ar4, initial, I);
  test2(tester, 9,  ar1, ar2, ar3, ar4, initial, I);
  test3(tester, 10, ar1, ar2, ar3, ar4, initial, I);
  test4(tester, 11, ar1, ar2, ar3, ar4, initial, I);

  tester.out() << "------------------------------------------------"
	       << std::endl;

  int retval = tester.results("array_test18");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test18.cpp,v $   $Author: richard $
// $Revision: 1.24 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
