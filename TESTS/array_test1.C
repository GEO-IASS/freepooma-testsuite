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
// Array test 1: slices.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/UniformGridLayout.h"
#include "Engine/BrickEngine.h"
#include "Engine/CompressibleBrick.h"
#include "Engine/MultiPatchEngine.h"
#include "Array/Array.h"

static bool OK = true;

template<class T>
inline void check(const T &ans, const T &correct, Pooma::Tester &tester)
{
  OK = (OK && (ans == correct));
  tester.check( OK );
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> D(6);
  Interval<3> I3(D,D,D);
  Array<3> a(I3), d(2,3,1);
  Array<2> b(2,3), e(2,3);
  Array<1> b2(2), b3(2);
  Array<2, bool> t(2,3);
  Array<1, bool> t1(2);
  int i0, i1, i2;
  
  Loc<3> blocks(2,2,2);
  UniformGridPartition<3> partition(blocks);   
  UniformGridLayout<3> layout(I3, partition,ReplicatedTag());
  
  Array<3, double, MultiPatch<UniformTag,Brick> > u(layout);
  Array<3, double, MultiPatch<UniformTag,CompressibleBrick> > c(layout);

  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  for (i2 = 0; i2 < 6; i2++)
    for (i1 = 0; i1 < 6; i1++)
      for (i0 = 0; i0 < 6; i0++)
	a(i0,i1,i2) = u(i0,i1,i2) = c(i0,i1,i2) = i2+10*(i1+10*i0);

  b(0,0) = 320; b(0,1) = 322; b(0,2) = 324;
  b(1,0) = 420; b(1,1) = 422; b(1,2) = 424;

  b2(0) = 420; b2(1) = 424;

  d = 0;
  e = 0;

  Interval<1> I(3,4),I2(1,2),I1(0,1);
  Range<1> R(0,4,2);
  
  tester.out() << "At start:" << std::endl;
  tester.out() << "a = " << a << std::endl;
  tester.out() << "b = " << b << std::endl;
  tester.out() << "c = " << c << std::endl;
  tester.out() << "d = " << d << std::endl;
  tester.out() << "e = " << e << std::endl;
  tester.out() << "u = " << u << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.out() << "I = " << I << ", R = " << R << std::endl;
  tester.out() << "I1 = " << I1 << ", I2 = " << I2 << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing t = (b == a(I,2,R)) ..." << std::endl;
  t = (b == a(I,2,R));
  check(all(t), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", t = " << t << std::endl;
  tester.out() << "a(I,2,R) = " << a(I,2,R) << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing t = (b == u(I,2,R)) ..." << std::endl;
  t = (b == u(I,2,R));
  check(all(t), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", t = " << t << std::endl;
  tester.out() << "u(I,2,R) = " << u(I,2,R) << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing t = (b == c(I,2,R)) ..." << std::endl;
  t = (b == c(I,2,R));
  check(all(t), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", t = " << t << std::endl;
  tester.out() << "c(I,2,R) = " << c(I,2,R) << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing b3 = a(I,2,R)(1,R2) - b2 ..." << std::endl;
  Range<1> R2(0,2,2);
  b3 = a(I,2,R)(1,R2) - b2;
  check(all(b3 == 0.0), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", b3 = " << b3 << std::endl;
  tester.out() << "a(I,2,R)(1,R2) = " << a(I,2,R)(1,R2) << std::endl;
  tester.out() << "a(I,2,R)(1,R2) - b3 = " << a(I,2,R)(1,R2) - b3 << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing b3 = u(I,2,R)(1,R2) - b2 ..." << std::endl;
  b3 = u(I,2,R)(1,R2) - b2;
  check(all(b3 == 0.0), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", b3 = " << b3 << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing b3 = c(I,2,R)(1,R2) - b2 ..." << std::endl;
  b3 = c(I,2,R)(1,R2) - b2;
  check(all(b3 == 0.0), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", b3 = " << b3 << std::endl;

  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing t = (b == d(AllDomain<2>(), 0)) ..." << std::endl;
  d(AllDomain<2>(),0) = a(I,2,R);
  t = (b == d(AllDomain<2>(),0));
  check(all(t), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", t = " << t << std::endl;
  
  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing t = (b == e(AllDomain<2>(), 0)) ..." << std::endl;
  e(AllDomain<2>()) = a(I,2,R);
  t = (b == e(AllDomain<2>()));
  check(all(t), true, tester);
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", t = " << t << std::endl;
    
  tester.out() << "-------------------------------------" << std::endl;
  tester.out() << "Testing t1 = (b(R4,1) == b(I1,I2)(R4,0)) ..." << std::endl;
  Range<1> R4(0,1,1);
  t1 = (b(R4,1) == b(I1,I2)(R4,0));
  tester.out() << "Finished: results = " << tester.ok();
  tester.out() << ", t1 = " << t1 << std::endl;
      
  tester.out() << "-------------------------------------" << std::endl;
  int ret = tester.results("array_test1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test1.cpp,v $   $Author: richard $
// $Revision: 1.34 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
