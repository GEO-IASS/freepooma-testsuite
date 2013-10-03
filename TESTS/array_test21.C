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
// Array test 21: Multi-patch engines with a GridLayout.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/GridLayout.h"
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
  int i0, i1, i2;
  
  Loc<3> blocks(2,2,2);
  UniformGridPartition<3> partition(blocks);   
  GridLayout<3> layout(I3, partition, ReplicatedTag());
  
  tester.out() << "Created GridLayout<3> = " << layout << std::endl;

  Array<3, double, MultiPatch<GridTag, Brick> > u(layout);
  Array<3, double, MultiPatch<GridTag, CompressibleBrick> > c(layout);

  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  for (i2 = 0; i2 < 6; i2++)
    for (i1 = 0; i1 < 6; i1++)
      for (i0 = 0; i0 < 6; i0++)
	a(i0,i1,i2) = u(i0,i1,i2) = c(i0,i1,i2) = i2+10*(i1+10*i0);

  b(0,0) = 320; b(0,1) = 322; b(0,2) = 324;
  b(1,0) = 420; b(1,1) = 422; b(1,2) = 424;

  b2(0) = 420; b2(1) = 424;

  tester.out() << "Created Array<3> u = " << u << std::endl;
  tester.out() << "Created Array<3> c = " << c << std::endl;

  Interval<1> I(3,4);
  Range<1> R(0,4,2);

  tester.out() << "u slice = ";
  tester.out() << u(I,2,R);
  tester.out() << std::endl;

  t = (b == a(I,2,R));
  check(all(t), true, tester);

  t = (b == u(I,2,R));
  check(all(t), true, tester);

  t = (b == c(I,2,R));
  check(all(t), true, tester);

  Range<1> R2(0,2,2);
  b3 = a(I,2,R)(1,R2) - b2;
  check(all(b3 == 0.0), true, tester);

  b3 = u(I,2,R)(1,R2) - b2;
  check(all(b3 == 0.0), true, tester);

  b3 = c(I,2,R)(1,R2) - b2;
  check(all(b3 == 0.0), true, tester);

  d(AllDomain<2>(),0) = a(I,2,R);
  t = (b == d(AllDomain<2>(),0));
  check(all(t), true, tester);
  
  e(AllDomain<2>()) = a(I,2,R);
  t = (b == e(AllDomain<2>()));
  check(all(t), true, tester);
    
  int ret = tester.results("array_test21");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test21.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
