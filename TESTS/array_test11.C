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
// Array test 11: negative strides.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Layout/UniformGridLayout.h"
#include "Engine/BrickEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "Array/Array.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Array<1, int> a(10), b(10);
  Array<2, int> c(10, 10);
  Loc<2> blocks(5,5);
  UniformGridLayout<2> layout(Interval<2>(10,10), blocks,ReplicatedTag());
  Array<2, int, MultiPatch<UniformTag,Brick> > u(layout);
  int i0, i1;

  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  for (i0 = 0; i0 < 9; i0++)
    a(i0) = b(i0) = i0;
  for (i1 = 0; i1 < 9; i1++)
    for (i0 = 0; i0 < 9; i0++)
      u(i0, i1) = c(i0,i1) = i1+10*i0;
    
  // Make some ranges with negative stride and use them.
  
  Range<1> R(7,3,-2), RR(2,0,-1);
  Range<1> Q(3,7,2);
  
  tester.check(all(a(R)(RR) == b(Q)));
  tester.check(all(c(1, R)(RR) == b(Q) + 10));
  tester.check(all(u(2, R)(RR) == b(Q) + 20));

  int ret = tester.results("array_test11");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test11.cpp,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
