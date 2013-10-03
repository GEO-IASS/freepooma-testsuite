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
// Array test 7: assignment.
//-----------------------------------------------------------------------------

// Bring in Pooma array machinery.

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/UniformGridLayout.h"
#include "Engine/BrickEngine.h"
#include "Engine/CompressibleBrick.h"
#include "Engine/MultiPatchEngine.h"
#include "Array/Array.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<3> I3(6,6,6);
  Array<3> a0(I3), b0(I3);
  Array<3, double, CompressibleBrick> a1(I3), b1(I3);
  
  Loc<3> blocks(2,2,2);
  UniformGridPartition<3> partition(blocks);   
  UniformGridLayout<3> layout(I3, partition, ReplicatedTag());
  
  Array<3, double, MultiPatch<UniformTag,Brick> > a2(layout), b2(layout);
  Array<3, double, MultiPatch<UniformTag,CompressibleBrick> > a3(layout), 
    b3(layout);

  b0 = 0.0;
  b1 = 1.0;
  b2 = 2.0;
  b3 = 3.0;
  
  tester.check("b3 #compressed", elementsCompressed(b3), 216L);
  tester.check("b3 fraction", compressedFraction(b3), 1.0);

  a0 = b0; tester.check(all(a0 == 0.0));
  a1 = b1; tester.check(all(a1 == 1.0));
  a2 = b2; tester.check(all(a2 == 2.0));
  a3 = b3; tester.check(all(a3 == 3.0));

  a0 = b1; tester.check(all(a0 == 1.0));
  a1 = b2; tester.check(all(a1 == 2.0));
  a2 = b3; tester.check(all(a2 == 3.0));
  a3 = b0; tester.check(all(a3 == 0.0));

  a0 = b2; tester.check(all(a0 == 2.0));
  a1 = b3; tester.check(all(a1 == 3.0));
  a2 = b0; tester.check(all(a2 == 0.0));
  a3 = b1; tester.check(all(a3 == 1.0));

  a0 = b3; tester.check(all(a0 == 3.0));
  a1 = b0; tester.check(all(a1 == 0.0));
  a2 = b1; tester.check(all(a2 == 1.0));
  a3 = b2; tester.check(all(a3 == 2.0));

  int ret = tester.results( "array_test7" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test7.cpp,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
