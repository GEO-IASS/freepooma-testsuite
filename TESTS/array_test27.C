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
// Array test 27: compressible operations: compress(), uncompress(),
//                elementsCompressed(), compressionFraction().
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
#include "Engine/RemoteEngine.h"
#include "Array/Array.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> I0(0,2), I1(3,5), I01(2,4);
  Interval<3> I3(6,6,6);
  Loc<3> blocks(2,2,2);
  UniformGridPartition<3> partition(blocks);   
  UniformGridLayout<3> replicated(I3, partition, ReplicatedTag()), 
    distributed(I3, partition, DistributedTag()); 
  Array<3, int, MultiPatch<UniformTag,CompressibleBrick> > a(replicated);
  Array<3, int, CompressibleBrick> b(I3);
  Array<3, int, MultiPatch<UniformTag,Remote<CompressibleBrick> > > c(distributed);

  a = 1;
  b = 1;
  c = 1;
  Pooma::blockAndEvaluate();

  tester.check("a #compressed", elementsCompressed(a), 216L);
  tester.check("a fraction", compressedFraction(a), 1.0);
    
  tester.check("b #compressed", elementsCompressed(b), 216L);
  tester.check("b fraction", compressedFraction(b), 1.0);
    
  tester.check("c #compressed", elementsCompressed(c), 216L);
  tester.check("c fraction", compressedFraction(c), 1.0);
    
  tester.check("bv #compressed", elementsCompressed(b(I1, I1, I0)), 27L);
  tester.check("bv fraction", compressedFraction(b(I1, I1, I0)), 1.0);

  a(4, 5, 1) = 2;
  b(4, 5, 1) = 2;
  c(4, 5, 1) = 2;

  tester.check("a #compressed", elementsCompressed(a), 189L);
  tester.check("a fraction", compressedFraction(a), 0.875, 1e-4);
  
  tester.check("b #compressed", elementsCompressed(b), 0L);
  tester.check("b fraction", compressedFraction(b), 0.0);

  tester.check("c #compressed", elementsCompressed(c), 189L);
  tester.check("c fraction", compressedFraction(c), 0.875, 1e-4);
  
  b(4, 5, 1) = 1;
  compress(b);
    
  tester.check("b #compressed", elementsCompressed(b), 216L);
  tester.check("b fraction", compressedFraction(b), 1.0);
  
  a(I1, I1, I0) = 2;
  c(I1, I1, I0) = 2;
  Pooma::blockAndEvaluate();

  tester.check("a #compressed", elementsCompressed(a), 216L);
  tester.check("a fraction", compressedFraction(a), 1.0);

  tester.check("c #compressed", elementsCompressed(c), 216L);
  tester.check("c fraction", compressedFraction(c), 1.0);
  
  uncompress(a);
  uncompress(b);
  uncompress(c);    

  tester.check("a #compressed", elementsCompressed(a), 0L);
  tester.check("a fraction", compressedFraction(a), 0.0);
     
  tester.check("b #compressed", elementsCompressed(b), 0L);
  tester.check("b fraction", compressedFraction(b), 0.0);
     
  tester.check("c #compressed", elementsCompressed(c), 0L);
  tester.check("c fraction", compressedFraction(c), 0.0);
  
  a(4, 5, 1) = 1;
  c(4, 5, 1) = 1;
  compress(a);
  compress(c);

  tester.check("a #compressed", elementsCompressed(a), 189L);
  tester.check("a fraction", compressedFraction(a), 0.875);
    
  tester.check("av #compressed", elementsCompressed(a(I01, I01, I01)), 23L);
  tester.check("av fraction", compressedFraction(a(I01, I01, I01)), 23./27);
    
  tester.check("c #compressed", elementsCompressed(c), 189L);
  tester.check("c fraction", compressedFraction(c), 0.875);
    
  tester.check("cv #compressed", elementsCompressed(c(I01, I01, I01)), 23L);
  tester.check("cv fraction", compressedFraction(c(I01, I01, I01)), 23./27);

  int ret = tester.results("array_test27");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test27.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
