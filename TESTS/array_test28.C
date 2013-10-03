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
// Array test 28: remote assignment.
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
#include "Tiny/Vector.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<3> I3(6,6,6);
  Array<3> a0(I3), b0(I3);
  Array<3, double, Remote<CompressibleBrick> > a1(I3), b1(I3);
  
  Loc<3> blocks(1,1,2);
  UniformGridPartition<3> partition(blocks);   
  UniformGridLayout<3> layout(I3, partition, DistributedTag());
  
  Array<3, double, MultiPatch<UniformTag,Remote<Brick> > > 
    a2(layout), b2(layout);
  Array<3, double, MultiPatch<UniformTag,Remote<CompressibleBrick> > > 
    a3(layout), b3(layout);

  b0 = 0.0;
  b1 = 1.0;
  b2 = 2.0;
  b3 = 3.0;

  a0 = b0; tester.check("Brick                      = Brick\n\t",
			all(a0 == 0.0));
  a1 = b1; tester.check("Remote<CBrick>             = Remote<CBrick>\n\t",
			all(a1 == 1.0));
  a2 = b2; tester.check("MultiPatch<Remote<Brick>>  = MultiPatch<Remote<Brick>>\n\t",
			all(a2 == 2.0));
  a3 = b3; tester.check("MultiPatch<Remote<CBrick>> = MultiPatch<Remote<CBrick>>\n\t",
			all(a3 == 3.0));

  a0 = b1; tester.check("Brick                      = Remote<CBrick>\n\t",
			all(a0 == 1.0));
  a1 = b2; tester.check("Remote<CBrick>             = MultiPatch<Remote<Brick>>\n\t",
			all(a1 == 2.0));
  a2 = b3; tester.check("MultiPatch<Remote<Brick>>  = MultiPatch<Remote<CBrick>>\n\t",
			all(a2 == 3.0));
  a3 = b0; tester.check("MultiPatch<Remote<CBrick>> = Brick\n\t",
			all(a3 == 0.0));

  a0 = b2; tester.check("Brick                      = MultiPatch<Remote<Brick>>\n\t",
			all(a0 == 2.0));
  a1 = b3; tester.check("Remote<CBrick>             = MultiPatch<Remote<CBrick>>\n\t",
			all(a1 == 3.0));
  a2 = b0; tester.check("MultiPatch<Remote<Brick>>  = Brick\n\t",
			all(a2 == 0.0));
  a3 = b1; tester.check("MultiPatch<Remote<CBrick>> = Remote<CBrick>\n\t",
			all(a3 == 1.0));

  a0 = b3; tester.check("Brick                      = MultiPatch<Remote<CBrick>>\n\t",
			all(a0 == 3.0));
  a1 = b0; tester.check("Remote<CBrick>             = Brick\n\t",
			all(a1 == 0.0));
  a2 = b1; tester.check("MultiPatch<Remote<Brick>>  = Remote<CBrick>\n\t",
			all(a2 == 1.0));
  a3 = b2; tester.check("MultiPatch<Remote<CBrick>> = MultiPatch<Remote<Brick>>\n\t",
			all(a3 == 2.0));

  Array<3, Vector<2, double>, Remote<Brick> > a4(I3);

  a4 = Vector<2, double>(1.0, 2.0);

  tester.check("a4.comp(1)", all(a4.comp(1) == 2.0));

  int ret = tester.results( "array_test28" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test28.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
