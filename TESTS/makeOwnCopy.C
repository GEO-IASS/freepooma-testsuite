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
// makeOwnCopy - test makeOwnCopy() function on several engines.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/UMPArrays.h"
#include "Engine/RemoteEngine.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  int i;
  
  // Create the total domain.
  
  Interval<1> domain(12);
  
  // Create the block sizes.
  
  Loc<1> blocks(3);

  // Create the layouts.
  
  UniformGridLayout<1> layout(domain, blocks, ReplicatedTag());
  
  // Make some UMP arrays and fill them.
  
  Array<1, double, MultiPatch<UniformTag,Brick> > a(layout);
  Array<1, double, MultiPatch<UniformTag,Brick> > b(a);

  b.makeOwnCopy();
  a = 2;
  b = 3;

  Pooma::blockAndEvaluate();

  tester.check("multipatch make own copy", sum(b - a - 1) == 0);

  tester.out() << a << b << std::endl;

  // dynamic array:

  Array<1, double, Dynamic > ad(domain);
  Array<1, double, Dynamic > bd(ad);

  bd.makeOwnCopy();
  ad = 2;
  bd = 3;

  Pooma::blockAndEvaluate();

  tester.check("dynamic make own copy", sum(bd - ad - 1) == 0);

  tester.out() << ad << bd << std::endl;

#if POOMA_MESSAGING

  // Create the layouts.
  
  UniformGridLayout<1> layout2(domain, blocks, DistributedTag());
  
  // Make some UMP arrays and fill them.
  
  Array<1, double, MultiPatch<UniformTag, Remote<Brick> > > a2(layout2);
  Array<1, double, MultiPatch<UniformTag, Remote<Brick> > > b2(a2);

  b2.makeOwnCopy();
  a2 = 2;
  b2 = 3;

  Pooma::blockAndEvaluate();

  tester.check("remote multipatch make own copy", sum(b2 - a2 - 1) == 0);

  tester.out() << a2 << b2 << std::endl;

  // remote dynamic array:

  Array<1, double, Remote<Dynamic> > ard(domain);
  Array<1, double, Remote<Dynamic> > brd(ard);

  brd.makeOwnCopy();
  ard = 2;
  brd = 3;

  Pooma::blockAndEvaluate();

  tester.check("remote dynamic make own copy", sum(brd - ard - 1) == 0);

  tester.out() << ard << brd << std::endl;

#endif // POOMA_MESSAGING
  
  int ret = tester.results("makeOwnCopy");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: makeOwnCopy.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
