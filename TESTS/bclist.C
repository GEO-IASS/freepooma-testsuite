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
// Particles test: ParticleBCList and ParticleBC
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Particles/ParticleBCList.h"
#include "Particles/PeriodicBC.h"
#include "Particles/ReflectBC.h"
#include "Particles/AbsorbBC.h"
#include "Particles/ReverseBC.h"
#include "Domain/Interval.h"
#include "Layout/DynamicLayout.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "DynamicArray/DynamicArray.h"

#include <iostream>

int main(int argc, char* argv[])
{
  // Initialize POOMA and Tester class.
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": ParticleBCList operations" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // First create some Attributes for BC's to act upon.

  tester.out() << "Creating DynamicArray objects for attributes ..."
               << std::endl;
  Interval<1> D(10);
  int blocks = 4;
  DynamicLayout layout(D,blocks);
#if POOMA_MESSAGING
  typedef MultiPatch< DynamicTag, Remote<Dynamic> > EngineTag_t;
#else
  typedef MultiPatch<DynamicTag,Dynamic> EngineTag_t;
#endif
  DynamicArray<int,EngineTag_t>    a1(layout);
  DynamicArray<long,EngineTag_t>   a2(layout);
  DynamicArray<float,EngineTag_t>  a3(layout);
  DynamicArray<int,EngineTag_t>    a4(layout);
  DynamicArray<double,EngineTag_t> a5(layout);

  // Construct a ParticleBCList to store our ParticleBC's

  tester.out() << "Constructing a ParticleBClist ..." << std::endl;
  ParticleBCList bclist;

  // Create some ParticleBC's

  tester.out() << "Creating some ParticleBC objects and adding to list ..."
               << std::endl;

  // For each BC, we construct the BCType with boundary values.
  // Then we add a ParticleBC with this type to our list, and we provide
  // the subject of the BC (and the object, if different).

  PeriodicBC<int> bc1(12, 18);
  bclist.addBC(a1, bc1);

  ReflectBC<long> bc2(100, 105);
  bclist.addBC(a2, bc2);

  AbsorbBC<float> bc3(0.15, 0.75);
  bclist.addBC(a3, bc3);

  ReverseBC<int> bc5(8, 42);
  bclist.addBC(a4, a5, bc5);

  tester.check(bclist.size() == 4);

  // Print out the list of BC's.

  tester.out() << "Printing contents of the ParticleBCList ... "
               << std::endl;
  tester.out() << bclist << std::endl;

  // Remove some of the ParticleBC's from the ParticleBCList

  tester.out() << "Removing every other ParticleBC from the list ... "
               << std::endl;
  
  int ibc, numbc = bclist.size();
  for (ibc = numbc-1; ibc >= 0; ibc-=2)
    bclist.removeBC(ibc);
  tester.out() << "There are now " << bclist.size() << " boundary conditions."
               << std::endl << std::endl;

  tester.check(bclist.size() == 2);

  // Return resulting error code and shut down POOMA.

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("ParticleBCList operations");
  Pooma::finalize();
  return retval;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: bclist.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
