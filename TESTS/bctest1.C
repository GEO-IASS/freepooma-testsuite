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
#include "Tulip/RemoteProxy.h"

#include <iostream>
#include <cmath>
#include <cfloat>


#if POOMA_MESSAGING
  typedef MultiPatch< DynamicTag, Remote<Dynamic> > EngineTag_t;
#else
  typedef MultiPatch<DynamicTag,Dynamic> EngineTag_t;
#endif

template <class T>
bool myDiff(const T& a, const T& b)
{
  return (std::abs(a-b) < FLT_EPSILON);
}

template <class T>
bool myDiff(const RemoteProxy<T>& a, const T& b)
{
  return (std::abs(a.value()-b) < FLT_EPSILON);
}


bool checkResults(DynamicArray<int,EngineTag_t> a1,
                  DynamicArray<long,EngineTag_t> a2,
                  DynamicArray<float,EngineTag_t> a3,
                  DynamicArray<int,EngineTag_t> a4,
                  DynamicArray<double,EngineTag_t> a5)
{
  // check against hard-coded results
  bool result = true;

  result = result && a1(0) == 16;
  result = result && a1(1) == 17;
  result = result && a1(2) == 12;
  result = result && a1(3) == 13;
  result = result && a1(4) == 14;
  result = result && a1(5) == 15;
  result = result && a1(6) == 16;
  result = result && a1(7) == 17;
  result = result && a1(8) == 18;
  result = result && a1(9) == 13;

  result = result && a2(0) == 100L;
  result = result && a2(1) == 101L;
  result = result && a2(2) == 102L;
  result = result && a2(3) == 103L;
  result = result && a2(4) == 104L;
  result = result && a2(5) == 105L;
  result = result && a2(6) == 104L;
  result = result && a2(7) == 103L;
  result = result && a2(8) == 102L;
  result = result && a2(9) == 101L;

  result = result && myDiff(a3(0), 0.15f);
  result = result && myDiff(a3(1), 0.15f);
  result = result && myDiff(a3(2), 0.2f);
  result = result && myDiff(a3(3), 0.3f);
  result = result && myDiff(a3(4), 0.4f);
  result = result && myDiff(a3(5), 0.5f);
  result = result && myDiff(a3(6), 0.6f);
  result = result && myDiff(a3(7), 0.7f);
  result = result && myDiff(a3(8), 0.75f);
  result = result && myDiff(a3(9), 0.75f);

  result = result && a4(0) == 16;
  result = result && a4(1) == 11;
  result = result && a4(2) == 10;
  result = result && a4(3) == 15;
  result = result && a4(4) == 20;
  result = result && a4(5) == 25;
  result = result && a4(6) == 30;
  result = result && a4(7) == 35;
  result = result && a4(8) == 40;
  result = result && a4(9) == 39;

  result = result && myDiff(a5(0), 1.5);
  result = result && myDiff(a5(1), 1.25);
  result = result && myDiff(a5(2), -1.0);
  result = result && myDiff(a5(3), -0.75);
  result = result && myDiff(a5(4), -0.5);
  result = result && myDiff(a5(5), -0.25);
  result = result && myDiff(a5(6), 0.0);
  result = result && myDiff(a5(7), 0.25);
  result = result && myDiff(a5(8), 0.5);
  result = result && myDiff(a5(9), -0.75);

  return result;
}


int main(int argc, char* argv[])
{
  // Initialize POOMA and Tester class.
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": ParticleBC operations" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // First create some Attributes for BC's to act upon.

  tester.out() << "Creating DynamicArray objects for attributes ..."
               << std::endl;
  Interval<1> D(10);
  int blocks = 4;
  DynamicLayout layout(D,blocks);
  DynamicArray<int,EngineTag_t>    a1(layout);
  DynamicArray<long,EngineTag_t>   a2(layout);
  DynamicArray<float,EngineTag_t>  a3(layout);
  DynamicArray<int,EngineTag_t>    a4(layout);
  DynamicArray<double,EngineTag_t> a5(layout);

  // Initialize the arrays with scalars.
  // Block since we're starting scalar code.

  Pooma::blockAndEvaluate();
  
  tester.out() << "Initializing DynamicArray objects ..."
               << std::endl;
  int i;
  for (i=0; i < D.size(); ++i) {
    a1(i) = 10 + i;
    a2(i) = 100 + i;
    a3(i) = 0.1 * i;
    a4(i) = 5 * i;
    a5(i) = 0.25 * i - 1.5;
  }
  tester.out() << "Initialization complete:" << std::endl;
  tester.out() << "  a1 = " << a1 << std::endl;
  tester.out() << "  a2 = " << a2 << std::endl;
  tester.out() << "  a3 = " << a3 << std::endl;
  tester.out() << "  a4 = " << a4 << std::endl;
  tester.out() << "  a5 = " << a5 << std::endl;

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

  // Print out the list of BC's.

  tester.out() << "Printing contents of the ParticleBCList ... "
               << std::endl;
  tester.out() << bclist << std::endl;

  // Apply each boundary condition and display the results

  tester.out() << "Applying the boundary conditions ..." << std::endl;
  int ibc, numbc = bclist.size();
  tester.out() << "There are now " << numbc << " boundary conditions."
               << std::endl << std::endl;
  for (ibc = 0; ibc < numbc; ++ibc) {
    bclist(ibc)->applyBoundaryCondition();
    Pooma::blockAndEvaluate();
    tester.out() << "Status after applying BC #" << ibc+1 << ": " << std::endl;
    tester.out() << "  a1 = " << a1 << std::endl;
    tester.out() << "  a2 = " << a2 << std::endl;
    tester.out() << "  a3 = " << a3 << std::endl;
    tester.out() << "  a4 = " << a4 << std::endl;
    tester.out() << "  a5 = " << a5 << std::endl;
  }

  bool success = checkResults(a1,a2,a3,a4,a5);
  tester.set(success);

  // Return resulting error code and shut down POOMA.

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("ParticleBC operations");
  Pooma::finalize();
  return retval;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: bctest1.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
