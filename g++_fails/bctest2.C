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
#include "Tiny/Vector.h"
#include "Tiny/Tensor.h"
#include "Tulip/RemoteProxy.h"

#include <iostream>


#if POOMA_MESSAGING
  typedef MultiPatch< DynamicTag, Remote<Dynamic> > EngineTag_t;
#else
  typedef MultiPatch<DynamicTag,Dynamic> EngineTag_t;
#endif

template <int Dim, class T>
bool myDiff(const Vector<Dim,T>& a, const Vector<Dim,T>& b)
{
  Vector<Dim,T> diff = a - b;
  return (norm(diff) < 1.0e-8);
}


bool checkResults(DynamicArray<Tensor<2,int>,EngineTag_t> a1,
                  DynamicArray<Vector<3,int>,EngineTag_t> a2,
                  DynamicArray<Vector<3,double>,EngineTag_t> a3)
{
  // check against hard-coded results
  bool result = true;

  result = result && (a1.read(0) == Tensor<2,int>(0,61,2,3));
  result = result && (a1.read(1) == Tensor<2,int>(10,71,12,13));
  result = result && (a1.read(2) == Tensor<2,int>(20,21,22,23));
  result = result && (a1.read(3) == Tensor<2,int>(30,31,32,33));
  result = result && (a1.read(4) == Tensor<2,int>(40,41,42,43));
  result = result && (a1.read(5) == Tensor<2,int>(50,51,52,53));
  result = result && (a1.read(6) == Tensor<2,int>(60,61,62,63));
  result = result && (a1.read(7) == Tensor<2,int>(70,71,72,73));
  result = result && (a1.read(8) == Tensor<2,int>(80,21,82,83));
  result = result && (a1.read(9) == Tensor<2,int>(90,31,92,93));

  result = result && (a2.read(0) == Vector<3,int>(0,2,12));
  result = result && (a2.read(1) == Vector<3,int>(5,7,9));
  result = result && (a2.read(2) == Vector<3,int>(10,12,14));
  result = result && (a2.read(3) == Vector<3,int>(15,17,19));
  result = result && (a2.read(4) == Vector<3,int>(20,22,24));
  result = result && (a2.read(5) == Vector<3,int>(25,27,29));
  result = result && (a2.read(6) == Vector<3,int>(30,32,34));
  result = result && (a2.read(7) == Vector<3,int>(35,37,37));
  result = result && (a2.read(8) == Vector<3,int>(40,42,32));
  result = result && (a2.read(9) == Vector<3,int>(45,47,27));

  result = result && myDiff(a3.read(0),Vector<3,double>(-1.5,-1.5,1.5));
  result = result && myDiff(a3.read(1),Vector<3,double>(-1.25,-1.25,-1.25));
  result = result && myDiff(a3.read(2),Vector<3,double>(-1.0,-1.0,-1.0));
  result = result && myDiff(a3.read(3),Vector<3,double>(-0.75,-0.75,-0.75));
  result = result && myDiff(a3.read(4),Vector<3,double>(-0.5,-0.5,-0.5));
  result = result && myDiff(a3.read(5),Vector<3,double>(-0.25,-0.25,-0.25));
  result = result && myDiff(a3.read(6),Vector<3,double>(0.0,0.0,0.0));
  result = result && myDiff(a3.read(7),Vector<3,double>(0.25,0.25,-0.25));
  result = result && myDiff(a3.read(8),Vector<3,double>(0.5,0.5,-0.5));
  result = result && myDiff(a3.read(9),Vector<3,double>(0.75,0.75,-0.75));

  return result;
}


int main(int argc, char* argv[])
{
  // Initialize POOMA and Tester class.
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": ParticleCompBC operations" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // First create some Attributes for BC's to act upon.

  tester.out() << "Creating DynamicArray objects for attributes ..."
               << std::endl;
  Interval<1> D(10);
  int blocks = 4;
  DynamicLayout layout(D,blocks);
  DynamicArray<Tensor<2,int>,EngineTag_t>    a1(layout);
  DynamicArray<Vector<3,int>,EngineTag_t>    a2(layout);
  DynamicArray<Vector<3,double>,EngineTag_t> a3(layout);

  // Initialize the arrays.
  // Block since we're starting scalar code.

  Pooma::blockAndEvaluate();
  
  tester.out() << "Initializing DynamicArray objects ..."
               << std::endl;
  int i;
  for (i=0; i < D.size(); ++i) {
    a1(i) = Tensor<2,int>(10*i,10*i+1,10*i+2,10*i+3);
    a2(i) = Vector<3,int>(5*i,5*i+2,5*i+4);
    a3(i) = Vector<3,double>(0.25*i-1.5);
  }
  tester.out() << "Initialization complete:" << std::endl;
  tester.out() << "  a1 = " << a1 << std::endl;
  tester.out() << "  a2 = " << a2 << std::endl;
  tester.out() << "  a3 = " << a3 << std::endl;

  // Construct a ParticleBCList to store our ParticleBC's

  tester.out() << "Constructing a ParticleBClist ..." << std::endl;
  ParticleBCList bclist;

  // Create some ParticleBC's

  tester.out() << "Creating some ParticleBC objects and adding to list ..."
               << std::endl;

  // For each BC, we construct the BCType with boundary values.
  // Then we add a ParticleBC with this type to our list, and we provide
  // the subject of the BC (and the object, if different).

  PeriodicBC<int> bc1(20, 80);
  ParticleCompBC< 2, PeriodicBC<int> > cbc1(bc1, 1, 0);
  bclist.addBC(a1, cbc1);

  ReverseBC<int> bc2(8, 38);
  ParticleCompBC< 1, ReverseBC<int> > cbc2(bc2, 2);
  bclist.addBC(a2, a3, cbc2);

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
  }

  bool success = checkResults(a1,a2,a3);
  tester.set(success);

  // Return resulting error code and shut down POOMA.

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("ParticleCompBC operations");
  Pooma::finalize();
  return retval;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: bctest2.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
