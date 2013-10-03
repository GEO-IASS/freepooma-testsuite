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
// Particles test: KillBC with expressions
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Particles/Particles.h"
#include "Particles/CommonParticleTraits.h"
#include "Particles/UniformLayout.h"
#include "Particles/KillBC.h"
#include "Domain/Interval.h"
#include "Domain/IndirectionList.h"
#include "Layout/DynamicLayout.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "DynamicArray/DynamicArray.h"

#include <iostream>


template <class PT>
class MyParticles : public Particles<PT>
{
public:
  // Useful typedefs

  typedef MyParticles<PT>                        This_t;
  typedef Particles<PT>                          Base_t;
  typedef typename PT::AttributeEngineTag_t      AttributeEngineTag_t;
  typedef typename PT::ParticleLayout_t          ParticleLayout_t;

  // Constructor: set up layouts, register attributes

  MyParticles(const ParticleLayout_t& pl)
    : Particles<PT>(pl)
    {
      this->addAttribute(a1);
      this->addAttribute(a2);
    }

  // List of attributes; we'll just make them public data members here,
  // you could also provide access via methods.

  DynamicArray< double, AttributeEngineTag_t >  a1;
  DynamicArray< double, AttributeEngineTag_t >  a2;
};


int main(int argc, char* argv[])
{
  // Initialize POOMA and Tester class.
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": KillBC with expressions" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // First create a Particles object with some Attributes for BC's to act upon.

  tester.out() << "Creating Particles object with DynamicArray attributes ..."
               << std::endl;
  UniformLayout pl(Pooma::contexts());
#if POOMA_MESSAGING
  MyParticles<MPRemoteDynamicUniform> P(pl);
#else
  MyParticles<MPDynamicUniform> P(pl);
#endif

  if (Pooma::context() == 0)
    P.create(10,0,false);
  P.sync(P.a1);

  // Initialize the arrays with scalars.
  // Block since we're starting scalar code.

  Pooma::blockAndEvaluate();
  
  tester.out() << "Initializing DynamicArray objects ..."
               << std::endl;
  int i;
  for (i=0; i < P.size(); ++i) {
    P.a1(i) = 0.1 * i;
    P.a2(i) = 0.25 * i - 1.5;
  }

  tester.out() << "Initialization complete:" << std::endl;
  tester.out() << "  a1 = " << P.a1 << std::endl;
  tester.out() << "  a2 = " << P.a2 << std::endl;
  tester.out() << "  a1*a1+a2*a2 = " << P.a1 * P.a1 + P.a2 * P.a2
               << std::endl;

  // Create a KillBC

  tester.out() << "Creating a Particle KillBC object ..."
               << std::endl;

  // For each BC, we construct the BCType with boundary values.
  // Then we add a ParticleBC with this type to our list, and we provide
  // the subject of the BC (and the object, if different).
  // For the KillBC, object must be the Particles object itself.

  KillBC<double> bc1(0.0, 0.8);
  P.addBoundaryCondition(P.a1 * P.a1 + P.a2 * P.a2, P, bc1);

  // Apply boundary condition and display the results

  tester.out() << "Applying the boundary conditions ..." << std::endl;
  tester.out() << "Before BC's, Particles = " << P << std::endl;
  P.applyBoundaryConditions();
  tester.out() << "After BC's, Particles = " << P << std::endl;
  P.performDestroy();
  Pooma::blockAndEvaluate();
  tester.out() << "Status after applying BC: " << std::endl;
  tester.out() << "  a1 = " << P.a1 << std::endl;
  tester.out() << "  a2 = " << P.a2 << std::endl;

  tester.check(P.size() == 5);

  // Let's also try a KillBC on a free-standing DynamicArray.

  tester.out() << "Creating a free-standing DynamicArray ..." << std::endl;
#if POOMA_MESSAGING
  DynamicArray< Vector<2,int>, MultiPatch< DynamicTag, Remote<Dynamic> > > a3;
#else
  DynamicArray< Vector<2,int>, MultiPatch<DynamicTag,Dynamic> > a3;
#endif

  Interval<1> empty;
  DynamicLayout layout(empty, Pooma::contexts());
  a3.initialize(layout);
  int npc = 20 / Pooma::contexts();
  int rem = 20 % Pooma::contexts();
  if (Pooma::context() < rem) npc++;
  a3.create(npc);
  a3.layout().sync();

  Pooma::blockAndEvaluate();
  for (i=0; i<a3.domain().size(); ++i)
    a3(i) = Vector<2,int>(i,2*i+1);

  tester.out() << "Initialization complete." << std::endl;
  tester.out() << "a3 = " << a3 << std::endl;

  // Now construct a KillBC for this DynamicArray and apply it

  tester.out() << "Creating a DynamicArray KillBC object ..." << std::endl;
  KillBC< Vector<2,int> > bc2(Vector<2,int>(2,2), Vector<2,int>(24,24));
  ParticleBCItem* killbc2 = bc2.create(a3);

  tester.out() << "Applying the boundary condition ..." << std::endl;
  killbc2->applyBoundaryCondition();
  a3.layout().sync();
  Pooma::blockAndEvaluate();
  tester.out() << "Status after applying BC:" << std::endl;
  tester.out() << "a3 = " << a3 << std::endl;

  tester.check(a3.domain().size() == 10);

  // Delete the ParticleBC that we created

  delete killbc2;

  // Return resulting error code and shut down POOMA.

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("KillBC with expression");
  Pooma::finalize();
  return retval;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: bctest3.cpp,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
