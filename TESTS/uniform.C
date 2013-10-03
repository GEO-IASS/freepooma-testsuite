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
// Particles test: UniformLayout with Particles subclass
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Particles/Particles.h"
#include "Particles/UniformLayout.h"
#include "Domain/Interval.h"
#include "Layout/DynamicLayout.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "DynamicArray/DynamicArray.h"
#include "Tiny/Vector.h"

//-----------------------------------------------------------------------------
// A traits class for a Particles object
//-----------------------------------------------------------------------------

template<class EngineTag>
struct PTraits
{
  // The type of engine to use in the attributes

  typedef EngineTag AttributeEngineTag_t;

  // The type of particle layout to use

  typedef UniformLayout ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// A Particles subclass, that defines a few attributes
//-----------------------------------------------------------------------------

template<class PT>
class Molecule : public Particles<PT>
{
public:
  // Useful typedefs to get from the base class

  typedef Particles<PT>                          Base_t;
  typedef typename Base_t::AttributeEngineTag_t  AttributeEngineTag_t;
  typedef typename Base_t::ParticleLayout_t      ParticleLayout_t;
  typedef double                                 AxisType_t;
  typedef Vector<2, AxisType_t>                  PointType_t;

  // Constructor: set up layouts, register attributes

  Molecule(const ParticleLayout_t &pl)
    : Particles<PT>(pl)
    {
      this->addAttribute(pos);
      this->addAttribute(mom);
      this->addAttribute(charge);
    }

  // List of attributes; we'll just make them public data members here,
  // you could also provide access via methods.

  DynamicArray< PointType_t, AttributeEngineTag_t >  pos;
  DynamicArray< PointType_t, AttributeEngineTag_t >  mom;
  DynamicArray< AxisType_t,  AttributeEngineTag_t >  charge;
};


//-----------------------------------------------------------------------------
// Typedefs for what we will compute
//-----------------------------------------------------------------------------

// Dimensionality of this problem

enum { PDim = 2 };

// Engine tag type for attributes

#if POOMA_MESSAGING
typedef MultiPatch< DynamicTag, Remote<Dynamic> > AttrEngineTag_t;
#else
typedef MultiPatch<DynamicTag, Dynamic> AttrEngineTag_t;
#endif

// The particle traits class we'll use

typedef PTraits<AttrEngineTag_t> PTraits_t;

// The particle layout type

typedef PTraits_t::ParticleLayout_t PLayout_t;


//-----------------------------------------------------------------------------
// The main routine for this test code
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Particles with uniform layout" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // Create a UniformLayout object for our use

  tester.out() << "Creating UniformLayout object ..." << std::endl;
  PLayout_t layout(5);		// create it with 5 patches per attribute

  // Create a Particles object, using our special subclass

  tester.out() << "Creating Molecule object ..." << std::endl;
  Molecule<PTraits_t> mol(layout);

  tester.out() << "Molecule created; initially, num attributes = ";
  tester.out() << mol.attributes() << ", num particles = " << mol.size();
  tester.out() << ", global patches = " << mol.attributeLayout().sizeGlobal();
  tester.out() << ", local patches = " << mol.attributeLayout().sizeLocal();
  tester.out() << std::endl;

  tester.check(mol.attributes() == 3);
  tester.check(mol.size() == 0);
  tester.check(mol.attributeLayout().sizeGlobal() == 5);

  // Create some particles, and then renumber.

  int createnum = 10;
  tester.out() << "Creating " << createnum << " particles "
               << "on context 0, patch 0 ..." << std::endl;
  if (Pooma::context() == 0)
    mol.create(createnum,0);
  else 
    mol.create(0);

  tester.out() << "Created (not yet initialized) ... attrib layout:\n";
  tester.out() << mol.attributeLayout() << std::endl;

  tester.check(mol.size() == 10);

  // Block before serial code.

  Pooma::blockAndEvaluate();

  // Initialize the positions.

  tester.out() << "Initializing values ..." << std::endl;
  int i;
  for (i = 0; i < createnum; ++i)
    mol.pos(i) = i;
  mol.mom = mol.pos * 100.0;
  mol.charge = 3.3;

  tester.out() << "Contents of particles:" << std::endl;
  tester.out() << mol << std::endl;

  // Sync the particles now that we've changed positions.
  // Note that in the case of a UniformLayout, the position values 
  // are irrelevant.  The layout will simply try to put an equal
  // number of particles in each patch.

  tester.out() << "Syncing particles ..." << std::endl;
  mol.sync(mol.pos);
  tester.out() << "After sync, contents of particles:" << std::endl;
  tester.out() << mol << std::endl;

  // Add more particles, and then resync

  tester.out() << "Adding " << createnum << " more particles "
               << "to last local patch of context " << Pooma::contexts()-1
               << " ..." << std::endl;
  if (Pooma::context() == Pooma::contexts()-1)
    mol.create(createnum);
  else
    mol.create(0);

  tester.check(mol.size() == 20);

  tester.out() << "Initializing attribute values for new particles ... "
               << std::endl;

  for (i = 0; i < createnum; ++i)
    mol.pos(i+createnum) = mol.pos(i);
  mol.mom = mol.pos * 50.0;
  mol.charge = 6.6;

  tester.out() << "Contents of particles:" << std::endl;
  tester.out() << mol << std::endl;

  tester.out() << "Syncing particles again ..." << std::endl;
  mol.sync(mol.pos);
  tester.out() << "After sync, contents of particles:" << std::endl;
  tester.out() << mol << std::endl;

  // Return resulting error code and exit

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("Particles with uniform layout");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: uniform.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
