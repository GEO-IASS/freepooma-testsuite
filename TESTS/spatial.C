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
// Particles test: Spatial layout with Particles subclass
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Particles/Particles.h"
#include "Particles/SpatialLayout.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/UniformGridLayout.h"
#include "Layout/DynamicLayout.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "Engine/BrickEngine.h"
#include "DynamicArray/DynamicArray.h"
#include "Field/Field.h"
#include "Field/FieldEngine/FieldEngine.h"
#include "Field/Mesh/UniformRectilinearMesh.h"

#include <iostream>
#include <cstdlib>


//-----------------------------------------------------------------------------
// A traits class for a Particles object
//-----------------------------------------------------------------------------

template <class EngineTag, class Mesh, class FL>
struct PTraits
{
  // The type of engine to use in the attributes

  typedef EngineTag AttributeEngineTag_t;

  // The type of particle layout to use

  typedef SpatialLayout< Mesh, FL > ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// A Particles subclass, that defines a few attributes
//-----------------------------------------------------------------------------

template <class PT>
class Molecule : public Particles<PT>
{
public:
  // Useful typedefs to get from the base class

  typedef Particles<PT>                          Base_t;
  typedef typename Base_t::AttributeEngineTag_t  AttributeEngineTag_t;
  typedef typename Base_t::ParticleLayout_t      ParticleLayout_t;
  typedef typename ParticleLayout_t::AxisType_t  AxisType_t;
  typedef typename ParticleLayout_t::PointType_t PointType_t;

  // Useful enums to get from the base class

  enum { dimensions = ParticleLayout_t::dimensions };

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

  DynamicArray<PointType_t,AttributeEngineTag_t>  pos;
  DynamicArray<PointType_t,AttributeEngineTag_t>  mom;
  DynamicArray<AxisType_t,AttributeEngineTag_t>   charge;
};


//-----------------------------------------------------------------------------
// Typedefs for what we will compute
//-----------------------------------------------------------------------------

// Dimensionality of this problem

static const int PDim = 2;

// Engine tag type for attributes

#if POOMA_MESSAGING
typedef MultiPatch< DynamicTag, Remote<Dynamic> > AttrEngineTag_t;
#else
typedef MultiPatch<DynamicTag, Dynamic> AttrEngineTag_t;
#endif

// Mesh type

typedef UniformRectilinearMesh<PDim> Mesh_t;

// Field type

#if POOMA_MESSAGING
typedef Field< Mesh_t, int, MultiPatch< UniformTag, Remote<Brick> > > Field_t;
#else
typedef Field< Mesh_t, int, MultiPatch<UniformTag,Brick> > Field_t;
#endif

// Field layout type

typedef Field_t::Layout_t FLayout_t;

// The particle traits class we'll use

typedef PTraits<AttrEngineTag_t,Mesh_t,FLayout_t> PTraits_t;

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

  tester.out() << argv[0] << ": Particles with spatial layout" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // Create a FieldLayout object.  We don't actually need a Field in
  // this example though, just the layout.

  tester.out() << "Creating FieldLayout object ..." << std::endl;
  Interval<PDim> meshDomain(12, 24);
  Loc<PDim> blocks(3,4);
  UniformGridPartition<PDim> gpar(blocks);
  DistributedMapper<PDim> cmap(gpar);
  FLayout_t flayout(meshDomain, gpar, cmap);

  // Create the UR mesh

  tester.out() << "Creating UniformRectilinearMesh object ..." << std::endl;
  Molecule<PTraits_t>::PointType_t meshOrigin(1.0, 2.0);
  Molecule<PTraits_t>::PointType_t meshSpacing(0.5, 0.5);
  Mesh_t mesh(flayout,meshOrigin,meshSpacing);

  // Create a spatial layout object for our use

  tester.out() << "Creating Particles SpatialLayout object ..." << std::endl;
  PLayout_t layout(mesh,flayout);

  // Create a Particles object, using our special subclass

  tester.out() << "Creating Molecule object ..." << std::endl;
  Molecule<PTraits_t> mol(layout);

  tester.out() << "Molecule created; initially, num attributes = ";
  tester.out() << mol.attributes() << ", num particles = " << mol.size();
  tester.out() << ", total patches = " << mol.attributeLayout().sizeGlobal();
  tester.out() << ", local patches = " << mol.attributeLayout().sizeLocal();
  tester.out() << std::endl;

  tester.check(mol.attributes() == 3);
  tester.check(mol.size() == 0);
  tester.check(mol.attributeLayout().sizeGlobal() == flayout.sizeGlobal());
  tester.check(mol.attributeLayout().sizeLocal() == flayout.sizeLocal());

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

  // Initialize the positions.

  tester.out() << "Initializing with random position values ..." << std::endl;
  std::srand(12345U);
  int i;
  for (i = 0; i < createnum; ++i) {
    double ranx = std::rand() / static_cast<double>(RAND_MAX);
    ranx = meshOrigin(0) +
           ranx * (meshDomain[0].length() - 1) * meshSpacing(0);
    double rany = std::rand() / static_cast<double>(RAND_MAX);
    rany = meshOrigin(1) +
           rany * (meshDomain[1].length() - 1) * meshSpacing(1);
    Molecule<PTraits_t>::PointType_t newpos(ranx, rany);
    mol.pos(i) = newpos;
  }
  mol.mom = mol.pos * 100.0;
  mol.charge = 3.3;

  tester.out() << "Contents of particles:" << std::endl;
  tester.out() << mol << std::endl;

  // Sync the particles now that we've changed positions.

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

  for (i = 0; i < createnum; ++i)
    mol.pos(i+createnum) = mol.pos(i);
  mol.mom = mol.pos * 50.0;
  mol.charge = 6.6;

  tester.out() << "Syncing particles again ..." << std::endl;
  mol.sync(mol.pos);
  tester.out() << "After sync, contents of particles:" << std::endl;
  tester.out() << mol << std::endl;

  // Return resulting error code and exit

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("Particles with spatial layout");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: spatial.cpp,v $   $Author: richard $
// $Revision: 1.25 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
