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
// Particles test: Gather/Scatter NGP Particle/Field interpolation
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/Particles.h"
#include "Pooma/DynamicArrays.h"
#include "Pooma/UMPArrays.h"
#include "Pooma/Fields.h"
#include "Pooma/Tiny.h"

#include <iostream>
#include <stdlib.h>


//-----------------------------------------------------------------------------
// A traits class for a Particles object
//-----------------------------------------------------------------------------

template <class EngineTag, class Mesh, class FL,
          class Interpolator>
struct PTraits
{
  // The type of engine to use in the attributes

  typedef EngineTag AttributeEngineTag_t;

  // The type of particle layout to use

  typedef SpatialLayout< Mesh, FL > ParticleLayout_t;

  // The type of interpolator to use
  
  typedef Interpolator Interpolator_t;
};


//-----------------------------------------------------------------------------
// A Particles subclass, that defines a few attributes
//-----------------------------------------------------------------------------

template <class PT>
class MyParticles : public Particles<PT>
{
public:
  // Useful typedefs

  typedef MyParticles<PT>                        This_t;
  typedef Particles<PT>                          Base_t;
  typedef typename PT::AttributeEngineTag_t      AttributeEngineTag_t;
  typedef typename PT::ParticleLayout_t          ParticleLayout_t;
  typedef typename ParticleLayout_t::AxisType_t  AxisType_t;
  typedef typename ParticleLayout_t::PointType_t PointType_t;
  typedef typename PT::Interpolator_t            Interpolator_t;
  typedef typename Interpolator_t::Cache_t       Cache_t;

  // Useful enums

  enum { dimensions = ParticleLayout_t::dimensions };

  // Constructor: set up layouts, register attributes

  MyParticles(const ParticleLayout_t& pl)
    : Particles<PT>(pl)
    {
      this->addAttribute(pos);
      this->addAttribute(efield);
      this->addAttribute(charge);
      this->addAttribute(cache);
    }

  // List of attributes; we'll just make them public data members here,
  // you could also provide access via methods.

  DynamicArray< PointType_t, AttributeEngineTag_t >  pos;
  DynamicArray< PointType_t, AttributeEngineTag_t >  efield;
  DynamicArray< AxisType_t,  AttributeEngineTag_t >  charge;
  DynamicArray< Cache_t,     AttributeEngineTag_t >  cache;
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
typedef MultiPatch<DynamicTag,Dynamic> AttrEngineTag_t;
#endif

// Mesh type

typedef UniformRectilinearMesh< PDim, double > Mesh_t;

// Field type

#if POOMA_MESSAGING
typedef Field< Mesh_t, double, MultiPatch< UniformTag, Remote<Brick> > > 
  DField_t;
typedef Field< Mesh_t, Vector<PDim,double>,
               MultiPatch< UniformTag, Remote<Brick> > > VecDField_t;
#else
typedef Field< Mesh_t, double, MultiPatch<UniformTag,Brick> > DField_t;
typedef Field< Mesh_t, Vector<PDim,double>,
               MultiPatch<UniformTag,Brick> > VecDField_t;
#endif

// Field layout type

typedef DField_t::Engine_t FEngine_t;
typedef FEngine_t::Layout_t FLayout_t;

// The interpolator types 

typedef Interpolator<PDim,double,NGP>   NGPInterpolator_t;
typedef Interpolator<PDim,double,CIC>   CICInterpolator_t;
typedef Interpolator<PDim,double,SUDS> SUDSInterpolator_t;

// The particle traits class we'll use

typedef PTraits<AttrEngineTag_t,Mesh_t,FLayout_t,NGPInterpolator_t>
  PTraits_t;

// The particle layout type

typedef PTraits_t::ParticleLayout_t PLayout_t;

// The particle type

typedef MyParticles<PTraits_t> Particles_t;


//-----------------------------------------------------------------------------
// The main routine for this test code
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Particle/Field interpolation" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // Create a cell-centered Mesh and Layout object

  tester.out() << "Creating URM object ..." << std::endl;
  
  Particles_t::PointType_t meshOrigin(-1.5, -2.0);
  Particles_t::PointType_t meshSpacing(0.5, 0.5);
  Interval<PDim>           meshDomain(8, 12);
  Mesh_t mesh(meshDomain, meshOrigin, meshSpacing);

  // Create a FieldLayout object. 

  tester.out() << "Creating Field layout object ..." << std::endl;

  GuardLayers<PDim> gl(1);
  Centering<PDim> cell = canonicalCentering<PDim>(CellType, Continuous);

  Loc<PDim> blocks(2,4);
  FLayout_t flayout(meshDomain, blocks, gl, DistributedTag());

  // Create a couple of Fields using this layout.
  // One is an electric field that the particles will gather.
  // The other is the particle charge density, which will be scattered.

  tester.out() << "Creating electric field and charge density field ..."
               << std::endl;
  VecDField_t electric(cell,flayout, mesh);
  DField_t    chargeDensity(cell,flayout, mesh);

  // Add boundary conditions that will manage the external guard layer values

  //FIXME: electric.addBoundaryConditions(AllLinearExtrapolateFaceBC());

  // Create a spatial layout object for the particles.

  tester.out() << "Creating SpatialLayout object ..." << std::endl;
  PLayout_t layout(mesh, flayout);

  // Create a Particles object, using our special subclass

  tester.out() << "Creating MyParticles object ..." << std::endl;
  Particles_t particles(layout);

  // Some quick checks on initialization.

  tester.out() << "Number of particle attributes = " << particles.attributes()
	       << std::endl
	       << "Number of particles = " << particles.size() << std::endl
	       << "Number of attribute patches = "
	       << particles.attributeLayout().sizeGlobal() << std::endl
	       << "Number of field patches = "
	       << flayout.sizeGlobal() << std::endl;
  tester.check("attributes() == 4", particles.attributes() == 4);
  tester.check("size() == 0", particles.size() == 0);
  tester.check("blocks",
    particles.attributeLayout().sizeGlobal() == flayout.sizeGlobal());

  // Create some particles, and then renumber.

  int createnum = 10;
  tester.out() << "Creating " << createnum << " particles "
               << "on context 0, patch 0 ..." << std::endl;
  if (Pooma::context() == 0)
    particles.create(createnum, 0);
  else
    particles.create(0);
  tester.out() << "Created (not yet initialized) ... attrib layout:\n";
  tester.out() << particles.attributeLayout() << std::endl;

  // Initialize the positions and other attributes

  tester.out() << "Initializing values ..." << std::endl;
  srand(12345U);
  Vector<PDim,PLayout_t::AxisType_t> initPos;
  int ip;
  for (ip = 0; ip < createnum; ++ip)
    {
      initPos(0) = 3.0 * rand() / static_cast<double>(RAND_MAX) - 1.5;
      initPos(1) = 4.0 * rand() / static_cast<double>(RAND_MAX) - 2.0;
      particles.pos(ip) = initPos;
    }
  particles.efield = Particles_t::PointType_t(0.0);
  particles.charge = 1.0;

  // Sync the particles now that we've changed positions.

  tester.out() << "Syncing particles ..." << std::endl;
  particles.sync(particles.pos);

  // Print out the particle attributes

  tester.out() << "Particle positions:" << std::endl
    << particles.pos << std::endl;
  tester.out() << "Particle electric field:" << std::endl
    << particles.efield << std::endl;
  tester.out() << "Particle charge:" << std::endl
    << particles.charge << std::endl;

  // Initialize the field values

  tester.out() << "Initializing Field values ..." << std::endl;
  Interval<PDim> dom = electric.physicalDomain();
  for (int i = dom[0].first(); i <= dom[0].last(); ++i)
    for (int j = dom[1].first(); j <= dom[1].last(); ++j)
      electric(i,j) = Particles_t::PointType_t(i+j,i-j);
  chargeDensity = 0.0;

  // Apply field boundary conditions

  electric.applyRelations(true);

  // Print initial field values

  tester.out() << "Electric field:" << std::endl
    << electric << std::endl;
  tester.out() << "Charge density field:" << std::endl
    << chargeDensity << std::endl;

  // Now gather the electric field and scatter the charge
  // using NGP interpolation
 
  tester.out() << "Gathering electric field to particle positions ..."
               << std::endl;
  gather(particles.efield, electric, particles.pos, NGP());
  tester.out()
    << "Scattering particle charge density into field and caching data ..."
    << std::endl;
  scatterCache(particles.charge, chargeDensity, particles.pos, particles.cache,
               NGP());

  // Check that the sum of the charge density is correct

  tester.out() << "Sum of charge density field = "
               << sum(chargeDensity)
	       << std::endl;
  tester.check("chargeDensity(NGP,attrib) == numparticles",
               fabs(sum(chargeDensity) -
		    static_cast<double>(createnum))<1.0e-5);

  // Print out the particle efield

  tester.out() << "Particle electric field:" << std::endl
    << particles.efield << std::endl;

  // Print charge density field

  tester.out() << "Charge density field:" << std::endl
    << chargeDensity << std::endl;

  // Now zero out the particles' electric field and the charge density field,
  // and recompute them with different gather/scatter calls

  tester.out()
    << "Clearing and recomputing electric field and charge density ... "
    << std::endl;
  particles.efield = Particles_t::PointType_t(0.0);
  chargeDensity = 0.0;

  tester.out()
    << "Gathering electric field using cached interpolation data ... "
    << std::endl;
  gatherCache(particles.efield, electric, particles.cache, NGP());
  tester.out()
    << "Scattering particle charge density using constant value ... "
    << std::endl;
  scatterValue(1.0, chargeDensity, particles.pos, NGP());

  // Check that the sum of the charge density is correct

  tester.out() << "Sum of charge density field = "
               << sum(chargeDensity)
	       << std::endl;
  tester.check("chargeDensity(NGP,value) == numparticles",
               fabs(sum(chargeDensity) -
		    static_cast<double>(createnum))<1.0e-5);

  // Print out the particle efield

  tester.out() << "Particle electric field:" << std::endl
    << particles.efield << std::endl;

  // Print charge density field

  tester.out() << "Charge density field:" << std::endl
    << chargeDensity << std::endl;

  // Now zero out the particles' electric field and the charge density field,
  // and recompute them with CIC interpolation

  tester.out()
    << "Clearing and recomputing electric field and charge density ... "
    << std::endl;
  particles.efield = Particles_t::PointType_t(0.0);
  chargeDensity = 0.0;

  // Now gather the electric field and scatter the charge
  // using CIC interpolation
 
  tester.out() << "Gathering electric field using CIC interpolation ..."
               << std::endl;
  gather(particles.efield, electric, particles.pos, CIC());
  tester.out()
    << "Scattering particle charge density using CIC interpolation ..."
    << std::endl;
  scatter(particles.charge, chargeDensity, particles.pos, CIC());

  // Check that the sum of the charge density is correct

  tester.out() << "Sum of charge density field = "
               << sum(chargeDensity(chargeDensity.totalDomain()))
	       << std::endl;
  tester.check("chargeDensity(CIC,attrib) == numparticles",
	       fabs(sum(chargeDensity(chargeDensity.totalDomain())) -
		    static_cast<double>(createnum))<1.0e-5);

  // Print out the particle efield

  tester.out() << "Particle electric field:" << std::endl
    << particles.efield << std::endl;

  // Print charge density field

  tester.out() << "Charge density field:" << std::endl
    << chargeDensity << std::endl;

  // Now zero out the particles' electric field and the charge density field,
  // and recompute them with SUDS interpolation

  tester.out()
    << "Clearing and recomputing electric field and charge density ... "
    << std::endl;
  particles.efield = Particles_t::PointType_t(0.0);
  chargeDensity = 0.0;

  // Now gather the electric field and scatter the charge
  // using SUDS interpolation
 
  tester.out() << "Gathering electric field using SUDS interpolation ..."
               << std::endl;
  gather(particles.efield, electric, particles.pos, SUDS());
  tester.out()
    << "Scattering particle charge density using SUDS interpolation ..."
    << std::endl;
  scatter(particles.charge, chargeDensity, particles.pos, SUDS());

  // Check that the sum of the charge density is correct

  tester.out() << "Sum of charge density field = "
               << sum(chargeDensity(chargeDensity.totalDomain()))
	       << std::endl;
  tester.check("chargeDensity(SUDS,attrib) == numparticles",
	       fabs(sum(chargeDensity(chargeDensity.totalDomain())) -
		    static_cast<double>(createnum))<1.0e-5);

  // Print out the particle efield

  tester.out() << "Particle electric field:" << std::endl
    << particles.efield << std::endl;

  // Print charge density field

  tester.out() << "Charge density field:" << std::endl
    << chargeDensity << std::endl;

  // Return resulting error code and exit

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("Particle/Field interpolation");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: interpolate.cpp,v $   $Author: richard $
// $Revision: 1.25 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
