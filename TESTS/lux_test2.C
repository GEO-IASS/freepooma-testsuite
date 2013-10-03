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
// Lux test 2: Display a 2D and a 3D particle set
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Lux.h"
#include "Pooma/Arrays.h"
#include "Pooma/Fields.h"
#include "Pooma/Particles.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"


// Traits class for Particles object
template <class EngineTag, class MeshType, class FL,
          class InterpolatorTag>
struct PTraits
{
  // The type of engine to use in the attributes
  typedef EngineTag AttributeEngineTag_t;

  // The type of particle layout to use
  typedef SpatialLayout<MeshType,FL> 
    ParticleLayout_t;

  // The type of interpolator to use
  typedef InterpolatorTag InterpolatorTag_t;
};

// Particles subclass with position and velocity
template <class PT>
class ChargedParticles : public Particles<PT>
{
public:
  // Typedefs
  typedef Particles<PT>                          Base_t;
  typedef typename Base_t::AttributeEngineTag_t  AttributeEngineTag_t;
  typedef typename Base_t::ParticleLayout_t      ParticleLayout_t;
  typedef typename ParticleLayout_t::AxisType_t  AxisType_t;
  typedef typename ParticleLayout_t::PointType_t PointType_t;
  typedef typename PT::InterpolatorTag_t         InterpolatorTag_t;

  // Dimensionality
  static const int dimensions = ParticleLayout_t::dimensions;

  // Constructor: set up layouts, register attributes
  ChargedParticles(const ParticleLayout_t &pl)
  : Particles<PT>(pl)
  {
    this->addAttribute(R);
    this->addAttribute(V);
    this->addAttribute(E);
    this->addAttribute(qm);
  }

  // Position and velocity attributes (as public members)
  DynamicArray<PointType_t,AttributeEngineTag_t> R;
  DynamicArray<PointType_t,AttributeEngineTag_t> V;
  DynamicArray<PointType_t,AttributeEngineTag_t> E;
  DynamicArray<double,     AttributeEngineTag_t> qm;
};


// Dimensionality of this problem
static const int PDim = 2;

// Engine tag type for attributes
typedef MultiPatch<GridTag,Brick> AttrEngineTag_t;

// Mesh type
typedef UniformRectilinearMesh<PDim,double> Mesh_t;

// Field types
typedef Field< Mesh_t, double,
               MultiPatch<UniformTag,Brick> > DField_t;
typedef Field< Mesh_t, Vector<PDim,double>,
               MultiPatch<UniformTag,Brick> > VecField_t;

// Field layout type, derived from Engine type
typedef DField_t::Engine_t Engine_t;
typedef Engine_t::Layout_t FLayout_t;

// Type of interpolator
typedef NGP InterpolatorTag_t;

// Particle traits class
typedef PTraits<AttrEngineTag_t,Mesh_t,FLayout_t,
                InterpolatorTag_t> PTraits_t;

// Type of particle layout
typedef PTraits_t::ParticleLayout_t PLayout_t;

// Type of actual particles
typedef ChargedParticles<PTraits_t> Particles_t;

// Grid sizes
const int nx = 200, ny = 200;

// Number of particles in simulation
const int NumPart = 400;

// Number of timesteps in simulation
const int NumSteps = 20;

// The value of pi
const double pi = acos(-1.0);

// Maximum value for particle q/m ratio
const double qmmax = 1.0;

// Timestep
const double dt = 1.0;


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Lux Particles PIC2d display test" << std::endl;
  tester.out() << "---------------------------------------------" << std::endl;

#if POOMA_LUX

  tester.out() << "Initializing particles ..." << std::endl;

  // Create mesh and geometry objects for cell-centered fields.
  Interval<PDim> meshDomain(nx+1,ny+1);
  Mesh_t mesh(meshDomain);
  Geometry_t geometry(mesh);

  // Create a second geometry object that includes a guard layer.
  GuardLayers<PDim> gl(1);
  Geometry_t geometryGL(mesh,gl);

  // Create field layout objects for our electrostatic potential
  // and our electric field.  Decomposition is 4 x 4.
  Loc<PDim> blocks(4,4);
  FLayout_t flayout(geometry.physicalDomain(),blocks);
  FLayout_t flayoutGL(geometryGL.physicalDomain(),blocks,gl);

  // Create and initialize electrostatic potential and electric field.
  DField_t phi(geometryGL,flayoutGL);
  VecField_t EFD(geometry,flayout);

  // potential phi = phi0 * sin(2*pi*x/Lx) * cos(4*pi*y/Ly)
  // Note that phi is a periodic Field
  // Electric field EFD = -grad(phi);
  Pooma::addAllPeriodicFaceBC(phi, 0.0);
  double phi0 = 0.01 * static_cast<double>(nx);
  phi = phi0 * sin(2.0*pi*phi.x().comp(0)/nx)
             * cos(4.0*pi*phi.x().comp(1)/ny);
  EFD = -grad<Centering_t>(phi);

  // Create a particle layout object for our use
  PLayout_t layout(geometry,flayout);

  // Create a Particles object and set periodic boundary conditions
  Particles_t P(layout);
  Particles_t::PointType_t lower(0.0,0.0), upper(nx,ny);
  PeriodicBC<Particles_t::PointType_t> bc(lower,upper);
  P.addBoundaryCondition(P.R,bc);

  // Create an equal number of particles on each processor
  // and recompute global domain.
  P.globalCreate(NumPart);

  // Random initialization for particle positions in nx by ny domain
  // Zero initialization for particle velocities
  // Random intialization for charge-to-mass ratio from -qmmax to qmmax
  P.V = Particles_t::PointType_t(0.0);
  srand(12345U);
  Particles_t::PointType_t initPos;
  for (int i = 0; i < NumPart; ++i)
  {
    initPos(0) = nx * rand() /
	         static_cast<Particles_t::AxisType_t>(RAND_MAX);
    initPos(1) = ny * rand() /
	         static_cast<Particles_t::AxisType_t>(RAND_MAX);
    P.R(i) = initPos;
    P.qm(i) = (2.0 * rand() / static_cast<double>(RAND_MAX) - 1.0) *
              qmmax;
  }

  // Redistribute particle data based on spatial layout
  P.swap(P.R);

  tester.out() << "PIC2d setup complete." << std::endl;
  tester.out() << "---------------------" << std::endl;

  // Create a Lux connection

  tester.out() << "Creating LuxConnection object ..." << std::endl;
  Connection<Lux> lux("test2");
  tester.out() << "Finished creating LuxConnection object." << std::endl;

  // Add attributes in to the display

  tester.out() << "Connecting qm for display ..." << std::endl;
  lux.connect("P-qm", P.R, P.qm, ConnectionBase::out);
  tester.out() << "Connecting velocity for display ..." << std::endl;
  lux.connect("P-velocity", P.R, P.V, ConnectionBase::out);
  tester.out() << "Connecting E.x and E.y for display ..." << std::endl;
  lux.connect("E-x", EFD.comp(0), ConnectionBase::out);
  lux.connect("E-y", EFD.comp(1), ConnectionBase::out);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  lux.ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Begin main timestep loop
  for (int it=1; it <= NumSteps; ++it)
  {
    // Advance particle positions
    tester.out() << "Advance particle positions ..." << std::endl;
    P.R = P.R + dt * P.V;

    // Invoke boundary conditions and update particle distribution
    tester.out() << "Synchronize particles ..." << std::endl;
    P.sync(P.R);
   
    // Gather the E field to the particle positions
    tester.out() << "Gather E field ..." << std::endl;
    gather( P.E, EFD, P.R, Particles_t::InterpolatorTag_t() );

    // Advance the particle velocities
    tester.out() << "Advance particle velocities ..." << std::endl;
    P.V = P.V + dt * P.qm * P.E;

    // Update display
    Pooma::blockAndEvaluate();
    tester.out() << "Updating for iters = " << it << std::endl;
    lux.ready();
  }

  // Display the final particle positions, velocities and qm values.
  tester.out() << "PIC2d timestep loop complete!" << std::endl;
  tester.out() << "-----------------------------" << std::endl;

#else // POOMA_LUX

  tester.out() << "Please configure with --lux to use this test code!"
	       << std::endl;

#endif // POOMA_LUX

  // Finish up and report results

  int retval = tester.results("Lux Particles PIC2d display test");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: lux_test2.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
