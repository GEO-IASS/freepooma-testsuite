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
// Helper items for the particle_test and particle_bench test codes.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Pooma/Particles.h"
#include "Pooma/DynamicArrays.h"
#include "Pooma/Domains.h"
#include "Pooma/Arrays.h"
#include "Pooma/Fields.h"
#include "Utilities/Clock.h"
#include "Utilities/Tester.h"
#include "Utilities/Inform.h"
#include "Tiny/Vector.h"

#include <stdlib.h>


// Static storage for a Tester class instance

Pooma::Tester *tester = 0;

// Static storage for a benchmark Inform output stream, and a macro
// for printing to it.

Inform *benchmsg = 0;
#define BCHOUT(x)     *benchmsg << x


//-----------------------------------------------------------------------------
// A traits class for a Particles object.  This example is templated
// on:
//
//   EngineTag == tag type for the attribute array engine (e.g., MP Dynamic)
//                This must be an engine type that supports shared layouts.
//
//   ParLayout == Layout type for the particles.  For example, SpatialLayout
//                templated on the proper other things like the geometry
//                and field layout types.
//
// This is used to construct the TestParticles class, defined after this.
//
//-----------------------------------------------------------------------------

template <class EngineTag, class ParLayout>
struct TestParTraits
{
  // The type of engine to use in the attributes.  Must be a shared-layout
  // engine type, that supports dynamic operations.

  typedef EngineTag AttributeEngineTag_t;

  // The type of particle layout to use.  This is where we get the number
  // of dimensions from.

  typedef ParLayout ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// A Particles subclass, that defines a few attributes for our tests.
//
// TestParticles is templated on a traits class type; for our tests,
// we will use the TestParTraits class from above.  This class will
// create and initialize three attributes, all provided as public data
// members:
//
//   pos - a Vector attribute, stores the position of each particle.
//   mom - a Vector attribute, could be anything (but say it's momentum).
//    da - a 'double' attribute, stores a scalar quantity.
//    di - a 'int' attribute, stores a scalar quantity.
//
// The only constructor for this class expects to be given a particle
// layout object, of which it makes a copy.
//
//-----------------------------------------------------------------------------

template <class PT>
class TestParticles : public Particles<PT>
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

  TestParticles(const ParticleLayout_t &pl)
    : Particles<PT>(pl)
    {
      addAllAttributes();
    }

  // Default constructor; if this is used, must call initialize() later

  TestParticles() { }

  // initialize this object, if it has not already been initialized.

  void initialize(const ParticleLayout_t &pl)
    {
      Particles<PT>::initialize(pl);
      addAllAttributes();
    }

  // List of attributes; we'll just make them public data members here,
  // you could also provide access via methods.

  DynamicArray<PointType_t,AttributeEngineTag_t> pos;
  DynamicArray<PointType_t,AttributeEngineTag_t> mom;
  DynamicArray<AxisType_t,AttributeEngineTag_t>  ad;
  DynamicArray<int,AttributeEngineTag_t>         ai;

private:
  // Add in the attributes that we have here

  void addAllAttributes()
    {
      this->addAttribute(pos);
      this->addAttribute(mom);
      this->addAttribute(ad);
      this->addAttribute(ai);
    }
};


//-----------------------------------------------------------------------------
//
// A routine to initialize the particle test
//
//-----------------------------------------------------------------------------

void startParticleTest(int argc, char *argv[], const char *msg)
{
  // Initialize POOMA and create a Tester

  Pooma::initialize(argc, argv);
  tester = new Pooma::Tester(argc, argv);
  benchmsg = new Inform(argv[0]);

  // Print an intro message

  tester->out() << argv[0] << ": " << msg << std::endl;
  tester->out() << "-------------------------------------------------------";
  tester->out() << std::endl;
}


//-----------------------------------------------------------------------------
//
// A routine to finish the particle test.  Returns the return code for
// the main routine.
//
//-----------------------------------------------------------------------------

int endParticleTest(const char *msg)
{
  // Print out the results

  tester->out() << "-------------------------------------------------------";
  tester->out() << std::endl;
  int retval = tester->results(msg);

  // Delete global objects

  delete tester;
  delete benchmsg;

  // Shut down POOMA

  Pooma::finalize();

  // Return the value to be used as the main routine return value.

  return retval;
}


//-----------------------------------------------------------------------------
//
// The main test routine.  This should be called with the Particle's
// object to test, and a Region<Dim,T> storing the bounding region
// in which the particle should be created and allowed to move.  The
// rest must have been set up by the user before calling this.
//
//-----------------------------------------------------------------------------

template <class PT, int Dim, class T>
void runParticleTest(TestParticles<PT> &P, const Region<Dim,T> &box)
{
  int i, n, d, createnum = 10;
  typedef typename TestParticles<PT>::PointType_t PointType_t;

  // Get the origin and lengths of the box

  PointType_t origin, len;
  DomainToVector(box, origin);
  DomainToVector(box.lengths(), len);

  // Initially, the particle object should be empty.

  tester->out() << "Starting test.  Initial Particles object:" << std::endl;
  tester->out() << P << std::endl;
  tester->out() << "Moving in box with origin = " << origin;
  tester->out() << ", lengths = " << len << std::endl;
  tester->out() << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  tester->check("A", P.size() == 0 && P.attributes() == 4);
  tester->check("Initialized", P.initialized());

  // Get context information

  int myContext = Pooma::context();
  int numContexts = Pooma::contexts();

  // Get the number of local patches in the attribute layout.

  int patches = P.attributeLayout().sizeLocal();
  tester->check("B", patches > 0);

  // Create particles in the last patch of the Particles object

  tester->out() << "Creating " << createnum << " particles in first and ";
  tester->out() << "last patches on each context, out of " << patches 
                << " local patches total" << std::endl;

  P.create(createnum, 0, false);
  P.create(createnum);		// this will also do a renumber
  for (i=0; i < patches; ++i)
    {
      int size = P.attributeLayout().ownedDomain(i).size();
      bool ok;
      if (patches > 1)
	ok = (size == (i == 0 || i == (patches-1) ? createnum : 0));
      else
	ok = (size == 2*createnum);
      tester->check("Patch size", ok);
    }
  //  tester->check("C", P.size() == 2 * createnum);

  // Initialize the particles, and then sync them and find out if they're in
  // the right place.

  tester->out() << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  tester->out() << "Initializing particles ..." << std::endl;
  createnum = P.size();
  for (i=0; i < createnum; ++i)
    {
      P.pos(i) = origin + len * ( T(i) / T(createnum) );
      P.mom(i) = P.pos(i) * 10.0;
      P.ad(i) = 0.01 * i;
      P.ai(i) = i + 1;
    }
  tester->out() << "Particles after initialization, before sync:" << std::endl;
  tester->out() << P << std::endl;
  tester->out() << "Syncing particles ..." << std::endl;
  P.sync(P.pos);
  //  P.sync(P.pos);
  tester->out() << "Particles after sync:" << std::endl;
  tester->out() << P << std::endl;
  tester->check("D", P.size() == createnum);

  // Destroy some of the particles ... first destroy all the even-numbered
  // ones (based on the ai attrib), then all with their ai value > createnum-4;

  tester->out() << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  tester->out() << "Destroying some particles ..." << std::endl;
  P.setDestroyMethod(ShiftUp());
  Array<1, int> dlist(createnum);
  for (i=0, n=0; i < createnum; ++i)
    {
      if ((P.ai(i) % 2) == 0)
	dlist(n++) = i;
    }
  tester->out() << "Destroying all particles w/ even ai values." << std::endl;
  P.destroy(IndirectionList<int>(dlist(Interval<1>(n))));
  tester->out() << "New Particles domain = " << P.ai.domain() << std::endl;
  for (i=0, n=0; i < P.size(); ++i)
    {
      if (P.ai(i) > (createnum - 4) && (P.ai(i) % 2) != 0)
	dlist(n++) = i;
    }
  tester->out() << "Destroying (cached) all particles w/ odd ai values > ";
  tester->out() << createnum - 4 << std::endl;
  P.deferredDestroy(IndirectionList<int>(dlist(Interval<1>(n))));
  tester->out() << "Carrying out destroy requests ..." << std::endl;
  P.sync(P.pos);
  tester->out() << "Particles after sync:" << std::endl;
  tester->out() << P << std::endl;
  for (i=0; i < P.size(); ++i)
    tester->check("E", (P.ai(i) % 2) != 0 && P.ai(i) <= (createnum - 4));

  // Create boundary conditions now:
  //   Absorb BC's for the integer attributes (should all get stuck at 0)
  //   Reverse BC's for the double attributes (should all become negative)
  //   Periodic BC's for the positions.

  tester->out() << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  tester->out() << "Adding Absorb BC for integer attribute ..." << std::endl;
  AbsorbBC<int> absorbBC(0, 0);
  P.addBoundaryCondition(P.ai, absorbBC);
  tester->out() << "Adding Reverse BC for double attribute ..." << std::endl;
  ReverseBC<double> reverseBC(10000, 20000);
  P.addBoundaryCondition(P.ad, reverseBC);
  tester->out() << "Adding Periodic BC for pos attribute ..." << std::endl;
  PeriodicBC<PointType_t> periodicBC(origin, origin + len);
  P.addBoundaryCondition(P.pos, periodicBC);

  tester->out() << "Doing sync to apply boundary conditions ..." << std::endl;
  P.sync(P.pos);
  tester->out() << "Particles after sync:" << std::endl;
  tester->out() << P << std::endl;
  for (i=0; i < P.size(); ++i)
    {
      tester->check(P.ai(i) == 0);
      tester->check(P.ad(i) <= 0.0);
      for (d=0; d < Dim; ++d)
	tester->check("F", P.pos(i)(d) >= origin(d) &&
		           P.pos(i)(d) <= (origin(d) + len(d)));
    }

  tester->out() << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  tester->out() << "Removing all the boundary conditions now." << std::endl;
  P.removeBoundaryConditions();
  tester->out() << "Doing final sync, should not change anything."<<std::endl;
  P.sync(P.pos);
  tester->out() << "Particles after sync:" << std::endl;
  tester->out() << P << std::endl;
  for (i=0; i < P.size(); ++i)
    {
      tester->check(P.ai(i) == 0);
      tester->check(P.ad(i) <= 0.0);
      for (d=0; d < Dim; ++d)
	tester->check("G", P.pos(i)(d) >= origin(d) &&
		           P.pos(i)(d) <= (origin(d) + len(d)));
    }
}


//-----------------------------------------------------------------------------
//
// A benchmark routine, that creates a large number of particles in
// the problem domain space, moves them around randomly, and times how
// long it takes to do expression evaluation and particle syncs.
// Results are printed to stdout; more detailed messages are printed to
// the tester stream (and so can be disabled by not giving the -v flag).
//
//-----------------------------------------------------------------------------

template <class PT, int Dim, class T>
void runParticleBenchmark(int argc, char *argv[],
			  TestParticles<PT> &P,
			  const Region<Dim,T> &box)
{
  int i, n, d;

  // Get the origin and lengths of the box

  typedef typename TestParticles<PT>::PointType_t PointType_t;
  PointType_t origin, len, lenfrac;
  DomainToVector(box, origin);
  DomainToVector(box.lengths(), len);

  // Default parameters for the benchmark.

#if POOMA_BOUNDS_CHECK
  int iters = 10;
#else
  int iters = 1000;
#endif
  int startnumparticles = 100;
  int endnumparticles = 10000;
  int multnumparticles = 10;
  double movefrac = 0.1;
  bool usesync = false;

  // Parse the command-line options and look for special ones indicating
  // how to run the benchmark.  Arguments are:
  //   -iters <iterations>
  //   -size <begin> <end> <multfactor>

  for (i=1; i < argc; ++i)
    {
      std::string arg(argv[i]);
      if (arg == "-iters" && i < (argc - 1))
	{
	  iters = atoi(argv[++i]);
	}
      else if (arg == "-frac" && i < (argc - 1))
	{
	  movefrac = atof(argv[++i]);
	}
      else if (arg == "-size" && i < (argc - 3))
	{
	  startnumparticles = atoi(argv[++i]);
	  endnumparticles = atoi(argv[++i]);
	  multnumparticles = atoi(argv[++i]);
	}
      else if (arg == "-sync")
	{
          usesync = true;
        }
    }

  // Compute the amount to add to position each time, based on the fraction

  lenfrac = len * movefrac;

  // Print summary of this benchmark.

  BCHOUT("Starting Particles benchmark." << std::endl);
  BCHOUT("-----------------------------------" << std::endl);
  BCHOUT("              Iterations: " << iters << std::endl);
  BCHOUT("  Starting particle size: " << startnumparticles << std::endl);
  BCHOUT("    Ending particle size: " << endnumparticles << std::endl);
  BCHOUT("Particle size multiplier: " << multnumparticles << std::endl);
  BCHOUT("Fraction moving off edge: " << movefrac << std::endl);
  BCHOUT("-----------------------------------" << std::endl);

  // Check for errors

  if (iters < 1)
    {
      BCHOUT("ERROR: Illegal iteration value.  Exiting." << std::endl);
      return;
    }

  if (startnumparticles < 1 || endnumparticles < 1 ||
      endnumparticles < startnumparticles || multnumparticles < 1)
    {
      BCHOUT("ERROR: Illegal particle size values.  Exiting." << std::endl);
      return;
    }

  // Add in periodic BC's for these particles

  tester->out() << "Setting up periodic BC's for particles, for origin = ";
  tester->out() << origin << " and size = " << len << std::endl;
  PeriodicBC<PointType_t> periodicBC(origin, origin + len);
  P.addBoundaryCondition(P.pos, periodicBC);

  // No errors; start looping over iterations and particle sizes.

  int numparticles = startnumparticles;
  while (numparticles <= endnumparticles)
    {
      tester->out() << "Starting work for iters = " << iters;
      tester->out() << ", numparticles = " << numparticles << std::endl;

      // Empty out any existing particles.

      if (P.size() > 0)
	{
	  tester->out() << "Removing existing "<< P.size() <<" particles.";
	  tester->out() << std::endl;
	  P.destroy(Interval<1>(P.size()));
	}
      tester->out() << "Finished clearing out old particles: P = ";
      tester->out() << P << std::endl;

      // Initialize the RND

      srand(12345U);

      // Create positions with random values

      tester->out() << "Creating and initializing " << numparticles;
      tester->out() << " particles in box with origin = " << origin;
      tester->out() << " and size = " << len << std::endl;
      P.globalCreate(numparticles); // will also renumber
      P.pos = PointType_t(0.0);
      P.mom = PointType_t(0.0);
      P.ad = 0;
      P.ai = 0;
      Pooma::blockAndEvaluate();
      tester->out() << "After create, attrib layout =\n";
      tester->out() << P.attributeLayout() << std::endl;
      for (n=0; n < numparticles; ++n)
	{
	  PointType_t initvec;
	  for (d=0; d < Dim; ++d)
	    {
	      initvec(d) = origin(d) +
		len(d) * 0.99 * rand() / static_cast<double>(RAND_MAX);
	      PAssert(initvec(d) >= origin(d) &&
		      initvec(d) < (origin(d) + len(d)));
	    }
	  P.pos(n) = initvec;
	}

      // Move these all to the proper location

      tester->out() << "Swapping particles after initialization ...";
      tester->out() << std::endl;
      P.swap(P.pos);

      // Do a series of iterations, moving the particles and performing
      // computations.

      double computetime = 0.0;
      double swaptime = 0.0;
      double bctime = 0.0;

      for (i=0; i < iters; ++i)
	{
	  tester->out() << "Performing iteration " << i << std::endl;

	  // Perform computations on the elements

	  tester->out() << "Timing computation of P.mom = P.pos * P.pos + len";
	  tester->out() << " and P.pos += lenfrac" << std::endl;
	  double compute = Pooma::Clock::value();
	  for (int comp=0; comp < 10; ++comp)
	    {
	      P.mom = P.pos * P.pos + len;
	      P.pos += (lenfrac * 0.1);
	    }
	  Pooma::blockAndEvaluate();
	  compute = Pooma::Clock::value() - compute;
	  computetime += compute;
	  tester->out() << "Computation took " << compute << std::endl;

          // Either do BC's and swap separately, or combine them

          if (usesync)
            {
              // Just time sync operation.

	      tester->out() << "Timing sync ..." << std::endl;
	      double swap = Pooma::Clock::value();
	      P.sync(P.pos);
	      swap = Pooma::Clock::value() - swap;
	      swaptime += swap;
	      tester->out() << "Syncing took " << swap << std::endl;
	    }
          else
            {
	      // Now apply BC's to particles, timing this

	      tester->out() << "Timing periodic BC's ..." << std::endl;
	      double bc = Pooma::Clock::value();
	      P.applyBoundaryConditions();
	      Pooma::blockAndEvaluate();
	      bc = Pooma::Clock::value() - bc;
	      bctime += bc;
	      tester->out() << "BC's took " << bc << std::endl;

	      // Move particles around, and compute time to do this

	      tester->out() << "Timing swap ..." << std::endl;
	      double swap = Pooma::Clock::value();
	      P.swap(P.pos);
	      swap = Pooma::Clock::value() - swap;
	      swaptime += swap;
	      tester->out() << "Swapping took " << swap << std::endl;
	    }
        }

      // Report the results

      double mflops = numparticles * 0.000001;
      mflops *= (10.0 * iters * 4.0 * Dim);
      mflops /= computetime;

      BCHOUT("For " << numparticles << " particles, " << iters);
      BCHOUT(" iterations:" << std::endl);
      BCHOUT("    Computation time: " << computetime << std::endl);
      BCHOUT("  Computation MFLOPS: " << mflops << std::endl);
      if (!usesync)
        {
          BCHOUT("             BC time: " << bctime << std::endl);
          BCHOUT("           Swap time: " << swaptime << std::endl);
        }
      else
        {
          BCHOUT("           Sync time: " << swaptime << std::endl);
        }
      BCHOUT("-----------------------------------" << std::endl);

      // Update the number of particles

      if (numparticles < endnumparticles)
	{	
	  numparticles *= multnumparticles;
	  if (numparticles > endnumparticles)
	    numparticles = endnumparticles;
	}
      else
	numparticles++;
    }
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: particle_tests.h,v $   $Author: richard $
// $Revision: 1.24 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
