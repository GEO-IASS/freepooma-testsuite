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
// Particles test: Particles create and destroy operations
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Particles/Particles.h"
#include "Particles/SpatialLayout.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Layout/UniformGridLayout.h"
#include "Layout/DynamicLayout.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "Engine/BrickEngine.h"
#include "DynamicArray/DynamicArray.h"
#include "Field/Mesh/UniformRectilinearMesh.h"
#include "Field/Field.h"
#include "Field/FieldEngine/FieldEngine.h"
#include "Tiny/Vector.h"

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

  typedef SpatialLayout<Mesh,FL> ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// A trivial Particles subclass that defines only a position attribute.
//-----------------------------------------------------------------------------

template <class PT>
class Point : public Particles<PT>
{
public:
  // Useful typedefs to get from the base class

  typedef Particles<PT>                         Base_t;
  typedef typename Base_t::AttributeEngineTag_t AttributeEngineTag_t;
  typedef typename Base_t::ParticleLayout_t     ParticleLayout_t;
  typedef typename ParticleLayout_t::AxisType_t AxisType_t;

  // Useful enums to get from the base class

  enum { dimensions = ParticleLayout_t::dimensions };

  // Constructor: set up layouts, register attributes

  Point(const ParticleLayout_t &pl)
    : Particles<PT>(pl)
    {
      this->addAttribute(Pos);
    }

  // List of attributes; we'll just make them public data members here,
  // you could also provide access via methods.

  DynamicArray< Vector<dimensions,AxisType_t>, AttributeEngineTag_t > Pos;
};


//-----------------------------------------------------------------------------
// Typedefs for what we will compute
//-----------------------------------------------------------------------------

// Dimensionality of this problem

static const int PDim = 3;

// Engine tag type for attributes

#if POOMA_MESSAGING
typedef MultiPatch< DynamicTag, Remote<Dynamic> > AttrEngineTag_t;
#else
typedef MultiPatch<DynamicTag,Dynamic> AttrEngineTag_t;
#endif

// Mesh type

typedef UniformRectilinearMesh<PDim> Mesh_t;

// Field type

#if POOMA_MESSAGING
typedef Field< Mesh_t, double, MultiPatch< UniformTag, Remote<Brick> > > 
  Field_t;
#else
typedef Field< Mesh_t, double, MultiPatch<UniformTag,Brick> > Field_t;
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

  tester.out() << argv[0] << ": Particles create/destroy operations\n";
  tester.out() << "------------------------------------------------\n";

  // Create a FieldLayout object

  tester.out() << "Creating Field layout object ..." << std::endl;
  int idim;
  Interval<PDim> meshDomain;
  for (idim=0; idim<PDim; ++idim)
    meshDomain[idim] = Interval<1>(6*idim+4);
  Loc<PDim> decomp;
  for (idim=0; idim<PDim; ++idim)
    decomp[idim] = Loc<1>(2);

  FLayout_t flayout(meshDomain, decomp, DistributedTag());
  tester.out() << "Layout created, Layout = " << flayout << std::endl;

  // Create a UniformRectilinearMesh object

  tester.out() << "Creating URM object ..." << std::endl;
  Vector<PDim,double> meshOrigin(0.0);
  Vector<PDim,double> meshSpacings(1.0);
  DomainLayout<PDim> meshLayout(meshDomain);
  Mesh_t mesh(meshLayout,meshOrigin,meshSpacings);

  // Create a spatial layout object for our use

  tester.out() << "Creating SpatialLayout object ..." << std::endl;
  PLayout_t layout(mesh, flayout);

  // Create a Particles object, using our special subclass

  tester.out() << "Creating Point object ..." << std::endl;
  Point<PTraits_t> point(layout);
  tester.out() << "Point created; initially, num attributes = ";
  tester.out() << point.attributes() << ", num particles = " << point.size();
  tester.out() << std::endl;
  tester.check(point.attributes() == 1);
  tester.check(point.size() == 0);

  // Create some particles and recompute the global domain ...

  tester.out() << "Creating 20 particles ..." << std::endl;
  point.globalCreate(20);
  tester.check(point.size() == 20);
  tester.out() << "Contents of Point object:" << std::endl;
  tester.out() << point << std::endl;

  // Block before serial code

  Pooma::blockAndEvaluate();

  // Initialize positions to random values within our domain ...

  tester.out() << "Initializing particle positions ..." << std::endl;
  Vector<PDim,Point<PTraits_t>::AxisType_t> initPos;
  srand(12345U);
  int ip;
  for (ip=0; ip<20; ++ip)
    {
      for (idim=0; idim<PDim; ++idim)
	initPos(idim) = (6*idim+4.0) *
	  (rand() / static_cast<Point<PTraits_t>::AxisType_t>(RAND_MAX));
      point.Pos(ip) = initPos;
    }
  
  // Print out Point object, sync, then print again ...

  typedef Point<PTraits_t>::AttributeLayout_t AttributeLayout_t;
  typedef AttributeLayout_t::iterator piterator;
  piterator p, pend;
  tester.out().setOutputContext(-1);

  pend = point.attributeLayout().endLocal();
  ip = 0;
  for (p = point.attributeLayout().beginLocal(); p != pend; ++p, ++ip)
    {
      tester.out() << "Size of Local Patch " << ip << " = ";
      tester.out() << p->domain().size() << std::endl;
    }

  tester.out().setOutputContext(0);
  tester.out() << "Pos attribute:\n" << point.Pos << std::endl;
  tester.out() << "Syncing particles ..." << std::endl;

  point.sync(point.Pos);
  tester.check(point.size() == 20);

  tester.out().setOutputContext(-1);
  pend = point.attributeLayout().endLocal();
  ip = 0;
  for (p = point.attributeLayout().beginLocal(); p != pend; ++p, ++ip)
    {
      tester.out() << "Size of Local Patch " << ip << " = ";
      tester.out() << p->domain().size() << std::endl;
    }

  tester.out().setOutputContext(0);
  tester.out() << "Pos attribute:" << std::endl << point.Pos << std::endl;

  // Now destroy some of the particles, renumber, and print out again ...

  tester.out() << "Destroying particles 5 thru 12 ..." << std::endl;
  point.destroy(Interval<1>(5,12));
  tester.check(point.size()==12);

  tester.out().setOutputContext(-1);
  pend = point.attributeLayout().endLocal();
  ip = 0;
  for (p = point.attributeLayout().beginLocal(); p != pend; ++p, ++ip)
    {
      tester.out() << "Size of Local Patch " << ip << " = ";
      tester.out() << p->domain().size() << std::endl;
    }

  tester.out().setOutputContext(0);
  tester.out() << "Pos attribute:" << std::endl << point.Pos << std::endl;

  // Now change to ShiftUp destroy method and do a deferred destroy ...

  tester.out() << "Doing deferred destroy of odd-numbered particles with ";
  tester.out() << "ShiftUp method ..." << std::endl;
  point.setDestroyMethod(ShiftUp());
  point.deferredDestroy(Range<1>(1,11,2));
  tester.check(point.size()==12);

  // Block before serial code

  Pooma::blockAndEvaluate();

  // assign new position values to the particles

  tester.out() << "Assigning new position values to particles ... "<<std::endl;
  for (ip=0; ip<12; ++ip)
    {
      for (idim=0; idim<PDim; ++idim)
	initPos(idim) = (6*idim+4.0) *
	  (rand() / static_cast<Point<PTraits_t>::AxisType_t>(RAND_MAX));
      point.Pos(ip) = initPos;
    }

  // Print out Pos attribute object, sync, then print again ...

  tester.out().setOutputContext(-1);
  pend = point.attributeLayout().endLocal();
  ip = 0;
  for (p = point.attributeLayout().beginLocal(); p != pend; ++p, ++ip)
    {
      tester.out() << "Size of Local Patch " << ip << " = ";
      tester.out() << p->domain().size() << std::endl;
    }

  tester.out().setOutputContext(0);
  tester.out() << "Pos attribute:" << std::endl << point.Pos << std::endl;

  tester.out() << "Syncing particles ..." << std::endl;
  point.sync(point.Pos);
  tester.check(point.size()==6);

  tester.out().setOutputContext(-1);
  pend = point.attributeLayout().endLocal();
  ip = 0;
  for (p = point.attributeLayout().beginLocal(); p != pend; ++p, ++ip)
    {
      tester.out() << "Size of Local Patch " << ip << " = ";
      tester.out() << p->domain().size() << std::endl;
    }

  tester.out().setOutputContext(0);
  tester.out() << "Pos attribute:" << std::endl << point.Pos << std::endl;

  // Return resulting error code and exit

  tester.out() << "------------------------------------------------";
  tester.out() << std::endl;
  int retval = tester.results("Particles create/destroy operations");
  Pooma::finalize();
  return retval;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: destroy.cpp,v $   $Author: richard $
// $Revision: 1.22 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
