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

#ifndef POOMA_PARTICLES_SPATIAL_LAYOUT_H
#define POOMA_PARTICLES_SPATIAL_LAYOUT_H

//-----------------------------------------------------------------------------
// Classes:
//   SpatialLayout<Mesh,FieldLayout>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Particles
 * @brief
 * SpatialLayout<Mesh,FieldLayout> is a particle layout class that is used to
 * determine where particles will be located in a parallel environment.
 *
 * SpatialLayout uses a domain decomposition algorithm to keep particles
 * distributed among patches relative to some Field's layout and mesh.
 * It is templated on the mesh type and layout of the system,
 * and must be used in conjunction with Fields with the same mesh type.
 * SpatialLayout is a PatchSwapLayout, and inherits the main "swap" routine
 * from that class.  It provides a "findPatchNumber" routine that calculates
 * patch numbers based on the spatial position of the particles.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/PatchSwapLayout.h"
#include "Partition/SpatialPartition.h"
#include "Utilities/PAssert.h"
#include "Domain/Loc.h"

#include <iosfwd>

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------


/**
 * SpatialLayout is one of the particle layout classes that determine
 * where particles should be located.  It derives from PatchSwapLayout,
 * since it is a layout that will use the infrastructure of PatchSwapLayout
 * to swap particles from one patch to another.  It does this with a
 * spatial decomposition algorithm: given a Mesh and a Layout from
 * a Field (which describes how the domain is decomposed into patches
 * assigned to memory contexts), particles are moved within patches
 * based on their spatial position.
 */

template <class M, class FL>
class SpatialLayout
  : public PatchSwapLayout< SpatialLayout<M,FL> >
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Utility typedef to refer to ourselves and our base class

  typedef SpatialLayout<M,FL> This_t;
  typedef PatchSwapLayout<This_t> Base_t;

  // The mesh type

  typedef M Mesh_t;

  // The field layout type

  typedef FL FieldLayout_t;

  // The domain type

  typedef typename FieldLayout_t::Domain_t Domain_t;

  // The type of data used to store a single coordinate value.
  // This is the type for a single component (e.g., double instead of
  // Vector<3,double>).  PointType_t is the type for a Vector to store
  // positions.

  typedef typename Mesh_t::T_t AxisType_t;
  typedef typename Mesh_t::PointType_t PointType_t;

  // Storage for hold particle amounts ... this should be the same
  // as the typedef in Particles.h

  typedef typename Base_t::Size_t Size_t;

  // The type of array used to store amounts to move to other patches

  typedef typename Base_t::AmountArray_t AmountArray_t;

  // The type of array used to store patch ID's and indices

  typedef typename Base_t::MoveArray_t MoveArray_t;

  // The number of dimensions in the mesh, which is the number of
  // indices needed to look up a mesh vertex or an Array value.

  enum { dimensions = Mesh_t::dimensions };

  //============================================================
  // Constructors
  //============================================================

  // Default constructor.  Initialize with assignment operator.

  SpatialLayout()
    {
      // The Mesh and Layout dimensions must be consistent

      CTAssert(dimensions == FieldLayout_t::dimensions);
    }

  // The main constructor that takes a mesh and a field layout object.
  // This class will make copies of these objects.

  SpatialLayout(const Mesh_t &mesh, const FieldLayout_t &layout)
    : mesh_m(mesh),
      fieldLayout_m(layout)
    {
      // The Mesh and Layout dimensions must be consistent
      
      CTAssert(dimensions == FieldLayout_t::dimensions);
    }

  // Copy constructor.

  SpatialLayout(const This_t &s)
    : mesh_m(s.mesh()),
      fieldLayout_m(s.layout())
    {
      // The Mesh and Layout dimensions must be consistent
      
      CTAssert(dimensions == FieldLayout_t::dimensions);
    }

  // Assignment operator

  inline This_t &operator=(const This_t &s)
    {
      mesh_m = s.mesh();
      fieldLayout_m = s.layout();

      return *this;
    }

  void initialize(const This_t &s)
    {
      mesh_m = s.mesh();
      fieldLayout_m = s.layout();
    }

  void initialize(const Mesh_t &mesh, const FieldLayout_t &layout)
    {
      mesh_m = mesh;
      fieldLayout_m = layout;
    }

  bool initialized() const
    {
      return fieldLayout_m.initialized();
    }

  //============================================================
  // Destructor
  //============================================================

  ~SpatialLayout()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  // Return the field layout we're using

  inline const FieldLayout_t &layout() const
    {
      return fieldLayout_m;
    }

  // Return the mesh we're using

  inline const Mesh_t &mesh() const
    {
      return mesh_m;
    }

  // Return the number of patches used by the field layout.  All
  // particle layout objects must provide these methods.

  inline int patchesGlobal() const
    {
      return layout().sizeGlobal();
    }

  inline int patchesLocal() const
    {
      return layout().sizeLocal();
    }

  inline int patchesRemote() const
    {
      return layout().sizeRemote();
    }


  //============================================================
  // Attribute layout initialization
  //============================================================

  // The following method is called by the user of this class to
  // initialize a layout object for attributes that will need to
  // be kept organized by this layout.  For spatial layout, this
  // makes sure the attribute layout has the same number of patches
  // as the array layout, located on the same contexts and with the
  // same processor affinity.  The attribute layout is initialized
  // with zero elements in each patch.

  template<class AL>
  void initializeAttributeLayout(AL &attriblayout)
    {
      // SpatialLayout is given an array layout for a multi-dimensional
      // array of data, that can have several patches.  We use that
      // with a SpatialPartition object to initialize the provided
      // attribute layout object properly.  SpatialPartition will add in
      // domains to attributelayout that are initially empty, have the
      // same total number as the array layout, and the same memory
      // affinity.  Later, the user will add in elements to that layout
      // via create() operations.

      typedef typename AL::Domain_t ALDomain_t;
      attriblayout.initialize(ALDomain_t(), SpatialPartition<FL>(layout()),
                              DefaultSPmapper(layout()));
    }


  //============================================================
  // Particle patch location calculation
  //============================================================

  // Calculate the patch ID that each particle should have in the
  // provided arrays.  The first argument is the global patch ID number
  // for the particles that are provided in this call, and the second
  // argument is an Array storing particle positions for that particular
  // patch.  This routine calculates the patch ID for each of the particles,
  // and stores that patch ID in the third argument (an array of the same
  // length that stores patch ID's).  The final argument is an Array with
  // length equal to the number of patches; this routine should increment
  // the value for the Nth patch in that array by the number of paricles that
  // this routine determines goes in that patch.  Finally, this returns
  // the total number of particles that are destined for patches OTHER than
  // the patch mentioned as the first argument.

  template <class Attr>
  Size_t findPatchNumber(int lid, int gid, const Attr &pos,
			 MoveArray_t &movepid, AmountArray_t &moveamount)
  {
    // Find the current size of this patch

    Size_t size = pos.domain().size();

    // Find the local domain for this patch in the mesh

    Domain_t ldom = layout().patchDomain(lid);

    // Create a debug output stream, used if POOMA_PATCHSWAPLAYOUT_DBG
    // is set properly.
    // NOTE: the POOMA_PATCHSWAPLAYOUT_DBG macro is defined in
    // PatchSwapLayout.h

    POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("SpatialLayout::findPatchNumber");)
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg.setOutputContext(-1);)
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finding patch numbers for "<< size)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " particles in local patch " << lid)
    POOMA_PATCHSWAPLAYOUT_DBG(       << ", global patch " << gid)
    POOMA_PATCHSWAPLAYOUT_DBG(       << std::endl;)

    Size_t totmove = 0;

    for (int i = 0; i < size; ++i)
    {
      // Convert the particle position to an index into the Field's
      // domain using the Mesh, and from this get the global
      // patch ID for where it goes.

      int npid = gid;
      Loc<dimensions> ploc = mesh_m.cellContaining(pos(i));
      if (!contains(ldom,ploc))
      {
	// Calculate the patch ID using the mesh and field layout

	npid = fieldLayout_m.globalID(ploc);
	PAssert(npid != gid);

	// Save the patch ID for this particle ... we know that npid
	// is not equal to pid since the particle is outside our domain

	moveamount(npid) += 1;
	++totmove;
      }

      // We must update the movepid array for all particles, even local ones.

      PAssert(npid >= 0 && npid < patchesGlobal());
      movepid(i) = npid;

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Examining particle " << i)
      POOMA_PATCHSWAPLAYOUT_DBG(       << " with pos = " << pos(i))
      POOMA_PATCHSWAPLAYOUT_DBG(       << ": Global PatchID = " << npid)
      POOMA_PATCHSWAPLAYOUT_DBG(       << std::endl;)
    }

    // Return the total number of particles to locate on other patches.

    return totmove;
  }


  //============================================================
  // I/O
  //============================================================

  template <class Out>
  void print(Out &o) const
    {
      o << "SpatialLayout:\n";
      o << "    Field Layout = " << fieldLayout_m << "\n";
      o << "    Mesh Origin = " << mesh_m.origin() << "\n";
      o << "    Mesh Spacings = " << mesh_m.spacings() << "\n";
    }

private:
  //============================================================
  // Private data storage
  //============================================================

  // The mesh for this spatial layout

  Mesh_t mesh_m;

  // The field layout to use for this spatial layout

  FieldLayout_t fieldLayout_m;
};


//-----------------------------------------------------------------------------
//
// A specialization of the Inform traits used to say that SpatialLayout has
// a print method.
//
//-----------------------------------------------------------------------------

template <class M, class FL>
std::ostream &operator<<(std::ostream &o, const SpatialLayout<M,FL> &sl)
{
  sl.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_SPATIAL_LAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SpatialLayout.h,v $   $Author: richard $
// $Revision: 1.36 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
