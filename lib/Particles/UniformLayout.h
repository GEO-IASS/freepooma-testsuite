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

#ifndef POOMA_PARTICLES_UNIFORM_LAYOUT_H
#define POOMA_PARTICLES_UNIFORM_LAYOUT_H

//-----------------------------------------------------------------------------
// Classes:
//   UniformLayout
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Particles
 * @brief
 * UniformLayout is a very simple particle layout class that is used to
 * determine where particles will be located in a parallel environment.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/PatchSwapLayout.h"
#include "Partition/UniformMapper.h"
#include "Partition/GridPartition.h"
#include "Utilities/PAssert.h"
#include "Pooma/Pooma.h"
#include "Domain/Loc.h"

#include <iosfwd>

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------

class UniformLayout;

/**
 * UniformLayout is a very simple particle layout class that is used to
 * determine where particles will be located in a parallel environment.
 * UniformLayout simply tries to keep the same number of particles on each
 * patch.
 * UniformLayout is a PatchSwapLayout, and inherits the main "swap" routine
 * from that class.  It provides a "findPatchNumber" routine that calculates
 * patch numbers based on getting an equal number of particles on each
 * patch.
 */

class UniformLayout : public PatchSwapLayout<UniformLayout>
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Utility typedef to refer to ourselves and our base class

  typedef UniformLayout             This_t;
  typedef PatchSwapLayout<This_t>   Base_t;

  // Storage for hold particle amounts ... this should be the same
  // as the typedef in Particles.h

  typedef Base_t::Size_t            Size_t;

  // The type of array used to store amounts to move to other patches

  typedef Base_t::AmountArray_t     AmountArray_t;

  // The type of array used to store patch ID's and indices

  typedef Base_t::MoveArray_t       MoveArray_t;


  //============================================================
  // Constructors
  //============================================================

  // The default constructor 

  UniformLayout()
    : numPatches_m(Pooma::contexts()),
      numLocalPatches_m(1)
    {
    }

  // The main constructor, which takes the number of patches.

  UniformLayout(int numPatches)
    : numPatches_m(numPatches),
      numLocalPatches_m(numPatches / Pooma::contexts())
    {
      int remainder = numPatches_m % Pooma::contexts();
      if (Pooma::context() < remainder)
	++numLocalPatches_m;
    }

  // Copy constructor.

  UniformLayout(const This_t &s)
    : numPatches_m(s.numPatches_m),
      numLocalPatches_m(s.numLocalPatches_m)
    {
    }

  // Assignment operator

  This_t& operator=(const This_t &s)
    {
      numPatches_m = s.numPatches_m;
      numLocalPatches_m = s.numLocalPatches_m;
      return *this;
    }

  // Initialization routines, taking same arguments as non-
  // default constructor.

  void initialize(int numPatches)
    {
      // determine local number of patches, assuming an even
      // distribution of patches across contexts.

      numPatches_m = numPatches;
      numLocalPatches_m = numPatches_m / Pooma::contexts();
      int remainder = numPatches_m % Pooma::contexts();
      if (Pooma::context() < remainder)
	++numLocalPatches_m;
    }

  void initialize(const This_t &s)
    {
      numPatches_m = s.numPatches_m;
      numLocalPatches_m = s.numLocalPatches_m;
    }

  //============================================================
  // Destructor
  //============================================================

  ~UniformLayout()
    {
    }

  //============================================================
  // Accessors
  //============================================================

  // Return whether this object has been initialized.  Since this is
  // so simple, it always is initialized.

  inline bool initialized() const
    {
      return true;
    }

  // Return the number of patches used.  All
  // particle layout objects must provide these methods.

  inline int patchesGlobal() const
    {
      return numPatches_m;
    }

  inline int patchesLocal() const
    {
      return numLocalPatches_m;
    }

  inline int patchesRemote() const
    {
      return (numPatches_m-numLocalPatches_m);
    }

  //============================================================
  // Attribute layout initialization
  //============================================================

  // The following method is called by the user of this class to
  // initialize a layout object for attributes that will need to
  // be kept organized by this layout.  For uniform layout, this
  // makes sure the attribute layout has the given number of patches
  // with roughly the same number of particles in each patch.
  // The attribute layout is initialized with zero elements in each patch.

  template <class AL>
  void initializeAttributeLayout(AL& attriblayout)
    {
      // UniformLayout is given only the number of patches.  We use that
      // with a GridPartition object to initialize the provided
      // attribute layout object properly.  GridPartition will add in
      // domains to attributelayout that are initially empty and have the
      // same total number as the uniform layout.  Later, the user will add
      // in elements to that layout via create() operations.
      // Note: We should use a GridPartition here rather than a
      // UniformGridPartition because we don't want the code to bomb
      // if the number of particles does not divide *exactly* evenly
      // across the given number of patches.
      // Also, explicitly use the UniformMapper here, which will 
      // distribute the patches across multiple contexts in an even manner.

      typedef typename AL::Domain_t ALDomain_t;
      Loc<1> blocks(numPatches_m);
      GridPartition<1> gpar(blocks);
      UniformMapper cmap(gpar);
      attriblayout.initialize(ALDomain_t(), gpar, cmap);
      // check that we got the expected number of patches
      PAssert(attriblayout.sizeGlobal() == numPatches_m);
      PAssert(attriblayout.sizeLocal() == numLocalPatches_m);
    }

  //============================================================
  // Particle patch location calculation
  //============================================================

  // Calculate the patch ID that each particle should have in the
  // provided arrays.  The first argument is that patch ID number
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

  template<class Attr>
  Size_t findPatchNumber(int lid, int gid, const Attr & /* pos */,
			 MoveArray_t &movepid, AmountArray_t &moveamount)
  {
    // Find the current size of this patch

    Size_t size = patchInfo(lid).size();

    // Create a debug output stream, used if POOMA_PATCHSWAPLAYOUT_DBG
    // is set properly.
    // NOTE: the POOMA_PATCHSWAPLAYOUT_DBG macro is defined in
    // PatchSwapLayout.h

    POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("UniformLayout::findPatchNumber");)
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg.setOutputContext(-1);)
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finding patch numbers for " << size)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " particles in local patch " << lid)
    POOMA_PATCHSWAPLAYOUT_DBG(       << ", global patch " << gid << std::endl;)

    // Compute the offset of this patch in the current distribution, and
    // the total number of particles.

    Size_t offset = 0, totalsize = 0;
    int myContext = Pooma::context();
    int numContexts = Pooma::contexts();
    for (int c = 0; c < numContexts; ++c)
    {
      totalsize += contextSize(c);
      if (c < myContext)
        offset += contextSize(c);
    }

    for (int p=0; p < lid; ++p)
    {
      Size_t psize = patchInfo(p).size();
      offset += psize;
    }

    // Compute the number of particles to put on each patch

    Size_t sizePerPatch = totalsize / patchesGlobal();

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Trying to put " << sizePerPatch)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " particles on each patch;")
    POOMA_PATCHSWAPLAYOUT_DBG(       << " global patch " << gid)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " has current offset = " << offset)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " out of " << totalsize)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " total particles." << std::endl;)

    // Now loop through the particles and find where they go.  Store
    // the global patch ID in movePatch.  If it does not move, store our own
    // global patch ID.  We basically shuffle the particles such that we have
    // sizePerPatch particles on each patch.

    Size_t totmove = 0;

    for (int i = 0; i < size; ++i)
    {
      int npid = (i+offset);
      if (sizePerPatch > 0)
        npid /= sizePerPatch;

      // check for a leftover particle

      if (npid >= patchesGlobal())
	npid = (i+offset) - (sizePerPatch*patchesGlobal());

      // Make sure this is kosher

      PAssert(npid >= 0 && npid < patchesGlobal());

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Particle " << i)
      POOMA_PATCHSWAPLAYOUT_DBG(       << " on global patch " << gid)
      POOMA_PATCHSWAPLAYOUT_DBG(       << " should go on global patch ")
      POOMA_PATCHSWAPLAYOUT_DBG(       << npid << std::endl;)

      // Save the patch ID for this particle

      movepid(i) = npid;
      if (npid != gid)
      {
	moveamount(npid) += 1;
	++totmove;
      }
    }

    // Return the total number of particles to locate on other patches.

    return totmove;
  }

    
  //============================================================
  // I/O
  //============================================================

  template <class Out>
  void print(Out& o) const
    {
      o << "UniformLayout:\n";
      o << "Number of global patches = " << numPatches_m << "\n";
      o << "Number of local patches = " << numLocalPatches_m << "\n"; 
    }

private:

  //============================================================
  // Private data storage
  //============================================================

  // The total number of patches

  int numPatches_m;

  // The number of patches local to this context

  int numLocalPatches_m;
};


//-----------------------------------------------------------------------------
//
// A specialization of the Inform traits used to say that SpatialLayout has
// a print method.
//
//-----------------------------------------------------------------------------

inline
std::ostream& operator<<(std::ostream& o, const UniformLayout& ul)
{
  ul.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_UNIFORM_LAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformLayout.h,v $   $Author: richard $
// $Revision: 1.25 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
