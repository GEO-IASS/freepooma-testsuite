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
// Classes:
//   PatchSwapInfo
//   PatchSwapFunctor<P>
//   PatchSwapLayout<L>
//   PSwapPack
//-----------------------------------------------------------------------------

#ifndef POOMA_PARTICLES_PATCH_SWAP_LAYOUT_H
#define POOMA_PARTICLES_PATCH_SWAP_LAYOUT_H

/** @file
 * @ingroup Particles
 * @brief
 * These classes define storage and functors for use in implementing
 * particle layout classes that manage particles by swapping them from
 * one patch to another.
 *
 * They work with multi-patch attributes by
 * computing on what patch each particle should be, and moving them around
 * when asked.  PatchSwapInfo stores per-patch info for use in the 
 * swapping operation, while PatchSwapFunctor's are the patch-functor invoked
 * to do the swapping for a single patch in a multi-threaded environment
 * for a particles object of type P.  PatchSwapLayout is a base class for
 * all swapping particle layouts, and defines the main "swap" routine.  It
 * is templated on the type of derived class that it works with, so
 * that it can provide that particular type to the swap functor.  Derived
 * classes should implement any specialized behavior and storage type 
 * type require, plus provide an implementation of a method
 *
 *    Size_t findPatchNumber(lid, gid, posarray, patcharray, amountarray)
 *
 * to calculate the global patch ID of all particles in a given array at the
 * present time.  This is used by the swap functor to determine where
 * things should go.  The different particle layout classes just
 * provide the algorithm that determines where things go.
 *
 * PSwapPack is a structure used to store the items 
 * needed to exchange particle data between patches on different contexts.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include "Pooma/Pooma.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Layout/DynamicEvents.h"
#include "Utilities/Inform.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
// Description of PatchSwapInfo:
//
// PatchSwapInfo store information for each patch in a multiple-patch
// particle layout that is needed to swap particles around.  It provides
// storage for each patch for the number of particles to swap from that
// patch to others, the list of swapped particle indices, etc.  The
// different particle layouts just need to apply their particular algorithms
// to determine the proper patch for each particle.
//
//-----------------------------------------------------------------------------

class PatchSwapInfo
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Utility typedef to refer to ourselves

  typedef PatchSwapInfo                      This_t;

  // Storage for hold particle amounts ... this should be the same
  // as the typedef in Particles.h

  typedef unsigned long                      Size_t;

  // The type of array used to store amounts to move to other patches

  typedef Array<1, int, Brick>               AmountArray_t;

  // The type of array used to store patch ID's and indices

  typedef Array<1, int, Brick>               MoveArray_t;


  //============================================================
  // Constructors
  //============================================================

  // Main constructor.

  PatchSwapInfo(int patchesLocal, int patchesRemote = 0)
    : send_m(0), copyPatches_m(0), destroySize_m(0), msgReceived_m(0)
    {
      initialize(patchesLocal,patchesRemote);
    }

  // Default constructor.  Must call initialize() before using object.

  PatchSwapInfo()
    : send_m(0), patchesLocal_m(0), patchesGlobal_m(0),
      copyPatches_m(0), destroySize_m(0), msgReceived_m(0)
    {
    }

  // Initialization

  void initialize(int patchesLocal, int patchesRemote)
    {
      PAssert(patchesLocal > 0);

      // Delete old particle send list storage

      if (send_m != 0)
	  delete [] send_m;

      // Save the number of local and global patches we are working with here.

      patchesLocal_m = patchesLocal;
      patchesGlobal_m = patchesLocal+patchesRemote;

      // Allocate storage for a set of arrays to store send lists

      send_m = new MoveArray_t[patchesGlobal_m];

      // Initialize amount array, which should be of length
      // patchesLocal + patchesRemote

      amount_m.initialize(patchesGlobal_m);
    }

  //============================================================
  // Destructor
  //============================================================

  ~PatchSwapInfo()
    {
      // delete particle send list storage

      if (send_m != 0)
	delete [] send_m;
    }


  //============================================================
  // Accessors and mutators
  //============================================================

  // Return the number of local patches this is set up for.

  int patchesLocal() const
    {
      return patchesLocal_m;
    }

  // Return the number of global patches this is set up for.

  int patchesGlobal() const
    {
      return patchesGlobal_m;
    }

  // Return the current patch size.

  Size_t size() const
    {
      return patchSize_m;
    }
  
  // Set the current patch size.

  void setSize(Size_t newsize)
    {
      patchSize_m = newsize;
    }

  // Return the number of elements to destroy after copy.

  Size_t destroySize() const
    {
      return destroySize_m;
    }
  
  // Set the number of elements to destroy after copy.

  void setDestroySize(Size_t newsize)
    {
      destroySize_m = newsize;
    }

  // Return the number of patches that have data that must be
  // copied into this patch.

  int copyPatches() const
    {
      return copyPatches_m;
    }

  // Set the number of patches that have data that must be
  // copied into this patch.

  void setCopyPatches(int p)
    {
      copyPatches_m = p;
    }

  // Return the number of particle swap messages that this patch has received.

  const int& msgReceived() const
    {
      return msgReceived_m;
    }

  int& msgReceived()
    {
      return msgReceived_m;
    }

  // Return the array used to store patch amounts

  inline AmountArray_t &amount()
    {
      return amount_m;
    }

  // Return the array used to store destroy lists

  inline MoveArray_t &destroyIndices()
    {
      return destroy_m;
    }

  // Return the array used to store send list for global patch p

  inline MoveArray_t &sendIndices(int p)
    {
      PAssert(send_m != 0);
      PAssert(p >= 0 && p < patchesGlobal_m);
      return send_m[p];
    }

private:

  // Patch info arrays

  AmountArray_t amount_m;
  MoveArray_t destroy_m;
  MoveArray_t *send_m;

  // Current number of patches (local and global)

  int patchesLocal_m;
  int patchesGlobal_m;

  // Current size of this local patch

  Size_t patchSize_m;

  // Most recent destroy count

  Size_t destroySize_m;

  // Most recent number of patches that will copy into this one

  int copyPatches_m;

  // Number of messages received from remote patches

  int msgReceived_m;
};


//-----------------------------------------------------------------------------
//
// Description of PatchSwapFunctor<P>:
//
// This simple functor is used to scan particles for an individual
// patch in a given particle objects's attribute array, to extend attribute
// storage, to copy particles, and to destroy old lists.
// It has one method "apply" that is invoked
// with a view on each patch, and a local patch ID number.  Each invocation
// will process that patch's data, and do the following steps, based on
// the operational "mode":
// MODE == syncScan:
//   Apply the BC's for the patch
//   Do deferred destroys for the patch
// MODE == swapScan:
//   Scan particles to find what patch they should go on
// MODE == swapSend:
//   Send out particle data that belongs on other patches that are remote
// MODE == swapExtend:
//   Extend the storage of all attributes to accomodate copied-in values.
// MODE == swapCopy:
//   Copy values from other local patches to the end of this patch.
// MODE == swapReceive:
//   Receive particle data from remote patches
// MODE == swapDestroy:
//   Destroy elements that were copied or sent out of this patch to others.
//
// After this, this functor will increment a semaphore, and finish.  The
// remaining functors will complete the swap operation, and they can only
// continue on after this step is finished.
//
// This functor stores a reference to a particle layout object used in the 
// swapping process, and the particles object being swapped.  The layout
// is accessed for information on other patches and to determine
// exactly which patch each particle should go on.  The Layout template
// parameter should be for a particle layout object that has methods
// 'patchInfo(int pid)' and others.
//
//-----------------------------------------------------------------------------

template<class P>
class PatchSwapFunctor
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumeration of possible modes of operation

  enum { syncScan, swapScan, swapSend, swapExtend, swapCopy,
         swapReceive, swapDestroy };

  // Utility typedef to refer to ourselves

  typedef PatchSwapFunctor<P>              This_t;

  // The Particles object that we use

  typedef P                                Particles_t;

  // The particle layout type that we use

  typedef typename P::ParticleLayout_t     Layout_t;

  // The data type for patch ID numbers

  typedef DynamicEvents::PatchID_t         PatchID_t;

  // Storage for hold particle amounts ... this is int for now due to
  // a bug in taking views of Array's with things other than int

  typedef int                              Size_t;

  // The type of array used to store amounts to move to other patches

  typedef PatchSwapInfo::AmountArray_t     AmountArray_t;

  // The type of array used to store patch ID's and indices

  typedef PatchSwapInfo::MoveArray_t       MoveArray_t;


  //============================================================
  // Constructors
  //============================================================

  // Construct with a reference to a particle layout and a particles
  // object.  We use these objects in the apply function.

  inline PatchSwapFunctor(Layout_t &layout,
			  Particles_t &particles,
			  int mode)
    : layout_m(layout),
      particles_m(particles),
      mode_m(mode)
    {
    }

  // Copy constructor.

  inline PatchSwapFunctor(const This_t &model)
    : layout_m(model.layout_m),
      particles_m(model.particles_m),
      mode_m(model.mode_m)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  inline ~PatchSwapFunctor()
    {
    }


  //============================================================
  // Apply method
  //============================================================

  // The apply method, invoked by the patch functor evaluator
  // for each patch in an array.  The first argument is a view
  // on one of the array's patches, the second argument is the
  // local ID number of that patch.  The "Patch" is part of the
  // position array for the particles.

  template<class ArrayPatch>
  inline void apply(const ArrayPatch& a, PatchID_t pid) const
    {
      if (mode_m == syncScan)
	performSync(a, pid);
      else if (mode_m == swapScan)
	performScan(a, pid);
      else if (mode_m == swapSend)
	performSend(a, pid);
      else if (mode_m == swapExtend)
	performExtend(a, pid);
      else if (mode_m == swapCopy)
	performCopy(a, pid);
      else if (mode_m == swapReceive)
	performReceive(a, pid);
      else if (mode_m == swapDestroy)
	performDestroy(a, pid);
      else
	PInsist(false, "Unknown mode in PatchSwapFunctor::apply");
    }

private:
  //============================================================
  // Private methods
  //============================================================

  // Carry out the actions from the apply method for the different
  // operation modes.

  template<class ArrayPatch>
  void performSync(const ArrayPatch& a, PatchID_t pid) const;

  template<class ArrayPatch>
  void performScan(const ArrayPatch& a, PatchID_t pid) const;

  template<class ArrayPatch>
  void performSend(const ArrayPatch& a, PatchID_t pid) const;

  template<class ArrayPatch>
  void performExtend(const ArrayPatch& a, PatchID_t pid) const;

  template<class ArrayPatch>
  void performCopy(const ArrayPatch& a, PatchID_t pid) const;

  template<class ArrayPatch>
  void performReceive(const ArrayPatch& a, PatchID_t pid) const;

  template<class ArrayPatch>
  void performDestroy(const ArrayPatch& a, PatchID_t pid) const;


  //============================================================
  // Private data
  //============================================================

  // A reference to the particle layout we should use in our apply function.

  Layout_t &layout_m;

  // A reference to the particles object we should be swapping for.

  Particles_t &particles_m;

  // Operating mode

  int mode_m;

};


//-----------------------------------------------------------------------------
//
// Description of PatchSwapLayout<P>:
//
// PatchSwapLayout<L> is a base class for all particle layout classes
// of type L that want to manage a particle distribution by swapping
// particles from one patch to another in a multi-patch attribute layout.
// The base classes provide the specialized algorithm to determine what
// patch each particle should go in; this class provides the generic
// algorithm to move them around once the destination patch ID is known
// for each particle.  It provides the standard routine
//
//   template<class P, class Attr>
//   void swap(P &particles, Attr &attribute)
//
// that swaps particles around based on a storage attribute.
//
// Derived class must, in addition to storing their particular data,
// provide the following methods:
//
// int patchesGlobal() const;
//
// int patchesLocal() const;
//
// int patchesRemote() const;
//
// template<class Attr>
// Size_t findPatchNumber(int lid, int gid, Attr &patch,
//                        MoveArray_t movepid, AmountArray_t moveamount)
// ==> This gets called by PatchSwapFunctor to find the proper global patch
// ID's of all the particles in the individual patch given in the second
// argument.  This routine should also set moveamount to the number of
// particles to move from this patch to other patches, and return the
// total number of particles to move.
//
// template<class AL>
// void initializeAttributeLayout(AL &attriblayout);
// ==> This method is called by the user of this class to
// initialize a layout object for attributes that will need to
// be kept organized by this layout.  E.g., for spatial layout, this
// makes sure the attribute layout has the same number of patches
// as the array layout, located on the same contexts and with the
// same processor affinity.  The attribute layout is initialized
// with zero elements in each patch.
//
//-----------------------------------------------------------------------------

template<class L>
class PatchSwapLayout
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Utility typedef to refer to ourselves

  typedef PatchSwapLayout<L>                 This_t;

  // The type of layout that is, generally, the derived class type

  typedef L                                  Layout_t;

  // Storage for hold particle amounts ... this should be the same
  // as the typedef in Particles.h

  typedef PatchSwapInfo::Size_t              Size_t;

  // The type of array used to store amounts to move to other patches

  typedef PatchSwapInfo::AmountArray_t       AmountArray_t;

  // The type of array used to store patch ID's and indices

  typedef PatchSwapInfo::MoveArray_t         MoveArray_t;


  //============================================================
  // Constructors
  //============================================================

  PatchSwapLayout()
    : patchInfo_m(0)
    {
      contextSizes_m.initialize(Pooma::contexts());
    }

  // The main constructor takes a reference to the Layout_t type
  // that we will use in the swap() routine.

  PatchSwapLayout(Layout_t &layout)
    : patchInfo_m(0)
    {
      contextSizes_m.initialize(Pooma::contexts());
    }

  // Copy constructor.

  PatchSwapLayout(const This_t &p)
    : patchInfo_m(0)
    {
      contextSizes_m.initialize(Pooma::contexts());
    }


  //============================================================
  // Destructor
  //============================================================

  // The destructor needs to clean up the patch-specific storage

  ~PatchSwapLayout()
    {
      if (patchInfo_m != 0)
        delete [] patchInfo_m;
    }


  //============================================================
  // Accessors
  //============================================================

  // Return the patch info for the Nth local patch.

  inline PatchSwapInfo &patchInfo(int pid)
    {
      PAssert(patchInfo_m != 0);
      return patchInfo_m[pid];
    }

  // Return the current size stored for the given context.

  inline Size_t contextSize(int c) const { return contextSizes_m(c); }

  //============================================================
  // The primary purpose for particle layout classes: swap
  //============================================================

  // Sync up a particle object, by letting all contexts know what the
  // other contexts have created and destroyed, plus move particles
  // around to their proper place.  sync() will also first apply boundary
  // conditions and perform deferred destroys. The no-argument version of
  // swap() should be used for particle layouts that do not need to know
  // the position or other information to do their job.
  // Since patch-swap layouts requires position information, if this version
  // is called it is an error.

  template<class P>
  void sync(P &)
    {
      PInsist(false, "You must call PatchSwapLayout::sync with positions.");
    }

  template<class P>
  void swap(P &)
    {
      PInsist(false, "You must call PatchSwapLayout::swap with positions.");
    }

  // Sync up this object, but also specify an object that can be used
  // as an attribute for the particle layout's distribution algorithm.
  // For example, this could be a position attribute when using a
  // spatial layout.  The object should have the same interface and
  // layout as the other particle attributes.  If more than one 
  // attribute is needed (e.g., separate components of the position),
  // they must be packaged up as a single attribute using, for example,
  // the global function "vector" that we haven't written yet.
  // sync() will do the same as swap, but also first apply boundary
  // conditions and do deferred destroys.

  template<class P, class A>
  void sync(P &particles, const A &pos)
    {
      performSync(particles, pos, true);
    }
  
  template<class P, class A>
  void swap(P &particles, const A &pos)
    {
      performSync(particles, pos, false);
    }

private:    
  //============================================================
  // Private methods
  //============================================================

  // Calculate the total size of each local patch in the particle's
  // attribute layout, and store them.

  template <class P>
  void findCurrentSizes(const P& particles);

  // Actually carry out the sync/swap operation.  If the final
  // argument is true, this does a full sync(), otherwise it just
  // does a swap.

  template<class P, class A>
  void performSync(P &particles, const A &pos, bool dosync);


  //============================================================
  // Private data storage
  //============================================================

  // Information about the patches, used in swapping

  PatchSwapInfo *patchInfo_m;

  // Current total size of all local patches on each context

  Array<1,Size_t,Brick> contextSizes_m;
};


#if POOMA_MESSAGING

#include "Tulip/Messaging.h"

//-----------------------------------------------------------------------------
//
// Description of PSwapPack
//
// This structure stores everything needed to send a set of particle
// data from one patch to another (remote) patch and to receive the 
// data on the other side.  We provide a specialization of the class
// template Serialize for this structure that contains the necessary
// pack() and unpack() functions.  PSwapPack stores the Particles 
// object, the local patch ID of the sending or receiving patch and
// the indirection list of local indices of particles to be sent.
// This indirection list is not used on the receiving side, because 
// we will create space for any new particles to be received.
//
//-----------------------------------------------------------------------------

template <class P>
struct PSwapPack
{
  // typedefs
  typedef IndirectionList<int> List_t;

  // constructors
  PSwapPack()
    : patchID_m(0), particles_m(0), 
      list_m(0), buffer_m(0) {}
  PSwapPack(int patchID, P& particles)
    : patchID_m(patchID), particles_m(&particles), 
      list_m(0), buffer_m(0) {}
  PSwapPack(int patchID, P& particles, List_t& list)
    : patchID_m(patchID), particles_m(&particles), 
      list_m(&list), buffer_m(0) {}

  // member data
  int patchID_m;
  P* particles_m;
  List_t* list_m;
  char* buffer_m;
};

namespace Cheetah {

template <class P>
class Serialize< CHEETAH, PSwapPack<P> >
{
public:
  static inline
  int size(const PSwapPack<P>& pack)
  {
    // first get pack size of the number of particles
    typedef typename P::Size_t Size_t;
    int nBytes = 0;
    nBytes += Serialize<CHEETAH,Size_t>::size(Size_t());

    // if there are particles to send, pack data for each attribute
    Size_t size = pack.list_m->size();
    if (size > 0)
    {
      int attribs = pack.particles_m->attributes();
      for (int i=0; i<attribs; ++i)
	nBytes += pack.particles_m->attribute(i)->packSize(size);
    }

    // need to add space for extra int to hold size of buffer in bytes
    nBytes += Serialize<CHEETAH,int>::size(int());
    return nBytes;
  }
  
  static inline 
  int pack(const PSwapPack<P>& pack, char *buffer)
  {
    // save original starting point of data buffer
    char* bufstart = buffer;

    // skip over space for storing buffer size; we'll fill that last
    int change, nBytes = 0;
    buffer += Serialize<CHEETAH,int>::size(nBytes);

    // first pack the number of particles
    typedef typename P::Size_t Size_t;
    Size_t size = pack.list_m->size();
    change = Serialize<CHEETAH,Size_t>::pack(size,buffer);
    buffer += change;
    nBytes += change;

    // if there are particles to send, pack data for each attribute
    if (size > 0)
    {
      int attribs = pack.particles_m->attributes();
      for (int i=0; i<attribs; ++i)
      {
	change = pack.particles_m->attribute(i)->
	  pack(pack.patchID_m, *(pack.list_m), buffer);
	buffer += change;
	nBytes += change;
      }
    }

    // now let's pack the number of bytes for data at the top of the buffer
    change = Serialize<CHEETAH,int>::pack(nBytes,bufstart);
    nBytes += change;
    return nBytes;
  }

  static inline 
  int unpack(PSwapPack<P>* &pack, char *buffer)
  {
    // allocate the PSwapPack object to store the received buffer
    pack = new PSwapPack<P>();

    // unpack the first item, which is the length of the message in bytes
    int *plen;
    int change, nBytes = 0;
    change = Serialize<CHEETAH,int>::unpack(plen,buffer);
    buffer += change;
    nBytes += change;

    // store pointer to message data in PSwapPack object
    pack->buffer_m = buffer;

    // return total size of unpacked message in bytes
    nBytes += *plen;
    return nBytes;
  }

  static inline 
  void cleanup(PSwapPack<P> *pack)
  {
    delete pack;
  }
};

} // namespace Cheetah

template <class P>
void pSwapUnpackFunc(PSwapPack<P> *pack, PSwapPack<P> &packbuf)
{
  // get reference to the data buffer from the PSwapPack buffer object
  char* &buffer = packbuf.buffer_m;

  // now we can unpack data into locally provided PSwapPack object
  // first unpack the number of particles to unpack
  typedef typename P::Size_t Size_t;
  Size_t *psize;
  int change;
  change = Cheetah::Serialize<Cheetah::CHEETAH,Size_t>::unpack(psize,buffer);
  buffer += change;

  // if there are particles to receive, unpack data for each attribute
  if (*psize > 0)
  {
    // first create new elements to hold received particle data
    // final argument indicates we do not want to do a global renumbering
    pack->particles_m->create(*psize, pack->patchID_m, false);

    // next, create Interval for *zero-based* domain of new particles
    Interval<1> patchDom = 
      pack->particles_m->attributeLayout().patchDomain(pack->patchID_m);
    int offset = patchDom.first();
    int last = patchDom.last();
    int first = last - *psize + 1;
    Interval<1> rdomain = Interval<1>(first-offset,last-offset);

    // now loop over attributes and unpack data into new elements
    int attribs = pack->particles_m->attributes();
    for (int i=0; i<attribs; ++i)
    {
      change = pack->particles_m->attribute(i)->
	unpack(pack->patchID_m, rdomain, buffer);
      buffer += change;
    }
  }

  // Note that we have received another message
  pack->particles_m->particleLayout().
    patchInfo(pack->patchID_m).msgReceived() += 1;
}

#endif // POOMA_MESSAGING

// Include out-of-line definitions

#include "Particles/PatchSwapLayout.cpp"


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_PATCH_SWAP_LAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PatchSwapLayout.h,v $   $Author: richard $
// $Revision: 1.23 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
