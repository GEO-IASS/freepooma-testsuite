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
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/PatchSwapLayout.h"
#include "Evaluator/PatchFunction.h"
#include "Tulip/RemoteProxy.h"

//-----------------------------------------------------------------------------
// A debugging macro; to enable printing debug messages during swap,
// set POOMA_PATCHSWAPLAYOUT_DBG(x) to x
//-----------------------------------------------------------------------------

#define POOMA_PATCHSWAPLAYOUT_DBG(x)


//============================================================
// Non-inline methods for PatchSwapFunctor
//============================================================


//-----------------------------------------------------------------------------
// The sync method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will apply
// particle BC's and do deferred destroys.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performSync(const ArrayPatch& a, PatchID_t lid) const
{
  // Get the current number of particles in this patch, and the
  // number of local patches.

  Size_t size = a.domain().size();
  int patchesLocal = layout_m.patchesLocal();

  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(sync)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local ID = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << " out of " << patchesLocal)
  POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl << "The working size is ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << size << std::endl;)

  // Apply boundary conditions and do deferred destroy for this
  // patch, if we need to.

  if (size > 0)
  {
    // Cache the number of particles that need to be destroyed from 
    // the deferred destroy list

    Size_t todestroy = particles_m.deferredDestroyAmount(lid);

    // Apply BC's if there will be any particles left to do BC's on

    if (todestroy < size)
    {
      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Applying boundary conditions to ")
      POOMA_PATCHSWAPLAYOUT_DBG(  << "local patch " << lid)
      POOMA_PATCHSWAPLAYOUT_DBG(  << " ..." << std::endl;)
      
      particles_m.applyBoundaryConditions(lid);
      todestroy = particles_m.deferredDestroyAmount(lid);
    }

    // Do deferred destroy's now

    if (todestroy > 0)
    {
      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Doing deferred destroys for ")
      POOMA_PATCHSWAPLAYOUT_DBG(  << "local patch " << lid)
      POOMA_PATCHSWAPLAYOUT_DBG(  << ", todestroy = " << todestroy)
      POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

      size -= todestroy;
      particles_m.performDestroy(lid, false);

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Destroy complete, new size = ")
      POOMA_PATCHSWAPLAYOUT_DBG(  << size << std::endl;)
    }
  }
}

//-----------------------------------------------------------------------------
// The scan method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will compute
// the patch ID of all particles.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performScan(const ArrayPatch& a, PatchID_t lid) const
{
  int i, p;

  // Get global patch ID from local patch ID

  PatchID_t gid = 
    particles_m.attributeLayout().nodeListLocal()[lid]->globalID();

  // Get the current number of particles in this patch, and the
  // total number of patches.

  Size_t size = a.domain().size();
  int patchesLocal = layout_m.patchesLocal();
  int patchesGlobal = layout_m.patchesGlobal();

  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(scan)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local ID = " << lid)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " out of " << patchesLocal << std::endl)
  POOMA_PATCHSWAPLAYOUT_DBG(  << "The working size is " << size << std::endl;)

  // Initialize the move amounts array

  AmountArray_t& moveAmount = layout_m.patchInfo(lid).amount();
  for (p = 0; p < patchesGlobal; ++p)
    moveAmount(p) = 0;

  // Check for zero size, nothing to do if we have no particles.

  Size_t totmove = 0;

  if (size == 0)
  {
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Empty patch, we can just return.")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)
  }
  else
  {
    // Set up an Array to store the number of particles to move to
    // each of the other patches, and a starting index for send ops.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Initializing to scan " << size)
    POOMA_PATCHSWAPLAYOUT_DBG(  << " particles " << "in local patch " << lid)
    POOMA_PATCHSWAPLAYOUT_DBG(  << " of " << patchesLocal << std::endl;)

    // Set up a list that will indicate to what patch each particle
    // should move.  This is the same size as our current particle count.

    MoveArray_t& movePatch = layout_m.patchInfo(lid).destroyIndices();
    if (movePatch.domain().size() < size)
      movePatch.initialize(size);

    // Now loop through the particles and find where they go.  Store
    // the global patch ID in movePatch.  Also store the number of 
    // particles to move to other patches in moveAmount.  The layout
    // routine should return the total number of particles that must
    // be moved from this patch to others.

    totmove = layout_m.findPatchNumber(lid, gid, a, movePatch, moveAmount);

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "At end of scan: amounts = ")
    POOMA_PATCHSWAPLAYOUT_DBG(  << moveAmount << std::endl;)

    if (totmove > 0)
    {
      // Compute indirection lists for other patches.  
      // First, get and resize send arrays.

      for (p = 0; p < patchesGlobal; ++p)
      {
	Size_t ma = moveAmount(p);
	if (ma > 0)
	{
	  MoveArray_t& sendPatch = layout_m.patchInfo(lid).sendIndices(p);
	  if (sendPatch.domain().size() < ma)
	    sendPatch.initialize(ma);
	  moveAmount(p) = 0;  // we'll restore this later
	}
      }

      // Go through the particles, examine the destination patch ID, and
      // copy the index to that patch's send list.  Also reuse movePatch
      // for the final destroy list.

      int destroyed = 0;
      for (i = 0; i < size; ++i)
      {
	p = movePatch(i);
	if (p != gid)
	{
	  layout_m.patchInfo(lid).sendIndices(p)(moveAmount(p)) = i;
	  moveAmount(p) += 1;
	  movePatch(destroyed++) = i;
	}
      }
      PAssert(destroyed == totmove);

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finished finding local send lists.")
      POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl << "Redone move amounts = ")
      POOMA_PATCHSWAPLAYOUT_DBG(  << moveAmount << std::endl;)
      POOMA_PATCHSWAPLAYOUT_DBG(for (p=0; p < patchesGlobal; ++p))
      POOMA_PATCHSWAPLAYOUT_DBG({)
      POOMA_PATCHSWAPLAYOUT_DBG(  if (p != gid && moveAmount(p) > 0))
      POOMA_PATCHSWAPLAYOUT_DBG(  {)
      POOMA_PATCHSWAPLAYOUT_DBG(    int amt = moveAmount(p);)
      POOMA_PATCHSWAPLAYOUT_DBG(    dbgmsg << "  Send list for local patch ")
      POOMA_PATCHSWAPLAYOUT_DBG(      << lid << " --> global patch " << p)
      POOMA_PATCHSWAPLAYOUT_DBG(      << " = " << layout_m.patchInfo(lid).)
      POOMA_PATCHSWAPLAYOUT_DBG(      sendIndices(p)(Interval<1>(amt)))
      POOMA_PATCHSWAPLAYOUT_DBG(      << std::endl;)
      POOMA_PATCHSWAPLAYOUT_DBG(  })
      POOMA_PATCHSWAPLAYOUT_DBG(})
    }
  }

  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finished with scan on local patch ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << " after finding " << totmove)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " particles to move to other patches.")
  POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

  layout_m.patchInfo(lid).setDestroySize(totmove);
}


//-----------------------------------------------------------------------------
// The extend method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will find out
// what particles should be copied to what patch, and extend the
// storage of the patch to accomodate what will be copied to it.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performExtend(const ArrayPatch& a, PatchID_t lid) const
{
  // Get global patch ID from local patch ID

  PatchID_t gid = 
    particles_m.attributeLayout().nodeListLocal()[lid]->globalID();

  // Get the current number of particles in this patch, and the
  // total number of patches.

  Size_t size = a.domain().size();
  int patchesLocal = layout_m.patchesLocal();

  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(extend)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local patch ID = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << " out of " << patchesLocal)
  POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl << "The working size is ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << size << std::endl;)

  // Find the number of particles that must be copied TO this patch from
  // other local patches, and the number of patches that that comes from.

  Size_t totmove = 0;
  int p, frompatch = 0;
  for (p = 0; p < patchesLocal; ++p)
  {
    Size_t extra = layout_m.patchInfo(p).amount()(gid);
    if (extra > 0)
    {
      totmove += extra;
      frompatch += 1;
    }
  }
  layout_m.patchInfo(lid).setCopyPatches(frompatch);

  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "There are " << totmove)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " particles to copy to this patch, from ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << frompatch << " other local patches.")
  POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

  // If all the particles go on this same patch, we can just return.

  if (totmove == 0)
  {
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Nothing to copy to this patch.")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)
  }
  else
  {
    // Create extra space in this patch to accomodate new elements.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Extending storage on local patch ")
    POOMA_PATCHSWAPLAYOUT_DBG(  << lid << " to have storage for " << totmove)
    POOMA_PATCHSWAPLAYOUT_DBG(  << " extra particles" << std::endl;)

    particles_m.attributeLayout().create(totmove, lid);
  }

  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finished with extend on local patch ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << std::endl;)
}


//-----------------------------------------------------------------------------
// The copy method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will copy lists
// of elements from other patches to this one, using the multiple-list
// dynamic copy operation.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performCopy(const ArrayPatch& a, PatchID_t lid) const
{
  // Get global patch ID from local patch ID

  PatchID_t gid =
    particles_m.attributeLayout().nodeListLocal()[lid]->globalID();

  // Get the number of patches that contain data to copy to here.

  int frompatches = layout_m.patchInfo(lid).copyPatches();
  int patchesLocal = layout_m.patchesLocal();

  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(copy)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local patch ID = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << ", copying from " << frompatches)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " patches." << std::endl;)

  // If there are any other patches, do the copy

  if (frompatches == 0)
  {
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Nothing to copy to this patch.")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)
  }
  else
  {
    // Form an array of IndirectionList's, one for each of the
    // other patches, that we'll use to copy from others to ourself.

    IndirectionList< IndirectionList<int> > clists(frompatches);
    IndirectionList<int> cpids(frompatches);
    int p, np;
    for (p = 0, np = 0; p < patchesLocal; ++p)
    {
      int extra = layout_m.patchInfo(p).amount()(gid);
      if (extra > 0) 
      {
	clists(np) =
          layout_m.patchInfo(p).sendIndices(gid)(Interval<1>(extra));
	cpids(np) = p;
	++np;
      }
    }

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Formed copy lists = " << clists)
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl << "Source patches = ")
    POOMA_PATCHSWAPLAYOUT_DBG(  << cpids << std::endl;)

    // Ask the attribute layout to use this "list of lists" to copy from
    // other patches to our own.  We also provide a list of patch ID's.
    // The final "false" argument says to overwrite existing storage
    // at the end of this patch ... that's fine, since we previously
    // created the storage.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking multiple-list copy ...")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

    particles_m.attributeLayout().copy(clists, cpids, lid, false);
  }

  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finished with copy on local patch ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << std::endl;)
}


//-----------------------------------------------------------------------------
// The destroy method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will destroy all
// the particles that were moved off to other patches.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performDestroy(const ArrayPatch& a, PatchID_t lid) const
{
  // Calculate how many particles are to be destroyed here.

  Size_t totdestroy = layout_m.patchInfo(lid).destroySize();

  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(destroy)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local patch ID = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << ", destroying " << totdestroy)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " particles." << std::endl;)

  // If there is something to destroy, do so now.

  if (totdestroy == 0)
  {
    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Nothing to destroy on this patch.")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)
  }
  else
  {
    // Get the array used to store the destoy list.

    MoveArray_t& destroyList = layout_m.patchInfo(lid).destroyIndices();

    // Invoke the per-patch destroy method.

    IndirectionList<int> destlist(destroyList(Interval<1>(totdestroy)));

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Destroying the following particles ")
    POOMA_PATCHSWAPLAYOUT_DBG(  << "from local patch " << lid << ": ")
    POOMA_PATCHSWAPLAYOUT_DBG(  << destlist << std::endl;)

    particles_m.destroy(destlist, lid, false);
  }

  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Finished with destroy on local patch ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << std::endl;)
}

//-----------------------------------------------------------------------------
// The send method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will process
// the sendIndices lists generated in the scan method and send
// particles to remote patches as needed.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performSend(const ArrayPatch& a, PatchID_t lid) const
{
  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(send)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local patch ID = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << "." << std::endl;)

#if POOMA_CHEETAH

  // get my global patch ID and total (these will be used for the message tag)

  int gid = particles_m.attributeLayout().nodeListLocal()[lid]->globalID();
  int gsize = particles_m.attributeLayout().sizeGlobal();

  // get number of remote patches to send to

  int remotePatches = particles_m.attributeLayout().sizeRemote();

  // loop over remote patches to send messages

  for (int i = 0; i < remotePatches; ++i)
  {
    // get global ID and context of next remote patch

    int togid = particles_m.attributeLayout().nodeListRemote()[i]->globalID();
    int toContext =
      particles_m.attributeLayout().nodeListRemote()[i]->context();

    // use my global ID and that of the recipient to create unique tag

    int tag = gid * gsize + togid;

    // create an indirection send list for this transaction

    IndirectionList<int> slist;
    int toSend = layout_m.patchInfo(lid).amount()(togid);
    if (toSend > 0)
      slist = layout_m.patchInfo(lid).sendIndices(togid)(Interval<1>(toSend));

    // create particle swap buffer object representing the list of
    // particles to be sent and the local patch ID

    PSwapPack<P> buf(lid, particles_m, slist);

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Send " << toSend << " particles")
    POOMA_PATCHSWAPLAYOUT_DBG(       << " from local patch " << lid)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " to global patch " << togid)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " on context " << toContext)
    POOMA_PATCHSWAPLAYOUT_DBG(       << " with tag " << tag << std::endl;)

    // Send this set of particle data

    Pooma::particleSwapHandler()->send(toContext, tag, buf);
  }

#elif POOMA_MPI
  PInsist(false, "Cross-context particles not supported for MPI");
#endif // POOMA_CHEETAH
}

//-----------------------------------------------------------------------------
// The receive method, invoked by the patch functor evaluator
// for each patch in an array.  The first argument is a view
// on one of the array's patches, the second argument is the
// local ID number of that patch.  The "Patch" is part of the
// position array for the particles.  This routine will receive
// particles from remote patches.
//-----------------------------------------------------------------------------

template <class P>
template <class ArrayPatch>
void 
PatchSwapFunctor<P>::performReceive(const ArrayPatch& a, PatchID_t lid) const
{
  // Create a debug output stream, that will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapFunctor(receive)",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "In apply, with patch of domain = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << a.domain() << " and local patch ID = ")
  POOMA_PATCHSWAPLAYOUT_DBG(  << lid << "." << std::endl;)

#if POOMA_CHEETAH

  // get number of remote patches to receive from

  int remotePatches = particles_m.attributeLayout().sizeRemote();

  // get my global patch ID and total (these will be used for the message tag)

  int gid = particles_m.attributeLayout().nodeListLocal()[lid]->globalID();
  int gsize = particles_m.attributeLayout().sizeGlobal();

  // reset message received counter to zero

  layout_m.patchInfo(lid).msgReceived() = 0;

  // loop over remote patches to receive messages

  for (int i = 0; i < remotePatches; ++i)
  {
    // get context and global ID of next remote patch

    int fromContext =
      particles_m.attributeLayout().nodeListRemote()[i]->context();
    int fromgid =
      particles_m.attributeLayout().nodeListRemote()[i]->globalID();

    // generate the unique tag for expected message

    int tag = fromgid * gsize + gid;

    // create particle swap buffer object representing the list of
    // particles to be received and the local patch ID

    PSwapPack<P> buf(lid, particles_m);

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Receive particles message")
    POOMA_PATCHSWAPLAYOUT_DBG(  << " on local patch " << lid)
    POOMA_PATCHSWAPLAYOUT_DBG(  << " from global patch " << fromgid)
    POOMA_PATCHSWAPLAYOUT_DBG(  << " on context " << fromContext)
    POOMA_PATCHSWAPLAYOUT_DBG(  << " with tag " << tag << std::endl;)

    // Receive this set of particle data

    typedef void (*pfunc_t)(PSwapPack<P>*, PSwapPack<P>&);
    pfunc_t pF = &pSwapUnpackFunc<P>;
    Pooma::particleSwapHandler()->request(fromContext, tag, pF, &buf);
  }

  while (layout_m.patchInfo(lid).msgReceived() < remotePatches)
    Pooma::poll();

#elif POOMA_MPI
  PInsist(false, "Cross-context particles not supported for MPI");
#endif // POOMA_CHEETAH
}


//============================================================
// Non-inline methods for PatchSwapLayout
//============================================================

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

template <class L>
template <class P, class A>
void 
PatchSwapLayout<L>::performSync(P& particles, const A& pos, bool dosync)
{
  // Get the number of local and global patches in the particles's attributes

  int patchesLocal = particles.attributeLayout().sizeLocal();
  int patchesGlobal = particles.attributeLayout().sizeGlobal();

  // Create a debug output stream, which will be used if the
  // POOMA_PATCHSWAPLAYOUT_DBG macro is enabled at the top of this file.

  POOMA_PATCHSWAPLAYOUT_DBG(Inform dbgmsg("PatchSwapLayout::swap",-1);)
  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Starting swap for " << patchesLocal)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " local patches and " << patchesGlobal)
  POOMA_PATCHSWAPLAYOUT_DBG(  << " global patches." << std::endl;)

  // Make sure the particle's attribute layout has the same number
  // of patches as the provided position attribute.

  PAssert(patchesLocal > 0);
  PAssert(patchesLocal == pos.layout().sizeLocal());
  PAssert(patchesGlobal == pos.layout().sizeGlobal());

  // If we only have one patch, we can just apply BC's and quit.

  if (patchesGlobal > 1)
  {
    // Create a swap functor type

    typedef PatchSwapFunctor<P>                           SFun_t;
    typedef PatchFunction< SFun_t, PatchParticle1<true> > Swap_t;

    if (dosync) 
    {
      // Before scanning the particles to determine where 
      // they belong, we need to take care of other sync tasks,
      // including applying particle BC's and processing any
      // deferred destroy requests.

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap sync functor ...")
      POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

      Swap_t swapSync( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::syncScan) );
      swapSync.block(pos);
    }

    // Get and store the current sizes of all patches.
    // Note: We need to make sure that the domains of all the 
    // global patches of particle data have been synchronized 
    // at this point by the attribute layout.

    findCurrentSizes(particles);

    // First swap step is to scan the positions and find the destination
    // patches and send lists.  These get stored in the patch info structures.
    // The block() syntax is used on the swap functor, so the main
    // thread will wait until all iterates spawned by the functor are done.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap scan functor ...")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

    Swap_t swapScan( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::swapScan) );
    swapScan.block(pos);

    // In case other threads are working on some of the attributes,
    // we must wait before going on to actually do the swap.  It would
    // be nice if we could remove this some time in the future.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Let things catch up before swap ...")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

    Pooma::blockAndEvaluate();

    if (Pooma::contexts() > 1) 
    {
      // At this point, we can fire off our sends of particle data
      // from local patches to remote patches, if any.

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap send functor ...")
      POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

      Swap_t swapSend( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::swapSend) );
      swapSend.block(pos);
    }

    // Next, we extend the existing storage, using a second functor.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap extend functor ...")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

    Swap_t swapExt( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::swapExtend) );
    swapExt.block(pos);

    // Now we copy particle data in from our neighboring local patches,
    // using the extended storage so that we don't have to remalloc.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap copy functor ...")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

    Swap_t swapCopy( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::swapCopy) );
    swapCopy.block(pos);

    if (Pooma::contexts() > 1)
    {
      // Local-to-local swapping is complete, so now receive particle 
      // data from remote patches and store.

      POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap receive functor ...")
      POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

      Swap_t swapReceive( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::swapReceive) );
      swapReceive.block(pos);
    }

    // Finally, destroy all outgoing particles

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Invoking swap destroy functor ...")
    POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

    Swap_t swapDest( SFun_t(static_cast<Layout_t&>(*this), particles, SFun_t::swapDestroy) );
    swapDest.block(pos);
  }
  else if (patchesGlobal == 1 && dosync)
  {
    // Just apply boundary conditions and do deferred destroys, if requested.

    POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Applying BC's and doing destroys ")
    POOMA_PATCHSWAPLAYOUT_DBG(  << "for one patch." << std::endl;)

    particles.applyBoundaryConditions(0);
    particles.performDestroy(0, false);
  }

  // At the end, we sync up with all contexts based on how many
  // particles have stayed here or moved elsewhere.

  POOMA_PATCHSWAPLAYOUT_DBG(dbgmsg << "Syncing up attribute layout ...")
  POOMA_PATCHSWAPLAYOUT_DBG(  << std::endl;)

  particles.renumber();
}


//-----------------------------------------------------------------------------
// Calculate the current size of each patch in the particle's attribute
// layout, and store them.
//-----------------------------------------------------------------------------

template <class L>
template <class P>
void 
PatchSwapLayout<L>::findCurrentSizes(const P& particles)
{
  // Find out how many local patches there are in the attribute layout

  int i;
  int patchesLocal = particles.attributeLayout().sizeLocal();

  // Create local patch info arrays if this has not yet been done.

  if (patchInfo_m == 0)
  {
    patchInfo_m = new PatchSwapInfo[patchesLocal];
    int patchesRemote = particles.attributeLayout().sizeRemote();
    for (i = 0; i < patchesLocal; ++i)
      patchInfo(i).initialize(patchesLocal,patchesRemote);
  }

  // Loop through the local patches and compute the total size.

  int size, mySize = 0;
  for (i = 0; i < patchesLocal; ++i) 
  {
    size = particles.attributeLayout().patchDomain(i).size();
    patchInfo(i).setSize(size);
    mySize += size;
  }


  // Now loop over all contexts and store the total size of each
  // one, using a RemoteProxy to broadcast the values.

  int contexts = Pooma::contexts();
  for (i = 0; i < contexts; ++i) 
  {
    size = mySize;
    RemoteProxy<int> proxySize(size,i);
    contextSizes_m(i) = proxySize;
  }
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PatchSwapLayout.cpp,v $   $Author: richard $
// $Revision: 1.18 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
