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
// MultiPatch-Engine non-inline template definitions.
//-----------------------------------------------------------------------------

#include "Engine/MultiPatchEngine.h"
#include "Evaluator/EngineTraits.h"
#include "Engine/CompressedFraction.h"
#include "Array/Array.h"
#include "Tulip/ReduceOverContexts.h"
#include "Tulip/SendReceive.h"
#include "Threads/PoomaCSem.h"
#include "Domain/IteratorPairDomain.h"

///////////////////////////////////////////////////////////////////////////////
//
// MultiPatch-Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,MultiPatch>::Engine()
//
// Default constructor ... you should use operator= to initialize this
// engine after using this constructor.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
Engine()
  : pDirty_m(0)
{
  // This object must be initialized later via operator=.  If it is not,
  // this is not a useful object.  Do not attach to the layout, since we
  // don't have one yet.
}


//-----------------------------------------------------------------------------
//
// Engine<Dim,T,MultiPatch>::Engine(const Layout_t &)
//
// Initialize with a layout - we take the total domain from this,
// and register as a user of the layout.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
Engine(const Layout_t &layout)
  : layout_m(layout),
    data_m(layout.sizeGlobal()),
    pDirty_m(new int)
{
  typedef typename Layout_t::Value_t Node_t;

  setDirty();

  // check for correct match of PatchTag and the mapper used to make the
  // layout.
  // THIS IS A HACK! we test on the context of the first patch, and if it
  // is -1, we have a Layout made with the LocalMapper.

#if POOMA_MESSAGING

  if( layout_m.nodeListGlobal().size() > 0)
  {    
    int fpID = layout_m.nodeListGlobal()[0]->context();

    bool compatible =
      DistributionTraits<PatchTag>::remote ? (fpID != -1) : (fpID == -1);
      
    PInsist(compatible,
	    "PatchTag is incompatible with the ContextMapper");
    }

#endif  
  
  // Initialize the patches. We iterate through the Nodes in the layout, which
  // are supposed to match up with the patches as they are stored in our
  // vector. For each patch, we pass Node information to enable initialization.
  
  int sz = data().size();

  typedef Pooma::CountingSemaphore CountingSemaphore_t;
  typedef typename Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
    template PatchAllocator<Node_t,CountingSemaphore_t> PatchAllocator_t;
  CountingSemaphore_t csem;
  
  csem.height(sz);
  
  typename Layout_t::const_iterator p = layout_m.beginGlobal();

  for (int i = 0; i < sz; ++i, ++p)
    {
      // Have threads perform data()[i] = PatchEngine_t(*p);
      
      PatchAllocator_t *spot = new PatchAllocator_t(data()[i], *p, csem);
      Pooma::addRunnable(spot);  // See spot run!
    }
  
  // Wait for all of the runnables to complete...
  
  csem.wait();
      
  // Attach ourself to the layout so we can receive messages.
  
  layout_m.attach(*this);
}


//-----------------------------------------------------------------------------
//
// Engine<Dim,T,MultiPatch>::Engine(const Engine_t &model)
//
// Copy constructor.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
Engine(const Engine_t &modelEngine)
  : layout_m(modelEngine.layout_m),
    data_m(modelEngine.data_m),
    pDirty_m(modelEngine.pDirty_m)
{
  // Attach ourself to the layout so we can receive messages.

  layout_m.attach(*this);  
}


//-----------------------------------------------------------------------------
//
// Engine<Dim,T,MultiPatch>::~Engine()
//
// Destructor ... detach from layout.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
~Engine()
{
  if (initialized())
    {
      // Detach ourselves from the layout.
  
      layout_m.detach(*this);
  
      // If this is the last reference to our data, delete the dirty flag.
  
      if (!data().isShared())
	delete pDirty_m;
    }
}


//-----------------------------------------------------------------------------
//
// Engine_t &Engine<Dim,T,MultiPatch>::operator=(const Engine_t &)
//
// Assignment operator.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> > &
Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
operator=(const Engine_t &model)
{
  // Check for self-assignment.
  
  if (&model == this)
    return *this;
 
  // Make sure the RHS is initialized.

  if (!model.initialized())
    return *this;

  // If we have been previously initialized, clean up ...

  if (initialized())
    {
      // If this is the last copy of our data, delete the old dirty flag.
  
      if (!data().isShared())
	delete pDirty_m;
    
      // Detach ourselves from the old layout.

      layout_m.detach(*this);
    }

  // Assign our data.

  data_m = model.data();
  pDirty_m = model.pDirty_m;

  // Copy and attach ourself to the layout so we can receive messages.

  layout_m = model.layout_m;
  layout_m.attach(*this);
  
  return *this;  
}


//-----------------------------------------------------------------------------
//
// Gets a private copy of this engine's data. Just loop through all of our
// engines and ask them to make their own copies.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> > &
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::
makeOwnCopy()
{
  if (data_m.isValid() && data_m.isShared()) {
    data_m.makeOwnCopy();
    pDirty_m = new int(*pDirty_m);
  }

  return *this;
}


//-----------------------------------------------------------------------------
//
// Fill the internal guard cells if needed, and clear the dirty flag.
// Current implementation is LOCAL ONLY!!!
//
//-----------------------------------------------------------------------------

/// Guard layer assign between non-remote engines, just use the
/// ET mechanisms

template <int Dim, class T, class Tag>
static inline
void simpleAssign(const Array<Dim, T, Tag>& lhs,
		  const Array<Dim, T, Tag>& rhs,
		  const Interval<Dim>& domain)
{
  lhs(domain) = rhs(domain);
}

/// Guard layer assign between remote engines, use Send/Receive directly
/// to avoid one extra copy of the data.

template <int Dim, class T, class Tag>
static inline
void simpleAssign(const Array<Dim, T, Remote<Tag> >& lhs,
		  const Array<Dim, T, Remote<Tag> >& rhs,
		  const Interval<Dim>& domain)
{
  if (lhs.engine().owningContext() == rhs.engine().owningContext())
    lhs(domain) = rhs(domain);
  else {
    typedef typename NewEngine<Engine<Dim, T, Tag>, Interval<Dim> >::Type_t ViewEngine_t;
    if (lhs.engine().engineIsLocal())
      Receive<ViewEngine_t>::receive(ViewEngine_t(lhs.engine().localEngine(), domain),
				     rhs.engine().owningContext());
    else if (rhs.engine().engineIsLocal())
      SendReceive::send(ViewEngine_t(rhs.engine().localEngine(), domain),
			lhs.engine().owningContext());
  }
}

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::
fillGuardsHandler(const GuardLayers<Dim>& g, const WrappedInt<true> &) const
{
  if (!isDirty()) return;

  int updated = 0;
  typename Layout_t::FillIterator_t p = layout_m.beginFillList();

  while (p != layout_m.endFillList())
    {
      int src  = p->ownedID_m;
      int dest = p->guardID_m;
      
      // Skip face, if not dirty.

      if (isDirty(p->face_m)) {

        // Check, if the p->domain_m is a guard which matches the
        // needed guard g.

	int d = p->face_m/2;
	int guardSizeNeeded = p->face_m & 1 ? g.upper(d) : g.lower(d);
        if (!(p->face_m != -1
	      && guardSizeNeeded == 0)) {

          // Create patch arrays that see the entire patch:
                  
          Array<Dim, T, PatchTag> lhs(data()[dest]), rhs(data()[src]);
      
          // Now do assignment from the subdomains.
#if POOMA_MPI
          simpleAssign(lhs, rhs, p->domain_m);
#else
          lhs(p->domain_m) = rhs(p->domain_m);
#endif

	  // Mark up-to-date.
	  updated |= 1<<p->face_m;

	}

      }

      ++p;
    }

  *pDirty_m &= ~updated;
}


//-----------------------------------------------------------------------------
//
// Set the internal guard cells to a particular value (default zero).
// Current implementation is LOCAL ONLY!!!
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::
setGuards(const T &val) const
{
  typename Layout_t::FillIterator_t p = layout_m.beginFillList();
   
  while (p != layout_m.endFillList())
    {
      int dest = p->guardID_m;
      
      // Create patch arrays that see the entire patch:
                  
      Array<Dim, T, PatchTag> lhs(data()[dest]);
      
      // Now do assignment from the subdomains.

      lhs(p->domain_m) = val;
      
      ++p;
    }

  setDirty();
}


//-----------------------------------------------------------------------------
//
// Accumulate from the internal guards into the owned domain.
// Current implementation is LOCAL ONLY!!!
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::
accumulateFromGuards() const
{
  typename Layout_t::FillIterator_t p = layout_m.beginFillList();
     
  while (p != layout_m.endFillList())
    {
      // This time we're going from the guards to the owned.
      
      int dest = p->ownedID_m;
      int src  = p->guardID_m;
      
      // Create patch arrays that see the entire patch:
                  
      Array<Dim, T, PatchTag> lhs(data()[dest]), rhs(data()[src]);
      
      // Now accumulate values from the guards.

      lhs(p->domain_m) += rhs(p->domain_m);
      
      ++p;
    }

  setDirty();
}


//-----------------------------------------------------------------------------
//
// void Engine<Dim,T,MultiPatch>::
// dynamicHandler(Observable_t, ObserverEvent, WrappedInt<true>)
//
// Handler for dynamic events for patch engines that have dynamic capabilities.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::
dynamicHandler(Observable_t &, 
               const ObserverEvent &event, 
               const WrappedInt<true> &)
{ 
  // From the event code, figure out what we should do.
  // (There has to be a better way to do this!!! [JAC])

  switch (event.event())
  {
    case DynamicEvents::create:
      {
        // Create new elements at the end of our block of data
        typedef const CreateEvent &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);
        performCreate(e.amount(), e.patch(), e.ID());
      }
      break;
    case DynamicEvents::destroyInterval:
      { 
        // Delete elements in our patch of data using an Interval  
      
        typedef const DestroyEvent< Interval<1> > &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);
	      
        switch (e.method())
        {
          case DynamicEvents::backfill:
            performDestroy(e.domain(), e.patch(), BackFill(), e.ID());
	    break;
	  case DynamicEvents::shiftup:
	    performDestroy(e.domain(), e.patch(), ShiftUp(), e.ID());
	    break;
	  default:
	    PInsist(false,
                    "Unsupported delete method MultiPatchEngine::destroy");
        }
      }
      break;
    case DynamicEvents::destroyRange:
      {
        // Delete elements in our patch of data using a Range
        
        typedef const DestroyEvent< Range<1> > &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);

        switch (e.method())
        {
          case DynamicEvents::backfill:
            performDestroy(e.domain(), e.patch(), BackFill(), e.ID());
	    break;
	  case DynamicEvents::shiftup:
	    performDestroy(e.domain(), e.patch(), ShiftUp(), e.ID());
	    break;
	  default:
	    PInsist(false,
                    "Unsupported delete method MultiPatchEngine::destroy");
        }
      }
      break;
    case DynamicEvents::destroyList:
      {
        // Delete elements in our patch of data using an IndirectionList
    
        typedef const DestroyEvent< IndirectionList<int> > &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);
    
        switch (e.method())
        {
          case DynamicEvents::backfill:
            performDestroy(e.domain(), e.patch(), BackFill(), e.ID());
	    break;
	  case DynamicEvents::shiftup:
	    performDestroy(e.domain(), e.patch(), ShiftUp(), e.ID());
	    break;
	  default:
	    PInsist(false,
                    "Unsupported delete method MultiPatchEngine::destroy");
        }
      }
      break;
    case DynamicEvents::destroyIterList:
      {
        // Delete elements in our patch of data using an IndirectionList

        using Pooma::IteratorPairDomain;
        typedef 
          const DestroyEvent< IteratorPairDomain<const int*> > &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);
    
        switch (e.method())
        {
          case DynamicEvents::backfill:
            performDestroy(e.domain(), e.patch(), BackFill(), e.ID());
	    break;
	  case DynamicEvents::shiftup:
	    performDestroy(e.domain(), e.patch(), ShiftUp(), e.ID());
	    break;
	  default:
	    PInsist(false,
                    "Unsupported delete method MultiPatchEngine::destroy");
        }
      }
      break;
    case DynamicEvents::copyInterval:
      {
        // Copy elements in a specific Interval
        typedef const CopyEvent< Interval<1> > &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);
        performCopy(e.domain(), e.fromPatch(), e.toPatch(), e.ID());
      }
      break;
    case DynamicEvents::copyRange:
      {
        // Copy elements in a specific Range
        typedef const CopyEvent< Range<1> > &EventRef_t;
	EventRef_t e = dynamic_cast<EventRef_t>(event);
	performCopy(e.domain(), e.fromPatch(), e.toPatch(), e.ID());
      }
      break;
    case DynamicEvents::copyList:
      {
        // Copy a list of elements
        typedef const CopyEvent< IndirectionList<int> > &EventRef_t;
	EventRef_t e = dynamic_cast<EventRef_t>(event);
	performCopy(e.domain(), e.fromPatch(), e.toPatch(), e.ID());
      }
      break;
    case DynamicEvents::copyPatchList:
      {
        // Copy from a list of lists.
        typedef const CopyPatchEvent &EventRef_t;
        EventRef_t e = dynamic_cast<EventRef_t>(event);
        performPatchCopy(e.domainLists(), e.fromPatch(), e.toPatch(), 
                         e.create(), e.ID());
      }
      break;
    case DynamicEvents::sync:
      {
        // loop across all patch engines, and change their (?)
        // domain layout domains. 

        for (int i=0; i < layout().sizeGlobal(); ++i)
	  data()[i].sync(layout().nodeListGlobal()[i]->domain());
      }
      break;
    default:
      PInsist(0,"Invalid dynamic op???");
  }
}

//-----------------------------------------------------------------------------
//
// void Engine<Dim,T,MultiPatch>::notify(Observable_t, ObserverEvent)
//
// Be notified of various events from the layout, including
// when the layout is deleted or when dynamic operations occur.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
notify(Observable_t &observed, const ObserverEvent &event)
{
  // Make sure this is an event for us.

  PAssert(observed.ID() == layout().ID());

  // If the event is for partitioning, do that

  if(event.event() == Layout_t::repartitionEvent)
    {
      typedef typename Layout_t::Value_t Node_t;
  
      // Reinitialize the patches. We need to make a new ref-counted
      // pointer since the number of nodes could have changed.
      // We then replace our data pointer with this new one.
      
      int sz = layout().sizeGlobal();
      
      PatchContainer_t newData(sz);
      
      typedef Pooma::CountingSemaphore CountingSemaphore_t;
      typedef typename Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
        template PatchAllocator<Node_t,CountingSemaphore_t> PatchAllocator_t;
      CountingSemaphore_t csem;
  
      csem.height(sz);
  
      typename Layout_t::const_iterator p = layout().beginGlobal();
      
      for (int i = 0; i < sz; ++i, ++p)
	{
	  // Have threads perform newData[i] = PatchEngine_t(*p);
	  
          PatchAllocator_t *spot = new PatchAllocator_t(newData[i], *p, csem);
          Pooma::addRunnable(spot);  // See spot run!
        }
        
      csem.wait();
      
      data_m = newData;
    }
  else
    {
      // The event is either dynamic, or unknown.
      // Dynamic events are deferred to dynamicHandler:

      if (DynamicEvents::isDynamic(event.event()))
        {
          dynamicHandler(observed, event, 
            WrappedInt<PatchEngine_t::dynamic>());
	} 
      else
	{
	  // an event we don't care about, do nothing
	}
    }
}


//-----------------------------------------------------------------------------
//
// void Engine<Dim,T,MultiPatch>::performCreate(CreateSize_t n, PatchID_t p)
//
// Carry out a request to perform a create operation in a particular
// patch.  The layout is responsible for figuring out what patch to do
// this in, so the patch number must be a valid index into our local
// patch list.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
performCreate(CreateSize_t num, PatchID_t localPatchID, DynamicID_t did)
{
  PAssert(num >= 0);
  PAssert(localPatchID >= 0 && localPatchID < layout().sizeLocal());

  // Check if this has been performed before.  If so, skip it.
  // Use the specialized routine "checkDynamicID(obj, ID)" to
  // check the ID, which should return true if the operation
  // should proceed (and also set the ID of the obj).

  int globalID = layout().nodeListLocal()[localPatchID]->globalID();

  if (! checkDynamicID(data()[globalID], did))
    return;

  // Ask the individual patch to do the create, since it has not yet.

  data()[globalID].create(num);
}

//-----------------------------------------------------------------------------
//
// void Engine<Dim,T,MultiPatch>::
// performDestroy(Dom &killlist, PatchID_t p, DelMethod method, DynamicID_t id)
//
// Carry out the work to perform a destroy operation on a particular
// patch.  The layout is responsible for figuring out what patch to do
// this in, so the patch number must be a valid index into our local
// patch list.  Also, the domain must be a "relative" domain, with zero-
// based values.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
template <class Dom, class DeleteMethod>
void Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
performDestroy(const Dom &killlist, PatchID_t localPatchID,
	       const DeleteMethod &method, DynamicID_t did)
{
  // Only the patch specific performDestroy is implemented, as the layout will
  // take care of breaking down a cross patch destroy call into a set
  // of patch specific destroy calls.

  PAssert(localPatchID >= 0 && localPatchID < layout().sizeLocal());

  // Check if this has been performed before.  If so, skip it.
  // Use the specialized routine "checkDynamicID(obj, ID)" to
  // check the ID, which should return true if the operation
  // should proceed (and also set the ID of the obj).

  int globalID = layout().nodeListLocal()[localPatchID]->globalID();

  if (! checkDynamicID(data()[globalID], did))
    return;

  // Ask the individual patch to do the destroy, since it has not yet.
  // Set the offsetFlag to true here, since kill list is zero-based

  data()[globalID].destroy(killlist, method, true);
}


//-----------------------------------------------------------------------------
//
// void Engine<Dim,T,MultiPatch>::performCopy(Dom &copylist, PatchID_t fpatch,
//                                            PatchID_t tpatch)
//
// Carry out the work to perform a copy of values from one patch
// to another (or perhaps to the same patch).
// The layout is responsible for figuring out what patch to do
// this in, so the patch number must be a valid index into our local
// patch list.  Also, the domain must be a "relative" domain, with zero-
// based values.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
template <class Dom>
void Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
performCopy(const Dom &copylist, PatchID_t frompatch, PatchID_t topatch,
	    DynamicID_t did)
{
  PAssert(Dim == 1);
  PAssert(frompatch >= 0 && frompatch < layout().sizeLocal());
  PAssert(topatch >= 0 && topatch < layout().sizeLocal());

  // Check if this has been performed before.  If so, skip it.
  // Use the specialized routine "checkDynamicID(obj, ID)" to
  // check the ID, which should return true if the operation
  // should proceed (and also set the ID of the obj).

  int from_gID =  layout().nodeListLocal()[frompatch]->globalID();
  int to_gID =  layout().nodeListLocal()[topatch]->globalID();


  bool chk1 = checkDynamicID(data()[from_gID], did);
  bool chk2 = chk1;
  if (frompatch != topatch)
    chk2 = checkDynamicID(data()[to_gID], did);
  PAssert(chk1 == chk2);
  if (!chk1 || !chk2)
    return;

  // We have to copy elements from one patch to another here (instead
  // of calling a routine in the single-patch engine) because the
  // data might span multiple patches.  The algorithm is the same
  // regardless of whether frompatch is the same as topatch.

  PAssert(copylist[0].max() < data()[from_gID].domain().size());

  // Create storage for copied elements, and note where we start
  // putting copied values in.

  int offs = data()[from_gID].domain()[0].first();
  int i = data()[to_gID].domain()[0].last() + 1;
  int num = copylist.size();
  data()[to_gID].create(num);

  // Copy over values from one patch to another.

  for (int n = 0; n < num; ++n, ++i)
  {
    localPatchEngine(data()[to_gID])(i) =
      localPatchEngine(data()[from_gID])(copylist[0](n) + offs);
  }
}


//-----------------------------------------------------------------------------
//
// void Engine<Dim,T,MultiPatch>::performPatchCopy(...)
//
// Do the actual work of a multiple-list copy.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
void Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> >::
performPatchCopy(const IndirectionList< IndirectionList<int> > &domlists,
		 const IndirectionList< int > &fromlist,
		 PatchID_t topatch,
		 bool docreate,
		 DynamicID_t did)
{
  PAssert(Dim == 1);
  PAssert(topatch >= 0 && topatch < layout().sizeLocal());

  // Check if this has been performed before.  If so, skip it.
  // Use the specialized routine "checkDynamicID(obj, ID)" to
  // check the ID, which should return true if the operation
  // should proceed (and also set the ID of the obj).

  int to_gID =  layout().nodeListLocal()[topatch]->globalID();

  if (! checkDynamicID(data()[to_gID], did))
    return;

  // We have to copy elements from one patch to another here (instead
  // of calling a routine in the single-patch engine) because the
  // data might span multiple patches.  The algorithm is the same
  // regardless of whether frompatch is the same as topatch.
  // Go through all the lists, and copy data to our end.  First make
  // sure we're not going to overflow anything.

  int i, p, fill, created = 0;
  int np = domlists.size();
  PAssert(fromlist.size() == np);
  for (p = 0; p < np; ++p)
    {
      int frompatch = fromlist(p);
      int from_gID = layout().nodeListLocal()[frompatch]->globalID();
      
      PAssert(frompatch >= 0 && frompatch < layout().sizeLocal());
      PAssert(domlists(p).first() >= 0);
      PAssert(domlists(p).last() <= data()[from_gID].domain()[0].last());
      created += domlists(p).size();
    }

  // First, create space at the end if necessary; otherwise, overwrite
  // storage at the end.

  if (docreate)
    {
      fill = data()[to_gID].domain()[0].last() + 1;
      data()[to_gID].create(created);
    }
  else
    {
      PAssert(created <= data()[to_gID].domain()[0].length());
      fill = data()[to_gID].domain()[0].last() + 1 - created;
    }

  // Now, copy elements from the given domain into the new storage, from
  // each of the indirection lists.

  for (p = 0; p < np; ++p)
    {
      int sz = domlists(p).size();
      int frompatch = fromlist(p);
      int from_gID = layout().nodeListLocal()[frompatch]->globalID();
      int offs = data()[from_gID].domain()[0].first();
      for (i = 0; i < sz; ++i, ++fill)
      {
	localPatchEngine(data()[to_gID])(fill) =
	  localPatchEngine(data()[from_gID])(domlists(p)(i) + offs);
      }
    }
}

//-----------------------------------------------------------------------------
//
// long elementsCompressed()
//
// Compute the number of elements that are currently compressed. Compute
// with the local patches and then do a cross-context reduction.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class LTag, class PatchTag>
long elementsCompressed(const Engine<Dim, T, MultiPatch<LTag, PatchTag> > 
  &engine)
{
  int size = engine.layout().sizeLocal();

  bool distributed = true;
  if (size > 0 && engine.layout().beginLocal()->context() == -1)
    distributed = false;
    
  long num = 0L;
  for (int i = 0 ; i < size ; ++i )
    num += elementsCompressed(engine.localPatch(i));
  
  if (distributed)
    {
      ReduceOverContexts<long, OpAddAssign> total(num);
      total.broadcast(num);
    }

  return num;
}

//-----------------------------------------------------------------------------
//
// bool compressed()
//
// Compute the number of elements that are currently compressed. Compute
// with the local patches and then do a cross-context reduction.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class LTag, class PatchTag>
bool compressed(const Engine<Dim, T, MultiPatch<LTag, PatchTag> > &engine)
{
  int size = engine.layout().sizeLocal();

  bool distributed = true;
  if (size > 0 && engine.layout().beginLocal()->context() == -1)
    distributed = false;
    
  int com = 1;
  for (int i = 0 ; i < size ; ++i )
    com &= (compressed(engine.localPatch(i)) ? 1 : 0);
  
  if (distributed)
    {
      ReduceOverContexts<int, OpBitwiseAndAssign> total(com);
      total.broadcast(com);
    }

  return com ? true : false;
}

//-----------------------------------------------------------------------------
//
// long elementsCompressed()
//
// Compute the number of elements that are currently compressed. Compute
// with the local patches and then do a cross-context reduction. This is
// a little tricky since we must iterate over Nodes here since patch indices
// don't really mean anything for views.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class LTag, class PatchTag, int Dim2>
long elementsCompressed(const 
  Engine<Dim, T, MultiPatchView<LTag, PatchTag, Dim2> > &engine)
{
  typedef Engine<Dim, T, MultiPatchView<LTag, PatchTag, Dim2> > Engine_t;
  typedef typename Engine_t::Layout_t Layout_t;
  
  int size = engine.layout().sizeLocal();
  bool distributed = true;
  if (size > 0 && engine.layout().beginLocal()->context() == -1)
    distributed = false;

  typename Layout_t::const_iterator i = engine.layout().beginLocal();
  long num = 0L;
  while (i != engine.layout().endLocal())
    {
      num += elementsCompressed(engine.globalPatch(*i));
      ++i;
    }

  if (distributed)
    {  
      ReduceOverContexts<long, OpAddAssign> total(num);
      total.broadcast(num);
    }
    
  return num;
}

//-----------------------------------------------------------------------------
//
// void compress()
//
// (Try to) compress all the local patches.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class LTag, class PatchTag>
void compress(Engine<Dim, T, MultiPatch<LTag, PatchTag> > &engine)
{
  // Iterate through patches and try to compress them all.
  
  for (int i = 0; i < engine.layout().sizeLocal(); ++i)
    compress(engine.localPatch(i));
}

//-----------------------------------------------------------------------------
//
// void uncompress()
//
// Manually uncompress all the local patches.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class LTag, class PatchTag>
void uncompress(Engine<Dim, T, MultiPatch<LTag, PatchTag> > &engine)
{
  // Iterate through patches and try to uncompress them all.
  
  for (int i = 0; i < engine.layout().sizeLocal(); ++i)
    uncompress(engine.localPatch(i));
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MultiPatchEngine.cpp,v $   $Author: richard $
// $Revision: 1.56 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
