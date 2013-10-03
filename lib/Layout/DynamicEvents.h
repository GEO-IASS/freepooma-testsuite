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

#ifndef POOMA_LAYOUT_DYNAMIC_EVENTS_H
#define POOMA_LAYOUT_DYNAMIC_EVENTS_H

//-----------------------------------------------------------------------------
// Classes: 
//   DynamicEvents
//   ShiftUp
//   BackFill
//   CreateEvent
//   DestroyEvent
//   CopyEvent
//   CopyPatchEvent
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Layout
 * @brief
 * DynamicEvents defines some simple enumerations used as codes to indicate
 * a type of "dynamic event".
 *
 * Dynamic Events are issued to engines by
 * objects like layouts in order to tell them to dynamically change their
 * size and contents.  This file also defines event objects for create,
 * destroy, and copy events, and simple "tag" classes used to indicate types
 * of destroy methods.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Loc.h"
#include "Domain/IndirectionList.h"
#include "Domain/IteratorPairDomain.h"
#include "Utilities/ObserverEvent.h"

using Pooma::IteratorPairDomain;


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

/**
 * A simple class that just defines enums for items involving
 * dynamic operations on layouts and data.  It defines enums for
 *  - type of operation (create, DESTROY, COPY, sync)
 *  - type of domain (INTERVAL, RANGE, LIST)
 *  - type of delete method (backfill, shiftup)
 */

struct DynamicEvents
{
  // Enumeration of events codes for layout notify events.  All start
  // with 1000 (a basically random number).

  enum EventCode { 
         create = 1000,
	 extend,
	 destroyInterval, destroyRange, destroyList, destroyIterList,
	 copyInterval, copyRange, copyList,
	 copyPatchList,
	 sync,
	 unknownEvent };

  // If the above list is extended, this needs changes as well.
  
  static bool isDynamic(int eventcode)
  {
    return eventcode >= DynamicEvents::create && 
           eventcode <  DynamicEvents::unknownEvent;
  }
  
  // Enumeration with types of delete methods

  enum { backfill = 100, shiftup, unknownMethod };

  // Typedefs for the type of data used to specify a patch ID and a
  // create amount.

  typedef int PatchID_t;
  typedef int CreateSize_t;
};


/**
 * BackFill: A tag class used to indicate that delete operation should
 *           proceed by 'backfilling', that is, moving data up from the
 *           end of a list to fill in holes.  More efficient than ShiftUp.
 *
 * All delete-method tag classes contain one enum, 'code', indicating their
 * integer code from DynamicEvents.
 */

struct BackFill
{
  enum { code = DynamicEvents::backfill };
  BackFill(){};
};

/**
 * ShiftUp: A tag class used to indicate that delete operations should
 *          proceed by 'shifting up', that is, moving the entire list up
 *          as a whole.  Less efficient than BackFill, but preserves relative
 *          element ordering.
 *
 * All delete-method tag classes contain one enum, 'code', indicating their
 * integer code from DynamicEvents.
 */

struct ShiftUp
{
  enum { code = DynamicEvents::shiftup };
  ShiftUp(){};
};


/**
 * DynamicEventType is a simple partially-specialized class used to
 * determine the event code type based on the input domain type for
 * destroy and copy operations.  It has a default of "LIST" for these
 * types, and is specialized for the known domain types.
 */

template<class T>
struct DynamicEventType
{
  enum { destroyCode = DynamicEvents::destroyList };
  enum { copyCode    = DynamicEvents::copyList };
  enum { dimensions  = 1 };
  typedef IndirectionList<int> Domain_t;
};

template <>
struct DynamicEventType<IteratorPairDomain<const int*> >
{
  enum { destroyCode = DynamicEvents::destroyIterList };
  enum { copyCode    = DynamicEvents::copyList };
  enum { dimensions  = 1 };
  typedef IteratorPairDomain<const int*> Domain_t;
};

template <>
struct DynamicEventType<IteratorPairDomain<int*> >
{
  enum { destroyCode = DynamicEvents::destroyIterList };
  enum { copyCode    = DynamicEvents::copyList };
  enum { dimensions  = 1 };
  typedef IteratorPairDomain<const int*> Domain_t;
};

template<>
struct DynamicEventType< IndirectionList< IndirectionList<int> > >
{
  enum { destroyCode = DynamicEvents::unknownEvent };
  enum { copyCode    = DynamicEvents::copyPatchList };
  enum { dimensions  = 1 };
  typedef IndirectionList< IndirectionList<int> > Domain_t;
};

template<int Dim>
struct DynamicEventType< Interval<Dim> >
{
  enum { destroyCode = DynamicEvents::destroyInterval };
  enum { copyCode    = DynamicEvents::copyInterval };
  enum { dimensions  = Dim };
  typedef Interval<Dim> Domain_t;
};

template<int Dim>
struct DynamicEventType< Range<Dim> >
{
  enum { destroyCode = DynamicEvents::destroyRange };
  enum { copyCode    = DynamicEvents::copyRange };
  enum { dimensions  = Dim };
  typedef Range<Dim> Domain_t;
};

template<int Dim>
struct DynamicEventType< Loc<Dim> >
{
  enum { destroyCode = DynamicEvents::destroyInterval };
  enum { copyCode    = DynamicEvents::copyInterval };
  enum { dimensions  = Dim };
  typedef Interval<Dim> Domain_t;
};

template<>
struct DynamicEventType<int>
{
  enum { destroyCode = DynamicEvents::destroyInterval };
  enum { copyCode    = DynamicEvents::copyInterval };
  enum { dimensions  = 1 };
  typedef Interval<1> Domain_t;
};


/**
 * CreateEvent: A class derived from ObserverEvent that stores information
 * on how many elements to create, in what patch, for an engine.  It stores
 * the create amount, create patch, and provides "amount()" and "patch()"
 * methods to query this info.
 */

class CreateEvent : public ObserverEvent
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef DynamicEvents::PatchID_t             PatchID_t;
  typedef DynamicEvents::CreateSize_t          CreateSize_t;


  //============================================================
  // Constructor
  //============================================================

  // There is only one way to initialize this object, by providing
  // the amount to create and the patch to create in.  If this patch number
  // is < 0, the create should be done in the last local patch.

  CreateEvent(CreateSize_t num, PatchID_t p)
    : ObserverEvent(DynamicEvents::create),
      amount_m(num), patch_m(p)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  // Nothing to see here, move along

  virtual ~CreateEvent()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  inline CreateSize_t amount() const
    {
      return amount_m;
    }

  inline PatchID_t patch() const
    {
      return patch_m;
    }

private:
  // The number of elements to create

  CreateSize_t amount_m;

  // The local patch to create in

  PatchID_t patch_m;
};



/**
 * DestroyEvent: A class derived from ObserverEvent that stores information
 * on how what elements to destroy in an engine.  It stores info on the
 * domain to destroy (domain()), the destroy method code (method()),
 * and the patch to destroy in (patch()).
 *
 * This is templated on the domain type, Dom.
 */

template<class Dom>
class DestroyEvent : public ObserverEvent
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef DynamicEvents::PatchID_t                  PatchID_t;
  typedef DynamicEvents::CreateSize_t               CreateSize_t;
  typedef typename DynamicEventType<Dom>::Domain_t  Domain_t;


  //============================================================
  // Constructor
  //============================================================

  // There is only one way to initialize this object, by providing
  // the domain to destroy and the patch to destroy from, plus the
  // destroy method code.

  template<class D>
  DestroyEvent(const D &d, PatchID_t p, int method)
    : ObserverEvent(DynamicEventType<Dom>::destroyCode),
      domain_m(d), patch_m(p), method_m(method)
    {
      CTAssert(DynamicEventType<Dom>::dimensions == 1);
    }

  //============================================================
  // Destructor
  //============================================================

  // Nothing to see here, move along.

  virtual ~DestroyEvent()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  inline const Domain_t &domain() const
    {
      return domain_m;
    }

  inline PatchID_t patch() const
    {
      return patch_m;
    }

  inline int method() const
    {
      return method_m;
    }

private:
  // The domain of the data to destroy.

  Domain_t domain_m;

  // The patch to destroy the data from.  If this is < 0, it means the
  // domain contains values within the total domain of the target.  If it
  // is >= 0, domain_m should contain zero-based values just for the
  // specified patch.

  PatchID_t patch_m;

  // The method code

  int method_m;
};


/**
 * CopyEvent: A class derived from ObserverEvent that stores information
 * on how what elements to copy in an engine.  It stores info on the
 * domain to copy (domain()), the patch to copy from (fromPatch()), and
 * the patch to copy to (toPatch()).
 *
 * This is templated on the domain type, Dom.
 */

template<class Dom>
class CopyEvent : public ObserverEvent
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef DynamicEvents::PatchID_t                  PatchID_t;
  typedef DynamicEvents::CreateSize_t               CreateSize_t;
  typedef typename DynamicEventType<Dom>::Domain_t  Domain_t;


  //============================================================
  // Constructor
  //============================================================

  // There is only one way to initialize this object, by providing
  // the domain to copy, the patch to copy from, and the patch to
  // copy to.

  CopyEvent(const Dom &d, PatchID_t fromp, PatchID_t top)
    : ObserverEvent(DynamicEventType<Dom>::copyCode),
      domain_m(d), from_m(fromp), to_m(top)
    {
      CTAssert(DynamicEventType<Dom>::dimensions == 1);
    }

  template<class D>
  CopyEvent(const D &d, PatchID_t fromp, PatchID_t top)
    : ObserverEvent(DynamicEventType<Dom>::copyCode),
      domain_m(d), from_m(fromp), to_m(top)
    {
      CTAssert(DynamicEventType<Dom>::dimensions == 1);
    }


  //============================================================
  // Destructor
  //============================================================

  // Nothing to see here, move along.

  virtual ~CopyEvent()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  inline const Domain_t &domain() const
    {
      return domain_m;
    }

  inline PatchID_t fromPatch() const
    {
      return from_m;
    }

  inline PatchID_t toPatch() const
    {
      return to_m;
    }

private:
  // The domain of the data to copy.

  Domain_t domain_m;

  // The patch to copy the data from.  If this is < 0, it means the
  // domain contains values within the total domain of the target.  If it
  // is >= 0, domain_m should contain zero-based values just for the
  // specified patch.

  PatchID_t from_m;

  // The patch to copy the data to.

  PatchID_t to_m;
};



/**
 * CopyPatchEvent: A class derived from ObserverEvent that stores information
 * on what elements to copy in an engine.  This is a special form of
 * copy that uses a lists of IndirectionLists, one for each of a set of
 * patches, and copies elements from the specified patches (also given
 * in an IndirectionList) into a requested destination patch.  This is
 * also special in that it can be used in a way that will not require
 * a "create", in case the user did this themselves earlier.  This is
 * not templated since it expects to work with int IndirectionLists's.
 *
 * The copy information is provided as an IndirectionList<IndirectionList<int>>
 * object, with one element per patch to copy from.  A second
 * IndirectionList<int> is required, of the same length as the first list,
 * with patch ID's for the patches to copy from.  The actual element indices
 * should be "relative" indices, and contained only within the relevent
 * source patch.
 */

class CopyPatchEvent : public ObserverEvent
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef DynamicEvents::PatchID_t                  PatchID_t;
  typedef DynamicEvents::CreateSize_t               CreateSize_t;
  typedef IndirectionList< IndirectionList<int> >   Domain_t;
  typedef IndirectionList<int>                      IDList_t;

  //============================================================
  // Constructor
  //============================================================

  // There is only one way to initialize this object, by providing
  // the list of IL's for a set of patches to copy, the patch ID's
  // for the lists (telling what the source of the data is), the
  // desination patch for the data, and a flag indicating if the
  // storage should be first created, or if the data should just
  // be copied into the last elements of the existing storage.

  CopyPatchEvent(const Domain_t &domlists, const IDList_t &fromlist,
		 PatchID_t top, bool create)
    : ObserverEvent(DynamicEvents::copyPatchList),
      lists_m(domlists), from_m(fromlist), to_m(top), create_m(create)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  // Nothing to see here, move along.

  virtual ~CopyPatchEvent()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  inline const Domain_t &domainLists() const
    {
      return lists_m;
    }

  inline const IDList_t &fromPatch() const
    {
      return from_m;
    }

  inline PatchID_t toPatch() const
    {
      return to_m;
    }

  inline bool create() const
    {
      return create_m;
    }

private:
  // The lists of domains of the data to copy.

  Domain_t lists_m;

  // The list of patch ID's, one for each list in the set of index lists.

  IDList_t from_m;

  // The patch to copy the data to.

  PatchID_t to_m;

  // Boolean flag; if this is true, we must also create storage instead
  // of just putting it at the end of the existing storage.

  bool create_m;
};


/**
 * SyncEvent: A class derived from ObserverEvent that stores information
 * about doing a sync.
 */

class SyncEvent : public ObserverEvent
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================


  //============================================================
  // Constructor
  //============================================================
 
  SyncEvent()
    : ObserverEvent(DynamicEvents::sync)
    {  
    }

  //============================================================
  // Destructor
  //============================================================

  virtual ~SyncEvent()
    {
    }

private:
  // data

};



// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_LAYOUT_DYNAMIC_EVENTS_H


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DynamicEvents.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
