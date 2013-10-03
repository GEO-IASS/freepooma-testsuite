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
// Dynamic-Engine non-inline template definitions.
//-----------------------------------------------------------------------------

#include "Engine/DynamicEngine.h"
#include "Domain/Contains.h"
#include "Domain/RangeIterator.h"
#include "Domain/IntervalIterator.h"
#include "Domain/IndirectionListIterator.h"
#include "Utilities/algorithms.h"
  
///////////////////////////////////////////////////////////////////////////////
//
// Dynamic-Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// Engine<1,T,Dynamic> Constructors
//
//-----------------------------------------------------------------------------

template <class T> Engine<1,T,Dynamic>::Engine() { }

template <class T>
Engine<1,T,Dynamic>::Engine(const Domain_t &dom)
  : domain_m(dom),
    data_m(dom.size()), 
    first_m(dom.first())
{ }

template <class T>
Engine<1,T,Dynamic>::Engine(const Node<Domain_t> &node)
  : domain_m(node.allocated()),
    data_m(node.allocated().size(), 
           node.affinity(), 
           typename DataBlockPtr<T>::WithAffinity_t()),
    first_m(node.allocated().first())
{ }

template <class T>
Engine<1,T,Dynamic>::Engine(const Layout_t &layout)
  : domain_m(layout.domain()),
    data_m(layout.domain().size()), 
    first_m(layout.domain().first())
{ }

template <class T>
Engine<1,T,Dynamic>::Engine(const Domain_t &dom, const T& model)
  : domain_m(dom),
    data_m(dom.size(), model), 
    first_m(dom.first())
{ }

template <class T>
Engine<1,T,Dynamic>::
Engine(const This_t &modelEngine)
  : domain_m(modelEngine.domain_m),
    data_m(modelEngine.data_m),
    first_m(modelEngine.first_m)
{
  PAssert(data_m.isAtBeginning());
}

//-----------------------------------------------------------------------------
//
// Engine<1,T,Dynamic> & operator=(const Engine<1,T,Dynamic> &)
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1,T,Dynamic> & 
Engine<1,T,Dynamic>::
operator=(const This_t &modelEngine)
{
  // Can skip the rest if we're trying to assign to ourselves

  if (this != &modelEngine)
    {
      domain_m = modelEngine.domain_m;
      data_m   = modelEngine.data_m;
      first_m  = modelEngine.first_m;

      PAssert(data_m.isAtBeginning());
    }
  return *this;
}

//-----------------------------------------------------------------------------
//
// ~Engine<1,T,Dynamic>()
//
// Out-of-line destructor to shorten compile times.
//
//-----------------------------------------------------------------------------

template <class T> Engine<1,T,Dynamic>::~Engine() { }

//-----------------------------------------------------------------------------
//
// Engine<1,T,Dynamic> & makeOwnCopy()
//
// Causes the Dynamic-Engine to obtain a private copy of the data
// that it refers to.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1,T,Dynamic> &
Engine<1,T,Dynamic>::
makeOwnCopy()
{
  if (data_m.isValid() && data_m.count() > 1) 
    {
      PAssert(data_m.isAtBeginning());
      data_m.makeOwnCopy();
    }

  return *this;
}

//-----------------------------------------------------------------------------
//
// Interval<1> create(CreateSize_t num)
//
// Create new elements by resizing the DataBlockPtr and updating the
// layout. Return Interval describing the new elements.
//
//-----------------------------------------------------------------------------

template <class T> 
Interval<1> 
Engine<1,T,Dynamic>::
create(CreateSize_t num)
{ 
  PAssert(num >= 0);

  // It would be nice to assert that no-one else is looking at the engine
  // when we perform dynamic operations, but all the particle swap operations
  // take place inside iterates which means the engine is a copy of another
  // engine, so the data is shared.
  //
  //  PAssert(!data_m.isShared());
  
  // Reallocate the storage

  data_m.resizeAndCopy(domain_m.size() + num); // initialize new elements!

  // Reset the domain (in the layout) to the new size

  const int newend = domain_m.last()  + num;

  domain_m = Interval<1>(domain().first(), newend);
				  
  PAssert(first_m == domain().first());
  
  // Return domain describing new elements
  
  const int newstart = domain_m.first() + num;

  return Interval<1>(newstart, newend);
}

//-----------------------------------------------------------------------------
//
// void destroy(...)
//
// Dynamic destroy interface. These bounce out to specialized functions
// for the appropriate fill methods.
//
//-----------------------------------------------------------------------------

template <class T>
template <class Dom>
void Engine<1,T,Dynamic>::destroy(const Dom &killList)
{
  //  PAssert(!data_m.isShared());
  performDestroy(killList, BackFill(), false);
}

template <class T>
template <class Iter>
void Engine<1,T,Dynamic>::destroy(Iter begin, Iter end)
{
  //  PAssert(!data_m.isShared());
  performDestroy(begin, end, BackFill(), false);
}

template <class T>
template <class Dom, class DeleteMethod>
void Engine<1,T,Dynamic>::
destroy(const Dom &killList, const DeleteMethod &method, bool offsetFlag)
{
  //  PAssert(!data_m.isShared());
  performDestroy(killList, method, offsetFlag);
}

template <class T>
template <class Iter, class DeleteMethod>
void Engine<1,T,Dynamic>::
destroy(Iter begin, Iter end, const DeleteMethod &method, bool offsetFlag)
{
  //  PAssert(!data_m.isShared());
  performDestroy(begin, end, method, offsetFlag);
}

// Utility class for generalizing the mechanism of getting the begin
// and end iterators for the performDestroy implementations.

template <class Domain>
struct DeleteDomainIteratorTraits
{
  typedef Domain                    Domain_t;
  typedef typename Domain::iterator Type_t;
  static Type_t begin(const Domain &d) { return d.begin(); }
  static Type_t end(const Domain &d)   { return d.end(); }
};

template <>
struct DeleteDomainIteratorTraits< Interval<1> >
{
  typedef Interval<1>      Domain_t;
  typedef IntervalIterator Type_t;
  static Type_t begin(const Domain_t &d) 
  { 
    return IntervalIterator(d); 
  }
  static Type_t end(const Domain_t &d)
  { 
    return IntervalIterator(d,d.size());
  }
};
 
template <>
struct DeleteDomainIteratorTraits< Range<1> >
{
  typedef Range<1>      Domain_t;
  typedef RangeIterator Type_t;
  static Type_t begin(const Domain_t &d) 
  { 
    return RangeIterator(d); 
  }
  static Type_t end(const Domain_t &d)
  { 
    return RangeIterator(d,d.size());
  }
};

template <>
struct DeleteDomainIteratorTraits< IndirectionList<int> >
{
  typedef IndirectionList<int>         Domain_t;
  typedef IndirectionListIterator<int> Type_t;
  static Type_t begin(const Domain_t &d) 
  { 
    return Type_t(d); 
  }
  static Type_t end(const Domain_t &d)
  { 
    return Type_t(d,d.size());
  }
};
 
//-----------------------------------------------------------------------------
//
// void performDestroy(const Domain &killList, BackFill, bool offsetFlag)
// void performDestroy(const Iter &killBegin, const Iter &killEnd,
//                     BackFill, bool offsetFlag)
//
// Setup and perform the destroy calculation using a backfill strategy.
// This option fills in the deleted elements with data from the end of 
// the list. Thus we must make N assignments, where N is the number of 
// elements killed. This is much less data shuffling than with shift 
// types of deletions, but it does not preserve the order of the data.
//
// Note that we could provide more efficient specializations for
// certain domains, notably Intervals. If these are frequently used,
// we should do so. For now, we'll go with the generic versions.
//
//-----------------------------------------------------------------------------

template <class T> 
template <class Iterator>
void 
Engine<1,T,Dynamic>::
performDestroy(const Iterator &killBegin,
               const Iterator &killEnd,
	       const BackFill &,
	       bool offsetFlag)
{
  // If offsetFlag is false (the default), the points in the killList
  // are interpreted as indices into the array; i.e. these points are 
  // a subset of our domain. This domain may not be zero-based, but
  // the underlying data block is, so we have to pass an offset value.
  //
  // If offsetFlag is true, the points in the killList are interpreted
  // to be offsets from the beginning of the array, so no additional 
  // offset is needed. 

  int koffset = (offsetFlag ? 0 : first_m);
  
  // Use the generic delete algorithm to do the work.
  
  int killed =
    Pooma::Algorithms::delete_backfill(data_m.begin(), data_m.end(),
                                       killBegin, killEnd, koffset);

  // Update the domain.
  // We have to handle the case of zero size specially since you can't
  // specify a zero-sized interval's "first" point. Thus we first create
  // a zero sized interval and then assign the appropriate interval to it
  // if there are elements left.
  
  Interval<1> newdom;
  
  if (killed < domain().size())
    newdom = Interval<1>(domain().first(), domain().last() - killed);

  domain_m = newdom;
  
  // Resize the data block to the new domain size.

  data_m.resize(domain().size(), typename DataBlockPtr<T>::NoInitTag());
}

// Now the version that takes a general domain. 
// This is implemented using the above version and the 
// DeleteDomainIteratorTraits class.

template <class T> 
template <class Domain>
void 
Engine<1,T,Dynamic>::
performDestroy(const Domain &killList,
	       const BackFill &,
	       bool offsetFlag)
{
  PAssert(killList.length() <= domain().length());
  if (offsetFlag)
    {
      PAssert( contains( Interval<1>(domain().length()),
        Interval<1>(killList.min(), killList.max()) ) );
    }
  else
    {
      PAssert( contains( Interval<1>(domain().length()),
        Interval<1>(killList.min()-first_m, killList.max()-first_m) ) );
    }
  
  typedef DeleteDomainIteratorTraits<Domain> Trait_t;

  performDestroy(Trait_t::begin(killList), 
                 Trait_t::end(killList),
                 BackFill(), 
                 offsetFlag);
}

//-----------------------------------------------------------------------------
//
// void performDestroy(const Domain &killList, ShiftUp)
// void performDestroy(const Iter &killBegin, const Iter &killEnd,
//                     ShiftUp, bool offsetFlag)
//
// Do the actual work of performing a destroy op, for a generic domain with
// a ShiftUp operation. 
// 
// See comments for Backfill version above.
//
//-----------------------------------------------------------------------------

template <class T> 
template <class Iterator>
void 
Engine<1,T,Dynamic>::
performDestroy(const Iterator &killBegin,
               const Iterator &killEnd,
               const ShiftUp &,
               bool offsetFlag)
{
  int koffset = (offsetFlag ? 0 : first_m);

  int killed =
    Pooma::Algorithms::delete_shiftup(data_m.begin(), data_m.end(),
                                      killBegin, killEnd, koffset);

  Interval<1> newdom;
  
  if (killed < domain().size())
    newdom = Interval<1>(domain().first(), domain().last() - killed);

  domain_m = newdom;
  
  data_m.resize(domain().size(), typename DataBlockPtr<T>::NoInitTag());
}

template <class T> 
template <class Domain>
void 
Engine<1,T,Dynamic>::
performDestroy(const Domain &killList, const ShiftUp &, bool offsetFlag)
{
  PAssert(killList.length() <= domain().length());
  if (offsetFlag)
    {
      PAssert( contains( Interval<1>(domain().length()),
        Interval<1>(killList.min(), killList.max()) ) );
    }
  else
    {
      PAssert( contains( Interval<1>(domain().length()),
        Interval<1>(killList.min()-first_m, killList.max()-first_m) ) );
    }

  typedef DeleteDomainIteratorTraits<Domain> Trait_t;

  performDestroy(Trait_t::begin(killList), 
                 Trait_t::end(killList),
                 ShiftUp(), 
                 offsetFlag);
}


//-----------------------------------------------------------------------------
//
// void sync(Domain_t)
//
// modify the domain (but not the size) of this engine.
//
//-----------------------------------------------------------------------------

template <class T>
void
Engine<1,T,Dynamic>::sync(const Domain_t & d)
{
  // Modify the block pointer's domain to reflect the new relative domain.

  first_m = d.first();

  // Update the domain

  domain_m = d;
}


///////////////////////////////////////////////////////////////////////////////
//
// DynamicView Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// Engine<1,T,DynamicView> Constructors
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1,T,DynamicView>::
Engine(const Engine<1,T,Dynamic> &engine, const Interval<1> &dom)
: domain_m(Interval<1>(dom.length())),
  data_m(engine.dataBlock(), dom.first() - engine.domain().first()),
  stride_m(1)
{    
  // The engine's data pointer should be at the beginning.
  PAssert(engine.dataBlock().isAtBeginning());
}

template <class T>
Engine<1,T,DynamicView>::
Engine(const Engine<1,T,Dynamic> &engine, const Range<1> &dom)
: domain_m(Interval<1>(dom.length())),
  data_m(engine.dataBlock(), dom.first() - engine.domain().first()),
  stride_m(dom.stride())
{    
  // The engine's data pointer should be at the beginning.
  PAssert(engine.dataBlock().isAtBeginning());
}

template <class T>
Engine<1,T,DynamicView>::
Engine(const This_t &engine, const Interval<1> &dom)
: domain_m(Interval<1>(dom.length())),
  data_m(engine.dataBlock(), engine.stride_m * dom.first()),
  stride_m(engine.stride_m)
{ }

template <class T>
Engine<1,T,DynamicView>::
Engine(const This_t &engine, const Range<1> &dom)
: domain_m(Interval<1>(dom.length())),
  data_m(engine.dataBlock(), engine.stride_m * dom.first()),
  stride_m(engine.stride_m*dom.stride())
{ }

template <class T>
Engine<1,T,DynamicView>::
Engine(const Engine_t &engine, const INode<1> &inode)
: domain_m(Interval<1>(inode.domain().length())),
  data_m(engine.dataBlock(), engine.stride_m * inode.domain().first()),
  stride_m(engine.stride_m)
{ }

//-----------------------------------------------------------------------------
//
// Engine<1,T,DynamicView>::
// Engine(const Engine &);
//
// Copy constructor.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1,T,DynamicView>::
Engine(const This_t &modelEngine)
  : domain_m(modelEngine.domain_m),
    data_m(modelEngine.data_m),
    stride_m(modelEngine.stride_m)
{ }

template <class T>
Engine<1,T,DynamicView>::
Engine(const This_t &modelEngine, const EngineConstructTag &)
  : domain_m(modelEngine.domain_m),
    data_m(modelEngine.data_m),
    stride_m(modelEngine.stride_m)
{ }

//-----------------------------------------------------------------------------
//
// ~Engine<1,T,DynamicView>()
//
// Moved destructor out-of-line to shorten compile time.  It's non-trivial 
// code since it contains a ref-counted-block
//
//-----------------------------------------------------------------------------

template <class T> Engine<1,T,DynamicView>:: ~Engine() { }

//-----------------------------------------------------------------------------
//
// Engine<1,T,DynamicView> & operator=(const Engine<1,T,DynamicView> &)
//
// Assignment operator.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1,T,DynamicView> &
Engine<1,T,DynamicView>::
operator=(const This_t &modelEngine)
{
  if (this != &modelEngine)
    {
      data_m       = modelEngine.data_m;
      domain_m     = modelEngine.domain_m;
      stride_m     = modelEngine.stride_m;
    }
  return *this;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DynamicEngine.cpp,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
