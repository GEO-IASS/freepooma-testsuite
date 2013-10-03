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
//   Dynamic                     - Dynamic-Engine specialization tag.
//   DynamicView                 - DynamicView-Engine specialization tag.
//   Engine<1,T,Dynamic>         - the "Dynamic-Engine" specialization.
//   Engine<1,T,DynamicView>     - the "DynamicView-Engine" specialization.
//
// Functions:
//   checkDynamicID
//
// Traits Classes (specialized for Dynamic & DynamicView engines):
//   NewEngine<Engine,SubDomain> 
//   NewEngineDomain<Engine,SubDomain>
//   ElementProperties<
//     Engine<1, T, Dynamic> > 
//   
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_DYNAMICENGINE_H
#define POOMA_ENGINE_DYNAMICENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 * Dynamic Engine
 *  - Dynamic & DynamicView
 *    tag classes used to select specializations of Engine
 *  - Engine<1,T,Dynamic>
 *    an Engine that manages a contiguous, local, resizable, 
 *    1-dimensional block of data. 
 *    See Engine.h for the general template.
 *  - Engine<1,T,DynamicView>
 *    an Engine that manages a view into a Dynamic-Engine.
 *  - NewEngine
 *  - NewEngineDomain
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/Loc.h"
#include "Domain/Touches.h"
#include "Engine/Engine.h"
#include "Layout/DomainLayout.h"
#include "Layout/Node.h"
#include "Layout/INode.h"
#include "Utilities/DataBlockPtr.h"
#include "Utilities/PAssert.h"

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

/**
 * These are tag classes used to select the "Dynamic" and "DynamicView" 
 * specializations of the Engine class template.
 */

struct Dynamic 
{ };

/**
 * These are tag classes used to select the "Dynamic" and "DynamicView" 
 * specializations of the Engine class template.
 */

struct DynamicView 
{ };


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

template <class T>
class Engine<1,T,DynamicView>;


/**
 *  Engine<1,T,Dynamic> (aka Dynamic-Engine) is an Engine that manages
 *  a contiguous, local, 1-dimensional, dynamically resizable, block of data.
 *
 *  Template Parameters:
 *   - T:  The type of object stored. 
 *         The only assumption made about T is that it have a copy
 *         constructor, and this is only required if the read method
 *         is invoked, which returns a T by value. All other properties
 *         of T are deferred to the ElementProperties class.
 *
 *  The Domain of this engine is an Interval<1>.
 *
 *  Subsetting Engine<1,T,Dynamic> returns an Engine<1,T,DynamicView>.
 *  See below.
 */

template <class T>
class Engine<1,T,Dynamic> 
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Dynamic                          Tag_t;
  typedef Engine<1,T,Tag_t>                This_t;
  typedef Engine<1,T,Tag_t>                Engine_t;
  typedef DomainLayout<1>                  Layout_t;
  typedef Layout_t::PatchID_t              PatchID_t;
  typedef Layout_t::CreateSize_t           CreateSize_t;
  typedef Interval<1>                      Domain_t;
  typedef T                                Element_t;
  typedef T&                               ElementRef_t;

  enum { brick         = true  };
  enum { dimensions    = 1     };
  enum { hasDataObject = true  };
  enum { dynamic       = true  };
  enum { zeroBased     = false };
  enum { multiPatch    = false };

  //============================================================
  // Constructors and Factory Methods
  //============================================================

  // Default constructor. Creates a Dynamic-Engine with no
  // data and an "empty" domain.  This is not really usable until it
  // has been initialized (via operator=) to a new engine with an
  // actual domain.

  Engine();

  // These constructors take an Interval<1> and create a new
  // Dynamic-Engine with data of type T on this Domain. This is where
  // storage gets allocated. The second version uses a model data
  // element to initialize storage.
  
  explicit 
  Engine(const Domain_t &domain);

  Engine(const Domain_t &domain, const T &elementModel);

  // You can build a dynamic-engine from a layout as well.

  explicit 
  Engine(const Layout_t &layout);

  // This constructor takes a Node object, extracts the domain, and
  // creates a new Dynamic-Engine with data of type T on this Domain. 
  // Use this if you want to specify the thread affinity of the patch.
  
  explicit 
  Engine(const Node<Domain_t> &node);

  // Copy constructor performs a SHALLOW copy.
  // But NOTE: the layouts will NOT be shared.

  Engine(const Engine_t &model);

  // Subsetting Constructors.
  
  Engine(const Engine_t &model, const Interval<1> &domain);
  
  //============================================================
  // Destructor
  //============================================================

  // All pointer members are smart pointers, so this is trivial.

  ~Engine(); 

  //============================================================
  // Assignment operators
  //============================================================

  // Assigment is SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);

  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  // Element access via Loc.

  inline Element_t read(const Loc<1> &l) const 
  { return data_m[l.first() - first_m]; }
  
  inline ElementRef_t operator()(const Loc<1> &l) const 
  { return data_m[l.first() - first_m]; }

  // Element access via ints for speed.

  inline Element_t read(int i) const { return data_m[i - first_m]; };
  inline ElementRef_t operator()(int i) const { return data_m[i - first_m]; };

  // Return the domain.

  inline const Domain_t &domain() const { return domain_m; }

  // Create a layout and return a copy. 
  
  inline Layout_t layout() const { return Layout_t(domain_m); }

  // Return whether the block controlled by this engine is shared.
  
  inline bool isShared() const { return data_m.isValid() && data_m.count() > 1; }

  // Get a private copy of data viewed by this Engine.

  Engine_t &makeOwnCopy();

  // Provide access to the data object. 

  inline Pooma::DataObject_t *dataObject() const { return data_m.dataObject(); }

  // Return access to our internal data block.  This is ref-counted,
  // so a copy is fine.  But you should really know what you're doing
  // if you call this method.

  inline const DataBlockPtr<T> & dataBlock() const { return data_m; }
  inline DataBlockPtr<T>         dataBlock()       { return data_m; }

  //============================================================
  // Dynamic interface methods.
  //============================================================

  // Create new elements by extending the current domain
  // on the local context by the requested number of elements.
  // Returns an Interval giving the domain of the newly created
  // elements.

  Interval<1> create(CreateSize_t num);

  // Delete the elements specified by the given domain. 
  // This backfills the deleted elements with elements from
  // the end of the list.
  
  template <class Dom>
  void destroy(const Dom &killList);

  // Same, but with iterators into some container holding the
  // points of the domain. These must be random-access iterators (a
  // requirement of the underlying delete algorithm).
  
  template <class Iter>
  void destroy(Iter begin, Iter end);

  // Delete the elements specified by the given domain, or by a
  // pair of iterators into some sort of collection, and the
  // appropriate fill method. If offsetFlag is set to true, the 
  // domain is interpreted as a set of offsets rather than a set
  // of points in our domain.

  // Available fill mechanisms are backfill and shift-up, selected by
  // passing either a BackFill or ShiftUp object. BackFill will move 
  // elements from the end up to fill the holes. ShiftUp will shift
  // elements up to fill in holes. The latter is considerably slower,
  // but maintains the relative ordering of the elements, which may
  // be important for some applications.

  template <class Dom, class DeleteMethod>
  void destroy(const Dom &killList, 
               const DeleteMethod &method,
               bool offsetFlag = false);
                      
  template <class Iter, class DeleteMethod>
  void destroy(Iter begin, 
               Iter end, 
               const DeleteMethod &method,
               bool offsetFlag = false);

  // sync() function is a no-op for a single-patch engine.
  // This version of sync() may be called via the DynamicArray interface.

  inline void sync() { }

  // Modify the domain (but not the size) of this engine.
  // This version of sync() may be called by MultiPatchEngine on its patches.

  void sync(const Domain_t & d);

#if POOMA_MESSAGING

  template <class Dom>
  inline int packSize(const Dom &) const
  {
    PInsist(false,"packSize() called on non-remote Dynamic Engine!!");
    return 0;
  }

  inline int pack(const IndirectionList<int> &, char *, bool = true) const
  {
    PInsist(false,"pack() called on non-remote Dynamic Engine!!");
    return 0;
  }

  inline int unpack(const Interval<1> &, char *, bool = true)
  {
    PInsist(false,"unpack() called on non-remote Dynamic Engine!!");
    return 0;
  }

#endif

private:

  //============================================================
  // Private dynamic-engine methods
  //============================================================

  // Implementations of the destroy operation.
  // These are specialized on the fill strategy.
  
  // The version 2.2 BrickEngine had specializations on the domain 
  // type, and there was a default generic version that would build an 
  // IndirectionList from the generic domain and then call the 
  // IndirectionList version.  I've scrapped this for now to simplify 
  // the code. The Interval versions would be more efficient than the 
  // current code, so if this is frequently requested, we could 
  // restore it (in fact, we could do it more efficiently than the
  // previous version). The generic version was really just a method to 
  // allow the killList to be stored in an array. Instead, I've added
  // versions that take iterators into the killList. This allows one
  // to use many more data structures for holding the kill list.
  
  // If offsetFlag is true (delete defaults it to false), the domain 
  // is interpreted as a list of offsets rather than a subset of 
  // the engine's domain.

  template <class Domain>
  void performDestroy(const Domain &domain, 
                      const BackFill &method,
                      bool offsetFlag);

  template <class Iterator>
  void performDestroy(const Iterator &killBegin,
                      const Iterator &killEnd,
                      const BackFill &method,
                      bool offsetFlag);

  template <class Domain>
  void performDestroy(const Domain &domain, 
                      const ShiftUp &method,
                      bool offsetFlag);

  template <class Iterator>
  void performDestroy(const Iterator &killBegin,
                      const Iterator &killEnd,
                      const ShiftUp &method,
                      bool offsetFlag);

  //============================================================
  // Private data
  //============================================================

  // Layout for this engine

  Domain_t domain_m;

  // Smart-pointer to Block-controller that manages the data
  // and the Smarts DataObject. 

  DataBlockPtr<T> data_m;

  // Index of the first point.
  
  int first_m;
};


/**
 *  A DynamicView-Engine is an Engine that manages a view of a
 *  Dynamic-Engine.  See Engine<1,T,Dynamic> for details.
 *
 *  Template Parameters:
 *   - T, the type of object stored. 
 *
 *  The Domain of this engine is an Interval<1>. For DynamicView-Engines, 
 *  these intervals will all be 0-based (i.e. [0..N]).
 *  
 *  Note that this is NOT the domain of the underlying data storage,
 *  but rather it is the domain as presented to the outside world.
 */

template <class T>
class Engine<1,T,DynamicView>
{
public:
  
  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef DynamicView                               Tag_t;
  typedef Engine<1,T,Tag_t>                         This_t;
  typedef Engine<1,T,Tag_t>                         Engine_t;
  typedef Interval<1>                               Domain_t;
  typedef T                                         Element_t;
  typedef T&                                        ElementRef_t;
  typedef DomainLayout<1>                           Layout_t;

  enum { dimensions    = 1     };
  enum { hasDataObject = true  };
  enum { dynamic       = false };
  enum { zeroBased     = true  };
  enum { multiPatch    = false };
  
  //============================================================
  // Constructors
  //============================================================

  // DynamicView-Engine is fundamentally a view - it never owns
  // its data, and thus there are no constructors that create
  // a DynamicView-Engine directly from a domain.

  // Copy constructor performs a SHALLOW copy:

  Engine(const Engine_t &);
  Engine(const Engine_t &, const EngineConstructTag &);

  // Subsetting Constructors.
  // A DynamicView-Engine is a strided view of a Dynamic-Engine.
  // Thus we write constructors to build DynamicViews from
  // Dynamic-Engines and all types 1-D Domains.

  // Build a DynamicView from Dynamic Engine and an Interval or Range.

  Engine(const Engine<1,T,Dynamic> &, const Interval<1> &);
  Engine(const Engine<1,T,Dynamic> &, const Range<1> &);

  // Build a DynamicView from another DynamicView and an Interval,
  // Range, or INode.
  
  Engine(const Engine_t &, const Interval<1> &);
  Engine(const Engine_t &, const Range<1> &);
  Engine(const Engine_t &, const INode<1> &);

  //============================================================
  // Destructor
  //============================================================

  ~Engine();

  //============================================================
  // Assignment operators
  //============================================================

  // Assigment is SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);

  //============================================================
  // Accessor functions:
  //============================================================

  // Element access via Loc.

  inline Element_t read(const Loc<1> &l) const 
  { return data_m[l.first()*stride_m]; }
  
  inline ElementRef_t operator()(const Loc<1> &l) const
  { return data_m[l.first()*stride_m]; }

  // Element access via ints for speed.

  inline Element_t read(int i)          const { return data_m[i*stride_m]; }
  inline ElementRef_t operator()(int i) const { return data_m[i*stride_m]; }

  // Return the domain:

  inline const Domain_t &domain() const { return domain_m; }

  // Return a DomainLayout built from our domain
  
  inline Layout_t layout() const { return Layout_t(domain_m); }

  // Return the stride.

  inline int stride() const { return stride_m; }
  
  // Provide access to the data object. 

  inline Pooma::DataObject_t *dataObject() const { return data_m.dataObject(); }

  // Return access to our internal data block.  This is ref-counted,
  // so a copy is fine.  But you should really know what you're doing
  // if you call this method.

  inline const DataBlockPtr<T> & dataBlock() const { return data_m; }
  inline DataBlockPtr<T>         dataBlock()       { return data_m; }

private:

  //============================================================
  // Data
  //============================================================

  // Domain for this engine:

  Domain_t domain_m;


  // Copy of the Block-controller that manages the data.

  DataBlockPtr<T> data_m;

  // Stride. 

  int stride_m;
};


//////////////////////////////////////////////////////////////////////
//
// Traits classes and utility functions...
//
//////////////////////////////////////////////////////////////////////

/**
 * Several specializations of NewEngine for combinations of 
 * Engines and Domains that produce DynamicView-Engines.
 */

template <class T>
struct NewEngine<Engine<1,T,Dynamic>, Interval<1> >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,Dynamic>, Range<1> >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,Dynamic>, Node<Interval<1> > >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,Dynamic>, INode<1> >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,DynamicView>, Interval<1> >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,DynamicView>, Range<1> >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,DynamicView>, Node<Interval<1> > >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

template <class T>
struct NewEngine<Engine<1,T,DynamicView>, INode<1> >
{
  typedef Engine<1,T,DynamicView> Type_t;
};

/**
 * Several specializations of NewEngineDomain for combinations of 
 * Dynamic/DynamicView-Engines and Node/INode..
 */

template <class T>
struct NewEngineDomain<Engine<1,T,Dynamic>, Node<Interval<1> > >
{
  typedef Interval<1> Type_t;
  typedef const Interval<1> &Return_t;
  static inline
  Return_t apply(const Engine<1,T,Dynamic> &,
		 const Node<Interval<1> > &node)
  {
    return node.domain();
  }
};

template <class T>
struct NewEngineDomain<Engine<1,T,Dynamic>, INode<1> >
{
  typedef Interval<1> Type_t;
  typedef const Interval<1> &Return_t;
  static inline
  Return_t apply(const Engine<1,T,Dynamic> &,
		 const INode<1> &inode)
  {
    return inode.domain();
  }
};

template <class T>
struct NewEngineDomain<Engine<1,T,DynamicView>, Node<Interval<1> > >
{
  typedef Interval<1> Type_t;
  typedef const Interval<1> &Return_t;
  static inline
  Return_t apply(const Engine<1,T,DynamicView> &,
		 const Node<Interval<1> > &node)
  {
    return node.domain();
  }
};

template <class T>
struct NewEngineDomain<Engine<1,T,DynamicView>, INode<1> >
{
  typedef Interval<1> Type_t;
  typedef const Interval<1> &Return_t;
  static inline
  Return_t apply(const Engine<1,T,DynamicView> &,
		 const INode<1> &inode)
  {
    return inode.domain();
  }
};

/**
 * Traits class telling RefCountedBlockPointer that this class has
 * shallow semantics and a makeOwnCopy method.
 */

template <class T>
struct ElementProperties<Engine<1, T, Dynamic> > 
  : public MakeOwnCopyProperties<Engine<1, T, Dynamic> >
{ };

/**
 * checkDynamicID(obj, ID) is a specializable function that is used
 * by some classes to check the dynamic ID value stored in the first
 * argument by some means.  If it is the same as the given ID, this
 * returns false.  If it is not the same, it should return true and
 * changethe state of obj to indicate that it has "seen" the given ID.
 *
 * For Dynamic, the dynamic ID is stored in the data block.
 */

template <class T>
inline bool checkDynamicID(Engine<1,T,Dynamic> &be, ObserverEvent::ID_t did)
{
  if (did == be.dataBlock().dynamicID())
    return false;

  be.dataBlock().setDynamicID(did);
  return true;
}



// Include .cpp file to get out-of-line functions.

#include "Engine/DynamicEngine.cpp"

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_ENGINE_DYNAMICENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DynamicEngine.h,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
