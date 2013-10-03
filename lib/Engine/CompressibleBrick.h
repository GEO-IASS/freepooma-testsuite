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


/** @file
 * @ingroup Engine
 * @brief
 * CompressibleBrick engine.
 *
 * Classes: 
 *  - CompressibleBrick, CompressibleBrick-Engine specialization tag.
 *  - CompressibleBrickView<int,bool>, CompressibleBrickView-Engine 
 *    specialization tag.
 *  - Engine<Dim,T,CompressibleBrick>, the "CompressibleBrick-Engine"
 *    specialization.
 *  - Engine<Dim,T,CompressibleBrickView>, the "CompressibleBrickView-Engine"
 *    specialization.
 *  - NewEngine<Engine,SubDomain>, specializations for 
 *    CompressibleBrickView-Engine.
 *  - ElementProperties<Engine <Dim, T, CompressibleBrick> >, specialization
 *    for CompressibleBrick-Engine.
 */

#ifndef POOMA_ENGINE_COMPRESSIBLEBRICK_H
#define POOMA_ENGINE_COMPRESSIBLEBRICK_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/BrickBase.h"

#include "Domain/Domain.h"
#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Domain/Range.h"
#include "Domain/SliceDomain.h"
#include "Engine/CompressibleBlock.h"
#include "Engine/Engine.h"
#include "Layout/DomainLayout.h"
#include "Layout/Node.h"
#include "Layout/INode.h"
#include "Threads/PoomaMutex.h"
#include "Utilities/Observer.h"
#include "Utilities/Observable.h"
#include "Utilities/DataBlockPtr.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//@{

/**
 * These are tag classes used to select the "CompressibleBrick" and 
 * "CompressibleBrickView" specializations of the Engine class template.
 */

struct CompressibleBrick { };

struct CompressibleBrickView { };

//@}

//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

template <int Dim, class T>
class Engine<Dim, T, CompressibleBrickView>;

template <int D1, int D2> class SliceInterval;
template <int D1, int D2> class SliceRange;

/**
 * Engine<Dim,T,CompressibleBrick>  (aka CompressibleBrick-Engine)
 *
 * Engine<Dim,T,CompressibleBrick> is an Engine that manages a contiguous, 
 * local, Dim-dimensional brick of data that can compress itself to a single
 * value if all the values are the same.
 *
 * Template Parameters:
 *  - Dim: An integer for the dimension of the CompressibleBrick.
 *  - T:   The type of object stored.  Often a numeric type.
 *         However, the Engine makes no assumptions beyond that
 *         T have a copy constructor and assignment operator.
 *
 *  The Domain of this engine is an Interval<Dim>, which is a 
 *  tensor product of Dim 1-D intervals.
 *
 *  Subsetting Engine<Dim,T,CompressibleBrick> returns an 
 *  Engine<Dim,T,CompressibleBrickView>. See below.
 */  

template <int Dim, class T>
class Engine<Dim, T, CompressibleBrick> 
  : public Pooma::BrickBase<Dim>, public Observer<T*>
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Engine<Dim,T,CompressibleBrick>    This_t;
  typedef Engine<Dim,T,CompressibleBrick>    Engine_t;
  typedef Pooma::BrickBase<Dim>              Base_t;
  typedef typename Base_t::Domain_t          Domain_t;
  typedef DomainLayout<Dim>                  Layout_t;
  typedef T                                  Element_t;
  typedef T&                                 ElementRef_t;
  typedef CompressibleBrick                  Tag_t;

  enum { brick         = true  };
  enum { dimensions    = Dim   };
  enum { hasDataObject = true  };
  enum { dynamic       = false };
  enum { zeroBased     = false };
  enum { multiPatch    = false };
  
  //============================================================
  // Constructors and Factory Methods
  //============================================================

  // Default constructor. Creates a CompressibleBrick-Engine with no
  // data and an "empty" domain.

  Engine() : data0_m(0) { }

  // Construct a CompressibleBrick-Engine representing a Dim-dimensional
  // brick (Fortran storage order) of elements of type T. The domain can
  // be specified directly or by passing a Node or Layout_t object. If
  // initializing with a Domain_t, one can optionally pass in a model 
  // element, which will be used to initialize storage. CompressibleBricks
  // are always born compressed, storing only a single value, which is 
  // initialized with ElementProperties::construct.
 
  explicit Engine(const Domain_t &domain);
  Engine(const Domain_t &domain, const T &elementModel);
  explicit Engine(const Layout_t &layout);  
  explicit Engine(const Node<Domain_t> &node);

  // Copy constructor performs a SHALLOW copy:
  // But NOTE: the layouts will NOT be shared.

  Engine(const Engine_t &model);
  
  // Subsetting Constructors.
  //   There are none - you cannot create a CompressibleBrick-Engine
  //   by taking a "view" of another Engine.  CompressibleBrick-
  //   Engines, by definition, view all of the data. See
  //   Engine<D,T,CompressibleBrickView> below.

  //============================================================
  // Destructor
  //============================================================

  // Stop observing the block controller on destruct.

  ~Engine();

  //============================================================
  // Assignment operators
  //============================================================

  // Assigment is SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);

  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  // Element access.
  // These do not lock the underlying block. The run-time system will
  // ensure that data-parallel statements use these safely. If the
  // user wishes to access individual elements in his code, the onus
  // is on him or her to ensure synchronization.

  // Read-write access via ints.
  // These functions are slow, since they have to check if the
  // block is compressed and uncompress it if it is.

  ElementRef_t operator()(int) const;
  ElementRef_t operator()(int, int) const;
  ElementRef_t operator()(int, int, int) const;
  ElementRef_t operator()(int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int, int, int) const;

  // Read-only access.
  // Guaranteed to be fast since there's no if-test.

  Element_t read(int) const;
  Element_t read(int, int) const;
  Element_t read(int, int, int) const;
  Element_t read(int, int, int, int) const;
  Element_t read(int, int, int, int, int) const;
  Element_t read(int, int, int, int, int, int) const;
  Element_t read(int, int, int, int, int, int, int) const;

  // Element access via Loc.

  ElementRef_t operator()(const Loc<Dim> &) const;
  Element_t read(const Loc<Dim> &) const;

  //---------------------------------------------------------------------------
  // Return the domain and base domain.

  inline const Domain_t &domain() const
  {
    return this->layout_m.domain();
  }

  // Get a private copy of data viewed by this Engine.

  Engine_t &makeOwnCopy();
  
  // Provide access to the data object.

  Pooma::DataObject_t *dataObject() const { return cblock_m.dataObject(); }

  // Return access to our internal data block (uncompressed).  
  // This is ref-counted, so a copy is fine.  But you should really 
  // know what you're doing if you call this method.
  // Taking this view renders the CBC in the "incompressible"
  // state, so any tryCompress will fail until all such views
  // go away.
  // (Brick has two of these (const & non-const) - why doesn't this???)
  
  DataBlockPtr<T> dataBlock() const { return cblock_m.view(); }

  // isShared???
  
  //============================================================
  // Compressibility-related Accessors & Mutators
  //============================================================
  
  // Return the CompressibleBlock.
  
  CompressibleBlock<T> cblock() const { return cblock_m; }

  // WARNING: These "Compress..." messages return mutable status
  // of the underlying CBC at a particular point in time. In a 
  // multithreaded environment, it is not guaranteed that the CBC 
  // will stay in that state since its state can be changed by 
  // other threads. 
  
  // Also, see comments in Evaluator/CompressibleEval.h. 
  
  // Check to see if this CompressibleBrick is compressed.
  
  bool compressed() const;

  // Return the number of compressed elements.
  
  long elementsCompressed() const;
  
  // Try to perform a compress.
  
  void tryCompress() { cblock_m.tryCompress(); }
  
  // Manually uncompress.
  
  void uncompress() { cblock_m.uncompress(); }
  
  // Read & read-write access to the compressed value.

  T compressedRead() const;
  T& compressedReadWrite() const;

  // Determine whether or not this view is over the whole domain.
  
  bool compressedBrickIsWholeView() const { return true; }

private:

  //============================================================
  // Data
  //============================================================

  // Smart-pointer to the controller.

  CompressibleBlock<T> cblock_m;

  // Plain pointer that points to the beginning of the actual data. 
  // If uncompressed, all index offsets are subtracted during
  // construction. If compressed, this simply points to the
  // compressed value.

  T *data0_m;

  // Mutex protection for the CompressibleBricks.
  // The strides_m and data0_m members can be changed asynchronously
  // by the notify method (when someone else compresses or uncompresses
  // a copy of the underlying CBC). Thus it is necessary to protect
  // most accesses to these members. Since Smarts will prevent an
  // block from being read to and written to at the same time, it
  // is not necessary to do the locking in the read() methods. 
  
  mutable Pooma::Mutex_t mutex_m;

  //============================================================
  // Private methods
  //============================================================

  // Mutex functions

  inline void lock()   const { mutex_m.lock(); }
  inline void unlock() const { mutex_m.unlock(); }

  // Notify function for the Observer<T*> base class:
  // Compressible bricks observe the CompressibleBlock
  // which notifies us when the data becomes compressed or uncompressed.
  // The notification comes with a pointer to the new data.
  
  // Note: The CBC is locked when this is called, so we don't have to
  // worry about contention for changing strides_m/data0_m, but we
  // do need to make sure no one tries to make a copy of these data
  // while they are being changed. Thus we lock our mutex.

  // Also note: A corallary is that if you're going to lock both
  // the engine and the CBC, always lock the CBC first. Otherwise
  // there is a potential deadlock!

  virtual void notify(T* &data, const ObserverEvent &event);

  // resetDataAndStrides is a utility function used by various
  // constructors.  It sets the strides based on the compression
  // status, and sets the pointer to the data.  Once the Compressible
  // Brick has been created, the strides and data pointer are updated
  // by the notify function.

  // NOTE: the cblock must be locked before this function is called.
  
  void resetDataAndStrides();

  // Helper function used in initialization
  
  void init();
};


/**
 * Engine<T,Dim,CompressibleBrickView> 
 * (aka CompressibleBrickView-Engine)
 *
 * A CompressibleBrickView-Engine is an Engine that manages a view of
 * a CompressibleBrick-Engine.  See Engine<Dim,T,CompressibleBrick>
 * for details.
 *
 * Template Parameters:
 *  - T, the type of object stored. 
 *  - Dim, the dimension of the CompressibleBrick-View.
 * 
 * The Domain of this engine is an Interval<Dim>, which is a tensor
 * product of Dim 1-D intervals. For CompressibleBrickView-Engines,
 * these intervals will all be 0-based (i.e. [0..N0]x[0..N1] etc.).
 *  
 * Note that this is NOT the domain of the underlying data storage,
 * but rather it is the domain as presented to the outside world.
 */

template <int Dim, class T>
class Engine<Dim,T,CompressibleBrickView> 
  : public Pooma::BrickViewBase<Dim>, public Observer<T*>
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Engine<Dim,T,CompressibleBrickView>           This_t;
  typedef Engine<Dim,T,CompressibleBrickView>           Engine_t;
  typedef Pooma::BrickViewBase<Dim>                     Base_t;
  typedef typename Base_t::Domain_t                     Domain_t;
  typedef DomainLayout<Dim>                             Layout_t;
  typedef T                                             Element_t;
  typedef T&                                            ElementRef_t;
  typedef CompressibleBrickView                         Tag_t;

  enum { dimensions    = Dim   };
  enum { hasDataObject = true  };
  enum { dynamic       = false };
  enum { zeroBased     = true  };
  enum { multiPatch    = false };

  //============================================================
  // Constructors and Factory Methods
  //============================================================

  // Default constructor. Creates a CompressibleBrickView-Engine with no
  // data and an "empty" domain.

  Engine() : data0_m(0) { }

  // Copy constructor performs a SHALLOW copy:

  Engine(const Engine_t &model);
  Engine(const Engine_t &model, const EngineConstructTag &);

  // Subsetting Constructors.
  // A CompressibleBrickView-Engine is a strided, brick-shaped view of
  // a CompressibleBrick-Engine. Thus we write constructors to build
  // CompressibleBrickViews from CompressibleBrick-Engines and all
  // types of brick-shaped Domains.

  // Build a CompressibleBrickView from a CompressibleBrick and a non-slice  
  // domain like an Interval<Dim> or Range<Dim>.

  template <class DT>
  Engine(const Engine<Dim,T,CompressibleBrick> &e, const Domain<Dim, DT> &dom)
  : Base_t(e, dom.unwrap()), cblock_m(e.cblock()),
    entire_m(e.domain() == dom.unwrap())
  {        
    init();
  }  

  // Build a CompressibleBrickView from a CompressibleBrick and a Node.

  template <class Domain>
  Engine(const Engine<Dim,T,CompressibleBrick> &e, const Node<Domain> &node)
    : Base_t(e, node.domain()), cblock_m(e.cblock()),
      entire_m(e.domain() == node.domain())
  {
    init();
  }  

  // Build a CompressibleBrickView from a CompressibleBrick and an INode.

  Engine(const Engine<Dim,T,CompressibleBrick> &e, const INode<Dim> &inode)
  : Base_t(e, inode.domain()), cblock_m(e.cblock()),
    entire_m(e.domain() == inode.domain())
  {
    init();
  }  
  
  // Build a CompressibleBrickView from CompressibleBrick and a 
  // SliceDomain<DT>.
  
  template <class DT, int Dim2>
  Engine(const Engine<Dim2,T,CompressibleBrick> &e, const SliceDomain<DT> &dom)
  : Base_t(e, dom.unwrap()), cblock_m(e.cblock()),
    entire_m(e.domain() == dom.totalDomain())
  {
    init();
  }

  // Build a CompressibleBrickView from another CompressibleBrickView and
  // a domain like an Interval<Dim> or Range<Dim>.

  template <class DT>
  Engine(const This_t &e, const Domain<Dim, DT> &dom)
  : Base_t(e, dom.unwrap()), cblock_m(e.cblock()),
    entire_m(e.entire_m && e.domain() == dom.unwrap())
  {
    init();
  }

  // Build a CompressibleBrickView from another CompressibleBrickView
  // and an INode.

  Engine(const This_t &e, const INode<Dim> &inode)
  : Base_t(e, inode.domain()), cblock_m(e.cblock()),
    entire_m(e.entire_m && e.domain() == inode.domain())
  {
    init();
  }
  
  // Build a CompressibleBrickView from another CompressibleBrickView and a 
  // SliceDomain<DT>

  template <int OrigDim, class DT>
  Engine(const Engine<OrigDim,T,CompressibleBrickView> &e, 
	 const SliceDomain<DT> &dom)
  : Base_t(e, dom.unwrap()), cblock_m(e.cblock()),
    entire_m(e.entire_m && e.domain() == dom.totalDomain())
  {
    init();
  }

  //============================================================
  // Destructor
  //============================================================

  // Stop observing the block controller on destruct.

  ~Engine();

  //============================================================
  // Assignment operators
  //============================================================

  // Assigment is SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);

  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  // Element access.
  // These do not lock the underlying block. The run-time system will
  // ensure that data-parallel statements use these safely. If the
  // user wishes to access individual elements in his code, the onus
  // is on him or her to ensure synchronization.

  // Read-write access via ints.
  // These functions are slow, since they have to check if the
  // block is compressed.

  ElementRef_t operator()(int) const;
  ElementRef_t operator()(int, int) const;
  ElementRef_t operator()(int, int, int) const;
  ElementRef_t operator()(int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int, int, int) const;

  // Read-only access.
  // Guaranteed to be fast since there's no if-test.

  Element_t read(int) const;
  Element_t read(int, int) const;
  Element_t read(int, int, int) const;
  Element_t read(int, int, int, int) const;
  Element_t read(int, int, int, int, int) const;
  Element_t read(int, int, int, int, int, int) const;
  Element_t read(int, int, int, int, int, int, int) const;

  // Element access via Loc.

  ElementRef_t operator()(const Loc<Dim> &) const;
  Element_t read(const Loc<Dim> &) const;
  
  //---------------------------------------------------------------------------
  // Return the domain and base domain.

  inline const Domain_t &domain() const
  {
    return this->domain_m;
  }

  // Return the block controller (uncompressed). 
  // See comments in CompressibleBrick above.
  
  DataBlockPtr<T> dataBlock() const { return cblock_m.view(); }

  // Provide access to the data object.

  inline 
  Pooma::DataObject_t *dataObject() const { return cblock_m.dataObject(); }

  //============================================================
  // Compressibility-related Accessors & Mutators
  //============================================================
  
  // Return the cblock.
  
  CompressibleBlock<T> cblock() const { return cblock_m; }
  
  // Check to see if this CompressibleBrickView is compressed.
  
  // See comments on thread safety in CompressibleBrick above.
  
  bool compressed() const { return cblock_m.compressed(); }

  // Read the compressed value for a CompressibleBrickView.

  T compressedRead() const;

  // Get a reference to the compressed value for a CompressibleBrickView.

  T& compressedReadWrite() const;

  // Determine whether or not this view is over the whole domain.
  
  bool compressedBrickIsWholeView() const { return entire_m; }

  // Return the number of compressed elements.
  
  long elementsCompressed() const;

private:

  // Mutex functions

  void lock()   const { mutex_m.lock(); }
  void unlock() const { mutex_m.unlock(); }

  // notify function for the Observer<T*> base class:
  // Compressible bricks observe the CompressibleBlockController
  // which notifies us when the data becomes compressed or uncompressed.
  // The notification comes with a pointer to the new data.

  // Note: The CBC is locked when this is called, so we don't have to
  // worry about contention for changind strides_m/data0_m, but we
  // do need to make sure no one tries to make a copy of these data
  // while they are being changed. Thus we lock our mutex.

  virtual void notify(T* &data, const ObserverEvent &event)
  {
    switch (event.event())
      {
      default:
      case CompressibleBlock<T>::notifyDestruct:
	// cblock has destructed. this should never happen
	// if the brick still exists.
	PAssert(false);
	break;
      case CompressibleBlock<T>::notifyUncompress: 
        lock();
        this->restoreStrides();
	data0_m = data + this->baseOffset();
	unlock();
	break;
      case CompressibleBlock<T>::notifyCompress:
        lock();
        this->zeroStrides();
	data0_m = data;
	unlock();
	break;
      }
  }

  // Helper function
  
  void init()
  {
    // resetDataAndStrides gets compression dependent data from
    // the CBC. We need to lock the CBC while getting this data,
    // and keep it locked until we are finished and have attached
    // as an observer.
    
    cblock_m.lock();

    resetDataAndStrides();
  
    PAssert(cblock_m.isControllerValidUnlocked());
    cblock_m.attach(this);

    cblock_m.unlock();
  }
  
  // resetDataAndStrides is a utility function used by various
  // constructors.  It sets the strides based on the compression
  // status, and sets the pointer to the data.  Once the Compressible
  // Brick has been created, the strides and data pointer are updated
  // by the notify function.

  // NOTE: the cblock must be locked before this function is called.

  void resetDataAndStrides()
  {
    if (cblock_m.compressed())
      {
        this->zeroStrides();
	data0_m = cblock_m.data();
      }
    else
      {
        this->restoreStrides();
	data0_m = cblock_m.data() + this->baseOffset();
      }
  }

  //============================================================
  // Data
  //============================================================

  // Smart-pointer to the controller.

  CompressibleBlock<T> cblock_m;

  // Plain pointer that points to the beginning of the actual data. 
  // If uncompressed, all index offsets are subtracted during
  // construction. If compressed, this simply points to the
  // compressed value.

  T *data0_m;
  
  // Flag that tells whether we're viewing the entire domain.
  
  bool entire_m;
  
  // Mutex protection for the CompressibleBricks.
  // This must be locked when changes are made via a notify
  // method, which can happen asynchronously, and when 
  // the items changed (strides_m and data0_m) are accessed.
  
  mutable Pooma::Mutex_t mutex_m;

};


/**
 * NewEngine<Engine,SubDomain>
 *
 * Several specializations of NewEngine for combinations of 
 * Engines and Domains that produce CompressibleBrickView-Engines.
 */

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrick>, Interval<Dim> >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrick>, Range<Dim> >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrick>,Node<Interval<Dim> > >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrick>,INode<Dim> >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrickView>, Interval<Dim> >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrickView>, Range<Dim> >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrickView>,
                 Node<Interval<Dim> > >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,CompressibleBrickView>,INode<Dim> >
{
  typedef Engine<Dim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,CompressibleBrick>,SliceInterval<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,CompressibleBrick>,SliceRange<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,CompressibleBrickView>, 
                 SliceInterval<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,CompressibleBrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,CompressibleBrickView>, 
                 SliceRange<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,CompressibleBrickView> Type_t;
};

/**
 * Traits class telling RefCountedBlockPointer that this class has
 * shallow semantics and a makeOwnCopy method.
 */

template <int Dim, class T>
struct ElementProperties<Engine<Dim, T, CompressibleBrick> > 
  : public MakeOwnCopyProperties<Engine<Dim, T, CompressibleBrick> >
{ };


///////////////////////////////////////////////////////////////////////////////
//
// Inline implementation of the functions for Engine<D,T,CompressibleBrick>
//
///////////////////////////////////////////////////////////////////////////////

// operator() and read() definitions

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(const Loc<Dim> &loc) const
{
  return data0_m[this->offsetC(loc)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1) const
{
  PAssert(Dim == 1);
  return data0_m[this->offsetC(i1)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1, int i2) const
{
  PAssert(Dim == 2);
  return data0_m[this->offsetC(i1,i2)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1, int i2, int i3) const
{
  PAssert(Dim == 3);
  return data0_m[this->offsetC(i1,i2,i3)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1, int i2, int i3, int i4) const
{
  PAssert(Dim == 4);
  return data0_m[this->offsetC(i1,i2,i3,i4)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1, int i2, int i3, int i4, int i5) const
{
  PAssert(Dim == 5);
  return data0_m[this->offsetC(i1,i2,i3,i4,i5)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1, int i2, int i3, int i4, int i5, int i6) const
{
  PAssert(Dim == 6);
  return data0_m[this->offsetC(i1,i2,i3,i4,i5,i6)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrick>::
read(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
{
  PAssert(Dim == 7);
  return data0_m[this->offsetC(i1,i2,i3,i4,i5,i6,i7)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(const Loc<Dim> &loc) const
{
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(loc)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1) const
{
  PAssert(Dim == 1);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1, int i2) const
{
  PAssert(Dim == 2);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1,i2)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1, int i2, int i3) const
{
  PAssert(Dim == 3);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1,i2,i3)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1, int i2, int i3, int i4) const
{
  PAssert(Dim == 4);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1,i2,i3,i4)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1, int i2, int i3, int i4, int i5) const
{
  PAssert(Dim == 5);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1,i2,i3,i4,i5)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1, int i2, int i3, int i4, int i5, int i6) const
{
  PAssert(Dim == 6);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1,i2,i3,i4,i5,i6)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrick>::
operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
{
  PAssert(Dim == 7);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offsetC(i1,i2,i3,i4,i5,i6,i7)];
}

template <int Dim, class T>
inline bool
Engine<Dim,T,CompressibleBrick>::
compressed() const
{
  PAssert(cblock_m.isControllerValidUnlocked());
  return cblock_m.compressed();
}

template <int Dim, class T>
inline T
Engine<Dim,T,CompressibleBrick>::
compressedRead() const 
{
  PAssert(cblock_m.isControllerValidUnlocked());
  PAssert(cblock_m.compressed());
  return *data0_m;
}

template <int Dim, class T>
inline T&
Engine<Dim,T,CompressibleBrick>::
compressedReadWrite() const
{
  PAssert(cblock_m.isControllerValidUnlocked());
  PAssert(cblock_m.compressed());
  return *data0_m;
}


//
// Free functions returning compressed status or doing compress/uncompress.
//

template <int Dim, class T>
inline bool compressed(const Engine<Dim, T, CompressibleBrick> &e)
{
  return e.compressed();
}

template <int Dim, class T>
inline long elementsCompressed(const Engine<Dim, T, CompressibleBrick> &e)
{
  return e.elementsCompressed();
}

template <int Dim, class T>
inline void compress(Engine<Dim, T, CompressibleBrick> &e)
{
  e.tryCompress();
}

template <int Dim, class T>
inline void uncompress(Engine<Dim, T, CompressibleBrick> &e)
{
  e.uncompress();
}

///////////////////////////////////////////////////////////////////////////////
//
// Inline implementation of the functions for Engine<D,T,CompressibleBrickView>
//
///////////////////////////////////////////////////////////////////////////////

// operator() and read() definitions

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(const Loc<Dim> &loc) const
{
  return data0_m[this->offset(loc)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1) const
{
  PAssert(Dim == 1);
  return data0_m[this->offset(i1)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1, int i2) const
{
  PAssert(Dim == 2);
  return data0_m[this->offset(i1,i2)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1, int i2, int i3) const
{
  PAssert(Dim == 3);
  return data0_m[this->offset(i1,i2,i3)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1, int i2, int i3, int i4) const
{
  PAssert(Dim == 4);
  return data0_m[this->offset(i1,i2,i3,i4)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1, int i2, int i3, int i4, int i5) const
{
  PAssert(Dim == 5);
  return data0_m[this->offset(i1,i2,i3,i4,i5)];
}

template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1, int i2, int i3, int i4, int i5, int i6) const
{
  PAssert(Dim == 6);
  return data0_m[this->offset(i1,i2,i3,i4,i5,i6)];
}


template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
read(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
{
  PAssert(Dim == 7);
  return data0_m[this->offset(i1,i2,i3,i4,i5,i6,i7)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(const Loc<Dim> &loc) const
{
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(loc)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1) const
{
  PAssert(Dim == 1);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1, int i2) const
{
  PAssert(Dim == 2);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1,i2)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1, int i2, int i3) const
{
  PAssert(Dim == 3);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1,i2,i3)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1, int i2, int i3, int i4) const
{
  PAssert(Dim == 4);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1,i2,i3,i4)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1, int i2, int i3, int i4, int i5) const
{
  PAssert(Dim == 5);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1,i2,i3,i4,i5)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1, int i2, int i3, int i4, int i5, int i6) const
{
  PAssert(Dim == 6);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1,i2,i3,i4,i5,i6)];
}

template <int Dim, class T>
inline T & Engine<Dim,T,CompressibleBrickView>::
operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
{
  PAssert(Dim == 7);
  if (cblock_m.compressed()) cblock_m.uncompress();
  return data0_m[this->offset(i1,i2,i3,i4,i5,i6,i7)];
}


template <int Dim, class T>
inline T Engine<Dim,T,CompressibleBrickView>::
compressedRead() const 
{
  PAssert(cblock_m.compressed());
  return *data0_m;
}

template <int Dim, class T>
inline T& Engine<Dim,T,CompressibleBrickView>::
compressedReadWrite() const 
{
  PAssert(cblock_m.compressed());
  return *data0_m;
}

//
// Free function returning compressed status.
//

template <int Dim, class T>
inline 
bool compressed(const Engine<Dim,T,CompressibleBrickView> &e)
{
  return e.compressed();
}

template <int Dim, class T>
inline 
long elementsCompressed(const Engine<Dim,T,CompressibleBrickView> &e)
{
  return e.elementsCompressed();
}

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

// Include .cpp file to get out-of-line functions.

#include "Engine/CompressibleBrick.cpp"

#endif // POOMA_ENGINE_COMPRESSIBLEBRICK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CompressibleBrick.h,v $   $Author: richard $
// $Revision: 1.75 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
