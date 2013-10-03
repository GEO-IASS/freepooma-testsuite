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
//   Brick                       - Brick-Engine specialization tag.
//   BrickView                   - BrickView-Engine specialization tag.
//   Engine<Dim,T,Brick>         - the "Brick-Engine" specialization.
//   Engine<Dim,T,BrickView>     - the "BrickView-Engine" specialization.
//   NewEngine<Engine,SubDomain> - specializations for BrickView-Engine.
//   ElementProperties<Engine
//     <Dim, T, Brick> >         - specialization for BrickEngine.
//   
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_BRICKENGINE_H
#define POOMA_ENGINE_BRICKENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 *   Brick & BrickView
 *    - tag classes used to select specializations of Engine
 *
 *   Engine<Dim,T,Brick>
 *    - an Engine that manages a contiguous, local, N-dimensional 
 *      brick of data. See Engine.h for the general template.
 *
 *   Engine<Dim,T,BrickView>
 *    - an Engine that manages a view into a Brick-Engine.
 *
 *   NewEngine
 *    - Specialized traits class for creating BrickView-Engines. 
 *      See Engine.h for the general version.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/BrickBase.h"
#include "Engine/Engine.h"
#include "Domain/Domain.h"
#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include "Domain/Touches.h"
#include "Engine/CompressibleBrick.h"
#include "Layout/DomainLayout.h"
#include "Layout/Node.h"
#include "Layout/INode.h"
#include "Utilities/DataBlockPtr.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

template <int Dim> class Range;

/**
 * These are tag classes used to select the "Brick" and "BrickView" 
 * specializations of the Engine class template.
 */

struct Brick {};

struct BrickView {};


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

template <int Dim, class T>
class Engine<Dim,T,BrickView>;

/**
 * Engine<Dim,T,Brick>  (aka Brick-Engine)
 *
 *  Engine<Dim,T,Brick> is an Engine that manages a contiguous, local,
 *  Dim-dimensional brick of data.
 *
 *  Template Parameters:
 *   - Dim, an integer for the dimension of the Brick.
 *   - T, the type of object stored.  Often a numeric type.
 *        However, the Engine makes no assumptions beyond that
 *        T have a copy constructor and assignment operator.
 *
 *  The Domain of this engine is an Interval<Dim>, which is a 
 *  tensor product of Dim 1-D intervals.
 *
 *  Subsetting Engine<Dim,T,Brick> returns an Engine<Dim,T,BrickView>.
 *  See below.
 */

template <int Dim, class T>
class Engine<Dim,T,Brick> : public Pooma::BrickBase<Dim>
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Engine<Dim,T,Brick>              This_t;
  typedef Engine<Dim,T,Brick>              Engine_t;
  typedef Pooma::BrickBase<Dim>            Base_t;
  typedef typename Base_t::Domain_t        Domain_t;
  typedef DomainLayout<Dim>                Layout_t;
  typedef T                                Element_t;
  typedef T&                               ElementRef_t;
  typedef Brick                            Tag_t;

  enum { brick         = true  };
  enum { dimensions    = Dim   };
  enum { hasDataObject = true  };
  enum { dynamic       = false };
  enum { zeroBased     = false };
  enum { multiPatch    = false };

  //============================================================
  // Constructors and Factory Methods
  //============================================================

  /// Default constructor. Creates a Brick-Engine with no data
  /// and an "empty" domain.  This is not really usable until it
  /// has been initialized (via operator=) to a new engine with an
  /// actual domain.

  Engine() { }

  //@{

  /// These constructors take an Interval<Dim> and create a new
  /// Brick-Engine with data of type T on this Domain. This is where
  /// storage gets allocated. The second version uses a model data
  /// element to initialize storage.
  
  explicit Engine(const Domain_t &domain);

  Engine(const Domain_t &domain, const T &elementModel);

  //@}

  /// You can build a brick from a layout as well.

  explicit Engine(const Layout_t &layout);

  /// This constructor takes a Node object, extracts the domain, and
  /// creates a new Brick-Engine with data of type T on this Domain. 
  
  explicit Engine(const Node<Domain_t> &node);

  /// This constructor takes a domain and a pointer to external memory,
  /// and constructs a Brick-Engine wrapper for this external
  /// memory. Used for layering Pooma II arrays on top of someone elses
  /// data.

  Engine(T * foreignData, const Domain_t &domain);
  
  /// Copy constructor performs a SHALLOW copy.
  /// But NOTE: the layouts will NOT be shared.

  Engine(const Engine_t &model);

  // Subsetting Constructors.
  //   There are none - you cannot create a Brick-Engine by taking
  //   a "view" of another Engine.  Brick-Engines, by definition,
  //   view all of the data. See Engine<D,T,BrickView> below.

  //============================================================
  // Destructor
  //============================================================

  /// All pointer members are smart pointers, so this is trivial.

  ~Engine() {}

  //============================================================
  // Assignment operators
  //============================================================

  /// Assigment is SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);

  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  //@{

  /// Element access via Loc.

  Element_t read(const Loc<Dim> &loc) const
  {
    return data_m[this->offset(loc)];
  }
  ElementRef_t operator()(const Loc<Dim> &loc) const
  {
    return data_m[this->offset(loc)];
  }

  //@}

  //@{

  /// Element access via ints for speed.

  Element_t read(int i1) const
  {
    CTAssert(Dim == 1);
    return data_m[this->offset(i1)];
  }
  Element_t read(int i1, int i2) const
  {
    CTAssert(Dim == 2);
    return data_m[this->offset(i1,i2)];
  }
  Element_t read(int i1, int i2, int i3) const
  {
    CTAssert(Dim == 3);
    return data_m[this->offset(i1,i2,i3)];
  }
  Element_t read(int i1, int i2, int i3, int i4) const
  {
    CTAssert(Dim == 4);
    return data_m[this->offset(i1,i2,i3,i4)];
  }
  Element_t read(int i1, int i2, int i3, int i4, int i5) const
  {
    CTAssert(Dim == 5);
    return data_m[this->offset(i1,i2,i3,i4,i5)];
  }
  Element_t read(int i1, int i2, int i3, int i4, int i5, int i6) const
  {
    CTAssert(Dim == 6);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6)];
  }
  Element_t read(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
  {
    CTAssert(Dim == 7);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6,i7)];
  }

  ElementRef_t operator()(int i1) const
  {
    CTAssert(Dim == 1);
    return data_m[this->offset(i1)];
  }
  ElementRef_t operator()(int i1, int i2) const
  {
    CTAssert(Dim == 2);
    return data_m[this->offset(i1,i2)];
  }
  ElementRef_t operator()(int i1, int i2, int i3) const
  {
    CTAssert(Dim == 3);
    return data_m[this->offset(i1,i2,i3)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4) const
  {
    CTAssert(Dim == 4);
    return data_m[this->offset(i1,i2,i3,i4)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5) const
  {
    CTAssert(Dim == 5);
    return data_m[this->offset(i1,i2,i3,i4,i5)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5, int i6) const
  {
    CTAssert(Dim == 6);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
  {
    CTAssert(Dim == 7);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6,i7)];
  }

  //@}

  /// Get a private copy of data viewed by this Engine.

  Engine_t &makeOwnCopy();

  /// Provide access to the data object. 

  inline
  Pooma::DataObject_t *dataObject() const { return dataBlock_m.dataObject(); }

  //@{

  /// Return access to our internal data block.  
  /// This is ref-counted, so a copy is fine.  But you should really 
  /// know what you're doing if you call this method.
  /// (cbrick's value version is const - why isn't this???)
  
  DataBlockPtr<T> dataBlock()               { return dataBlock_m; }
  const DataBlockPtr<T> & dataBlock() const { return dataBlock_m; }

  //@}

  /// Return whether the block controlled by this engine is shared.
  /// (Is this used??? CompressibleBrick didn't have one!!!)
  
  bool isShared() const { return dataBlock_m.isValid() && dataBlock_m.count() > 1; }

private:

  //============================================================
  // Private data
  //============================================================

  // The structural data (domain, strides, etc.) are in BrickBase.
  
  /// Smart-pointer to Block-controller that manages the data
  /// and the Smarts DataObject. 

  DataBlockPtr<T> dataBlock_m;

  T *data_m;
};



/**
 * Engine<T,Dim,BrickView> (aka BrickView-Engine)
 *
 *  A BrickView-Engine is an Engine that manages a view of a
 *  Brick-Engine.  See Engine<Dim,T,Brick> for details.
 *
 *  Template Parameters:
 *   - T, the type of object stored. 
 *   - Dim, an integer for the dimension of the BrickView.
 *
 *  The Domain of this engine is an Interval<Dim>, which is a tensor
 *  product of Dim 1-D intervals. For BrickView-Engines, these
 *  intervals will all be 0-based (i.e. [0..N0]x[0..N1] etc.).
 *  
 *  Note that this is NOT the domain of the underlying data storage,
 *  but rather it is the domain as presented to the outside world.
 */

template <int Dim, class T>
class Engine<Dim,T,BrickView> 
 : public Pooma::BrickViewBase<Dim>
{
public:
  
  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Engine<Dim,T,BrickView>                   This_t;
  typedef Engine<Dim,T,BrickView>                   Engine_t;
  typedef Pooma::BrickViewBase<Dim>                 Base_t;
  typedef typename Base_t::Domain_t                 Domain_t;
  typedef DomainLayout<Dim>                         Layout_t;
  typedef T                                         Element_t;
  typedef T&                                        ElementRef_t;
  typedef BrickView                                 Tag_t;

  enum { dimensions    = Dim   };
  enum { hasDataObject = true  };
  enum { dynamic       = false };
  enum { zeroBased     = true  };
  enum { multiPatch    = false };
  
  //============================================================
  // Constructors
  //============================================================

  // BrickView-Engine is fundamentally a view - it never owns
  // its data, and thus there are no constructors that create
  // a BrickView-Engine directly from a domain.

  // Default constructor is required for containers.

  Engine();

  // Copy constructor performs a SHALLOW copy:

  Engine(const Engine_t &);
  Engine(const Engine_t &, const EngineConstructTag &);

  // Subsetting Constructors.
  // A BrickView-Engine is a strided, brick-shaped view of a
  // Brick-Engine. Thus we write constructors to build BrickViews from
  // Brick-Engines and all types of brick-shaped Domains.

  // Build a BrickView from Brick and a domain like an Interval<Dim> 
  // or Range<Dim>.

  // Why is this templated on ETag???
  
  template <class ETag, class DT>
  Engine(const Engine<Dim,T,ETag> &e, const Domain<Dim, DT> &dom)
  : Base_t(e, dom.unwrap()), dataBlock_m(e.dataBlock(), e.offset(dom.unwrap())),
    data_m(dataBlock_m.currentPointer())
  { 
    // The engine's data pointer should be at the beginning.
    PAssert(e.dataBlock().isAtBeginning());
  }

  // Build a BrickView from Brick and a SliceRange<Dim2,Dim>.

  template<int Dim2>
  Engine(const Engine<Dim2,T,Brick> &e, const SliceRange<Dim2,Dim> &dom)
    : Base_t(e, dom), dataBlock_m(e.dataBlock(), e.offset(dom.totalDomain())),
    data_m(dataBlock_m.currentPointer())
  {    
    // The engine's data pointer should be at the beginning.
    PAssert(e.dataBlock().isAtBeginning());
  }
  
  template<int Dim2>
  Engine(const Engine<Dim2,T,Brick> &e, const SliceInterval<Dim2,Dim> &dom)
    : Base_t(e, dom), dataBlock_m(e.dataBlock(), e.offset(dom.totalDomain())),
    data_m(dataBlock_m.currentPointer())
  {    
    // The engine's data pointer should be at the beginning.
    PAssert(e.dataBlock().isAtBeginning());
  }
  
  // what is this #if???
#if 0
  template <int ODim, class ETag>
  Engine(const Engine<ODim,T,ETag> &e, const SliceRange<ODim,Dim> &dom)
    : Base_t(e, dom), dataBlock_m(e.dataBlock(), e.offset(dom.totalDomain())),
    data_m(dataBlock_m.currentPointer())
  {    
    // The engine's data pointer should be at the beginning.
    PAssert(e.dataBlock().isAtBeginning());
  }
  // Build a BrickView from Brick and a SliceInterval<Dim2,Dim>.

  template <int ODim, class ETag>
  Engine(const Engine<ODim,T,ETag> &e, const SliceInterval<ODim,Dim> &dom)
    : Base_t(e, dom), dataBlock_m(e.dataBlock(), e.offset(dom.totalDomain())),
    data_m(dataBlock_m.currentPointer())
  {    
    // The engine's data pointer should be at the beginning.
    PAssert(e.dataBlock().isAtBeginning());
  }
#endif


  // Build a BrickView from another BrickView and a domain like
  // an Interval<Dim> or Range<Dim>.

  template <class DT>
  Engine(const This_t &e, const Domain<Dim, DT> &d)
    : Base_t(e, d.unwrap()), dataBlock_m(e.dataBlock(), e.offset(d.unwrap())),
    data_m(dataBlock_m.currentPointer())
  { }

  // Build a BrickView from another BrickView and an INode.

  Engine(const This_t &e, const INode<Dim> &inode)
    : Base_t(e,inode.domain()), dataBlock_m(e.dataBlock(), e.offset(inode.domain())),
    data_m(dataBlock_m.currentPointer())
  { }
  
  // Build a BrickView from another BrickView and a SliceRange<ODim,Dim>.

  template <int ODim>
  Engine(const Engine<ODim,T,BrickView> &e, 
         const SliceRange<ODim,Dim> &dom)
    : Base_t(e, dom), dataBlock_m(e.dataBlock(), e.offset(dom.totalDomain())),
    data_m(dataBlock_m.currentPointer())
  { }

  // Build a BrickView from another BrickView and a SliceRange<ODim,Dim>.

  template <int ODim>
  Engine(const Engine<ODim,T,BrickView> &e, 
    const SliceInterval<ODim,Dim> &dom)
    : Base_t(e, dom), dataBlock_m(e.dataBlock(), e.offset(dom.totalDomain())),
    data_m(dataBlock_m.currentPointer())
  { }

  // Build a BrickView-Engine from a compressible brick.
  
  explicit Engine(const Engine<Dim,T,CompressibleBrick> &);
  explicit Engine(const Engine<Dim,T,CompressibleBrickView> &);
  
  //============================================================
  // Destructor
  //============================================================

  ~Engine() {}

  //============================================================
  // Assignment operators
  //============================================================

  // Assigment is SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);

  //============================================================
  // Accessor functions:
  //============================================================

  //@{

  /// Element access via Loc.

  Element_t read(const Loc<Dim> &loc) const
  {
    return data_m[this->offset(loc)];
  }
  ElementRef_t operator()(const Loc<Dim> &loc) const
  {
    return data_m[this->offset(loc)];
  }

  //@}

  //@{

  /// Element access via ints for speed.

  Element_t read(int i1) const
  {
    CTAssert(Dim == 1);
    return data_m[this->offset(i1)];
  }
  Element_t read(int i1, int i2) const
  {
    CTAssert(Dim == 2);
    return data_m[this->offset(i1,i2)];
  }
  Element_t read(int i1, int i2, int i3) const
  {
    CTAssert(Dim == 3);
    return data_m[this->offset(i1,i2,i3)];
  }
  Element_t read(int i1, int i2, int i3, int i4) const
  {
    CTAssert(Dim == 4);
    return data_m[this->offset(i1,i2,i3,i4)];
  }
  Element_t read(int i1, int i2, int i3, int i4, int i5) const
  {
    CTAssert(Dim == 5);
    return data_m[this->offset(i1,i2,i3,i4,i5)];
  }
  Element_t read(int i1, int i2, int i3, int i4, int i5, int i6) const
  {
    CTAssert(Dim == 6);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6)];
  }
  Element_t read(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
  {
    CTAssert(Dim == 7);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6,i7)];
  }

  ElementRef_t operator()(int i1) const
  {
    CTAssert(Dim == 1);
    return data_m[this->offset(i1)];
  }
  ElementRef_t operator()(int i1, int i2) const
  {
    CTAssert(Dim == 2);
    return data_m[this->offset(i1,i2)];
  }
  ElementRef_t operator()(int i1, int i2, int i3) const
  {
    CTAssert(Dim == 3);
    return data_m[this->offset(i1,i2,i3)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4) const
  {
    CTAssert(Dim == 4);
    return data_m[this->offset(i1,i2,i3,i4)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5) const
  {
    CTAssert(Dim == 5);
    return data_m[this->offset(i1,i2,i3,i4,i5)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5, int i6) const
  {
    CTAssert(Dim == 6);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6)];
  }
  ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
  {
    CTAssert(Dim == 7);
    return data_m[this->offset(i1,i2,i3,i4,i5,i6,i7)];
  }

  //@}

  // Return the DataBlockPtr. See comments in BrickEngine above.
  
  DataBlockPtr<T> dataBlock()              { return dataBlock_m; }
  const DataBlockPtr<T> &dataBlock() const { return dataBlock_m; }

  // Provide access to the data object.  This should really be
  // a reference instead of a pointer.
  
  inline 
  Pooma::DataObject_t *dataObject() const { return dataBlock_m.dataObject(); }

private:

  //============================================================
  // Data
  //============================================================

  // Structural data is in BrickViewBase.

  // We just have a pointer to our data block.
  
  DataBlockPtr<T> dataBlock_m;

  T *data_m;
};



/**
 * NewEngine<Engine,SubDomain>
 *
 * Several specializations of NewEngine for combinations of 
 * Engines and Domains that produce BrickView-Engines.
 */

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,Brick>, Interval<Dim> >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,Brick>, Range<Dim> >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,Brick>, Node<Interval<Dim> > >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,Brick>, INode<Dim> >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,BrickView>, Interval<Dim> >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,BrickView>, Range<Dim> >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,BrickView>, 
                 Node<Interval<Dim> > >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim,T,BrickView>, INode<Dim> >
{
  typedef Engine<Dim,T,BrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,Brick>, SliceInterval<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,BrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,Brick>, SliceRange<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,BrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,BrickView>, 
                 SliceInterval<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,BrickView> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,BrickView>, 
                 SliceRange<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,BrickView> Type_t;
};



/**
 * NewEngineDomain<Engine,SubDomain>
 *
 * Several specializations of NewEngineDomain for combinations of 
 * Brick/BrickView-Engines and Node/INode..
 */

template <int Dim, class T>
struct NewEngineDomain<Engine<Dim,T,Brick>, Node<Interval<Dim> > >
{
  typedef Interval<Dim> Type_t;
  typedef const Interval<Dim> &Return_t;
  static inline
  Return_t apply(const Engine<Dim,T,Brick> &,
		 const Node<Interval<Dim> > &node)
  {
    return node.domain();
  }
};

template <int Dim, class T>
struct NewEngineDomain<Engine<Dim,T,Brick>, INode<Dim> >
{
  typedef Interval<Dim> Type_t;
  typedef const Interval<Dim> &Return_t;
  static inline
  Return_t apply(const Engine<Dim,T,Brick> &,
		 const INode<Dim> &inode)
  {
    return inode.domain();
  }
};

template <int Dim, class T>
struct NewEngineDomain<
  Engine<Dim,T,BrickView>, 
  Node<Interval<Dim> > >
{
  typedef Interval<Dim> Type_t;
  typedef const Interval<Dim> &Return_t;
  static inline
  Return_t apply(const Engine<Dim,T,BrickView> &,
		 const Node<Interval<Dim> > &node)
  {
    return node.domain();
  }
};

template <int Dim, class T>
struct NewEngineDomain<Engine<Dim,T,BrickView>, INode<Dim> >
{
  typedef Interval<Dim> Type_t;
  typedef const Interval<Dim> &Return_t;
  static inline
  Return_t apply(const Engine<Dim,T,BrickView> &,
		 const INode<Dim> &inode)
  {
    return inode.domain();
  }
};



/**
 * Traits class telling RefCountedBlockPointer that this class has
 * shallow semantics and a makeOwnCopy method.
 */

template <int Dim, class T>
struct ElementProperties<Engine<Dim, T, Brick> > 
  : public MakeOwnCopyProperties<Engine<Dim, T, Brick> >
{ };



// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

// Include .cpp file to get out-of-line functions.

#include "Engine/BrickEngine.cpp"

#endif // POOMA_ENGINE_BRICKENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BrickEngine.h,v $   $Author: richi $
// $Revision: 1.136 $   $Date: 2004/11/29 16:46:29 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
