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
//   LayoutBase<Dim>
//   LayoutBaseData<Dim>
//   LayoutBaseView<Dim, Dim2>
//   LayoutBaseViewData<Dim, Dim2>
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_LAYOUTBASE_H
#define POOMA_LAYOUT_LAYOUTBASE_H

/** @file
 * @ingroup Layout
 * @brief
 * LayoutBase<Dim> and related classes providing domain access.
 *
 *   LayoutBase<Dim>
 *     - Layout class that breaks Dim-dimensional domain into equal
 *       sized sub-domains arranged in a Dim-dimensional grid.
 *   LayoutBaseView<Dim, Dim2>
 *     - view of a LayoutBase
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Layout/MultiPatchLayoutTraits.h"
#include "Layout/INode.h"
#include "Layout/TouchesConstruct.h"
#include "Layout/GuardLayers.h"
#include "Domain/Interval.h"
#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/NewDomain.h"
#include "Domain/SliceRange.h"
#include "Partition/UniformGridPartition.h"
#include "Partition/ContextMapper.h"
#include "Utilities/DerefIterator.h"
#include "Utilities/ViewIndexer.h"
#include "Utilities/Observable.h"
#include "Utilities/Observer.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/RefCounted.h"
#include "Utilities/Unique.h"
#include "Utilities/PAssert.h"

#include <vector>
#include <iosfwd>


/**
 * Tag class specifying domain replication on all nodes,
 * implying a use of LocalMapper for mapping patches to the
 * only context.
 */

struct ReplicatedTag {};

/**
 * Tag class specifying domain distribution between all nodes,
 * implying a use of DistributedMapper for mapping patches
 * to contexts.
 */

struct DistributedTag {};


// Open Pooma namespace:
// namespace Pooma {

/**
 * Data container for LayoutBase class.
 */

template <int Dim>
class LayoutBaseData 
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef Interval<Dim>           Domain_t;
  typedef Interval<Dim>           BaseDomain_t;
  typedef int                     Context_t;
  typedef Unique::Value_t         ID_t;
  typedef Node<Domain_t>          Value_t;
  typedef std::vector<Value_t *>  List_t;           // for convenience
  typedef GuardLayers<Dim>        GuardLayers_t;    // for convenience
    
  // Guard-cell filling information

  struct GCFillInfo
  {
    GCFillInfo(const Domain_t &dom, int ownedID, int guardID, int face=-1)
    : domain_m(dom), ownedID_m(ownedID), guardID_m(guardID), face_m(face) { }
    
    // Get a CW warning about this not having a default constructor
    // when we instantiate the vector<GCFillInfo> below. This never
    // gets called, so it seems bogus, but necessary.

    GCFillInfo() { PInsist(0,"Shouldn't get here!"); }
    
    Domain_t domain_m;    // guard layer domain
    int      ownedID_m;   // node ID for which domain_m is owned 
    int      guardID_m;   // node ID for which domain_m is in the guards
    int      face_m;      // destination face of the guard layer (or -1, if unknown)
    
    Domain_t & domain() { return domain_m;}
    int & ownedID() { return ownedID_m;}
    int & guardID() { return guardID_m;}

  };
  
  typedef GCFillInfo GCFillInfo_t;
  
  typedef typename std::vector<GCFillInfo>::const_iterator FillIterator_t;

  // constructor

  LayoutBaseData()
    :
    ID_m(Unique::get()),
    domain_m(Interval<Dim>()), 
    innerdomain_m(Interval<Dim>()),
    hasInternalGuards_m(false), 
    hasExternalGuards_m(false), 
    internalGuards_m(0),
    externalGuards_m(0)
  {
  }

  LayoutBaseData(bool hasIG, bool hasEG,
		 GuardLayers_t eg, GuardLayers_t ig,
		 Domain_t d, Domain_t id)
    :
    ID_m(Unique::get()),
    domain_m(d), 
    innerdomain_m(id),
    hasInternalGuards_m(hasIG), 
    hasExternalGuards_m(hasEG), 
    internalGuards_m(ig),
    externalGuards_m(eg) 
  {
  }

  // destructor

  ~LayoutBaseData()
  {
    // All the destruction is done in the derived class destructor. 
  }

  /// Shortcut for allocatedDomain(int).

  inline const Domain_t & domain(int i) const
  {
    PAssert(i >= 0 && i < all_m.size());
    return all_m[i]->allocated();
  }

  /// The domain for patch i without internal guards.

  inline const Domain_t & ownedDomain(int i) const
  {
    PAssert(i >= 0 && i < all_m.size());
    return all_m[i]->domain();
  }

  /// The domain for patch i with internal guards.

  inline const Domain_t & allocatedDomain(int i) const
  {
    PAssert(i >= 0 && i < all_m.size());
    return all_m[i]->allocated();
  }

  inline const GuardLayers_t& internalGuards() const
  {
    return internalGuards_m;
  }

  inline const GuardLayers_t& externalGuards() const
  {
    return externalGuards_m;
  }

  inline List_t &nodeListGlobal() 
  {
    return all_m;
  }

  inline List_t &nodeListLocal()
  {
    return local_m;
  }

  inline List_t &nodeListRemote()
  {
    return remote_m;
  }

  /// Check if we've been initialized. Used in touches.

  inline bool initialized() const { return all_m.size() > 0; }
  
  /// first is the inner domain starting point

  inline int first(int d) const { return firsti_m[d]; }

  /// firsts is the domain (external guard region included) starting point
   
  inline int firsts(int d) const { return firste_m[d]; }

  /// number of blocks along each axis. 

  inline const Loc<Dim>& blocks() const  { return blocks_m; }

  ///@name Guard-cell related functions.
  /// Iterators into the fill list. These are MultiPatch's interface to
  /// the information needed to fill guard cells, which is cached here
  /// in the layout.
  //@{
   
  FillIterator_t beginFillList() const
  {
    return gcFillList_m.begin();
  }
   
  FillIterator_t endFillList() const
  {
    return gcFillList_m.end();
  }

  //@}

  //============================================================
  // Data
  //============================================================

  /// Our ID value, which is simply a unique value.

  ID_t ID_m;

  /// The global domain of this LayoutBase, including external
  /// guards.

  Domain_t domain_m;
  
  /// The global domain, excluding external guards.
  
  Domain_t innerdomain_m;

  /// The list of all, local, and remote subdomains.

  List_t all_m;
  List_t local_m;
  List_t remote_m;

  /// Do we have internal guards?
  
  bool hasInternalGuards_m;
  
  /// Do we have external guards?
  
  bool hasExternalGuards_m;
  
  /// The internal guard cell sizes:
  
  GuardLayers_t internalGuards_m;
  
  /// The external guard cell sizes:

  GuardLayers_t externalGuards_m;
  
  /// Cached guard-cell filling info.
  
  std::vector<GCFillInfo> gcFillList_m;

  /// The initial points for the domain, including external guards.
  
  int firste_m[Dim];

  /// the initial points for the inner domain, excluding external guards.
  int firsti_m[Dim];

  /// for UGL and GL, the number of blocks along each axis

  Loc<Dim> blocks_m;
};


/**
 * Base class for all layouts.
 */

template <int Dim, class LBD>
class LayoutBase
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.
  
  typedef LayoutBaseData<Dim>                  LayoutData_t; // for convenience
  typedef typename LayoutData_t::Domain_t      Domain_t;
  typedef typename LayoutData_t::BaseDomain_t  BaseDomain_t;
  typedef typename LayoutData_t::Context_t     Context_t;
  typedef typename LayoutData_t::ID_t          ID_t;
  typedef typename LayoutData_t::Value_t       Value_t;
  typedef typename LayoutData_t::List_t        List_t; // for convenience
  typedef LayoutBase<Dim,LBD>                  This_t; // for convenience
  typedef Observable<This_t>                   Observable_t;

  /// Iterator through nodes. Basically the same as the vector iterator
  /// except it dereferences automatically.
  
  typedef DerefIterator<Value_t>               iterator;
  typedef ConstDerefIterator<Value_t>          const_iterator;
  
  /// Iterator through guard-cell-fill requests. 
  
  typedef typename LayoutData_t::GCFillInfo_t GCFillInfo_t;
  
  typedef typename 
  std::vector<GCFillInfo_t>::const_iterator FillIterator_t;
  
  typedef GuardLayers<Dim> GuardLayers_t;

  enum { supportsGuards = true };
  
  //-------------------------------------------------------------------
  // Constructors
  //-------------------------------------------------------------------

  LayoutBase(LBD * Ldata)
    : pdata_m(Ldata)
  {
  }

   LayoutBase(RefCountedPtr<LBD> pdata)
    : pdata_m(pdata)
  {
  }

  //-------------------------------------------------------------------
  // Destructor
  //-------------------------------------------------------------------

  ~LayoutBase()
  {
  }


  //-------------------------------------------------------------------
  // Accessors
  //-------------------------------------------------------------------

  inline ID_t ID() const
  {
    return pdata_m->ID_m;
  }

  inline ID_t baseID() const
  {
    return pdata_m->ID_m;
  }

  /// Return whether or not this layout has been initialized.
  
  inline bool initialized() const
  {
    return (sizeGlobal() > 0);
  }
 
  /// Domain translation utility. 
  
  template <class DT>
  BaseDomain_t &localToBase(const Domain<Dim, DT> &dlocal, 
    BaseDomain_t &base) const 
  {
    return pdata_m->indexer_m.localToBase(dlocal,base);
  }

  ///@name Domain accessors
  //@{

  /// The global domain, including external guards.
  inline const Domain_t &domain() const
  {
    return pdata_m->domain_m;
  }

  /// The global domain, excluding external guards.
  inline const Domain_t &innerDomain() const
  {
    return pdata_m->innerdomain_m;
  }

  /// The global domain, including external guards.
  inline const Domain_t &baseDomain() const
  {
    return pdata_m->domain_m;
  }

  /// The global domain, including internal guards, i'th patch
  inline const Domain_t &domain(int i) const
  {
    return pdata_m->domain(i);
  }

  /// The global domain, excluding internal guards, i'th patch
  inline const Domain_t &ownedDomain(int i) const
  {
    return pdata_m->ownedDomain(i);
  }

  /// The global domain, including internal guards, i'th patch
  inline const Domain_t &allocatedDomain(int i) const
  {
    return pdata_m->allocatedDomain(i);
  }

  //@}

  ///@name Node list accessors.
  //@{

  inline const List_t & nodeListGlobal() const
  {
    return pdata_m->nodeListGlobal();
  }
  
  inline const List_t &nodeListLocal() const
  {
    return pdata_m->nodeListLocal();
  }

  inline const List_t &nodeListRemote() const
  {
    return pdata_m->nodeListRemote();
  }

  //@}

  ///@name Guard layer access.
  //@{

  /// Internal (between cells) guard layers.
  inline GuardLayers_t internalGuards() const
  {
    return pdata_m->internalGuards();
  }

  /// External (beyond physical domain) guard layers.
  inline GuardLayers_t externalGuards() const
  {
    return pdata_m->externalGuards();
  }

  //@}

  /// d'th component of the lower left of the inner domain.
  inline int first(int d) const { return pdata_m->first(d); }


  inline Loc<Dim>  blocks() const { return pdata_m->blocks(); }

  /// Return the domain of a local patch

  inline const Domain_t patchDomain(int lid) const
  {
    return nodeListLocal()[lid]->domain();
  }

  /// Convert a local patch ID to a global patch ID

  inline int localToGlobalPatchID(int lid) const
  {
    return nodeListLocal()[lid]->globalID();
  }

  //-----------------------------------------------------------------------------
  //
  /// @name globalID accessors
  /// Return the globalID of the node containing the given element,
  /// expressed either as a single Loc or as a set of int's.
  //
  //-----------------------------------------------------------------------------
  //@{
  // These just bounce the question to the LayoutData object.

  int globalID(const Loc<Dim> &loc) const
  { return pdata_m->globalID(loc); }

  int globalID(int i0) const
  { return pdata_m->globalID(i0); }

  int globalID(int i0, int i1) const
  { return pdata_m->globalID(i0,i1); }

  int globalID(int i0, int i1, int i2) const
  { return pdata_m->globalID(i0,i1,i2); }

  int globalID(int i0, int i1, int i2, int i3) const
  { return pdata_m->globalID(i0,i1,i2,i3); }

  int globalID(int i0, int i1, int i2, int i3, int i4) const
  { return pdata_m->globalID(i0,i1,i2,i3,i4); }

  int globalID(int i0, int i1, int i2, int i3, int i4, int i5) const
  { return pdata_m->globalID(i0,i1,i2,i3,i4,i5); }

  int globalID(int i0, int i1, int i2,int i3, int i4, int i5, int i6) const
  { return pdata_m->globalID(i0,i1,i2,i3,i4,i5,i6); }
  //@}

  //============================================================
  /// @name Partition methods
  /// Repartition the layout using a new Partitioner scheme.  The
  /// initial domain lists are cleared out, the partitioner is invoked,
  /// and then all the observers are notified.  This can only be done
  /// with a GridParition partitioner.
  //@{
  template <class Partitioner>
  bool repartition(const Partitioner &gp,const ContextMapper<Dim> &cmap)
  { 
    return pdata_m->repartition(gp,cmap); 
  }


  template <class Partitioner>
  bool repartition(const Partitioner &gp)
  { 
    typename Partitioner::DefaultMapper_t cmap(gp);
    return pdata_m->repartition(gp,cmap); 
  }
  //@}

  //============================================================
  // Compare 
  //============================================================

  /// Compare to another Layout.  The layouts are the same if:
  ///   -# They have the same base ID value.
  ///   -# They have the same base domain.
  
  template <class L>
  inline bool operator==(const L &layout) const
  {
    return (baseID()     == layout.baseID() && 
            baseDomain() == layout.baseDomain());
  }

  /// Compare for inequality.
  
  template <class L>
  inline bool operator!=(const L &layout) const 
  {
    return !(*this == layout);
  }

  //============================================================
  // Iterators

  /// @name Return begin and end iterators for the list of all subdomains
  //@{
  
  inline iterator beginGlobal() 
  { 
    return iterator(pdata_m->all_m.begin()); 
  }
  inline iterator endGlobal() 
  { 
    return iterator(pdata_m->all_m.end()); 
  }
  inline const_iterator beginGlobal() const 
  { 
    return const_iterator(pdata_m->all_m.begin()); 
  }
  inline const_iterator endGlobal() const 
  { 
    return const_iterator(pdata_m->all_m.end()); 
  }
  inline int sizeGlobal() const 
  { 
    return pdata_m->all_m.size(); 
  }
  //@}

  /// @name Return begin and end iterators for the list of local subdomains
  //@{
  inline iterator beginLocal() 
  { 
    return iterator(pdata_m->local_m.begin()); 
  }
  inline iterator endLocal() 
  { 
    return iterator(pdata_m->local_m.end());
  }
  inline const_iterator beginLocal() const
  { 
    return const_iterator(pdata_m->local_m.begin()); 
  }
  inline const_iterator endLocal() const 
  { 
    return const_iterator(pdata_m->local_m.end());
  }
  inline int sizeLocal() const 
  { 
    return pdata_m->local_m.size(); 
  }
  //@}

  /// @name Return begin and end iterators for the list of remote subdomains
  //@{
  
  inline iterator beginRemote() 
  { 
    return iterator(pdata_m->remote_m.begin()); 
  }
  inline iterator endRemote() 
  { 
    return iterator(pdata_m->remote_m.end()); 
  }
  inline const_iterator beginRemote() const 
  { 
    return const_iterator(pdata_m->remote_m.begin()); 
  }
  inline const_iterator endRemote() const 
  {
    return const_iterator(pdata_m->remote_m.end()); 
  }
  inline int sizeRemote() const 
  { 
    return pdata_m->remote_m.size(); 
  }
  //@}

  /// @name Iterators through local guard-cell-fill requests.
  //@{
  
  FillIterator_t beginFillList() const
  {
    return pdata_m->beginFillList();
  }
   
  FillIterator_t endFillList() const
  {
    return pdata_m->endFillList();
  }
  //@}

  //============================================================
  /// @name Touch methods
  //============================================================
  //@{

  /// Find all subdomains that touch on a given domain, and insert the
  /// intersection of these subdomains into the given output iterator.
  /// Return the number of touching elements. This version of touches
  /// can build either pointers or objects.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d,o,ctag);
  }
  
  /// Version that does the calculation based on allocated instead
  /// of owned domains.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesAlloc(const OtherDomain &d, OutIter o, 
                   const ConstructTag &ctag) const
  {
    return pdata_m->touchesAlloc(d, o, ctag);
  }  

  /// Find local subdomains that touch on a given domain, and insert
  /// the intersection of these subdomains into the given output
  /// iterator.  Return the number of touching elements. This version
  /// of touches can build either pointers or objects.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesLocal(const OtherDomain &d, OutIter o, 
			  const ConstructTag &ctag) const 
  {
    return pdata_m->touchesLocal(d, o, ctag);
  }
  
  /// Allocated version:
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesAllocLocal(const OtherDomain &d, OutIter o, 
			       const ConstructTag &ctag) const 
  {
    return pdata_m->touchesAllocLocal(d, o, ctag);
  }  

  /// Find remote subdomains that touch on a given domain, and insert
  /// the intersection of these subdomains into the given output
  /// iterator.  Return the number of touching elements. This version
  /// of touches can build either pointers or objects.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesRemote(const OtherDomain & d, OutIter o, 
			   const ConstructTag &ctag) const 
  {
    return pdata_m->touchesRemote(d,o,ctag);
  }

  /// Allocated version:
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesAllocRemote(const OtherDomain &d, OutIter o, 
			        const ConstructTag & ctag) const 
  {
    return pdata_m->touchesAllocRemote(d,o,ctag);
  }

  /// Find all subdomains that touch on a given domain,
  /// and insert the intersection of these subdomains into the given
  /// output iterator.  Return the number of touching elements.  These
  /// versions of touches can build only objects objects.
  
  template <class OtherDomain, class OutIter>
  inline int touches(const OtherDomain &d, OutIter o) const 
  {
    return touches(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesLocal(const OtherDomain &d, OutIter o) const 
  {
    return touchesLocal(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesRemote(const OtherDomain &d, OutIter o) const 
  {
    return touchesRemote(d, o, TouchesConstructNodeObj());
  }

  /// Versions of the above that do the allocated, instead of the
  /// owned, domains.
  
  template <class OtherDomain, class OutIter>
  int touchesAlloc(const OtherDomain &d, OutIter o) const
  {
    return touchesAlloc(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  int touchesAllocLocal(const OtherDomain &d, OutIter o) const
  {
    return touchesAllocLocal(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  int touchesAllocRemote(const OtherDomain &d, OutIter o) const
  {
    return touchesAllocRemote(d, o, TouchesConstructNodeObj());
  }
  //@}

#if !POOMA_NO_TEMPLATE_FRIENDS

  template <int Dim1, int Dim2, class lbd>
  friend class LayoutBaseView;

#endif
  
  //============================================================
  // Data
  //============================================================

  friend class LayoutBaseData<Dim>;
  
  /// LayoutBase stores its data in a RefCounted class to
  /// simplify memory management.
  
  RefCountedPtr<LBD> pdata_m;
};


/**
 * This is the actual data for the LayoutBaseView class, which is
 * simply a wrapper that holds a reference counted instance of this
 * data class.
 */

template <int Dim, int Dim2, class L>
class LayoutBaseViewData 
{
public:

  typedef L                                    Layout_t;

  typedef Interval<Dim>                        Domain_t;
  typedef Range<Dim2>                          BaseDomain_t;
  typedef int                                  Context_t;
  typedef Unique::Value_t                      ID_t;

  typedef typename Layout_t::Domain_t          AllocatedDomain_t;
  typedef ViewIndexer<Dim,Dim2>                Indexer_t;

  typedef Node<Domain_t,AllocatedDomain_t>     Value_t;
  typedef std::vector<Value_t *>               List_t;        // convenience
  typedef GuardLayers<Dim>                     GuardLayers_t; // convenience
  
  //============================================================
  // Constructors
  //============================================================
  
  LayoutBaseViewData() 
   : id_m(Unique::get())
  { }
  
  template<class DT>
  LayoutBaseViewData(const L & layout, const Domain<Dim,DT> & dom)
    : id_m(Unique::get()), layout_m(layout), 
      internalGuards_m(layout.internalGuards()),
      externalGuards_m(layout.externalGuards()),
      indexer_m(dom),
      subdomainsComputed_m(false)
  {
    // We cannot logically be a slice here.

    CTAssert(Dim == Dim2);

    // The layout passed in must be initialized.  

    PAssert(layout_m.initialized());

    // The domain we're passing in must be contained in the base
    // layout.

    PAssert(contains(layout_m.domain(), dom.unwrap()));
  }
  
  template <class DT>
  LayoutBaseViewData(const L &layout, const SliceDomain<DT> &dom)
  : id_m(Unique::get()), layout_m(layout), indexer_m(dom), 
    subdomainsComputed_m(false) 
  {
    // We are a slice and our dimensions must be consistent with us
    // and the layout we're being spawned by.

    CTAssert(Dim == DT::sliceDimensions);
    CTAssert(Dim2 == DT::dimensions);

    // The layout passed in must be initialized.  

    PAssert(layout_m.initialized());

    // The total domain we're passing in must be contained in the
    // base layout.

    PAssert(contains(layout_m.domain(), dom.totalDomain()));

    // Set up guard cell specifications on the reduced dimensionality domain.

    int dt, d;
    for (d = 0, dt = 0; dt < Dim2; ++dt)
      {
        if (!dom.ignorable(dt))
          {
            internalGuards_m.lower(d) = layout_m.internalGuards().lower(dt);
            internalGuards_m.upper(d) = layout_m.internalGuards().upper(dt);
            externalGuards_m.lower(d) = layout_m.externalGuards().lower(dt);
            externalGuards_m.upper(d) = layout_m.externalGuards().upper(dt);
            PAssert(d < Dim);
            ++d;
          }
      }
  }
  // construct a View of a View..
   template <class DT,class LV>
   LayoutBaseViewData(const L &layout,
		      const LV & viewLayout,
		      const Indexer_t & indexer, 
		      const Domain<Dim, DT> &dom,
		      GuardLayers_t ig,
		      GuardLayers_t eg)
  :
    id_m(Unique::get()),
    layout_m(layout), 
    internalGuards_m(ig),
    externalGuards_m(eg),
    indexer_m(indexer, dom),
    subdomainsComputed_m(false)
  {
    // The layout passed in must be initialized. 

    PAssert(layout_m.initialized());

    // The domain we're passing in must be contained in the base
    // layout.

    PAssert(contains(viewLayout.domain(), dom.unwrap()));
  }

// construct a slice of a View..

   template <class DT,class LV>
   LayoutBaseViewData(const L &layout,
		      const LV &viewLayout,
		      const Indexer_t indexer, 
		      const SliceDomain<DT> &dom)
     : id_m(Unique::get()), layout_m(layout), 
     indexer_m(indexer),
     subdomainsComputed_m(false)
  {
    // Our dimensionality must be the same as the slice's reduced
    // dimensionality.

    CTAssert((int)DT::sliceDimensions == Dim);

    // The slice's dimensionality must match that of the previous
    // view.

    CTAssert((int)DT::dimensions == LV::dimensions);

    // The layout passed in must be initialized.  

    PAssert(layout_m.initialized());

    // The total domain we're passing in must be contained in the
    // base layout.

    PAssert(contains(viewLayout.domain(), dom.totalDomain()));

    // Set up guard cell specifications on the reduced dimensionality domain.
    // NEEDS TESTED!!!

    int dt, d;
    for (d = 0, dt = 0; dt < LV::dimensions ; ++dt)
      {
        if (!dom.ignorable(dt))
          {
            internalGuards_m.lower(d) = viewLayout.internalGuards().lower(dt);
            internalGuards_m.upper(d) = viewLayout.internalGuards().upper(dt);
            externalGuards_m.lower(d) = viewLayout.externalGuards().lower(dt);
            externalGuards_m.upper(d) = viewLayout.externalGuards().upper(dt);
            PAssert(d < Dim);
            ++d;
          }
      }
  }

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, 
              const ConstructTag &ctag) const 
  {
    // Transform the local domain to base coordinates.

    BaseDomain_t bd = Pooma::NoInit();
    indexer_m.localToBase(d, bd);

    // Run the touches function for our underlying layout.

    std::vector<Node<BaseDomain_t,AllocatedDomain_t> > tnodes;
    int count = layout_m.touches(bd, std::back_inserter(tnodes));

    // Now, run through the nodes we've got and patch up the 
    // domains.

    Range<Dim> ld = Pooma::NoInit(); 

    for (int i = 0; i < count; i++)
      {
        *o++ = 
          touchesConstruct(indexer_m.baseToLocal(tnodes[i].domain(), ld),
                           tnodes[i].allocated(), // Do not convert this!
                           tnodes[i].affinity(), tnodes[i].context(),
                           tnodes[i].globalID(), tnodes[i].localID(), ctag);
      }

    // Done; return number we touched

    return count;
  }

  //============================================================
  // Utility functions
  //============================================================

  void computeSubdomains() const
  {
    // We don't need to do anything if we've already done this work.

    if (subdomainsComputed_m)
      return;

    // We need to find the nodes that intersect with our base
    // domain.  To do this, run the touches function for our
    // underlying layout.

    std::vector<Node<BaseDomain_t,AllocatedDomain_t> > tnodes;
    int count = layout_m.touches(indexer_m.baseDomain(), 
                                 std::back_inserter(tnodes));

    // Now, run through the nodes we've got and patch up the domains.

    Domain_t ld = Pooma::NoInit();

    for (int i = 0; i < count; ++i)
      {
        Value_t *pt =
          touchesConstruct(indexer_m.baseToLocal(tnodes[i].domain(), ld),
                           tnodes[i].allocated(),  // Do not convert this!
                           tnodes[i].affinity(),tnodes[i].context(),
                           tnodes[i].globalID(),tnodes[i].localID(),
                           TouchesConstructNodePtr());
        all_m.push_back(pt);
        if (pt->context() == Pooma::context()
	    ||pt->context() == -1 )
          local_m.push_back(pt);
        else
          remote_m.push_back(pt);
      }
    subdomainsComputed_m = true;
  }

  //============================================================
  // Data
  //============================================================

  // Our unique ID number.

  ID_t id_m;

  // A copy of the ultimate Layout object that we are viewing.

  L layout_m;

  //  LayoutBase<Dim2,lbd> layout_m;

  // Guard layer specification (different from layout_m due to potential
  // reduced dimensionality).

  GuardLayers_t internalGuards_m;
  GuardLayers_t externalGuards_m;

  // A ViewIndexer object, used to calculate indexes in the original
  // domain.

  Indexer_t indexer_m;

  // The list of all, local, and remote subdomains. These are evaluated in
  // a lazy manner so that we don't waste time if they're not actually used. 
  // Thus we declare them mutable since the actual evaluation may be forced  
  // by a const member function.

  mutable List_t all_m;
  mutable List_t local_m;
  mutable List_t remote_m;

  // Have we filled our subdomain lists?

  mutable bool subdomainsComputed_m;

};

/**
 * LayoutBaseView is a base class for all Layout classes that provides 
 * a view of an existing LayoutBase object. Dim is the logical dimension of
 * the layout. Dim2 is the dimension of the LayoutBase
 * contained within.
 *
 * To construct a LayoutBaseView, you need an existing
 * Layout or a LayoutView and the subdomain that
 * is being viewed. This class does not have a useful default
 * constructor since it is based on an existing LayoutBase.
 *
 * int globalID(const Loc<Dim> &pos) : return the globalID
 *    of the node that contains the point.
 */

template <int Dim, int Dim2, class lvd>
class LayoutBaseView
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.

  enum { dimensions = Dim };

  // General public typedefs.
  
  //  typedef LayoutBaseViewData<Dim, Dim2, lvd>      LayoutData_t; 

  typedef lvd                                       LayoutData_t; 

  typedef typename LayoutData_t::Domain_t           Domain_t;
  typedef typename LayoutData_t::BaseDomain_t       BaseDomain_t;
  typedef typename LayoutData_t::Context_t          Context_t;
  typedef typename LayoutData_t::ID_t               ID_t;
  
  typedef typename LayoutData_t::Layout_t           Layout_t;
  typedef typename LayoutData_t::AllocatedDomain_t  AllocatedDomain_t;
  typedef typename LayoutData_t::Value_t            Value_t;
  
  typedef typename LayoutData_t::List_t             List_t;
  typedef typename LayoutData_t::Indexer_t          Indexer_t;
  typedef typename LayoutData_t::GuardLayers_t      GuardLayers_t;
    
  typedef LayoutBaseView<Dim, Dim2, lvd>          This_t;      // convenience
  typedef LayoutBaseView<Dim, Dim2, lvd>          ViewLayout_t;// convenience

  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.
  
  typedef DerefIterator<Value_t>                    iterator;
  typedef ConstDerefIterator<Value_t>               const_iterator;

  //============================================================
  // Constructors
  //============================================================

  LayoutBaseView(LayoutData_t * lvdp)
    : pdata_m(lvdp)
  {}

  LayoutBaseView(const  RefCountedPtr<LayoutData_t> & pdata)
    : pdata_m(pdata)
  {}


  //============================================================
  // Accessors
  //============================================================

  //@{

  /// Return ID value. Our ID is unique. However, use the one from the
  /// layout we're viewing as the baseID.
  
  inline ID_t ID() const { return pdata_m->id_m; }
  inline ID_t baseID() const { return pdata_m->layout_m.baseID(); }

  //@}

  /// Return that this layout has been initialized. There is no other
  /// way to legally construct an object of this type.
  
  inline bool initialized() const { return true; }

  /// Return the global domain in our coordinate system and in that of
  /// our original domain.
  
  inline const Domain_t &domain() const 
  { 
    return pdata_m->indexer_m.domain(); 
  }
  
  inline const Domain_t &innerDomain() const 
  { 
    return pdata_m->indexer_m.innerDomain(); 
  }
  
  inline const BaseDomain_t &baseDomain() const 
  {
    return pdata_m->indexer_m.baseDomain(); 
  }
  
  inline const Layout_t &baseLayout() const
  {
    return pdata_m->layout_m;
  }

  template <class DT>
  BaseDomain_t &localToBase(const Domain<Dim, DT> &dlocal, 
    BaseDomain_t &base) const 
  {
    return pdata_m->indexer_m.localToBase(dlocal,base);
  }

  template <class DT>
  SliceRange<Dim2, Dim> &localToBase(const Domain<Dim, DT> &dlocal, 
    SliceRange<Dim2, Dim> &base) const 
  {
    return pdata_m->indexer_m.localToBase(dlocal,base);
  }

  inline GuardLayers_t internalGuards() const
  {
    return pdata_m->internalGuards_m;
  }

  inline GuardLayers_t externalGuards() const
  {
    return pdata_m->externalGuards_m;
  }

  /// Return the first index in the specified direction.
  /// (Always zero since this is a zero-based layout.)
  
  inline int first(int) const { return 0; }

  /// Compare to another Layout.  The layouts are the same if:
  ///   1. They have the same base ID value.
  ///   2. They have the same base domain.
  
  template <class L>
  inline bool operator==(const L &layout) const 
  {
    return (baseID() == layout.baseID() && 
            baseDomain() == layout.baseDomain());
  }

  /// Compare for inequality.

  template <class L>
  inline bool operator!=(const L &layout) const
  {
    return !(*this == layout);
  }

  //@{

  //-----------------------------------------------------------------------------
  //
  /// Return the globalID to the given element, expressed either as
  /// a single Loc or as a set of int's.
  ///
  /// These simply call the indexer's translate function and then pass the
  /// resulting Loc to the original layout to get the patch ID.
  //
  //-----------------------------------------------------------------------------

  inline int
  globalID(const Loc<Dim> &loc, Loc<Dim2> &oloc) const
  {
    pdata_m->indexer_m.translate(loc,oloc);  
    return pdata_m->layout_m.globalID(oloc);
  }

  inline int
  globalID(int i0, Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  inline int
  globalID(int i0, int i1, Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,i1,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  inline int
  globalID(int i0, int i1, int i2, Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,i1,i2,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  inline int
  globalID(int i0, int i1, int i2, int i3, 
           Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,i1,i2,i3,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  inline int
  globalID(int i0, int i1, int i2, int i3, 
           int i4, Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,i1,i2,i3,i4,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  inline int
  globalID(int i0, int i1, int i2, int i3, 
           int i4, int i5, Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,i1,i2,i3,i4,i5,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  inline int
  globalID(int i0, int i1, int i2, int i3, 
           int i4, int i5, int i6, Loc<Dim2> &loc) const
  {
    pdata_m->indexer_m.translate(i0,i1,i2,i3,i4,i5,i6,loc);  
    return pdata_m->layout_m.globalID(loc);
  }

  //@}

  //============================================================
  // Touch methods
  //============================================================

  /// Find all subdomains that touch on a given domain, and insert the
  /// intersection of these subdomains into the given output iterator.
  /// Return the number of touching elements. This version of touches
  /// can build either pointers or objects.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, 
    const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d,o,ctag);
  }

  /// Find local subdomains that touch on a given domain, and insert
  /// the intersection of these subdomains into the given output
  /// iterator.  Return the number of touching elements. This version
  /// of touches can build either pointers or objects. 

  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesLocal(const OtherDomain &d, OutIter o, 
    const ConstructTag &ctag) const {
    return pdata_m->touches(d, o, ctag);
  }

  /// Find remote subdomains that touch on a given domain, and insert
  /// the intersection of these subdomains into the given output
  /// iterator.  Return the number of touching elements. This version
  /// of touches can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesRemote(const OtherDomain &, OutIter, 
    const ConstructTag &) const {
    return 0;
  }

  //@{

  /// Find all/local/remote subdomains that touch on a given domain,
  /// and insert the intersection of these subdomains into the given
  /// output iterator.  Return the number of touching elements.  These
  /// versions of touches can build only objects objects.

  template <class OtherDomain, class OutIter>
  inline int touches(const OtherDomain &d, OutIter o) const {
    return touches(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesLocal(const OtherDomain &d, OutIter o) const {
    return touchesLocal(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesRemote(const OtherDomain &d, OutIter o) const {
    return touchesRemote(d, o, TouchesConstructNodeObj());
  }      

  //@}

  //============================================================
  // Iterators
  //============================================================

  //@{

  /// Return begin and end iterators for the list of all subdomains;
  /// take care to compute subdomains first
  
  inline iterator beginGlobal() { 
    computeSubdomains(); 
    return iterator(pdata_m->all_m.begin()); 
  }
  inline iterator endGlobal() { 
    computeSubdomains();
    return iterator(pdata_m->all_m.end()); 
  }
  inline const_iterator beginGlobal() const { 
    computeSubdomains(); 
    return const_iterator(pdata_m->all_m.begin()); 
  }
  inline const_iterator endGlobal() const { 
    computeSubdomains(); 
    return const_iterator(pdata_m->all_m.end()); 
  }
  inline int sizeGlobal() const { 
    computeSubdomains(); 
    return pdata_m->all_m.size(); 
  }

  //@}
  //@{

  /// Return begin and end iterators for the list of all subdomains;
  /// take care to compute subdomains first
  
  inline iterator beginLocal() { 
    computeSubdomains(); 
    return iterator(pdata_m->local_m.begin()); 
  }
  inline iterator endLocal() { 
    computeSubdomains(); 
    return iterator(pdata_m->local_m.end()); 
  }
  inline const_iterator beginLocal() const { 
    computeSubdomains(); 
    return const_iterator(pdata_m->local_m.begin()); 
  }
  inline const_iterator endLocal() const { 
    computeSubdomains(); 
    return const_iterator(pdata_m->local_m.end()); 
  }
  inline int sizeLocal() const { 
    computeSubdomains(); 
    return pdata_m->local_m.size(); 
  }

  //@}
  //@{

  /// Return begin and end iterators for the list of all subdomains;
  /// take care to compute subdomains first
  
  inline iterator beginRemote() {
    computeSubdomains(); 
    return iterator(pdata_m->remote_m.begin()); 
  }
  inline iterator endRemote() { 
    computeSubdomains(); 
    return iterator(pdata_m->remote_m.end()); 
  }
  inline const_iterator beginRemote() const { 
    computeSubdomains(); 
    return const_iterator(pdata_m->remote_m.begin()); 
  }
  inline const_iterator endRemote() const { 
    computeSubdomains(); 
    return const_iterator(pdata_m->remote_m.end());
  }
  inline int sizeRemote() const { 
    computeSubdomains(); 
    return pdata_m->remote_m.size(); 
  }

  //@}
  
#if !POOMA_NO_TEMPLATE_FRIENDS

  template <int OtherDim, int OtherDim2, class OtherLayoutData>
  friend class LayoutBaseView;

  template <int OtherDim, int OtherDim2, class OtherLayout>
  friend class LayoutBaseViewData;

#endif

  //============================================================
  // Private utility functions
  //============================================================

  /// Fill our subdomain lists.
  
  void computeSubdomains() const { pdata_m->computeSubdomains(); }
   
  //============================================================
  // Data
  //============================================================

  /// The data is stored in a RefCounted class to simplify memory
  /// management.  This is probably not as important for ViewLayout
  /// classes as for Layout classes, but we do it for
  /// consistency. Currently ViewLayouts are not observed directly by
  /// anyone. Of course, the Layout that we have a copy of is observed.

  RefCountedPtr<LayoutData_t> pdata_m;

};

// } // namespace POOMA

#endif // POOMA_LAYOUT_LAYOUTBASE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LayoutBase.h,v $   $Author: richi $
// $Revision: 1.30 $   $Date: 2004/11/10 22:17:08 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
