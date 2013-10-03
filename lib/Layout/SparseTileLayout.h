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

#ifndef POOMA_LAYOUT_SPARSETILELAYOUT_H
#define POOMA_LAYOUT_SPARSETILELAYOUT_H

//-----------------------------------------------------------------------------
// Classes: 
//   SparseTileLayout<Dim>
//   SparseTileLayoutView<Dim, Dim2>
//   SparseTileTag
//   MultiPatchLayoutTraits<SparseTileTag,Dim>
//-----------------------------------------------------------------------------



/** @file
 * @ingroup Layout
 * @brief
 *  SparseTileLayout<Dim>
 *        - Layout class that tiles a Dim-dimensional bounding box with
 *          non-overlaping sub-domains. The tiling does not have to be
 *          complete. 
 *   SparseTileLayoutView<Dim, Dim2>
 *     - view of a SparseTileLayout
 *   SparseTileTag
 *     - tag used to specialize MultiPatchLayoutTraits
 *   MultiPatchLayoutTraits<SparseTileTag,Dim>
 *     - traits class used by MultiPatch-engine to determine layout type.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

// The kitchen sink, this should be trimmed. 

#include "Layout/MultiPatchLayoutTraits.h"
#include "Layout/INode.h"
#include "Layout/TouchesConstruct.h"
#include "Layout/GuardLayers.h"
#include "Layout/DynamicEvents.h"
#include "Partition/ContextMapper.h"
#include "Partition/TilePartition.h"
#include "Domain/Interval.h"
#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/NewDomain.h"
#include "Domain/SliceRange.h"
#include "Domain/DomainMap.h"
#include "Utilities/DerefIterator.h"
#include "Utilities/ViewIndexer.h"
#include "Utilities/Observable.h"
#include "Utilities/Observer.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/RefCounted.h"
#include "Utilities/Unique.h"
#include "Utilities/PAssert.h"

#include "Layout/LayoutBase.h"

#include <vector>
#include <map>
#include <iosfwd>



///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//============================================================
// Forward declarations
//============================================================


template <int Dim> class SparseTileLayoutData;
template <int Dim> class SparseTileLayout;

template <int Dim, int Dim2> class SparseTileLayoutViewData;
template <int Dim, int Dim2> class SparseTileLayoutView;


/**
 * SparseTileTag class.
 */

struct SparseTileTag { };


/**
 * Specialization of MultiPatchLayoutTraits for SparseTileLayout.
 */

template <int Dim>
struct MultiPatchLayoutTraits<SparseTileTag,Dim>
{
  typedef SparseTileLayout<Dim> Layout_t;
  
  template <int ViewDim>
  struct View
  {
    typedef SparseTileLayoutView<ViewDim,Dim> Layout_t;
  };
};

//
//
// The touches operations of the of SparseTileLayout
// is done by looking up the patches in a DomainMap

template<int Dim> 
class SparseTileLayoutData
 : public LayoutBaseData<Dim>, 
   public RefCounted,
   public Observable< SparseTileLayoutData<Dim> >
{
public:
  
   // General public typedefs.

  typedef SparseTileLayoutData<Dim>            This_t;
  typedef Observable<This_t>                   Observable_t;
  typedef Interval<Dim>                        Domain_t;
  typedef Interval<Dim>                        BaseDomain_t;
  typedef Interval<Dim>                        AllocatedDomain_t;
  typedef int                                  Context_t;
  typedef Unique::Value_t                      ID_t;
  typedef Node<Domain_t,AllocatedDomain_t>     Value_t;
  typedef std::vector<Value_t *>               List_t;
  typedef std::map<int,Value_t>                Map_t;
  typedef GuardLayers<Dim>                     GuardLayers_t;    
  typedef std::pair<int,int>                   pidx_t;
  typedef typename DynamicEvents::PatchID_t    PatchID_t;
  typedef typename DynamicEvents::CreateSize_t CreateSize_t;
  typedef BaseDomain_t                         SubPatch_t;
  typedef std::vector<SubPatch_t>              PatchList_t;

 
  typedef typename LayoutBaseData<Dim>::GCFillInfo_t GCFillInfo_t;

  typedef typename std::vector<GCFillInfo_t>::const_iterator FillIterator_t;

  //=======================================================================

  struct GCBorderFillInfo
  {
    GCBorderFillInfo(const Domain_t &dom, int patchID)
      : domain_m(dom), patchID_m(patchID)
      {
      }
    
    // Get a CW warning about this not having a default constructor
    // when we instantiate the vector<GCBorderFillInfo> below. This never
    // gets called, so it seems bogus, but necessary.

    GCBorderFillInfo()
      {
	PInsist(0,"Shouldn't get here!");
      }
  
    Domain_t domain_m;    // guard layer domain
 
    int patchID_m;   // node ID for which domain_m is in the guards
   
    Domain_t domain() const { return domain_m;}
 
    int patchID() const { return patchID_m;}
  };
  
  typedef GCBorderFillInfo GCBorderFillInfo_t;

  typedef typename std::vector<GCBorderFillInfo>::const_iterator
                                               BorderFillIterator_t;

  enum { dimensions = Dim };
  enum { repartitionEvent = 1 };
  enum { dynamic = false };


  //============================================================
  // Constructors
  //============================================================

  // Default constructor: initially no blocks, etc.
  
  SparseTileLayoutData();

  // the "normal" constructors

  SparseTileLayoutData(const Domain_t &,
		       const PatchList_t  &,
		       const ContextMapper<Dim> &);
  
  SparseTileLayoutData(const Domain_t &boundingbox,
		       const GuardLayers_t & globalGL,
		       const PatchList_t  & PatchList,
		       const ContextMapper<Dim> &);

  SparseTileLayoutData(const Domain_t &boundingbox,
		       const GuardLayers_t & internalGL,
		       const GuardLayers_t & externalGL,
		       const PatchList_t  & PatchList,
		       const ContextMapper<Dim> &);
 
 // constructors that don't take mapper tags

  SparseTileLayoutData(const Domain_t &boundingbox);
  
  SparseTileLayoutData(const Domain_t &boundingbox,
		       const GuardLayers_t & globalGL);

  SparseTileLayoutData(const Domain_t &boundingbox,
		       const GuardLayers_t & internalGL,
		       const GuardLayers_t & externalGL);


  // Constructor based on a Partitioner

template <class Partitioner>
SparseTileLayoutData(const Domain_t &bbox, 
		     const Partitioner & gpar,
		     const ContextMapper<Dim> &cmap);


  //============================================================
  // Destructor
  //============================================================

  // This is only called when all references to this data go away, in
  // which case we need to delete our nodes. The Observable destructor
  // will broadcast messages up to all observers of the Layout.

  ~SparseTileLayoutData() ;

  //============================================================
  // Mutators
  //============================================================
  void initialize(const Domain_t &bbox);

  void initialize(const Domain_t &bbox,
		  const GuardLayers_t & globalGL);
  
  void initialize(const Domain_t &bbox,
		  const GuardLayers_t & internalGL,
		  const GuardLayers_t & externalGL);


  void initialize(const Domain_t &bbox,
		  const PatchList_t &plist,
		  const ContextMapper<Dim> &cmap);

  void initialize(const Domain_t &bbox,
		  const GuardLayers_t & globalGL,
		  const PatchList_t &plist,
		  const ContextMapper<Dim> &cmap);

  void initialize(const Domain_t &bbox,
		  const GuardLayers_t & internalGL,
		  const GuardLayers_t & externalGL,
		  const PatchList_t &plist,
		  const ContextMapper<Dim> &cmap);
 
  template <class Partitioner>
  void initialize(const Domain_t &bbox,
		  const Partitioner &gpar,
		  const ContextMapper<Dim> &cmap);
   
  void syncPatch();
  
  void calcMaps();
  
  void calcAllocMaps() ;

  //============================================================
  // Accessors
  //============================================================


  BorderFillIterator_t beginBorderFillList() const 
  { 
    return gcBorderFillList_m.begin(); 
  }

  BorderFillIterator_t endBorderFillList() const 
  { 
    return gcBorderFillList_m.end(); 
  }

  // Accessors for getting the global ID of the patch containing
  // a particular element.

  int globalID(const Loc<Dim> &loc) const;
  int globalID(int) const;
  int globalID(int,int) const;
  int globalID(int,int,int) const;
  int globalID(int,int,int,int) const;
  int globalID(int,int,int,int,int) const;
  int globalID(int,int,int,int,int,int) const;
  int globalID(int,int,int,int,int,int,int) const;


  //============================================================
  // touches operations
  //============================================================

  // Find all subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.
  // Return the number of touching elements. This version of touches
  // can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, 
	      OutIter o, 
	      const ConstructTag &ctag) const;

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesAlloc(const OtherDomain &d, 
		   OutIter o, 
		   const ConstructTag &ctag) const;


  template<class Out>
  void print(Out & o) const;


private:
  //============================================================
  // Private Methods
  //============================================================

  // This function calculates the cached guard-cell filling information.
  
  void calcGCFillList();

  // This function recalculates what the total domain of each patch
  // should be, since this can change due to dynamic operations.

  void calcDomains();

  // This function recalculates the domain maps, since this can
  // change due to dynamic operations.

  void calcMaps() const;
  void calcAllocMaps() const;

  //============================================================
  // Private Data
  //============================================================

  
  // Cached border guard-cell filling info.
  
  std::vector<GCBorderFillInfo> gcBorderFillList_m;

  // This DomainMap is used for touches operations on the non GC'ed
  // patches.

  mutable DomainMap<Interval<Dim>,pidx_t> map_m;

  // This DomainMap is used for touches operations on the  GC'ed
  // patches.

  mutable DomainMap<Interval<Dim>,pidx_t> mapAloc_m;
};


template <int Dim>
class SparseTileLayout : public LayoutBase<Dim,SparseTileLayoutData<Dim> >,
                         public Observable<SparseTileLayout<Dim> >,
                         public Observer<SparseTileLayoutData<Dim> >
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.

  enum { dimensions = Dim };
  enum { repartitionEvent = 1 };
  enum { dynamic = true };

  // General public typedefs.

  typedef SparseTileLayout<Dim>                This_t;
  typedef Observable<This_t>                   Observable_t;
  typedef SparseTileLayoutData<Dim>            LayoutData_t;
  typedef typename LayoutData_t::Domain_t      Domain_t;
  typedef typename LayoutData_t::BaseDomain_t  BaseDomain_t;
  typedef typename LayoutData_t::Context_t     Context_t;
  typedef typename LayoutData_t::ID_t          ID_t;
  typedef typename LayoutData_t::Value_t       Value_t;
  typedef typename LayoutData_t::List_t        List_t;
  typedef DynamicEvents::PatchID_t             PatchID_t;
  typedef DynamicEvents::CreateSize_t          CreateSize_t;
  typedef GuardLayers<Dim>                     GuardLayers_t;
  typedef typename LayoutData_t::SubPatch_t    SubPatch_t;
  typedef typename LayoutData_t::PatchList_t   PatchList_t;

 

  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.  
  
  typedef DerefIterator<Value_t>               iterator;
  typedef ConstDerefIterator<Value_t>          const_iterator;
   
  // Iterator through guard-cell-fill requests. 
  
  typedef typename LayoutData_t::GCFillInfo_t           GCFillInfo_t;
  
  typedef typename LayoutData_t::FillIterator_t         FillIterator_t;
  
  typedef typename LayoutData_t::BorderFillIterator_t   BorderFillIterator_t;

  //============================================================
  // Constructors
  //============================================================

  // The default constructor does not initialize the layout.  In this
  // case, layout initialization must be completed with the
  // 'initialize' method before the layout can be used.  A default
  // layout has an empty global domain, and empty subdomain lists.
  
  // constructors with no mapper tags

  SparseTileLayout();

  SparseTileLayout(const Domain_t &boundingbox);
  
  SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & globalGL);
  
  SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & internalGL,
		   const GuardLayers_t & externalGL);

  // ReplicatedTag mapper constructors.

  SparseTileLayout(Domain_t & boundingbox,
		   const PatchList_t &patchlist,
		   const ReplicatedTag &);

  SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & globalGL,
		   const PatchList_t  & PatchList,
		   const ReplicatedTag &);
  
  SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & internalGL,
		   const GuardLayers_t & externalGL,
		   const PatchList_t  & PatchList,
		   const ReplicatedTag &);
  
  template<class Partitioner>
  SparseTileLayout(const Domain_t &bbox,
		   const Partitioner &gpar,
		   const ReplicatedTag &);

  // DistributedTag mapper 

  SparseTileLayout(Domain_t & boundingbox,
		   const PatchList_t &patchlist,
		   const DistributedTag &);

  SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & globalGL,
		   const PatchList_t  & PatchList,
		   const DistributedTag &);
  
  SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & internalGL,
		   const GuardLayers_t & externalGL,
		   const PatchList_t  & PatchList,
		   const DistributedTag &);
  
  template<class Partitioner>
  SparseTileLayout(const Domain_t &bbox,
		   const Partitioner &gpar,
		   const DistributedTag &);

  // fully specified constructor

  template<class Partitioner>
  SparseTileLayout(const Domain_t &bbox,
		   const Partitioner &gpar,
		   const ContextMapper<Dim> &cmap);


  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
  SparseTileLayout(const This_t &);

  This_t & 
  operator=(const This_t &model)
    
  {
    if (this != &model)
      {
	this->pdata_m->detach(*this);
	this->pdata_m = model.pdata_m;
	this->pdata_m->attach(*this);
      }
    
    return *this;
  }
  

  //============================================================
  // Destructor
  //============================================================
  // The actual data will be cleaned up by the LayoutData_t destructor
  // if all references to the data go away.
  // If any Observers remain, they will be notified by the Observable 
  // destructor. 

  ~SparseTileLayout()
  {
    this->pdata_m->detach(*this);
  }
  

  //============================================================
  // Initialize methods
  //============================================================


  // Initialize a layout with nothing else but a global domain.

  void initialize(const Domain_t & a);
  
  void initialize(const Domain_t &,
         	    const GuardLayers_t &);

  void initialize(const Domain_t &,
         	    const GuardLayers_t &,
		   const PatchList_t & );


  template<class Partitioner>
  void initialize(const Domain_t &bbox,
		    const Partitioner &gpar);

  //============================================================
  // Data lookup
  //============================================================

  // Iterators through border guard-cell-fill requests.
  
  BorderFillIterator_t beginBorderFillList() const
    {
      return this->pdata_m->beginBorderFillList();
    }
   
  BorderFillIterator_t endBorderFillList() const
    {
      return this->pdata_m->endBorderFillList();
    }


  void syncPatch();


  //============================================================
  // Observer methods
  //============================================================

  // Respond to events generated by the LayoutData_t.
  // These are just passed on to our observers.

  virtual void notify(LayoutData_t &d, const ObserverEvent &event)
    {
      // We should only get this message from our LayoutData_t object
      PAssert(&d == this->pdata_m.rawPointer());
      Observable_t::notify(event);
    }


  //============================================================
  // Output
  //============================================================

  // Print a SparseTileLayout on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const {
//     ostr << "SparseTileLayout " << ID() << " on global domain " 
//       << domain() << ":" << '\n';
//     ostr << "   Total subdomains: " << size() << '\n';
//     ostr << "   Local subdomains: " << sizeLocal() << '\n';
//     ostr << "  Remote subdomains: " << sizeRemote() << '\n';
//     typename SparseTileLayout<Dim>::const_iterator a;
//     for (a = begin(); a != end(); ++a)
//       ostr << "  Global subdomain = " << *a << '\n';
//     for (a = beginLocal(); a != endLocal(); ++a)
//       ostr << "   Local subdomain = " << *a << '\n';
//     for (a = beginRemote(); a != endRemote(); ++a)
//       ostr << "  Remote subdomain = " << *a << '\n';
    this->pdata_m->print(ostr);
  }

#if !POOMA_NO_TEMPLATE_FRIENDS

  //private:

  template <int Dim1, int Dim2>
  friend class SparseTileLayoutView;

#endif

  //============================================================
  // Data
  //============================================================

  // SparseTileLayout stores its data in a LayoutBase::RefCounted class to
  // simplify memory management.
  
  friend class SparseTileLayoutData<Dim>;
};

/**
 * The data object held by a SparseTileLayoutView object.
 */

template <int Dim, int Dim2>
class SparseTileLayoutViewData 
: public LayoutBaseViewData<Dim, Dim2, SparseTileLayout<Dim2> >,
  public RefCounted
{
public:

  typedef SparseTileLayout<Dim2>                  Layout_t;
  typedef SparseTileLayoutView<Dim, Dim2>         ViewLayout_t;

  typedef Interval<Dim>                     Domain_t;
  typedef Range<Dim2>                       BaseDomain_t;
  typedef int                               Context_t;
  typedef Unique::Value_t                   ID_t;

  typedef typename Layout_t::Domain_t       AllocatedDomain_t;
  typedef ViewIndexer<Dim,Dim2>             Indexer_t;

  typedef Node<Domain_t,AllocatedDomain_t>  Value_t;
  typedef std::vector<Value_t *>            List_t;        // for convenience
  typedef GuardLayers<Dim>                  GuardLayers_t; // for convenience

  typedef SparseTileLayoutViewData<Dim,Dim2>      LayoutData_t;

  //============================================================
  // Constructors
  //============================================================

  SparseTileLayoutViewData() { }
  
  template <class DT>
  inline SparseTileLayoutViewData(const Layout_t &layout, const Domain<Dim, DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,SparseTileLayout<Dim2> >(layout,dom)
  { 
    // We cannot logically be a slice here.

    CTAssert(Dim == Dim2);

    // The layout passed in must be initialized.  

    PAssert(this->layout_m.initialized());

    // The domain we're passing in must be contained in the base
    // layout.

    PAssert(contains(this->layout_m.domain(), dom.unwrap()));
  }

  template <class DT>
  inline SparseTileLayoutViewData(const Layout_t &layout, const SliceDomain<DT> &dom)
  :LayoutBaseViewData<Dim,Dim2,SparseTileLayout<Dim2> >(layout,dom)
  {  
    // We are a slice and our dimensions must be consistent with us
    // and the layout we're being spawned by.

    CTAssert(Dim == DT::sliceDimensions);
    CTAssert(Dim2 == DT::dimensions);

    // The layout passed in must be initialized.  

    PAssert(this->layout_m.initialized());

    // The total domain we're passing in must be contained in the
    // base layout.

    PAssert(contains(this->layout_m.domain(), dom.totalDomain()));

    // Set up the guard cell specifications on reduced dimensionality domain.

    int dt, d;
    for (d = 0, dt = 0; dt < Dim2; ++dt)
      {
	if (!dom.ignorable(dt))
	  {
	    this->internalGuards_m.lower(d) = this->layout_m.internalGuards().lower(dt);
	    this->internalGuards_m.upper(d) = this->layout_m.internalGuards().upper(dt);
	    this->externalGuards_m.lower(d) = this->layout_m.externalGuards().lower(dt);
	    this->externalGuards_m.upper(d) = this->layout_m.externalGuards().upper(dt);
	    PAssert(d < Dim);
	    ++d;
	  }
      }
  }

  template <class DT>
  SparseTileLayoutViewData(const ViewLayout_t &layout, const Domain<Dim, DT> &dom)
  :  LayoutBaseViewData<Dim,Dim2,SparseTileLayout<Dim2> >(
					      layout.pdata_m->layout_m,
					      layout,
					      layout.pdata_m->indexer_m, 
					      dom,
					      layout.internalGuards(),
					      layout.externalGuards())
  {
    // The layout passed in must be initialized. 

    PAssert(this->layout_m.initialized());

    // The domain we're passing in must be contained in the base layout.

    PAssert(contains(layout.domain(), dom.unwrap()));
  }

  template <int OrigDim, class DT>
  SparseTileLayoutViewData(const SparseTileLayoutView<OrigDim, Dim2> &layout, 
			   const SliceDomain<DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,SparseTileLayout<Dim2> >(
			     layout.pdata_m->layout_m,
			     layout,
			     Indexer_t(layout.pdata_m->indexer_m,dom),
			     dom)
  {  
    // Our dimensionality must be the same as the slice's  reduced
    // dimensionality.

    CTAssert(DT::sliceDimensions == Dim);

    // The slice's dimensionality must match that of the previous view.

    CTAssert(DT::dimensions == OrigDim);

    // The layout passed in must be initialized.  

    PAssert(this->layout_m.initialized());

    // The total domain we're passing in must be contained in the
    // base layout.

    PAssert(contains(layout.domain(), dom.totalDomain()));

    // Set up the guard cell specifications on reduced dimensionality domain.
    // NEEDS TESTED!!!

    int dt, d;
    for (d = 0, dt = 0; dt < OrigDim; ++dt)
      {
	if (!dom.ignorable(dt))
	  {
	    this->internalGuards_m.lower(d) = layout.internalGuards().lower(dt);
	    this->internalGuards_m.upper(d) = layout.internalGuards().upper(dt);
	    this->externalGuards_m.lower(d) = layout.externalGuards().lower(dt);
	    this->externalGuards_m.upper(d) = layout.externalGuards().upper(dt);
	    PAssert(d < Dim);
	    ++d;
	  }
      }
  }

  // Destructor

  ~SparseTileLayoutViewData() 
  {
    typename List_t::iterator a;
    for (a = this->all_m.begin(); a != this->all_m.end(); ++a)
      delete (*a);
  }

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, 
	      const ConstructTag &ctag) const 
  {
    // transform the local domain to base coordinates

    BaseDomain_t bd = Pooma::NoInit();
    this->indexer_m.localToBase(d, bd);

    // run the touches function for our underlying layout

    std::vector<Node<BaseDomain_t,AllocatedDomain_t> > tnodes;
    int count = this->layout_m.touches(bd, std::back_inserter(tnodes));

    // now, run through the nodes we've got and patch the 
    // domains up

    Range<Dim> ld = Pooma::NoInit();

    for (int i = 0; i < count; i++)
      {
	*o++ = 
	  touchesConstruct(this->indexer_m.baseToLocal(tnodes[i].domain(), ld),
			   tnodes[i].allocated(), // Don't try to convert this!
			   tnodes[i].affinity(),tnodes[i].context(),
			   tnodes[i].globalID(), tnodes[i].localID(), ctag);
      }

    // done; return number we touched

    return count;
  }

  //============================================================
  // Utility functions
  //============================================================

  void computeSubdomains() const
  {
    // We don't need to do anything if we've already done this work.

    if (this->subdomainsComputed_m)
      return;

    // We need to find the nodes that intersect with our base domain.
    // To do this, run the touches function for our underlying layout.

    std::vector<Node<BaseDomain_t,AllocatedDomain_t> > tnodes;
    int count = this->layout_m.touches(this->indexer_m.baseDomain(), 
				       std::back_inserter(tnodes));

    // Now, run through the nodes we've got and patch the domains up.

    Domain_t ld = Pooma::NoInit();

    for (int i = 0; i < count; i++)
      {
	Value_t *pt =
	  touchesConstruct(this->indexer_m.baseToLocal(tnodes[i].domain(), ld),
			   tnodes[i].allocated(), // Don't try to convert this!
			   tnodes[i].affinity(),tnodes[i].context(),
			   tnodes[i].globalID(),tnodes[i].localID(),
			   TouchesConstructNodePtr());
	this->all_m.push_back(pt);
      }

 
    // Leave after setting the flag indicating we've computed these
    // subdomains.

    this->subdomainsComputed_m = true;

  }

  // All data members are inherited from LayoutBaseViewData
  
};


template <int Dim, int Dim2>
class SparseTileLayoutView
: public LayoutBaseView<Dim, Dim2, SparseTileLayoutViewData<Dim,Dim2> >
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.
  enum { dimensions = Dim };

  // General public typedefs.

  typedef SparseTileLayoutViewData<Dim, Dim2>            LayoutData_t; 

  typedef typename LayoutData_t::Domain_t          Domain_t;
  typedef typename LayoutData_t::BaseDomain_t      BaseDomain_t;
  typedef typename LayoutData_t::Context_t         Context_t;
  typedef typename LayoutData_t::ID_t              ID_t;
  
  typedef typename LayoutData_t::Layout_t          Layout_t;
  typedef typename LayoutData_t::AllocatedDomain_t AllocatedDomain_t;
  typedef typename LayoutData_t::Value_t           Value_t;
  
  typedef typename LayoutData_t::List_t            List_t;
  typedef typename LayoutData_t::Indexer_t         Indexer_t;
  typedef typename LayoutData_t::GuardLayers_t     GuardLayers_t;

  typedef SparseTileLayoutView<Dim, Dim2>                This_t;       // convenience
  typedef SparseTileLayoutView<Dim, Dim2>                ViewLayout_t; // convenience
  typedef LayoutBaseView<Dim,Dim2,LayoutData_t>    Base_t;
  
  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.
  
  typedef DerefIterator<Value_t>                   iterator;
  typedef ConstDerefIterator<Value_t>              const_iterator;


  //============================================================
  // Constructors
  //============================================================

  // Default constructor
  
  SparseTileLayoutView()
  : Base_t(new LayoutData_t())
  { }
  
  // Constructor building a SparseTileLayoutView from a
  // SparseTileLayout and a non-slice domain like an Interval<Dim> or
  // Range<Dim>.
  
  template <class DT>
  SparseTileLayoutView(const Layout_t &layout, const Domain<Dim2, DT> &dom)
  : LayoutBaseView<Dim,Dim2,SparseTileLayoutViewData<Dim,Dim2> >
    (new SparseTileLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Constructor building a SparseTileLayoutView from a
  // SparseTileLayout and a slice domain like a
  // SliceInterval<Dim2,Dim> or SliceRange<Dim2,Dim>.
  
  template <class DT>
  SparseTileLayoutView(const Layout_t &layout, const SliceDomain<DT> &dom)
  : LayoutBaseView<Dim,Dim2,SparseTileLayoutViewData<Dim,Dim2> >
    (new SparseTileLayoutViewData<Dim,Dim2>(layout,dom))
  { }
  
  // Constructor building a SparseTileLayoutView from another
  // SparseTileLayoutView and a non-slice domain like an
  // Interval<Dim> or Range<Dim>.
  
  template <class DT>
  SparseTileLayoutView(const ViewLayout_t &layout, const Domain<Dim, DT> &dom)
  :  LayoutBaseView<Dim,Dim2,SparseTileLayoutViewData<Dim,Dim2> >
    (new SparseTileLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Constructor building a SparseTileLayoutView from another
  // SparseTileLayoutView and a slice domain like a
  // SliceInterval<Dim2,Dim> or SliceRange<Dim2,Dim>.
  
  template <int OldViewDim, class DT>
  SparseTileLayoutView(const SparseTileLayoutView<OldViewDim, Dim2> &layout, 
                        const SliceDomain<DT> &dom)
  :  LayoutBaseView<Dim,Dim2,SparseTileLayoutViewData<Dim,Dim2> >
    (new SparseTileLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
  inline SparseTileLayoutView(const This_t &model) 
    : LayoutBaseView<Dim,Dim2,SparseTileLayoutViewData<Dim,Dim2> >
  (model.pdata_m)
  { }
  
  inline This_t &operator=(const This_t &model) 
  {
    if (this != &model)
      {
        this->pdata_m = model.pdata_m;
      }
    return *this;
  }


  //============================================================
  // Destructor
  //============================================================

  // The actual data will be cleaned up by the SparseTileLayoutData
  // destructor if all references to the data go away, so there is
  // nothing to do here.
   
  inline ~SparseTileLayoutView() 
  { }


  //============================================================
  // Accessors
  //============================================================

  // Return ID value. Our ID is unique. However, use the one from the 
  
  //============================================================
  // Iterators
  //============================================================


  //============================================================
  // Output
  //============================================================
    
  // Print a SparseTileLayoutView on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const 
  {
    ostr << "SparseTileLayoutView " << this->ID() << " on global domain " 
      << this->domain() << ":" << '\n';
    ostr << "   Base ID:          " << this->baseID() << '\n';
    ostr << "   Base domain:      " << this->baseDomain() << '\n';
    ostr << "   Total subdomains: " << this->sizeGlobal() << '\n';
    ostr << "   Local subdomains: " << this->sizeLocal() << '\n';
    ostr << "  Remote subdomains: " << this->sizeRemote() << '\n';
    const_iterator a;
    for (a = this->beginGlobal(); a != this->endGlobal(); ++a)
      ostr << "  Global subdomain = " << *a << '\n';
    for (a = this->beginLocal(); a != this->endLocal(); ++a)
      ostr << "   Local subdomain = " << *a << '\n';
    for (a = this->beginRemote(); a != this->endRemote(); ++a)
      ostr << "  Remote subdomain = " << *a << '\n';
  }

#if !POOMA_NO_TEMPLATE_FRIENDS

  //private:

  template <int OtherDim, int OtherDim2>
  friend class SparseTileLayoutView;

  template <int OtherDim, int OtherDim2>
  friend class SparseTileLayoutViewData;

#endif

  //============================================================
  // Private utility functions
  //============================================================

  // Fill our subdomain lists.
  
  void computeSubdomains() const { this->pdata_m->computeSubdomains(); }
 
};


//============================================================
// SparseTileLayoutView inline method definitions
//============================================================

//============================================================
// NewDomain1 traits classes for SparseTileLayout and SparseTileLayoutView
//============================================================

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a SparseTileLayout.
//
//-----------------------------------------------------------------------------

template <int Dim>
struct NewDomain1<SparseTileLayout<Dim> >
{
  typedef SparseTileLayout<Dim> &Type_t;

  inline static Type_t combine(const SparseTileLayout<Dim> &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a SparseTileLayoutView.
//
//-----------------------------------------------------------------------------

template <int Dim, int Dim2>
struct NewDomain1<SparseTileLayoutView<Dim, Dim2> >
{
  typedef SparseTileLayoutView<Dim, Dim2> &Type_t;

  inline static Type_t combine(const SparseTileLayoutView<Dim, Dim2> &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// ostream inserters for SparseTileLayout and SparseTileLayoutView:
//
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &ostr, 
			 const SparseTileLayout<Dim> &layout)
{
  layout.print(ostr);
  return ostr;
}

template <int Dim, int Dim2>
std::ostream &operator<<(std::ostream &ostr, 
			 const SparseTileLayoutView<Dim, Dim2> &layout)
{
  layout.print(ostr);
  return ostr;
}

//-----------------------------------------------------------------------------
//
// IsValidLocations specialization for S.T.L.. Used in PrintArray
//
//-----------------------------------------------------------------------------

template<int Dim> struct IsValid;

template<class LayoutTag, class PatchTag> struct MultiPatch;
template<class LayoutTag, class PatchTag, int Dim2> struct MultiPatchView;
template<class Expr> struct ExpressionTag;
template<class Eng, class Tag> struct EngineFunctor;

template<class Object,class Dom,class PatchTag>
inline bool isValidLocation(const Object  & e,
			    Dom & domain,
			    MultiPatch<SparseTileTag,PatchTag> &) 
{
  typedef typename Object::Domain_t domain_t;
  typedef Node<domain_t,domain_t> node_t;
  std::vector<node_t> v;
  int count = e.engine().layout().touches(domain,std::back_inserter(v));
  return (count!=0);
}

template<class Object,class Dom,class PatchTag,int Dim2>
inline bool isValidLocation(const Object  & e,
			    Dom & domain,
			    MultiPatchView<SparseTileTag,
			                   PatchTag,
			                   Dim2> &) 
{
  typedef typename Object::Domain_t domain_t;
  typedef Node<domain_t,domain_t> node_t;
  std::vector<node_t> v;
  int count = e.engine().layout().touches(domain,std::back_inserter(v));
  return (count!=0);
}

template<class Object,class Dom,class expr>
inline bool isValidLocation(const Object  & e,
			    Dom & domain,
			    ExpressionTag<expr> &) 
{
  typedef typename Object::Domain_t domain_t;
  typedef Node<domain_t,domain_t> node_t;
  std::vector<node_t> v;
  typedef typename Object::Engine_t Engine_t;

  IsValid<Engine_t::dimensions> l(domain);
  EngineFunctor<Engine_t, IsValid<Engine_t::dimensions> > ef;

  return ef.apply(e.engine(),l);
  //  int count = e.engine().layout().touches(domain,std::back_inserter(v));
}

// } // namespace POOMA

#include "Layout/SparseTileLayout.cpp"

#endif // POOMA_LAYOUT_SPARSETILELAYOUT_H
