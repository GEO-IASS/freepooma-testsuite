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

#ifndef POOMA_LAYOUT_GRIDLAYOUT_H
#define POOMA_LAYOUT_GRIDLAYOUT_H

//-----------------------------------------------------------------------------
// Classes: 
//   GridLayout<Dim>
//   GridLayoutView<Dim, Dim2>
//   GridTag
//   MultiPatchLayoutTraits<GridTag,Dim>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Layout
 * @brief
 * Layout with domain breaked into sub-domains specified by a grid.
 *
 *   GridLayout<Dim>
 *     - Layout class that breaks Dim-dimensional domain into
 *       sub-domains arranged in a Dim-dimensional grid, where the
 *       sub-domains are defined by a set of intervals along each axis. 
 *
 *   GridLayoutView<Dim, Dim2>
 *     - view of a GridLayout
 *   GridTag
 *     - tag used to specialize MultiPatchLayoutTraits
 *   MultiPatchLayoutTraits<GridTag,Dim>
 *     - traits class used by MultiPatch-engine to determine layout type.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Layout/MultiPatchLayoutTraits.h"
#include "Layout/INode.h"
#include "Layout/TouchesConstruct.h"
#include "Layout/GuardLayers.h"
#include "Layout/DynamicEvents.h"
#include "Partition/UniformGridPartition.h"
#include "Partition/GridPartition.h"
#include "Partition/ContextMapper.h"
#include "Domain/Interval.h"
#include "Domain/Grid.h"
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
#include <iosfwd>


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//============================================================
// Forward declarations
//============================================================

template <int Dim> class GridLayoutData;
template <int Dim> class GridLayout;

template <int Dim, int Dim2> class GridLayoutViewData;
template <int Dim, int Dim2> class GridLayoutView;


/**
 * GridTag class.
 */

struct GridTag { };


/**
 * Specialization of MultiPatchLayoutTraits for GridLayout.
 */

template <int Dim>
struct MultiPatchLayoutTraits<GridTag,Dim>
{
  typedef GridLayout<Dim> Layout_t;
  
  template <int ViewDim>
  struct View
  {
    typedef GridLayoutView<ViewDim,Dim> Layout_t;
  };
};


/**
 * Holds the data for a GridLayout. That class has a ref counted
 * instance of this class
 */

template<int Dim>
class GridLayoutData 
  : public LayoutBaseData<Dim>,
    public RefCounted,
    public Observable< GridLayoutData<Dim> >
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef GridLayoutData<Dim>                  This_t;
  typedef Observable<This_t>                   Observable_t;
  typedef Interval<Dim>                        Domain_t;
  typedef Interval<Dim>                        BaseDomain_t;
  typedef int                                  Context_t;
  typedef Unique::Value_t                      ID_t;
  typedef Node<Domain_t>                       Value_t;
  typedef std::vector<Value_t *>               List_t;
  typedef GuardLayers<Dim>                     GuardLayers_t;    
  typedef int                                  AxisIndex_t;
  typedef typename DynamicEvents::PatchID_t    PatchID_t;
  typedef typename DynamicEvents::CreateSize_t CreateSize_t;

  typedef typename LayoutBaseData<Dim>::GCFillInfo_t GCFillInfo_t;

  typedef typename std::vector<GCFillInfo_t>::const_iterator FillIterator_t;

  // Enumerations.

  enum { dimensions = Dim };
  enum { repartitionEvent = 1 };
  enum { dynamic = false };


  //============================================================
  // Constructors
  //============================================================

 
  /// Default constructor for GridLayoutData, this sets up this object
  /// to look like an "empty" layout, with no patches and an empty domain
  /// with no guard cells.  The "initialize" method can be used to change
  /// it to a different state.

  GridLayoutData();


  /// If we specify a partitioner, then we'll get the guard-cell sizes
  /// from it.

  template <class Partitioner>
  GridLayoutData(const Domain_t &gdom, 
		 const Partitioner &gpar,
		 const ContextMapper<Dim> & cmap);
  
  /// If we specify a partitioner, then we'll get the guard-cell sizes
  /// from it.

  template <class Partitioner>
  GridLayoutData(const Grid<Dim> &gdom, 
		 const Partitioner &gpar,
                 const ContextMapper<Dim> &cmap);


  //============================================================
  // Destructor
  //============================================================

  /// This is only called when all references to this data go away, in
  /// which case we need to delete our nodes. The Observable destructor
  /// will broadcast messages up to all observers of the Layout.

  ~GridLayoutData();
  
  //============================================================
  // Mutators
  //============================================================

  /// Initialize this object by invoking the partitioner and setting
  /// up the domains and guard cell filling info.  This can be called
  /// after using the default constructor.

  template <class Partitioner>
  void initialize(const Domain_t &gdom, 
		  const Partitioner &gpar,
		  const ContextMapper<Dim> &cmap);

  /// Initialize this object by invoking the partitioner and setting
  /// up the domains and guard cell filling info.  This can be called
  /// after using the default constructor.

  template <class Partitioner>
  inline void initialize(const Grid<Dim> &gdom, 
			 const Partitioner &gpar,
                         const ContextMapper<Dim> &cmap);

  /// Used by the I/O or data management system to initialize the layout based
  /// on detailed state information previously stored.

  void initialize(const Domain_t& idom,
		  const List_t& nodes,
		  const Loc<Dim>& blocks,
		  bool hasIG, bool hasEG,
		  const GuardLayers_t& ig,
		  const GuardLayers_t& eg);

  //============================================================
  ///@name Accessors
  //============================================================
  //@{

  inline const Loc<Dim> &blocks() const
    {
      return this->blocks_m;
    }

  inline bool dirty() const
    {
      return dirtyLayout_m;
    }

  /// Accessors for getting the global ID of the patch containing
  /// a particular element.

  int globalID(const Loc<Dim> &loc) const;
  int globalID(int) const;
  int globalID(int,int) const;
  int globalID(int,int,int) const;
  int globalID(int,int,int,int) const;
  int globalID(int,int,int,int,int) const;
  int globalID(int,int,int,int,int,int) const;
  int globalID(int,int,int,int,int,int,int) const;

  //@}

  /// Iterators into the fill list. These are MultiPatch's interface to
  /// the information needed to fill guard cells, which is cached here
  /// in the layout.
   
  FillIterator_t beginFillList() const
    {
      return this->gcFillList_m.begin();
    }
   
  FillIterator_t endFillList() const
    {
      return this->gcFillList_m.end();
    }

  //============================================================
  // touches operations
  //============================================================

  /// Find all subdomains that touch on a given domain, and insert the
  /// intersection of these subdomains into the given output iterator.
  /// Return the number of touching elements. This version of touches
  /// can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &fulld, OutIter o,
	      const ConstructTag &ctag) const;

  /// Find all subdomains that touch on a given domain, and insert the
  /// intersection of these subdomains into the given output iterator.
  /// Return the number of touching elements. This version of touches
  /// can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesAlloc(const OtherDomain &fulld, OutIter o,
		   const ConstructTag &ctag) const;

  void sync();

  template<class Out>
  void print(Out & ostr);

private:
  //============================================================
  // Private Methods
  //============================================================

  /// This function calculates the cached guard-cell filling information.
  
  void calcGCFillList();

  /// This function recalculates what the total domain of each patch
  /// should be, since this can change due to dynamic operations.

  void calcDomains();

  /// This function recalculates the domain maps, since this can
  /// change due to dynamic operations.

  void calcMaps() const;
  void calcAllocMaps() const;

  /// compute the serial index of the given block

  inline int blockIndex(const Loc<Dim> &loc) const
    {
      int pos = loc[0].first();
      for (int i=1; i < Dim; ++i)
	pos += loc[i].first() * blockStrides_m[i];
      return pos;
    }

  //============================================================
  // Private Data
  //============================================================

//   // copy of the Grid<Dim> object that can be used to describe
//   // the layout. 
// **** Not currently stored, however it seems to me that we want to
// store this information, as with it, and the info in LayoutBase
// et al, the layout can be easily reconstructed. 
//    Grid<Dim> grid_m;


  /// Is the current layout information out-of-date (due to dynamic ops
  /// that have yet to have sync called?)

  bool dirtyLayout_m;

  /// Strides used to convert block i,j,k indices to serialized indices.

  int blockStrides_m[Dim];

  /// This domain map is created for each axis. For each axis, the domain
  /// is tiled exactly, with no overlap by each DomanMapNode. 

  mutable DomainMap<Interval<1>,AxisIndex_t> map_m[Dim];

  /// This DM is initialized with the real domain, extended at either end by
  /// the external guard cell size. Internally, the DMN's overlap by the amount
  /// of the internal guard cell layers. 
  ///
  /// This will only yield ordered touches
  /// (where the touch operation returns iterators that deref to an increasing
  /// fortran style index) if the domain on which the grid is built is uniformly
  /// increasing or decreasing (i.e. Grid(Interval<1>(0,5),Range<1>(8,4,-2)) may
  /// not yield ordered results.)

  mutable DomainMap<Interval<1>,AxisIndex_t> mapAloc_m[Dim]; 
};


/**
 * GridLayout is a Layout class that breaks an N-dimensional
 * Interval into sub-domains arranged in an N-dimensional
 * grid, where the sub-domain sizes are specified by a Grid domain object
 *
 * This is an alternative to the more general tile Layout class that should
 * perform faster since subdomains can be found using a set of 1-dimensional
 * domainMap's, rather than by a more general search.
 *
 * To construct a GridLayout, you can do any of the following:
 *   -# provide a global domain, and let the GridLayout perform
 *      its default partitioning by just using one single block;
 *   -# provide a global domain, a Loc with the number of blocks to
 *      use along each dimension, and an optional context number;
 *   -# provide a global domain and a GridPartition  or UniformGridPartition 
 *      object.
 *   -# provide a Grid domain object
 *   -# provide a Grid domain object and a GuardLayers objects to create
 *      the grid with external and/or internal guard cells. 
 */      

template <int Dim>
class GridLayout : public LayoutBase<Dim,GridLayoutData<Dim> >,
                   public Observable<GridLayout<Dim> >,
                   public Observer<GridLayoutData<Dim> >
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef GridLayout<Dim>                      This_t;
  typedef Observable<This_t>                   Observable_t;
  typedef GridLayoutData<Dim>                  LayoutData_t;
  typedef typename LayoutData_t::Domain_t      Domain_t;
  typedef typename LayoutData_t::BaseDomain_t  BaseDomain_t;
  typedef typename LayoutData_t::Context_t     Context_t;
  typedef typename LayoutData_t::ID_t          ID_t;
  typedef typename LayoutData_t::Value_t       Value_t;
  typedef typename LayoutData_t::List_t        List_t;
  typedef DynamicEvents::PatchID_t             PatchID_t;
  typedef DynamicEvents::CreateSize_t          CreateSize_t;
  typedef GuardLayers<Dim>                     GuardLayers_t;

  /// Iterator through nodes. Basically the same as the vector iterator
  /// except it dereferences automatically.  
  
  typedef DerefIterator<Value_t>               iterator;
  typedef ConstDerefIterator<Value_t>          const_iterator;
   
  /// Iterator through guard-cell-fill requests. 
  
  typedef typename LayoutData_t::GCFillInfo_t GCFillInfo_t;
  
  typedef typename 
    std::vector<GCFillInfo_t>::const_iterator FillIterator_t;
  
  // Enumerations.

  enum { dimensions = Dim };
  enum { repartitionEvent = 1 };
  enum { dynamic = true };

  
  //============================================================
  // Constructors
  //============================================================

  /// The default constructor does not initialize the layout.  In this
  /// case, layout initialization must be completed with the
  /// 'initialize' method before the layout can be used.  A default
  /// layout has an empty global domain, and empty subdomain lists.
  
  GridLayout();

  /// Construct a layout with nothing else but a global domain.  In
  /// this case, a default partitioner will be used, the GridPartition
  /// object, which will just make a grid with one block and no guard cells.
  
  GridLayout(const Domain_t &,const DistributedTag &);

  GridLayout(const Domain_t &,const ReplicatedTag &);

  /// Same, but also specify guard cell info

  GridLayout(const Domain_t &,
             const GuardLayers_t &,const DistributedTag &);

  GridLayout(const Domain_t &,
             const GuardLayers_t &,const ReplicatedTag &);

  /// Domain + block count constructors.
  /// In this case, the user specifies the global domain and the number
  /// of blocks in each dimension, which will cause the domain to be
  /// partitioned in a uniform manner. One or two GuardLayers objects may
  /// be specified. If only one is specified, it is used for both
  /// internal and external guards.
  GridLayout(const Domain_t &,
	     const Loc<Dim> &,const DistributedTag &);

  GridLayout(const Domain_t &,
	     const Loc<Dim> &, 
	     const GuardLayers_t &,const DistributedTag &);

  GridLayout(const Domain_t &,
	     const Loc<Dim> &, 
	     const GuardLayers_t &,
	     const GuardLayers_t &,const DistributedTag &);

  GridLayout(const Domain_t &,
	     const Loc<Dim> &,const ReplicatedTag &);

  GridLayout(const Domain_t &,
	     const Loc<Dim> &, 
	     const GuardLayers_t &,const ReplicatedTag &);

  GridLayout(const Domain_t &,
	     const Loc<Dim> &, 
	     const GuardLayers_t &,
	     const GuardLayers_t &,const ReplicatedTag &);

  /// Domain + Grid Domain constructors. 
  /// Like the block count constructors, but a Grid object is specified,
  /// allowing a general non-uniform grid to be created.  Again, one or two
  /// GuardLayers objects may be specified. If only one is specified, it is
  /// used for both internal and external guards.
  
  GridLayout(const Grid<Dim> &,
	     const DistributedTag &);

  /// Internal guard layers same as external guard layers, specified by gcs:
  GridLayout(const Grid<Dim> &,
	     const GuardLayers_t &,
	     const DistributedTag &);

  /// Independently-specified internal and external guard layers:
  GridLayout(const Grid<Dim> &,
	     const GuardLayers_t &,
	     const GuardLayers_t &,
	     const DistributedTag &);

  GridLayout(const Grid<Dim> &,
	     const ReplicatedTag &);

  /// Internal guard layers same as external guard layers, specified by gcs:
  GridLayout(const Grid<Dim> &,
	     const GuardLayers_t &,
	     const ReplicatedTag &);

  /// Independently-specified internal and external guard layers:
  GridLayout(const Grid<Dim> &,
	     const GuardLayers_t &,
	     const GuardLayers_t &,
	     const ReplicatedTag &);

  template <class Partitioner>
  GridLayout(const Domain_t &,
	     const Partitioner &,
	     const DistributedTag &);

  template <class Partitioner>
  GridLayout(const Domain_t &,
	     const Partitioner &,
	     const ReplicatedTag &);
  
  /// Domain + partition + mapper constructor. 

  template <class Partitioner>
  GridLayout(const Domain_t &,
	     const Partitioner &,
	     const ContextMapper<Dim> &);

  /// Copy constructor & assignment operator
  /// Shallow copies with reference counting.

  GridLayout(const This_t &);

  This_t &operator=(const This_t &);
  

  //============================================================
  // Destructor
  //============================================================

  /// The actual data will be cleaned up by the LayoutData_t destructor
  /// if all references to the data go away.
  /// If any Observers remain, they will be notified by the Observable 
  /// destructor. 
  
  inline ~GridLayout() 
    { 
      this->pdata_m->detach(*this);
    }


  //============================================================
  // Initialize methods
  //============================================================

  /// Initialize a layout with nothing else but a global domain.  In
  /// this case, a default partitioner will be used, the GridPartition
  /// object, which will try to make a grid with one block.
  
  void initialize(const Domain_t &,
		  const DistributedTag &);

  void initialize(const Domain_t &,
		  const GuardLayers_t &,
		  const DistributedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &,
		  const DistributedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &,
		  const GuardLayers_t &,
		  const DistributedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &, 
		  const GuardLayers_t &,
		  const GuardLayers_t &,
		  const DistributedTag &);

  void initialize(const Grid<Dim> &,
		  const DistributedTag &);

  void initialize(const Grid<Dim> &,
		  const GuardLayers_t &,
		  const DistributedTag &);

  void initialize(const Grid<Dim> &,
		  const GuardLayers_t &,
		  const GuardLayers_t &,
		  const DistributedTag &);

  template <class Partitioner>
  void initialize(const Domain_t &,
		  const Partitioner &,
		  const DistributedTag &);

  // ReplicatedTag
  void initialize(const Domain_t &,
		  const ReplicatedTag &);

  void initialize(const Domain_t &,
		  const GuardLayers_t &,
		  const ReplicatedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &,
		  const ReplicatedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &,
		  const GuardLayers_t &,
		  const ReplicatedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &, 
		  const GuardLayers_t &,
		  const GuardLayers_t &,
		  const ReplicatedTag &);

  void initialize(const Grid<Dim> &,
		  const ReplicatedTag &);

  void initialize(const Grid<Dim> &,
		  const GuardLayers_t &,
		  const ReplicatedTag &);

  void initialize(const Grid<Dim> &,
		  const GuardLayers_t &,
		  const GuardLayers_t &,
		  const ReplicatedTag &);

  template <class Partitioner>
  void initialize(const Domain_t &,
		  const Partitioner &,
		  const ReplicatedTag &);

  /// Domain + partition + mapper initializer.

  template <class Partitioner>
  void initialize(const Domain_t &,
		  const Partitioner &,
		  const ContextMapper<Dim> &);

  /// Used by the I/O or data management system to initialize the layout based
  /// on detailed state information previously stored. Since this is specialized
  /// for the I/O system, no trailing tag is used. 

  void initialize(const Domain_t& idom,
		  const List_t& nodes,
		  const Loc<Dim>& blocks,
		  bool hasIG, bool hasEG,
		  const GuardLayers_t& ig,
		  const GuardLayers_t& eg);

  //============================================================
  // GridLayout-specific accessors
  //============================================================

  /// Return the number of blocks in each dimension, as a Loc.
  
  inline const Loc<Dim> &blocks() const 
    {
      return this->pdata_m->blocks();
    }

 
  //============================================================
  /// Repartition the layout using a new Partitioner scheme.  The
  /// initial domain lists are cleared out, the partitioner is invoked,
  /// and then all the observers are notified.  This can only be done
  /// with a partitioner that generates grid-like blocks.

  template <class Partitioner>
  inline bool repartition(const Partitioner &gp,const ContextMapper<Dim> &cm)
    { 
      this->pdata_m->initialize(this->domain(),gp,cm);
      // Can't this just be Observable_t::notify(repartitionEvent)???
      this->pdata_m->notify(repartitionEvent); 
      return true;
    }

  inline void sync()
    { 
      this->pdata_m->sync();
    }


  //============================================================
  // Observer methods
  //============================================================

  /// Respond to events generated by the LayoutData_t.
  /// These are just passed on to our observers.

  virtual void notify(LayoutData_t &d, const ObserverEvent &event)
    {
      // We should only get this message from our LayoutData_t object

      PAssert(&d == this->pdata_m.rawPointer());
      Observable_t::notify(event);
    }


  //============================================================
  // Output
  //============================================================

  /// Print a GridLayout on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const;

#if !POOMA_NO_TEMPLATE_FRIENDS

  //private:

  template <int Dim1, int Dim2>
  friend class GridLayoutView;

#endif

  //============================================================
  // Data
  //============================================================

  /// GridLayout stores its data in a RefCounted class to
  /// simplify memory management.
  
  friend class GridLayoutData<Dim>;

};


/**
 * The data object held by a GridLayoutView object.
 */

template <int Dim, int Dim2>
class GridLayoutViewData : 
                public LayoutBaseViewData<Dim, Dim2, GridLayout<Dim2> >,
                public RefCounted
{
public:

  typedef GridLayout<Dim2>                  Layout_t;
  typedef GridLayoutView<Dim, Dim2>         ViewLayout_t;

  typedef Interval<Dim>                     Domain_t;
  typedef Range<Dim2>                       BaseDomain_t;
  typedef int                               Context_t;
  typedef Unique::Value_t                   ID_t;

  typedef typename Layout_t::Domain_t       AllocatedDomain_t;
  typedef ViewIndexer<Dim,Dim2>             Indexer_t;

  typedef Node<Domain_t,AllocatedDomain_t>  Value_t;
  typedef std::vector<Value_t *>            List_t;        // for convenience
  typedef GuardLayers<Dim>                  GuardLayers_t; // for convenience

  typedef GridLayoutViewData<Dim,Dim2>      LayoutData_t;
  
  // Enumerations.

  enum { dim = Dim };
  enum { dim2 = Dim2 };
  

  //============================================================
  // Constructors
  //============================================================

  GridLayoutViewData() { }
  
  template <class DT>
  inline GridLayoutViewData(const Layout_t &layout, const Domain<Dim, DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,GridLayout<Dim2> >(layout,dom)
    // layout_m(layout), indexer_m(dom), 
    //id_m(Unique::get()), subdomainsComputed_m(false),
    //internalGuards_m(layout.internalGuards()),
    //externalGuards_m(layout.externalGuards())
  { 
  }

  template <class DT>
  inline GridLayoutViewData(const Layout_t &layout, 
                            const SliceDomain<DT> &dom)
  :LayoutBaseViewData<Dim,Dim2,GridLayout<Dim2> >(layout,dom)
  { 
  }

  template <class DT>
  GridLayoutViewData(const ViewLayout_t &layout, 
                     const Domain<Dim, DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,GridLayout<Dim2> >(
					      layout.pdata_m->layout_m,
					      layout,
					      layout.pdata_m->indexer_m, 
					      dom,
					      layout.internalGuards(),
					      layout.externalGuards())
  {
  }

  template <int OrigDim, class DT>
  GridLayoutViewData(const GridLayoutView<OrigDim, Dim2> &layout, 
                     const SliceDomain<DT> &dom)
    : LayoutBaseViewData<Dim,Dim2,GridLayout<Dim2> >(
			     layout.pdata_m->layout_m,
			     layout,
			     Indexer_t(layout.pdata_m->indexer_m,dom),
			     dom)
  { 
  }

  // Destructor

  ~GridLayoutViewData() 
  {
    typename List_t::iterator a;
    for (a = this->all_m.begin(); a != this->all_m.end(); ++a)
      delete (*a);
  }

};

/**
 * GridLayoutView is a Layout class that provides a view of an
 * existing GridLayout object. Dim is the logical dimension of
 * the layout. Dim2 is the dimension of the GridLayout
 * contained within.
 *
 * To construct a GridLayoutView, you need an existing
 * GridLayout or a GridLayoutView and the subdomain that
 * is being viewed. This class does not have a useful default
 * constructor since it is based on an existing GridLayout.
 *
 * Once created, GridLayoutView has the same interface as
 * Layout (see Layout.h). It also provides this extra interface:
 *
 * int globalID(const Loc<Dim> &pos) : return the globalID
 *    of the node that contains the point.
 */

template <int Dim, int Dim2>
class GridLayoutView  
: public LayoutBaseView<Dim, Dim2, GridLayoutViewData<Dim,Dim2> >
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.
  enum { dimensions = Dim };
  enum { dim = Dim };
  enum { dim2 = Dim2 };

  // General public typedefs.

  typedef GridLayoutViewData<Dim, Dim2>            LayoutData_t; 

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

  typedef GridLayoutView<Dim, Dim2>                This_t;       // convenience
  typedef GridLayoutView<Dim, Dim2>                ViewLayout_t; // convenience
  typedef LayoutBaseView<Dim, Dim2, LayoutData_t>  Base_t;
  
  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.
  
  typedef DerefIterator<Value_t>                   iterator;
  typedef ConstDerefIterator<Value_t>              const_iterator;


  //============================================================
  // Constructors
  //============================================================

  // Default constructor.
  
  GridLayoutView()
    : Base_t(new LayoutData_t())
  { }
  
  // Constructor building a GridLayoutView from a
  // GridLayout and a non-slice domain like an Interval<Dim> or
  // Range<Dim>.
  
  template <class DT>
  GridLayoutView(const Layout_t &layout, const Domain<Dim2, DT> &dom)
    : LayoutBaseView<Dim,Dim2,GridLayoutViewData<Dim,Dim2> >
  (new GridLayoutViewData<Dim,Dim2>(layout,dom)) 
  { }

  // Constructor building a GridLayoutView from a
  // GridLayout and a slice domain like a
  // SliceInterval<Dim2,Dim> or SliceRange<Dim2,Dim>.
  
  template <class DT>
  GridLayoutView(const Layout_t &layout, const SliceDomain<DT> &dom)
    : LayoutBaseView<Dim,Dim2,GridLayoutViewData<Dim,Dim2> >
  (new GridLayoutViewData<Dim,Dim2>(layout,dom))  
  { }
  
  // Constructor building a GridLayoutView from another
  // GridLayoutView and a non-slice domain like an
  // Interval<Dim> or Range<Dim>.
  
  template <class DT>
  GridLayoutView(const ViewLayout_t &layout, const Domain<Dim, DT> &dom)
 : LayoutBaseView<Dim,Dim2,GridLayoutViewData<Dim,Dim2> >
  (new GridLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Constructor building a GridLayoutView from another
  // GridLayoutView and a slice domain like a
  // SliceInterval<Dim2,Dim> or SliceRange<Dim2,Dim>.
  
  template <int OldViewDim, class DT>
  GridLayoutView(const GridLayoutView<OldViewDim, Dim2> &layout, 
                        const SliceDomain<DT> &dom)
  : LayoutBaseView<Dim,Dim2,GridLayoutViewData<Dim,Dim2> > 
  (new GridLayoutViewData<Dim,Dim2>(layout,dom)) 
  { }

  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
  inline GridLayoutView(const This_t &model) 
    : LayoutBaseView<Dim,Dim2,GridLayoutViewData<Dim,Dim2> >(model.pdata_m)
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

  // The actual data will be cleaned up by the GridLayoutData
  // destructor if all references to the data go away, so there is
  // nothing to do here.
   
  inline ~GridLayoutView() 
  { }

  //============================================================
  // Output
  //============================================================
    
  // Print a GridLayoutView on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const;

#if !POOMA_NO_TEMPLATE_FRIENDS

  //private:

  template <int OtherDim, int OtherDim2>
  friend class GridLayoutView;

  template <int OtherDim, int OtherDim2>
  friend class GridLayoutViewData;

#endif

  //============================================================
  // Private utility functions
  //============================================================

  // Fill our subdomain lists.
  
  void computeSubdomains() const { this->pdata_m->computeSubdomains(); }
   
  //============================================================
  // Data
  //============================================================

  // The data is stored in a RefCounted class to simplify memory
  // management.  This is probably not as important for ViewLayout
  // classes as for Layout classes, but we do it for
  // consistency. Currently ViewLayouts are not observed directly by
  // anyone. Of course, the Layout that we have a copy of is observed.

};

//============================================================
// NewDomain1 traits classes for GridLayout and GridLayoutView
//============================================================

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a GridLayout.
//
//-----------------------------------------------------------------------------

template <int Dim>
struct NewDomain1<GridLayout<Dim> >
{
  typedef GridLayout<Dim> &Type_t;

  inline static Type_t combine(const GridLayout<Dim> &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a GridLayoutView.
//
//-----------------------------------------------------------------------------

template <int Dim, int Dim2>
struct NewDomain1<GridLayoutView<Dim, Dim2> >
{
  typedef GridLayoutView<Dim, Dim2> &Type_t;

  inline static Type_t combine(const GridLayoutView<Dim, Dim2> &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// ostream inserters for GridLayout and GridLayoutView:
//
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &ostr, 
			 const GridLayout<Dim> &layout)
{
  layout.print(ostr);
  return ostr;
}

template <int Dim, int Dim2>
std::ostream &operator<<(std::ostream &ostr, 
			 const GridLayoutView<Dim, Dim2> &layout)
{
  layout.print(ostr);
  return ostr;
}


// } // namespace POOMA

#include "Layout/GridLayout.cpp"

#endif // POOMA_LAYOUT_GRIDLAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: GridLayout.h,v $   $Author: richi $
// $Revision: 1.112 $   $Date: 2004/11/10 22:02:10 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
