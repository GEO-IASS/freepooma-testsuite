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
//   UniformGridLayout<Dim>
//   UniformGridLayoutView<Dim, Dim2>
//   UniformTag
//   MultiPatchLayoutTraits<UniformTag,Dim>
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_UNIFORMGRIDLAYOUT_H
#define POOMA_LAYOUT_UNIFORMGRIDLAYOUT_H

/** @file
 * @ingroup Layout
 * @brief
 *   UniformGridLayout<Dim>
 *     - Layout class that breaks Dim-dimensional domain into equal
 *       sized sub-domains arranged in a Dim-dimensional grid.
 *   UniformGridLayoutView<Dim, Dim2>
 *     - view of a UniformGridLayout
 *   UniformTag
 *     - tag used to specialize MultiPatchLayoutTraits
 *   MultiPatchLayoutTraits<UniformTag,Dim>
 *     - traits class used by MultiPatch-engine to determine layout type.
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

#include "Layout/LayoutBase.h"

#include <vector>
#include <iosfwd>


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//============================================================
// Forward declaration
//============================================================

template <int Dim> class UniformGridLayout;
template <int Dim, int Dim2> class UniformGridLayoutView;

//-----------------------------------------------------------------------------
// Full Description:
// UniformTag
//
// Tag class.
//-----------------------------------------------------------------------------

struct UniformTag { };


/**
 * Specialization of MultiPatchLayoutTraits for UniformGridLayout.
 */

template <int Dim>
struct MultiPatchLayoutTraits<UniformTag,Dim>
{
  typedef UniformGridLayout<Dim> Layout_t;
  
  template <int ViewDim>
  struct View
  {
    typedef UniformGridLayoutView<ViewDim,Dim> Layout_t;
  };
};


/**
 * This is the actual data for the UniformGridLayout class, which is
 * simply a wrapper that holds a reference counted instance of this
 * data class.
 */

template <int Dim>
class UniformGridLayoutData 
 : public LayoutBaseData<Dim>, 
   public RefCounted, 
   public Observable<UniformGridLayoutData<Dim> >
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
  
  typedef typename LayoutBaseData<Dim>::GCFillInfo GCFillInfo_t;
  
  typedef typename std::vector<GCFillInfo_t>::const_iterator FillIterator_t;
  
  // Enumerations.

  enum { dimensions = Dim };
  enum { repartitionEvent = 1 };
  enum { dynamic = false };


  //============================================================
  // Constructors
  //============================================================

  // Default version creates an empty domain.
  
  UniformGridLayoutData();
  
  // All initialization is done via the following constructor.
  // Originally we had several constructors, but once we went
  // to not storing a partitioner object, it was easier to just
  // have the Layout construct the partitioner appropriately
  // and pass it to the LayoutData constructor.
  
  template <class Partitioner>
  UniformGridLayoutData(const Domain_t &gdom, 
			const Partitioner &gpar,
			const ContextMapper<Dim> & cmap );


  //============================================================
  // Special I/O initializer
  //============================================================
  // Used by the I/O or data management system to initialize the layout based
  // on detailed state information previously stored.

  void initialize(const Domain_t& idom,
		  const List_t& nodes,
		  const Loc<Dim>& blocks,
		  bool hasIG, bool hasEG,
		  const GuardLayers_t& ig,
		  const GuardLayers_t& eg);


  //============================================================
  // Destructor
  //============================================================

  // This is only called when all references to this data go away, in
  // which case we need to delete our nodes. The Observable destructor
  // will broadcast messages up to all observers of the Layout.

  ~UniformGridLayoutData()
  {
    typename List_t::iterator a;
    for (a = this->all_m.begin(); a != this->all_m.end(); ++a)
      delete (*a);
  }


  // Find all subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.
  // Return the number of touching elements. This version of touches
  // can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, const ConstructTag &ctag) const;

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesLocal(const OtherDomain &d, 
		   OutIter o, 
		   const ConstructTag &ctag) const;

  // Find all Remote subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.
  // Return the number of touching elements. This version of touches
  // can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesRemote(const OtherDomain &d, 
		    OutIter o, 
		    const ConstructTag &ctag) const;

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesAlloc(const OtherDomain &d, OutIter o, 
                   const ConstructTag &ctag) const;

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesAllocLocal(const OtherDomain &d, OutIter o, 
			const ConstructTag &ctag) const;

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touchesAllocRemote(const OtherDomain &d, OutIter o, 
			 const ConstructTag &ctag) const;

  //private:

  friend class UniformGridLayout<Dim>;

  //============================================================
  // Accessors
  //============================================================


  // These are deferred to UniformGridLayoutData in order to avoid a
  // bunch of indirections in their implementation.

  // Accessors for getting the global ID of the patch containing a
  // particular element.

  int globalID(const Loc<Dim> &loc) const;
  int globalID(int) const;
  int globalID(int,int) const;
  int globalID(int,int,int) const;
  int globalID(int,int,int,int) const;
  int globalID(int,int,int,int,int) const;
  int globalID(int,int,int,int,int,int) const;
  int globalID(int,int,int,int,int,int,int) const;
  
  // Mutators used for partitioning.

  template <class Partitioner>
  void partition(const Partitioner &, const ContextMapper<Dim>& cmap);

  template <class Partitioner>
  bool repartition(const Partitioner &,const ContextMapper<Dim>&);
  
  // This function calculates the cached guard-cell filling information.
  
  void calcGCFillList();
  
  
  //============================================================
  // Data
  //============================================================

  // The stride array for indexing into the 1-D list of blocks;
  // i.e. in 2D, block (i,j) is element (i + j * blockstride_m[1]) in
  // the list.
  
  int blockstride_m[Dim];

  // The patch size. Could be stored in an interval, but use an array
  // since this info is zero-based.
  
  int blocksizes_m[Dim];

  // The domain of the "brick" of patches stored in the all_m list.
  
  Interval<Dim> allDomain_m;
};


/**
 * UniformGridLayout is a Layout class that breaks an N-dimensional
 * Interval into equal sized sub-domains arranged in an N-dimensional
 * grid.
 *
 * This is an alternative to the more general Layout class that should
 * perform somewhat faster since subdomains can be found
 * arithmetically, rather than via a search.  It is only able to
 * represent grid-like layouts, however.
 *
 * To construct a UniformGridLayout, you can do any of the following:
 *   -# provide a global domain, and let the UniformGridLayout perform
 *      its default partitioning by just using one single block;
 *   -# provide a global domain, a Loc with the number of blocks to
 *      use along each dimension
 *   -# provide a global domain and a UniformGridPartition object.
 *      
 * Alternatively, you can use the default UniformGridLayout
 * constructor, and call the 'initialize' method later with the same
 * possible set of arguments.
 *
 * You can also specify internal and external guard layers for the
 * domains. See the comments in UniformGridLayout below.
 */

template <int Dim>
class UniformGridLayout : public LayoutBase<Dim,UniformGridLayoutData<Dim> >,
                          public Observable<UniformGridLayout<Dim> >,
                          public Observer<UniformGridLayoutData<Dim> >
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.
  
  typedef UniformGridLayoutData<Dim>           LayoutData_t; // for convenience
  typedef typename LayoutData_t::Domain_t      Domain_t;
  typedef typename LayoutData_t::BaseDomain_t  BaseDomain_t;
  typedef typename LayoutData_t::Context_t     Context_t;
  typedef typename LayoutData_t::ID_t          ID_t;
  typedef typename LayoutData_t::Value_t       Value_t;
  typedef typename LayoutData_t::List_t        List_t; // for convenience
  typedef UniformGridLayout<Dim>               This_t; // for convenience
  typedef Observable<This_t>                   Observable_t;

  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.
  
  typedef DerefIterator<Value_t>               iterator;
  typedef ConstDerefIterator<Value_t>          const_iterator;
  
  // Iterator through guard-cell-fill requests. 
  
  typedef typename LayoutData_t::GCFillInfo_t GCFillInfo_t;
  
  typedef typename 
    std::vector<GCFillInfo_t>::const_iterator FillIterator_t;
  
  typedef GuardLayers<Dim> GuardLayers_t;

  // Enumerations.

  enum { dimensions = Dim };
  enum { repartitionEvent = LayoutData_t::repartitionEvent };
  enum { dynamic = false };


  //============================================================
  // Constructors
  //============================================================

  /**
  // The default constructor does not initialize the layout.  In this
  // case, layout initialization must be completed with the
  // 'initialize' method before the layout can be used.  A default
  // layout has an empty global domain, and empty subdomain lists.
  
  // This is also the only constructor that doesn't demand either
  // ReplicatedTag or DistributedTag
  */
  
  UniformGridLayout();

  // Distributed versions


  UniformGridLayout(const Domain_t &,
		    const DistributedTag &);

  UniformGridLayout(const Domain_t &, 
                    const GuardLayers_t &,
		    const DistributedTag &);

  UniformGridLayout(const Domain_t &, 
                    const Loc<Dim> &,
		    const DistributedTag &);

  UniformGridLayout(const Domain_t &, 
                    const Loc<Dim> &, 
                    const GuardLayers_t &,
		    const DistributedTag &);

  UniformGridLayout(const Domain_t &, 
                    const Loc<Dim> &, 
                    const GuardLayers_t &, 
                    const GuardLayers_t &,
		    const DistributedTag &);
  //========================= Replicated versions ===================

  UniformGridLayout(const Domain_t &,
		    const ReplicatedTag &);

  UniformGridLayout(const Domain_t &, 
                    const GuardLayers_t &,
		    const ReplicatedTag &);

  UniformGridLayout(const Domain_t &, 
                    const Loc<Dim> &,
		    const ReplicatedTag &);

  UniformGridLayout(const Domain_t &, 
                    const Loc<Dim> &, 
                    const GuardLayers_t &,
		    const ReplicatedTag &);

  UniformGridLayout(const Domain_t &, 
                    const Loc<Dim> &, 
                    const GuardLayers_t &, 
                    const GuardLayers_t &,
		    const ReplicatedTag &);

  // Domain + partition constructor. 
  // In this case, the partitioner must be a UniformGridPartition
  // object. The partitioner must know about guard cells, so we don't
  // specify them separately here.
  
  template <class Partitioner>
  UniformGridLayout(const Domain_t &,
		    const Partitioner &,
		    const ContextMapper<Dim> & );

  template <class Partitioner>
  UniformGridLayout(const Domain_t &, 
                    const Partitioner &,
		    const DistributedTag &);

  template <class Partitioner>
  UniformGridLayout(const Domain_t &, 
                    const Partitioner &,
		    const ReplicatedTag &);

  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
  UniformGridLayout(const This_t &);
  
  This_t &operator=(const This_t &);

  // I think we need a version of this here ... just forward to
  // the RefCountedPtr.
  
  // makeOwnCopy???


  //============================================================
  // Destructor
  //============================================================

  // The actual data will be cleaned up by the LayoutData_t destructor
  // if all references to the data go away.  If any Observers remain,
  // they will be notified by the Observable destructor.
  
  inline ~UniformGridLayout() 
  { 
    this->pdata_m->detach(*this);
  }


  //============================================================
  // Initialize methods
  //============================================================

  // Initialize a layout with nothing else but a global domain.  In
  // this case, a default partitioner will be used, the
  // UniformGridPartition object, which will try to make a grid with
  // one block.  
  
  void initialize(const Domain_t &,
		  const DistributedTag &);
  
  void initialize(const Domain_t &,
		  const ReplicatedTag &);
  
  void initialize(const Domain_t &, 
		  const GuardLayers_t &,
		  const DistributedTag &);
  
  void initialize(const Domain_t &, 
		  const GuardLayers_t &,
		  const ReplicatedTag &);
  
  void initialize(const Domain_t &, 
		  const Loc<Dim> &,
		  const DistributedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &,
		  const ReplicatedTag &);

  void initialize(const Domain_t &, 
		  const Loc<Dim> &, 
		  const GuardLayers_t &,
		  const DistributedTag & );

  void initialize(const Domain_t &, 
		  const Loc<Dim> &, 
		  const GuardLayers_t &,
		  const ReplicatedTag & );

  void initialize(const Domain_t &, 
                  const Loc<Dim> &, 
                  const GuardLayers_t &, 
                  const GuardLayers_t &,
		  const DistributedTag &);

  void initialize(const Domain_t &, 
                  const Loc<Dim> &, 
                  const GuardLayers_t &, 
                  const GuardLayers_t &,
		  const ReplicatedTag &);

  // Domain + partition initializer.
  // In this case, the partitioner must be a UniformGridPartition
  // object.
  
  template <class Partitioner>
  void initialize(const Domain_t &, 
		  const Partitioner &,
		  const DistributedTag &);

  template <class Partitioner>
  void initialize(const Domain_t &, 
		  const Partitioner &,
		  const ReplicatedTag &);
  
  // Domain + partition + mapper  initializer.
  // In this case, the partitioner must be a UniformGridPartition
  // object.
  
  template <class Partitioner>
  void initialize(const Domain_t &, 
		  const Partitioner &,
		  const ContextMapper<Dim> &);
 
  // Used by the I/O or data management system to initialize the layout based
  // on detailed state information previously stored. Since this is specialized
  // for the I/O system, no trailing tag is used. 

  void initialize(const Domain_t& idom,
		  const List_t& nodes,
		  const Loc<Dim>& blocks,
		  bool hasIG, bool hasEG,
		  const GuardLayers_t& ig,
		  const GuardLayers_t& eg);
  //============================================================
  // Accessors
  //============================================================




  //============================================================
  // Observer methods
  //============================================================

  // Respond to events generated by the LayoutData_t.
  // These are just passed on to our observers.
  // this needs to be here, and not in the base class

  virtual void notify(LayoutData_t &d, const ObserverEvent &event)
  {
    // We should only get this message from our LayoutData_t object
    PAssert(&d == this->pdata_m.rawPointer());
    Observable_t::notify(event);
  }
    

  //============================================================
  // Output
  //============================================================
    
  // Print a UniformGridLayout on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const {
    ostr << "UniformGridLayout " << this->ID() << " on global domain " 
      << this->domain() << ":" << '\n';
    ostr << "   Total subdomains: " << this->sizeGlobal() << '\n';
    ostr << "   Local subdomains: " << this->sizeLocal() << '\n';
    ostr << "  Remote subdomains: " << this->sizeRemote() << '\n';
    ostr << "        Grid blocks: " << this->blocks() << '\n';
    typename UniformGridLayout<Dim>::const_iterator a;
    for (a = this->beginGlobal(); a != this->endGlobal(); ++a)
      ostr << "  Global subdomain = " << *a << '\n';
    for (a = this->beginLocal(); a != this->endLocal(); ++a)
      ostr << "   Local subdomain = " << *a << '\n';
    for (a = this->beginRemote(); a != this->endRemote(); ++a)
      ostr << "  Remote subdomain = " << *a << '\n';
  }


#if !POOMA_NO_TEMPLATE_FRIENDS
  // can't seem to get specialized template friends working...
  //private:

  template <int Dim1, int Dim2>
  friend class UniformGridLayoutView;

#endif
  
  //============================================================
  // Data
  //===========================================================

  // UniformGridLayout stores its data in a LayoutBase::RefCounted class to
  // simplify memory management.
  
  friend class UniformGridLayoutData<Dim>;

};



/**
 * This is the actual data for the UniformGridLayoutView class, which is
 * simply a wrapper that holds a reference counted instance of this
 * data class.
 */

template <int Dim, int Dim2>
class UniformGridLayoutViewData 
  : public LayoutBaseViewData<Dim, Dim2, UniformGridLayout<Dim2> >,
    public RefCounted
{
public:

  typedef UniformGridLayout<Dim2>              Layout_t;
  typedef UniformGridLayoutView<Dim, Dim2>     ViewLayout_t;
  typedef LayoutBaseViewData<Dim,Dim2,Layout_t> Base_t;
  typedef Interval<Dim>                        Domain_t;
  typedef Range<Dim2>                          BaseDomain_t;
  typedef int                                  Context_t;
  typedef Unique::Value_t                      ID_t;

  typedef typename Layout_t::Domain_t          AllocatedDomain_t;
  typedef ViewIndexer<Dim,Dim2>                Indexer_t;

  typedef Node<Domain_t,AllocatedDomain_t>     Value_t;
  typedef std::vector<Value_t *>               List_t;        // convenience
  typedef GuardLayers<Dim>                     GuardLayers_t; // convenience


  // Enumerations.

  enum { dim = Dim };
  enum { dim2 = Dim2};

  //============================================================
  // Constructors
  //============================================================

  UniformGridLayoutViewData() { };
  
  template <class DT>
  inline 
  UniformGridLayoutViewData(const Layout_t &layout, const Domain<Dim, DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,UniformGridLayout<Dim2> >(layout,dom)
  { 
  }

  template <class DT>
  inline 
  UniformGridLayoutViewData(const Layout_t &layout, const SliceDomain<DT> &dom)
  :LayoutBaseViewData<Dim,Dim2,UniformGridLayout<Dim2> >(layout,dom)
  {  
  }

  template <class DT>
  UniformGridLayoutViewData(const ViewLayout_t &layout, 
                            const Domain<Dim, DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,UniformGridLayout<Dim2> >(
					      layout.pdata_m->layout_m,
					      layout,
					      layout.pdata_m->indexer_m, 
					      dom,
					      layout.internalGuards(),
					      layout.externalGuards())
  {
  }

  template <int OrigDim, class DT>
  UniformGridLayoutViewData(const UniformGridLayoutView<OrigDim,Dim2> &layout, 
                            const SliceDomain<DT> &dom)
  : LayoutBaseViewData<Dim,Dim2,UniformGridLayout<Dim2> >(
			     layout.pdata_m->layout_m,
			     layout,
			     Indexer_t(layout.pdata_m->indexer_m,dom),
			     dom)
  { 
  }

  // Destructor

  ~UniformGridLayoutViewData() 
  {
    typename List_t::iterator a;
    for (a = this->all_m.begin(); a != this->all_m.end(); ++a)
      delete (*a);
  }

};


/**
 * UniformGridLayoutView is a Layout class that provides a view of an
 * existing UniformGridLayout object. Dim is the logical dimension of
 * the layout. Dim2 is the dimension of the UniformGridLayout
 * contained within.
 *
 * To construct a UniformGridLayoutView, you need an existing
 * UniformGridLayout or a UniformGridLayoutView and the subdomain that
 * is being viewed. This class does not have a useful default
 * constructor since it is based on an existing UniformGridLayout.
 *
 * Once created, UniformGridLayoutView has the same interface as
 * Layout (see Layout.h). It also provides this extra interface:
 *
 * int globalID(const Loc<Dim> &pos) : return the globalID
 *    of the node that contains the point.
 */

template <int Dim, int Dim2>
class UniformGridLayoutView
  : public LayoutBaseView<Dim, Dim2, UniformGridLayoutViewData<Dim,Dim2> >
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.

  enum { dimensions = Dim };

  enum { dim = Dim };
  enum { dim2 = Dim2};
  // General public typedefs.
  
  typedef UniformGridLayoutViewData<Dim, Dim2>      LayoutData_t; 

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
    
  typedef UniformGridLayoutView<Dim, Dim2>          This_t;      // convenience
  typedef UniformGridLayoutView<Dim, Dim2>          ViewLayout_t;// convenience
  typedef LayoutBaseView<Dim,Dim2,LayoutData_t>     Base_t;

  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.
  
  typedef DerefIterator<Value_t>                    iterator;
  typedef ConstDerefIterator<Value_t>               const_iterator;


  //============================================================
  // Constructors
  //============================================================

  // Default constructor. Final initialization should be done with
  // the assignment operator. 
  
  UniformGridLayoutView()
  : Base_t(new LayoutData_t())
  { }
  
  // Constructor building a UniformGridLayoutView from a
  // UniformGridLayout and a non-slice domain like an Interval<Dim> or
  // Range<Dim>.
  
  template <class DT>
  UniformGridLayoutView(const Layout_t &layout, const Domain<Dim2, DT> &dom)
  : LayoutBaseView<Dim,Dim2,UniformGridLayoutViewData<Dim,Dim2> >
    (new UniformGridLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Constructor building a UniformGridLayoutView from a
  // UniformGridLayout and a slice domain like a
  // SliceInterval<Dim2,Dim> or SliceRange<Dim2,Dim>.
  
  template <class DT>
  UniformGridLayoutView(const Layout_t &layout, const SliceDomain<DT> &dom)
  : LayoutBaseView<Dim,Dim2,UniformGridLayoutViewData<Dim,Dim2> >
    (new UniformGridLayoutViewData<Dim,Dim2>(layout,dom))
  { }
  
  // Constructor building a UniformGridLayoutView from another
  // UniformGridLayoutView and a non-slice domain like an
  // Interval<Dim> or Range<Dim>.
  
  template <class DT>
  UniformGridLayoutView(const ViewLayout_t &layout, const Domain<Dim, DT> &dom)
  : LayoutBaseView<Dim,Dim2,UniformGridLayoutViewData<Dim,Dim2> >
    (new UniformGridLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Constructor building a UniformGridLayoutView from another
  // UniformGridLayoutView and a slice domain like a
  // SliceInterval<Dim2,Dim> or SliceRange<Dim2,Dim>.
  
  template <int OldViewDim, class DT>
  UniformGridLayoutView(const UniformGridLayoutView<OldViewDim, Dim2> &layout, 
                        const SliceDomain<DT> &dom)
  : LayoutBaseView<Dim,Dim2,UniformGridLayoutViewData<Dim,Dim2> >
    (new UniformGridLayoutViewData<Dim,Dim2>(layout,dom))
  { }

  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
  inline UniformGridLayoutView(const This_t &model) 
    : LayoutBaseView<Dim,
                     Dim2,
                     UniformGridLayoutViewData<Dim,Dim2> >(model.pdata_m)
  { }
  
  inline This_t &operator=(const This_t &model) 
  {
    if (this != &model)
      {
        this->pdata_m = model.pdata_m;
      }
    return *this;
  }

  // I think we need a version of this here ... just forward to
  // the RefCountedPtr.
  
  // makeOwnCopy???



  //============================================================
  // Destructor
  //============================================================

  // The actual data will be cleaned up by the UniformGridLayoutData
  // destructor if all references to the data go away, so there is
  // nothing to do here.
   
  inline ~UniformGridLayoutView() 
  { }

  //============================================================
  // Output
  //============================================================
    
  // Print a UniformGridLayoutView on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const 
  {
    ostr << "UniformGridLayoutView " << this->ID() << " on global domain " 
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

  //protected:

  template <int OtherDim, int OtherDim2>
  friend class UniformGridLayoutView;
  
  template <int OtherDim, int OtherDim2>
  friend class UniformGridLayoutViewData;

#endif

  //============================================================
  // Data
  //============================================================

  // The data is stored in a LayoutBaseView::RefCounted class to simplify 
  // memory management. (LayoutViews have shallow copy semantics and they
  // own copies of the underlying Layout as well as their own sets of
  // Node lists, which must be cleaned up when all views are destroyed.)

};


//=============================================================================
// UniformGridLayout & UniformGridLayoutData inline method definitions
//=============================================================================

//-----------------------------------------------------------------------------
//
// Constructors and Initialize methods
//
//-----------------------------------------------------------------------------

// See comments in class definition above.

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout()
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t()),
  Observable<This_t>(*this)
{ 
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom,
		  const DistributedTag& t)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(),
		      DistributedMapper<Dim>(UniformGridPartition<Dim>()))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom,
		  const ReplicatedTag & t)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(),
		      LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const GuardLayers_t &gcs,
		  const DistributedTag &)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(gcs),
		      DistributedMapper<Dim>(UniformGridPartition<Dim>(gcs)))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const GuardLayers_t &gcs,
		  const ReplicatedTag & )
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(gcs),
		      LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const Loc<Dim> &blocks,
		  const DistributedTag & )
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(blocks),
		      DistributedMapper<Dim>(
		        UniformGridPartition<Dim>(blocks)))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const Loc<Dim> &blocks,
		  const ReplicatedTag & t)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(blocks),
		      LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
                  const Loc<Dim> &blocks, 
                  const GuardLayers_t &igcs,
		  const DistributedTag &)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
   (new LayoutData_t(gdom,
		     UniformGridPartition<Dim>(blocks,igcs),
		     DistributedMapper<Dim>(
		      UniformGridPartition<Dim>(blocks,igcs)))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
                  const Loc<Dim> &blocks, 
                  const GuardLayers_t &igcs,
		  const ReplicatedTag &)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
     (new LayoutData_t(gdom,
		       UniformGridPartition<Dim>(blocks,igcs),
		       LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
                  const Loc<Dim> &blocks, 
                  const GuardLayers_t &igcs, 
                  const GuardLayers_t &egcs,
		  const DistributedTag &)
     
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(blocks,igcs,egcs),
		      DistributedMapper<Dim>(
                       UniformGridPartition<Dim>(blocks,igcs,egcs)))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
                  const Loc<Dim> &blocks, 
                  const GuardLayers_t &igcs, 
                  const GuardLayers_t &egcs,
		  const ReplicatedTag &t)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,
		      UniformGridPartition<Dim>(blocks,igcs,egcs),
		      LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
template <class Partitioner>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const Partitioner &gpar,
		  const DistributedTag & )
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
   (new LayoutData_t(gdom,gpar,DistributedMapper<Dim>(gpar))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
template <class Partitioner>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const Partitioner &gpar,
		  const ReplicatedTag &)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
   (new LayoutData_t(gdom,gpar,LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
template <class Partitioner>
inline UniformGridLayout<Dim>::
UniformGridLayout(const Domain_t &gdom, 
		  const Partitioner &gpar,
		  const ContextMapper<Dim> & cmap)
: LayoutBase<Dim,UniformGridLayoutData<Dim> >
    (new LayoutData_t(gdom,gpar,cmap)),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
inline UniformGridLayout<Dim>::
UniformGridLayout(const This_t &model) 
: LayoutBase<Dim,UniformGridLayoutData<Dim> >(model.pdata_m),
  Observable<This_t>(*this)
{ 
   this->pdata_m->attach(*this);
}
  
template <int Dim>
inline UniformGridLayout<Dim> & UniformGridLayout<Dim>::
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

// Initialize methods...

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
	   const DistributedTag &)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.
  
  this->pdata_m->domain_m = gdom;
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(),
		     DistributedMapper<Dim>(UniformGridPartition<Dim>()));
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
	   const ReplicatedTag &)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.
  
  this->pdata_m->domain_m = gdom;
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(),
		     LocalMapper<Dim>());
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom, 
	   const GuardLayers_t &gcs,
	   const DistributedTag &)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(gcs),
		     DistributedMapper<Dim>(UniformGridPartition<Dim>(gcs) ));
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom, 
	   const GuardLayers_t &gcs,
	   const ReplicatedTag &)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(gcs),
		     LocalMapper<Dim>());
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
	   const Loc<Dim> &blocks,
	   const DistributedTag &)
{
  PAssert(!this->initialized());
  
  // Initialize our global domain, and then do the partitioning.
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(blocks),
		    DistributedMapper<Dim>(UniformGridPartition<Dim>(blocks)));
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
	   const Loc<Dim> &blocks,
	   const ReplicatedTag &)
{
  PAssert(!this->initialized());
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(blocks),
		     LocalMapper<Dim>());
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Loc<Dim> &blocks, 
           const GuardLayers_t &gcs,
	   const DistributedTag &)
{
  PAssert(!this->initialized()); 
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(blocks, gcs),
		     DistributedMapper<Dim>(
		       UniformGridPartition<Dim>(blocks, gcs)));
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Loc<Dim> &blocks, 
           const GuardLayers_t &gcs,
	   const ReplicatedTag &)
{
  PAssert(!this->initialized()); 
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(blocks, gcs),
		     LocalMapper<Dim>());
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Loc<Dim> &blocks, 
           const GuardLayers_t &igcs,
           const GuardLayers_t &egcs,
	   const DistributedTag &)
{
  PAssert(!this->initialized());
  
  // Initialize our global domain, and then do the partitioning.
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->partition(UniformGridPartition<Dim>(blocks, igcs, egcs),
		     DistributedMapper<Dim>(
		       UniformGridPartition<Dim>(blocks, igcs, egcs)));
}

template <int Dim>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Loc<Dim> &blocks, 
           const GuardLayers_t &igcs,
           const GuardLayers_t &egcs,
	   const ReplicatedTag &)
{
  PAssert(!this->initialized());
  
  // Initialize our global domain, and then do the partitioning.
  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->blocks_m = blocks; 
  this->pdata_m->partition(UniformGridPartition<Dim>(blocks, igcs, egcs),
		     LocalMapper<Dim>());
}


template <int Dim>
template <class Partitioner>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Partitioner &p,
	   const DistributedTag &)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.

  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->blocks_m = p.blocks(); 
  this->pdata_m->partition(p,DistributedMapper<Dim>(p));
}

template <int Dim>
template <class Partitioner>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Partitioner &p,
	   const ReplicatedTag &)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.

  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->blocks_m = p.blocks(); 
  this->pdata_m->partition(p,LocalMapper<Dim>());
}
template <int Dim>
template <class Partitioner>
inline void
UniformGridLayout<Dim>::
initialize(const Domain_t &gdom,
           const Partitioner &p,
	   const ContextMapper<Dim> &cmap)
{
  PAssert(!this->initialized());

  // Initialize our global domain, and then do the partitioning.

  this->pdata_m->innerdomain_m = gdom;
  this->pdata_m->domain_m = gdom;
  this->pdata_m->blocks_m = p.blocks();
  this->pdata_m->partition(p,cmap); 
}

// This initializer is intented to be used by the I/O system

template <int Dim>
void UniformGridLayout<Dim>::initialize(const Domain_t& idom,
				 const List_t& nodes,
				 const Loc<Dim>& blocks,
				 bool hasIG, bool hasEG,
				 const GuardLayers_t& ig,
				 const GuardLayers_t& eg)
{
  this->pdata_m->initialize(idom,nodes,blocks,hasIG,hasEG,ig,eg);
}

// Here are the implementations for globalID:

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(const Loc<Dim> &loc) const
{
  // Make sure the point is in our domain.
  PAssert(contains(this->domain_m, loc));
  int currloc;
  
  if (!this->hasExternalGuards_m) 
    {
      currloc = (loc[0].first() - this->firsti_m[0]) / blocksizes_m[0];
      for (int d = 1; d < Dim; ++d)
        currloc += blockstride_m[d] * 
          ((loc[d].first() - this->firsti_m[d]) / blocksizes_m[d]);
    }
  else
    {
      currloc = 0;
      for (int d = 0; d < Dim; ++d)
        {
          int l = loc[d].first();
      
          // If l < this->firsti_m[0], currloc is unchanged.
          
          if (l >= this->firsti_m[d]) 
            {
              if (l <= this->innerdomain_m[d].last())
                {
                  // The usual expression in this direction.
                
                  currloc += blockstride_m[d] * 
                    ((l - this->firsti_m[d]) / blocksizes_m[d]);
                }
              else
                {
                  // Must be in the last block in this direction.
                  
                  currloc += blockstride_m[d] * allDomain_m[d].last();
                }
            }
        }
    }
    
  // Return the globalID for the currloc's node

  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0) const
{
  PAssert(Dim == 1);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  
  // Compute fortran-order index from position in block grid
  // See the Loc<Dim> version for comments. 
  
  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
    }
  else
    {
      currloc = 0;
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
          currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
          currloc = allDomain_m[0].last();
      }
    }
    
  // Return the globalID for the currloc's node.
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0, int i1) const
{
  PAssert(Dim == 2);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  PAssert(i1 >= this->domain_m[1].first() && i1 <= this->domain_m[1].last());

  // Compute fortran-order index from position in block grid

  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0]
              + blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
    }
  else
    {
      currloc = 0;
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
          currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
          currloc = allDomain_m[0].last();
      }
      if (i1 >= this->firsti_m[1]) {
        if (i1 <= this->innerdomain_m[1].last())
          currloc += blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
        else
          currloc += blockstride_m[1] * allDomain_m[1].last();
      }
    }
    
  // Return the globalID for the currloc's node
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0, int i1, int i2) const
{
  PAssert(Dim == 3);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  PAssert(i1 >= this->domain_m[1].first() && i1 <= this->domain_m[1].last());
  PAssert(i2 >= this->domain_m[2].first() && i2 <= this->domain_m[2].last());

  // Compute fortran-order index from position in block grid

  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0]
              + blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1])
              + blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2]);
    }
  else
    {
      currloc = 0; 
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
          currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
          currloc = allDomain_m[0].last();
      }
      if (i1 >= this->firsti_m[1]) {
        if (i1 <= this->innerdomain_m[1].last())
          currloc += blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
        else
          currloc += blockstride_m[1] * allDomain_m[1].last();
      }
      if (i2 >= this->firsti_m[2]) {
        if (i2 <= this->innerdomain_m[2].last())
          currloc += blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2]);
        else
          currloc += blockstride_m[2] * allDomain_m[2].last();
      }
    }
    
  // Return the globalID for the currloc's node
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3) const
{
  PAssert(Dim == 4);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  PAssert(i1 >= this->domain_m[1].first() && i1 <= this->domain_m[1].last());
  PAssert(i2 >= this->domain_m[2].first() && i2 <= this->domain_m[2].last());
  PAssert(i3 >= this->domain_m[3].first() && i3 <= this->domain_m[3].last());

  // Compute fortran-order index from position in block grid

  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0]
              + blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1])
              + blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2])
              + blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3]);
    }
  else
    {
      currloc = 0;
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
           currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
           currloc = allDomain_m[0].last();
      }
      if (i1 >= this->firsti_m[1]) {
        if (i1 <= this->innerdomain_m[1].last())
          currloc += blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
        else
          currloc += blockstride_m[1] * allDomain_m[1].last();
      }
      if (i2 >= this->firsti_m[2]) {
        if (i2 <= this->innerdomain_m[2].last())
          currloc += blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2]);
        else
          currloc += blockstride_m[2] * allDomain_m[2].last();
      }
      if (i3 >= this->firsti_m[3]) {
        if (i3 <= this->innerdomain_m[3].last())
          currloc += blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3]);
        else
          currloc += blockstride_m[3] * allDomain_m[3].last();
      }
    }
    
  // Return the globalID for the currloc's node
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
                                     int i4) const
{
  PAssert(Dim == 5);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  PAssert(i1 >= this->domain_m[1].first() && i1 <= this->domain_m[1].last());
  PAssert(i2 >= this->domain_m[2].first() && i2 <= this->domain_m[2].last());
  PAssert(i3 >= this->domain_m[3].first() && i3 <= this->domain_m[3].last());
  PAssert(i4 >= this->domain_m[4].first() && i4 <= this->domain_m[4].last());

  // Compute fortran-order index from position in block grid

  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0]
              + blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1])
              + blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2])
              + blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3])
              + blockstride_m[4] * ((i4 - this->firsti_m[4]) / blocksizes_m[4]);
    }
  else
    {
      currloc = 0;
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
           currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
           currloc = allDomain_m[0].last();
      }
      if (i1 >= this->firsti_m[1]) {
        if (i1 <= this->innerdomain_m[1].last())
          currloc += blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
        else
          currloc += blockstride_m[1] * allDomain_m[1].last();
      }
      if (i2 >= this->firsti_m[2]) {
        if (i2 <= this->innerdomain_m[2].last())
          currloc += blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2]);
        else
          currloc += blockstride_m[2] * allDomain_m[2].last();
      }
      if (i3 >= this->firsti_m[3]) {
        if (i3 <= this->innerdomain_m[3].last())
          currloc += blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3]);
        else
          currloc += blockstride_m[3] * allDomain_m[3].last();
      }
      if (i4 >= this->firsti_m[4]) {
        if (i4 <= this->innerdomain_m[4].last())
          currloc += blockstride_m[4] * ((i4 - this->firsti_m[4]) / blocksizes_m[4]);
        else
          currloc += blockstride_m[4] * allDomain_m[4].last();
      }
    }
    
  // Return the globalID for the currloc's node
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
                                     int i4, int i5) const
{
  PAssert(Dim == 6);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  PAssert(i1 >= this->domain_m[1].first() && i1 <= this->domain_m[1].last());
  PAssert(i2 >= this->domain_m[2].first() && i2 <= this->domain_m[2].last());
  PAssert(i3 >= this->domain_m[3].first() && i3 <= this->domain_m[3].last());
  PAssert(i4 >= this->domain_m[4].first() && i4 <= this->domain_m[4].last());
  PAssert(i5 >= this->domain_m[5].first() && i5 <= this->domain_m[5].last());

  // Compute fortran-order index from position in block grid

  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0]
              + blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1])
              + blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2])
              + blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3])
              + blockstride_m[4] * ((i4 - this->firsti_m[4]) / blocksizes_m[4])
              + blockstride_m[5] * ((i5 - this->firsti_m[5]) / blocksizes_m[5]);
    }
  else
    {
      currloc = 0;
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
           currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
           currloc = allDomain_m[0].last();
      }
      if (i1 >= this->firsti_m[1]) {
        if (i1 <= this->innerdomain_m[1].last())
          currloc += blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
        else
          currloc += blockstride_m[1] * allDomain_m[1].last();
      }
      if (i2 >= this->firsti_m[2]) {
        if (i2 <= this->innerdomain_m[2].last())
          currloc += blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2]);
        else
          currloc += blockstride_m[2] * allDomain_m[2].last();
      }
      if (i3 >= this->firsti_m[3]) {
        if (i3 <= this->innerdomain_m[3].last())
          currloc += blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3]);
        else
          currloc += blockstride_m[3] * allDomain_m[3].last();
      }
      if (i4 >= this->firsti_m[4]) {
        if (i4 <= this->innerdomain_m[4].last())
          currloc += blockstride_m[4] * ((i4 - this->firsti_m[4]) / blocksizes_m[4]);
        else
          currloc += blockstride_m[4] * allDomain_m[4].last();
      }
      if (i5 >= this->firsti_m[5]) {
        if (i5 <= this->innerdomain_m[5].last())
          currloc += blockstride_m[5] * ((i5 - this->firsti_m[5]) / blocksizes_m[5]);
        else
          currloc += blockstride_m[5] * allDomain_m[5].last();
      }
    }
    
  // Return the globalID for the currloc's node
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

template <int Dim>
inline int
UniformGridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
                                     int i4, int i5, int i6) const
{
  PAssert(Dim == 7);
  PAssert(i0 >= this->domain_m[0].first() && i0 <= this->domain_m[0].last());
  PAssert(i1 >= this->domain_m[1].first() && i1 <= this->domain_m[1].last());
  PAssert(i2 >= this->domain_m[2].first() && i2 <= this->domain_m[2].last());
  PAssert(i3 >= this->domain_m[3].first() && i3 <= this->domain_m[3].last());
  PAssert(i4 >= this->domain_m[4].first() && i4 <= this->domain_m[4].last());
  PAssert(i5 >= this->domain_m[5].first() && i5 <= this->domain_m[5].last());
  PAssert(i6 >= this->domain_m[6].first() && i6 <= this->domain_m[6].last());

  // Compute fortran-order index from position in block grid

  int currloc;
  if (!this->hasExternalGuards_m)
    {
      currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0]
              + blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1])
              + blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2])
              + blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3])
              + blockstride_m[4] * ((i4 - this->firsti_m[4]) / blocksizes_m[4])
              + blockstride_m[5] * ((i5 - this->firsti_m[5]) / blocksizes_m[5])
              + blockstride_m[6] * ((i6 - this->firsti_m[6]) / blocksizes_m[6]);
    }
  else
    {
      currloc = 0;
      if (i0 >= this->firsti_m[0]) {
        if (i0 <= this->innerdomain_m[0].last())
           currloc = (i0 - this->firsti_m[0]) / blocksizes_m[0];
        else
           currloc = allDomain_m[0].last();
      }
      if (i1 >= this->firsti_m[1]) {
        if (i1 <= this->innerdomain_m[1].last())
          currloc += blockstride_m[1] * ((i1 - this->firsti_m[1]) / blocksizes_m[1]);
        else
          currloc += blockstride_m[1] * allDomain_m[1].last();
      }
      if (i2 >= this->firsti_m[2]) {
        if (i2 <= this->innerdomain_m[2].last())
          currloc += blockstride_m[2] * ((i2 - this->firsti_m[2]) / blocksizes_m[2]);
        else
          currloc += blockstride_m[2] * allDomain_m[2].last();
      }
      if (i3 >= this->firsti_m[3]) {
        if (i3 <= this->innerdomain_m[3].last())
          currloc += blockstride_m[3] * ((i3 - this->firsti_m[3]) / blocksizes_m[3]);
        else
          currloc += blockstride_m[3] * allDomain_m[3].last();
      }
      if (i4 >= this->firsti_m[4]) {
        if (i4 <= this->innerdomain_m[4].last())
          currloc += blockstride_m[4] * ((i4 - this->firsti_m[4]) / blocksizes_m[4]);
        else
          currloc += blockstride_m[4] * allDomain_m[4].last();
      }
      if (i5 >= this->firsti_m[5]) {
        if (i5 <= this->innerdomain_m[5].last())
          currloc += blockstride_m[5] * ((i5 - this->firsti_m[5]) / blocksizes_m[5]);
        else
          currloc += blockstride_m[5] * allDomain_m[5].last();
      }
      if (i6 >= this->firsti_m[6]) {
        if (i6 <= this->innerdomain_m[6].last())
          currloc += blockstride_m[6] * ((i6 - this->firsti_m[6]) / blocksizes_m[6]);
        else
          currloc += blockstride_m[6] * allDomain_m[6].last();
      }
    }
    
  // Return the globalID for the currloc's node
  
  PAssert(currloc >= 0 && currloc < this->all_m.size());
  return currloc;
}

//============================================================
// NewDomain1 traits classes for UniformGridLayout and
// UniformGridLayoutView
//============================================================

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a UniformGridLayout.
//
//-----------------------------------------------------------------------------

template <int Dim>
struct NewDomain1<UniformGridLayout<Dim> >
{
  typedef UniformGridLayout<Dim> &Type_t;

  inline static Type_t combine(const UniformGridLayout<Dim> &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a UniformGridLayoutView.
//
//-----------------------------------------------------------------------------

template <int Dim, int Dim2>
struct NewDomain1<UniformGridLayoutView<Dim, Dim2> >
{
  typedef UniformGridLayoutView<Dim, Dim2> &Type_t;

  inline static Type_t combine(const UniformGridLayoutView<Dim, Dim2> &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// ostream inserters for UniformGridLayout and UniformGridLayoutView:
//
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &ostr, 
                         const UniformGridLayout<Dim> &layout)
{
  layout.print(ostr);
  return ostr;
}

template <int Dim, int Dim2>
std::ostream &operator<<(std::ostream &ostr, 
                         const UniformGridLayoutView<Dim, Dim2> &layout)
{
  layout.print(ostr);
  return ostr;
}


// } // namespace POOMA

#include "Layout/UniformGridLayout.cpp"

#endif // POOMA_LAYOUT_UNIFORMGRIDLAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformGridLayout.h,v $   $Author: richard $
// $Revision: 1.88 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
