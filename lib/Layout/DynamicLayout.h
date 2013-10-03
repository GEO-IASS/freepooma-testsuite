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

#ifndef POOMA_LAYOUT_DYNAMICLAYOUT_H
#define POOMA_LAYOUT_DYNAMICLAYOUT_H

//-----------------------------------------------------------------------------
// Classes: 
//   DynamicLayout
//   DynamicLayoutView
//   DynamicTag
//   MultiPatchLayoutTraits<DynamicTag,1>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Layout
 * @brief
 *   DynamicLayout
 *     - Layout class that breaks dynamically sized 1-dimensional domain 
 *       into contiguous sub-domains arranged in a 1-dimensional grid.
 *
 *   DynamicLayoutView
 *     - view of a DynamicLayout
 *   DynamicTag
 *     - tag used to specialize MultiPatchLayoutTraits
 *   MultiPatchLayoutTraits<DynamicTag,1>
 *     - traits class used by MultiPatch-engine to determine layout type.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Layout/MultiPatchLayoutTraits.h"
#include "Layout/TouchesConstruct.h"
#include "Layout/DynamicEvents.h"
#include "Partition/ContextMapper.h"
#include "Partition/UniformMapper.h"
#include "Domain/Interval.h"
#include "Domain/Grid.h"
#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/NewDomain.h"
#include "Domain/DomainMap.h"
#include "Utilities/DerefIterator.h"
#include "Utilities/Observable.h"
#include "Utilities/Observer.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/RefCounted.h"
#include "Utilities/Unique.h"
#include "Utilities/PAssert.h"
#include "Threads/PoomaSmarts.h" // Smarts::concurrency()
#include "Pooma/Pooma.h" // Pooma::context() (ugh!!! brings in too much stuff)

#include <vector>
#include <iosfwd>

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//============================================================
// Forward declarations
//============================================================

class DynamicLayoutData;
class DynamicLayout;

class DynamicLayoutViewData;
class DynamicLayoutView;

// template <int Dim> class ContextMapper;

//-----------------------------------------------------------------------------
// Full Description:
// DynamicTag
//
// DynamicTag class.
//-----------------------------------------------------------------------------

struct DynamicTag { };

/**
 * Specialization of MultiPatchLayoutTraits for DynamicLayout.
 */

template <>
struct MultiPatchLayoutTraits<DynamicTag,1>
{
  typedef DynamicLayout Layout_t;
  
  template <int ViewDim>
  struct View
  {
    typedef DynamicLayoutView Layout_t;
  };
};

/**
 * Holds the data for a DynamicLayout. That class has a ref counted
 * instance of this class
 */

class DynamicLayoutData 
  : public RefCounted,
    public Observable<DynamicLayoutData>
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef DynamicLayoutData                    This_t;
  typedef Observable<This_t>                   Observable_t;
  typedef Interval<1>                          Domain_t;
  typedef Interval<1>                          BaseDomain_t;
  typedef int                                  Context_t;
  typedef Unique::Value_t                      ID_t;
  typedef Node<Domain_t>                       Value_t;
  typedef std::vector<Value_t *>               List_t;
  typedef int                                  AxisIndex_t;
  typedef DynamicEvents::PatchID_t             PatchID_t;
  typedef DynamicEvents::CreateSize_t          CreateSize_t;

  // Enumerations.

  enum { dimensions = 1 };
  enum { repartitionEvent = 1 };
  enum { dynamic = true };


  //============================================================
  // Constructors
  //============================================================

  // Default constructor: initially no blocks, etc.

  DynamicLayoutData();

  // The partitioner provides the information for actually constructing
  // the layout's data.
  
  template <class Partitioner>
  DynamicLayoutData(const Domain_t &gdom, const Partitioner &gpar,
		    const ContextMapper<1> &cmap)
    : Observable_t(*this),
    ID_m(Unique::get()), 
    dirtyLayout_m(true)
  {
    initialize(gdom, gpar, cmap);
  }

  
  //============================================================
  // Destructor
  //============================================================

  // This is only called when all references to this data go away, in
  // which case we need to delete our nodes. The Observable destructor
  // will broadcast messages up to all observers of the Layout.

  ~DynamicLayoutData();

  //============================================================
  // Mutators
  //============================================================

  // Initialize this object by invoking the partitioner and setting
  // up the domains.  
  // These can be called after using the default constructor.

  template <class Partitioner>
  void initialize(const Domain_t &gdom, const Partitioner &gpar,
                  const ContextMapper<1> &cmap)
  {
    // This will work with grid (and simpler) partitioners.

    CTAssert(Partitioner::gridded);

    // Delete existing nodes and clear all the lists.

    if (all_m.size() > 0)
      {
        for (int i = 0; i < all_m.size(); ++i) delete all_m[i];
        all_m.clear();
        local_m.clear();
        remote_m.clear();
      }

    domain_m = gdom;

    PAssert(!gpar.hasInternalGuards());
    PAssert(!gpar.hasExternalGuards());
  
    // Invoke the partitioner, which adds the subdomains directly to 
    // the all_m list.

    gpar.partition(domain_m, all_m, cmap);
    PAssert(gpar.blocks().first() == all_m.size());
  
    // Create the local and remote lists. 

    Context_t thisContext = Pooma::context();

    for (int i = 0; i < all_m.size(); ++i)
      {
        Context_t ci = all_m[i]->context(); 
        if (ci == thisContext || ci == -1)
          local_m.push_back(all_m[i]);
        else
          remote_m.push_back(all_m[i]);
      }
        
    calcMaps();

    dirtyLayout_m = false;

  }

  // Used by the I/O or data management system to initialize the layout based
  // on detailed state information previously stored.

void initialize(const Domain_t& gdom,
		const List_t& nodes)
{
  int i;

  // delete existing nodes and clear all the lists

  if (all_m.size() > 0)
    {
      for (i=0; i < all_m.size(); ++i)
	delete all_m[i];
      all_m.clear();
      local_m.clear();
      remote_m.clear();
    }

  domain_m = gdom;

  // Assign the given list of nodes to the total list.
  all_m= nodes;

  // Iterate through the complete list of nodes provided and assign to the
  // appropriate subcategories.

  List_t::iterator start = all_m.begin();
  List_t::iterator end   = all_m.end();
  
  for ( ; start!=end ;++start )
    {
      if( (*start)->context() == Pooma::context() ||
	  (*start)->context() == -1 )
	local_m.push_back(*start);
      else
	remote_m.push_back(*start);
    }

  // Calculate the domain maps.
  calcMaps();

  // Set the dirty layout flag.
  dirtyLayout_m = false;

}

  template <class Partitioner>
  inline void initialize(const Grid<1> &gdom, const Partitioner &gpar)
  {
    Domain_t idom(gdom.first(), gdom.last() - 1);
    initialize(idom, gpar);
  }

  //============================================================
  // Accessors
  //============================================================

  inline ID_t  ID() const                    { return ID_m; }
  inline bool  initialized() const           { return (all_m.size() > 0); }
  inline int   blocks() const                { return all_m.size(); }
  inline bool  dirty() const                 { return dirtyLayout_m; }
  inline const Domain_t &domain() const      { return domain_m; }
  inline const Domain_t &ownedDomain() const { return domain_m; }

  // Accessors to get node domains: (do we need all of these???)
  
  inline const Domain_t &domain(int i) const 
  {
    PAssert(i >= 0 && i < local_m.size());
    return local_m[i]->domain();
  }

  inline const Domain_t &ownedDomain(int i) const      { return domain(i); }
  inline const Domain_t &allocatedDomain(int i) const  { return domain(i); }
  inline const Domain_t &patchDomain(int i) const      { return domain(i); }
  inline const Domain_t &patchDomainOwned(int i) const { return domain(i); }

  // Accessors to get node lists:
  
  inline List_t &nodeListGlobal() { return all_m; }
  inline List_t &nodeListLocal()  { return local_m; }
  inline List_t &nodeListRemote() { return remote_m; }

  // Accessors for getting the global ID of the patch containing
  // a particular element.

  int globalID(const Loc<1> &) const;
  int globalID(int) const;

  //============================================================
  // touches operations
  //============================================================

  // Find all subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.
  // Return the number of touching elements. This version of touches
  // can build either pointers or objects.

  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &fulld, OutIter o, 
              const ConstructTag &ctag) const
  {
    // Make sure we have a valid layout

    PAssert(initialized());

    // Figure the type of the domain resulting from the intersection

    typedef typename 
      IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;

    // We only need to do touches for the overlapping domain.  If there is
    // nothing left, we can just return.

    OutDomain_t d = intersect(domain_m, fulld);
    if (d.empty())
      return 0;

    // Create an object of this output domain type for use later.

    OutDomain_t outDomain = Pooma::NoInit();

    // Find the begin/end iterator pair of the touching domains.

    typedef DomainMap<Interval<1>,AxisIndex_t>::Touch_t Touch_t;
    typedef DomainMap<Interval<1>,AxisIndex_t>::touch_iterator MapIterator_t;
    // HACK ALERT!!! DomainMap<Interval,int>::touches only takes an Interval
    // as an argument. This seems like a DomainMap deficiency.
    Touch_t touchIter = map_m.touch(Interval<1>(d));
    MapIterator_t begin = touchIter.first;
    MapIterator_t end = touchIter.second;

    // Go through all the blocks and output the values.

    int count = 0;
  
    for (MapIterator_t iter = begin; iter != end; ++iter)
      { 
        AxisIndex_t i = *iter;

	// Make sure that block is OK ... this is a sanity check.

	outDomain = intersect(fulld,all_m[i]->domain());
	PAssert(!outDomain.empty());

	// Output this block.

	*o++ = touchesConstruct(outDomain,
		  		all_m[i]->allocated(),
				all_m[i]->affinity(),
				all_m[i]->context(),
				all_m[i]->globalID(),
				all_m[i]->localID(),
				ctag);
	++count;
      }
    return count;
  }

  //============================================================
  // Dynamic engine methods
  //============================================================

  // Create new elements by extending the current domain
  // of the specified local patch by the requested number of elements.
  // 'local' means on this same context.  The patch is referred to
  // by local index, from 0 ... # local patches - 1.
  // The default is to create elements in the last local patch.
  // All observers are notified of the change, then we change our
  // domain value.

  void create(CreateSize_t num, PatchID_t patch = (-1));

  // Destroy the elements in given patch using the provided
  // domain as offsets into that patch, and using the specified
  // delete method.
  // The domain values in this case should be zero-based.

  template <class Dom, class DestroyMethod>
  void destroy(const Dom &dom, PatchID_t fromPatch,
	       const DestroyMethod &)
  {
    PAssert(fromPatch < static_cast<int>(local_m.size()));
    PAssert(fromPatch >= 0);

    typedef typename DynamicEventType<Dom>::Domain_t DynDomain_t;

    // Perform destroy operation for the specified patch and the
    // "local" domain (all domain values are in the range 0 ...
    // patchsize - 1)

    Value_t *lp = local_m[fromPatch];
  
    PAssert(contains(lp->domain() - lp->domain().first(),
		     Interval<1>(dom.first(), dom.last())));
    
    // Let all registered engine's know that they must destroy these

    notify(DestroyEvent<DynDomain_t>(dom, fromPatch, DestroyMethod::code));

    // Modify the domain for this local patch.  When sync is called,
    // everything else will get updated.

    deleteElements(lp->domain(), dom.size());

    // Note that we will need to rebuild things, and return.

    dirtyLayout_m = true;
  }

  // Destroy the elements specified by the global domain gdom.

  template <class Dom, class DestroyMethod>
  void destroy(const Dom &gdom, const DestroyMethod &)
  {
    PAssert(contains(domain_m, Interval<1>(gdom.first(),gdom.last())));

    typedef DynamicEventType<IndirectionList<int> >::Domain_t ILDynDomain_t;

    // Find pieces of this total destroy domain in each subdomain,
    // and destroy them.

    int is = 0, ie = 0;
    PatchID_t ip = 0;

    // Skip to the first non-empty local patch

    bool found = false;
    while (!found && ip < local_m.size()) {
      found = !local_m[ip]->domain().empty();
      if (!found) ip++;
    }
    if (!found) return;

    // Some portion of the destroy domain may precede all of the 
    // domain controlled by this context, so skip that part.

    while (gdom(is) < local_m[ip]->domain().first() && is < gdom.size())
      is++;
    ie = is;

    while (ip < local_m.size() && ie < gdom.size())
      {
	for ( ; ip < local_m.size() && ie < gdom.size(); ++ie)
	  if (gdom(ie) > local_m[ip]->domain().last()) break;

	if (ie == is)
	  {
	    ++ip;
	    continue;
	  }

	Value_t *lp = local_m[ip];
          
	IndirectionList<int> iltemp(ie-is);

	// Patch specific dynamic ops are zero based within the patch...

	for (int j = is; j < ie; ++j)
	  iltemp(j-is) = gdom(j) - lp->domain().first();

	notify(DestroyEvent<ILDynDomain_t>(iltemp, ip, DestroyMethod::code));

	// Modify the domain for this local patch.  When sync is called,
	// everything else will get updated.

	deleteElements(lp->domain(), iltemp.size());
            
	// Move on to next local patch

	++ip;
	is = ie;
      }

    // Note that we will need to rebuild things, and return.

    dirtyLayout_m = true;
  }


  // Copy all elements of domain n to the end of patch p.  If p < 0, 
  // copy to the end of the last local patch.
  // This is for a domain in the global domain space.

  template <class Dom>
  void copy(const Dom &dom, PatchID_t toPatch = (-1))
  {
    PAssert(contains(domain(), Interval<1>(dom.first(),dom.last()) ));

    typedef DynamicEventType<IndirectionList<int> >::Domain_t ILDynDomain_t;

    int is = 0, ie = 0;
    PatchID_t ip = 0;

    // Adjust the toPatch, if necessary

    if (toPatch < 0) toPatch = local_m.size() - 1;
    PAssert(toPatch < local_m.size());

    // Go through the patches, and copy the interesecting domains
    // to the specified patch.

    while (ip < local_m.size() && ie < dom.size())
      {
        Value_t *lp = local_m[ip];
        for ( ; ip < local_m.size() && ie < dom.size(); ++ie)
	  if (dom(ie) > lp->domain().last()) break;

        if (ie == is)
	  {
	    ++ip;
	    continue;
	  }
      
        IndirectionList<int> iltemp(ie-is);
      
        // Patch specific dynamic ops are zero based within the patch...

        for (int j = is; j < ie; ++j)
	  iltemp(j-is) = dom(j) - lp->domain().first();

        // Let all engine's know they must copy data between these two patches.

        notify(CopyEvent<ILDynDomain_t>(iltemp, ip, toPatch));

        // Modify the toPatch domain

        addElements(local_m[toPatch]->domain(), iltemp.size());
        
        // Move on to next local patch
      
        ++ip;
        is = ie;
      }
  
    // Note that we will need to rebuild things, and return.

    dirtyLayout_m = true;
  }


  // Copy all elements of domain dom to the end of patch p.  This version
  // also specifies the patch to copy values from.
  // In this case, the domain values should be zero-based.

  template <class Dom>
  void copy(const Dom &dom, PatchID_t fromPatch, PatchID_t toPatch)
  {
    PAssert(fromPatch < local_m.size());

    // If the toPatch number is < 0, change it to the last local patch.

    if (toPatch < 0) toPatch = local_m.size() - 1;
    PAssert(toPatch < local_m.size());

    // This is a local copy, so check the domain in the "local"
    // domain space (dom contains zero-based domain values).

    PAssert(dom.max() < local_m[fromPatch]->domain().size());

    // Let all users know of the copy request.

    typedef typename DynamicEventType<Dom>::Domain_t  DynDomain_t;
    notify(CopyEvent<DynDomain_t>(dom, fromPatch, toPatch));

    // Modify the domain for this local patch.  When sync is called,
    // everything else will get updated.

    addElements(local_m[toPatch]->domain(), dom.size());
    
    // Note that we will need to rebuild things.

    dirtyLayout_m = true;
  }


  // Perform a "multiple patch" copy, using a list of IndirectionList's
  // for a set of source patches, and an IndirectionList giving the
  // patch ID for the source patches.  Copy data into the destination
  // patch.  The source and desination patches must be specified, this
  // is only for "zero-based" index lists.  If the last argument is
  // true, storage is created at the end, otherwise elements are
  // just copied to the end of the existing storage.

  void copy(const IndirectionList< IndirectionList<int> > &domlists,
	    const IndirectionList< int > &fromlist,
	    PatchID_t toPatch,
	    bool docreate);

  // Sync up the layout with any other contexts, taking into account
  // that other contexts may have performed create/destroy operations.
  // This will reset all the local domains to be properly contiguous,
  // and let all engine's using this layout reset their domains.

  void sync();

  template <class Out>
  void print(Out &ostr)
  {
    ostr << " dirtyLayout_m = " << dirtyLayout_m << std::endl;
  }

private:

  //============================================================
  // Private Methods
  //============================================================

  // This function recalculates what the total domain of each patch
  // should be, since this can change due to dynamic operations.

  void calcDomains();

  // This function recalculates the domain maps, since this can
  // change due to dynamic operations.

  void calcMaps() const;

  // This function does a cross-context recalculation of the
  // patch domains. 

  void syncGlobalDomains();

  // Add or subtract elements from the given domain
  
  inline void addElements(Domain_t &domain, int num)
  {
    PAssert(num > 0);
    if (domain.size() > 0)
       domain = Domain_t(domain.first(), domain.last() + num);
    else
       domain = Domain_t(num);
  }
  
  inline void deleteElements(Domain_t &domain, int num)
  {
    PAssert(num <= domain.size());
    if (num < domain.size())
      domain = Domain_t(domain.first(), domain.last() - num);
    else
      domain = Domain_t();
  }
  
  //============================================================
  // Private Data
  //============================================================

  // Our ID value, which is simply a unique value.

  ID_t ID_m;

  // The global domain of this DynamicLayout.

  Domain_t domain_m;
  
  // The list of all, local, and remote subdomains.

  List_t all_m;
  List_t local_m;
  List_t remote_m;

  // Is the current layout information out-of-date (due to dynamic ops
  // that have yet to have sync called?)

  bool dirtyLayout_m;

  // Domain map storing the subdomains in a tree for fast lookup.

  mutable DomainMap<Interval<1>,AxisIndex_t> map_m;

};


/**
 * DynamicLayout is a Layout class that breaks an N-dimensional
 * Interval into sub-domains arranged in an N-dimensional
 * grid, where the sub-domain sizes are specified by a Grid domain object
 *
 * This is an alternative to the more general tile Layout class that should
 * perform faster since subdomains can be found using a set of 1-dimensional
 * domainMap's, rather than by a more general search.
 *
 * To construct a DynamicLayout, you can do any of the following:
 *   -# provide a global domain, and let the DynamicLayout perform
 *      its default partitioning by just using one single block;
 *   -# provide a global domain, a Loc with the number of blocks to
 *      use along each dimension, and an optional context number;
 *   -# provide a global domain and a GridPartition or UniformGridPartition 
 *      object.
 *   -# provide a Grid domain object
 */      

class DynamicLayout : public Observable<DynamicLayout>,
                      public Observer<DynamicLayoutData>
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.

  enum { dimensions       = 1     };
  enum { repartitionEvent = 1     };
  enum { dynamic          = true  };
  enum { supportsGuards   = false };

  // General public typedefs.

  typedef DynamicLayout                        This_t;
  typedef Observable<This_t>                   Observable_t;
  typedef DynamicLayoutData                    LayoutData_t;
  typedef LayoutData_t::Domain_t               Domain_t;
  typedef LayoutData_t::BaseDomain_t           BaseDomain_t;
  typedef LayoutData_t::Context_t              Context_t;
  typedef LayoutData_t::ID_t                   ID_t;
  typedef LayoutData_t::Value_t                Value_t;
  typedef LayoutData_t::List_t                 List_t;
  typedef DynamicEvents::PatchID_t             PatchID_t;
  typedef DynamicEvents::CreateSize_t          CreateSize_t;

  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.  
  
  typedef DerefIterator<Value_t>               iterator;
  typedef ConstDerefIterator<Value_t>          const_iterator;
  
  //============================================================
  // Constructors
  //============================================================

  // The default constructor does not initialize the layout.  In this
  // case, layout initialization must be completed with the
  // 'initialize' method before the layout can be used.  A default
  // layout has an empty global domain, and empty subdomain lists.
  
  DynamicLayout();

  // Construct a layout with nothing else but a global domain.  In
  // this case, a default partitioner will be used, the GridPartition
  // object, which will just make a grid with one block.
  
  DynamicLayout(const Domain_t &);

  // Domain + block count constructors.
  // In this case, the user specifies the global domain and the number
  // of blocks in each dimension, which will cause the domain to be
  // partitioned in a near-uniform manner. 
  
  DynamicLayout(const Domain_t &,
	        int);

  // Grid Domain constructors. 
  // Like the block count constructors, but a Grid object is specified,
  // allowing a general non-uniform grid to be created.  
  
  DynamicLayout(const Grid<1> &);

  // Domain + partition constructor. 

  template <class Partitioner>
  DynamicLayout(const Domain_t &gdom, 
                const Partitioner &gpar)
    : Observable<This_t>(*this), 
      pdata_m( new LayoutData_t( gdom, gpar, UniformMapper(gpar) ) )
  {
    pdata_m->attach(*this);
  }

  // Domain + partition + mapper constructor. 

  template <class Partitioner>
  DynamicLayout(const Domain_t &gdom, 
                const Partitioner &gpar,
                const ContextMapper<1> &cmap)
    : Observable<This_t>(*this), 
      pdata_m( new LayoutData_t(gdom, gpar, cmap) )
  {
    pdata_m->attach(*this);
  }

  // Copy constructor & assignment operator
  // Shallow copies with reference counting.

  DynamicLayout(const This_t &);

  This_t &operator=(const This_t &);
  

  //============================================================
  // Destructor
  //============================================================

  // The actual data will be cleaned up by the LayoutData_t destructor
  // if all references to the data go away.
  // If any Observers remain, they will be notified by the Observable 
  // destructor. 
  
  inline ~DynamicLayout() { pdata_m->detach(*this); }

  //============================================================
  // Initialize methods
  //============================================================

  // Initialize a layout with nothing else but a global domain.  In
  // this case, a default partitioner will be used, the GridPartition
  // object, which will try to make a grid with one block.
  
  void initialize(const Domain_t &);

  // Domain + block count initializer.

  void initialize(const Domain_t &, 
		  int);

  // Domain + Grid object based initializers.

  void initialize(const Domain_t &,
		  const Grid<1> &);

  // Grid object based initializers.

  void initialize(const Grid<1> &);

  // Domain + partition initializer.

  template <class Partitioner>
  void initialize(const Domain_t &gdom, 
                  const Partitioner &gpar)
  {
    UniformMapper cmap(gpar);
    pdata_m->initialize(gdom, gpar, cmap);
  }

  // Domain + partition + mapper initializer.

  template <class Partitioner>
  void initialize(const Domain_t &gdom, 
                  const Partitioner &gpar,
                  const ContextMapper<1> &cmap)
  {
    pdata_m->initialize(gdom, gpar, cmap);
  }

  // Used by the I/O or data management system to initialize the layout based
  // on detailed state information previously stored.

  void initialize(const Domain_t& gdom,
	     const List_t& nodes){
    pdata_m->initialize(gdom, nodes);
  }

  //============================================================
  // Accessors
  //============================================================

  // Return ID value and base ID (which is same as ID).
  
  inline ID_t ID() const
  {
    return pdata_m->ID();
  }

  inline ID_t baseID() const
  {
    return pdata_m->ID();
  }

  // Return whether or not this layout has been initialized.
  
  inline bool initialized() const
  {
    return (sizeGlobal() > 0);
  }

  // Return starting point of domain along given dimension.
  // The only valid dimension is 1 in this case.

  inline int first(int) const
  {
    return pdata_m->domain().first();
  }

  // Return the global domain and the base domain.
  
  inline const Domain_t &domain() const
  {
    return pdata_m->domain();
  }

  inline const Domain_t &ownedDomain() const
  {
    return pdata_m->ownedDomain();
  }

  inline const Domain_t &domain(int i) const
  {
    return pdata_m->domain(i);
  }

  inline const Domain_t &ownedDomain(int i) const
  {
    return pdata_m->ownedDomain(i);
  }

  inline const Domain_t &allocatedDomain(int i) const
  {
    return pdata_m->allocatedDomain(i);
  }

  inline const Domain_t &baseDomain() const
  {
    return pdata_m->domain();
  }

  inline const Domain_t &patchDomain(int i) const
  {
    return pdata_m->patchDomain(i);
  }

  inline const List_t & nodeListGlobal() const
  {
    return pdata_m->nodeListGlobal();
  }

  inline const List_t & nodeListLocal() const
  {
    return pdata_m->nodeListLocal();
  }

  inline const List_t & nodeListRemote() const
  {
    return pdata_m->nodeListRemote();
  }

  // Compare to another Layout.  The layouts are the same if:
  //   1. They have the same base ID value.
  //   2. They have the same base domain.
  
  template <class Layout>
  inline bool operator==(const Layout &layout) const
  {
    return (baseID()==layout.baseID() && baseDomain()==layout.baseDomain());
  }

  // Compare for inequality.
  
  template <class Layout>
  inline bool operator!=(const Layout &layout) const 
  {
    return !(*this == layout);
  }


  //============================================================
  // DynamicLayout-specific accessors
  //============================================================

  // Return the number of blocks.
  
  inline int blocks() const 
  {
    return pdata_m->blocks();
  }


  //============================================================
  // Data lookup
  //============================================================

  // Return the globalID of the node containing the given element.
  
  inline int globalID(const Loc<1> &loc) const
  {
    return pdata_m->globalID(loc);
  }

  inline int globalID(int a1) const
  {
    return pdata_m->globalID(a1);
  }

  //============================================================
  // Iterators
  //============================================================

  // Return begin and end iterators for the list of all subdomains
  
  inline iterator beginGlobal()
  { 
    return iterator(pdata_m->nodeListGlobal().begin()); 
  }

  inline iterator endGlobal()
  { 
    return iterator(pdata_m->nodeListGlobal().end());
  }

  inline const_iterator beginGlobal() const
  {
    return const_iterator(pdata_m->nodeListGlobal().begin());
  }

  inline const_iterator endGlobal() const
  { 
    return const_iterator(pdata_m->nodeListGlobal().end());
  }

  inline int sizeGlobal() const
  { 
    return pdata_m->nodeListGlobal().size(); 
  }

  // Return begin and end iterators for the list of all subdomains
  
  inline iterator beginLocal()
  {
    return iterator(pdata_m->nodeListLocal().begin()); 
  }

  inline iterator endLocal()
  {
    return iterator(pdata_m->nodeListLocal().end());
  }

  inline const_iterator beginLocal() const
  { 
    return const_iterator(pdata_m->nodeListLocal().begin()); 
  }

  inline const_iterator endLocal() const
  {
    return const_iterator(pdata_m->nodeListLocal().end());
  }

  inline int sizeLocal() const
  {
    return pdata_m->nodeListLocal().size();
  }

  // Return begin and end iterators for the list of all subdomains
  
  inline iterator beginRemote()
  {
    return iterator(pdata_m->nodeListRemote().begin()); 
  }

  inline iterator endRemote()
  {
    return iterator(pdata_m->nodeListRemote().end()); 
  }

  inline const_iterator beginRemote() const
  {
    return const_iterator(pdata_m->nodeListRemote().begin()); 
  }

  inline const_iterator endRemote() const
  {
    return const_iterator(pdata_m->nodeListRemote().end()); 
  }

  inline int sizeRemote() const
  {
    return pdata_m->nodeListRemote().size(); 
  }

  //============================================================
  // Repartition the layout using a new Partitioner scheme.  The
  // initial domain lists are cleared out, the partitioner is invoked,
  // and then all the observers are notified.  This can only be done
  // with a partitioner that generates grid-like blocks.

  template <class Partitioner>
  inline bool repartition(const Partitioner &gp)
  {
    pdata_m->initialize(domain(),gp);
    // Can't this be Observable_t::notify???
    pdata_m->notify(repartitionEvent);
    return true;
  }

  //============================================================
  // Dynamic operations
  //============================================================

  void create(CreateSize_t num, PatchID_t patch=(-1) )
  {
    pdata_m->create(num, patch);
  }

  // This is for the global domain.

  template <class Dom, class DeleteMethod>
  void destroy(const Dom &killlist, const DeleteMethod &method)
  {
    pdata_m->destroy(killlist, method);
  }

  template <class Dom, class DeleteMethod>
  void destroy(const Dom &killlist, PatchID_t patch, 
               const DeleteMethod &method)
  {
    pdata_m->destroy(killlist, patch, method);
  }

  template <class Dom>
  void copy(const Dom &copylist, PatchID_t toPatch=(-1))
  {
    pdata_m->copy(copylist, toPatch);
  }

  template <class Dom>
  void copy(const Dom &copylist, PatchID_t fromPatch, PatchID_t toPatch)
  {
    // If this is a multi-patch copy, and fromPatch < 0, we're copying
    // with global domain values.  The other version of copy will
    // break the total copy domain up into pieces for each patch, with
    // relative domain values and patch indices >= 0.

    if (fromPatch < 0)
      pdata_m->copy(copylist, toPatch);
    else
      pdata_m->copy(copylist, fromPatch, toPatch);
  }

  void copy(const IndirectionList<IndirectionList<int> > &domlists,
	    const IndirectionList<int> &fromlist,
	    PatchID_t toPatch,
	    bool docreate)
  {
    pdata_m->copy(domlists, fromlist, toPatch, docreate);
  }

  void sync()
  { 
    pdata_m->sync();
  }


  //============================================================
  // Touch methods
  //============================================================

  // Find all subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.
  // Return the number of touching elements. This version of touches
  // can build either pointers or objects.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, 
              const ConstructTag &ctag) const 
  {
   return pdata_m->touches(d, o, ctag);
  }

  // Dynamic layouts have no guards, so touchesAlloc just calls 
  // the underlying touches.
  
  template <class OtherDomain, class  OutIter, class ConstructTag>
  int touchesAlloc(const OtherDomain &d, OutIter o, 
                   const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d, o, ctag);
  }

  // Find local subdomains that touch on a given domain, and insert
  // the intersection of these subdomains into the given output
  // iterator.  Return the number of touching elements. This version
  // of touches can build either pointers or objects. 
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesLocal(const OtherDomain &d, OutIter o, 
		          const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d, o, ctag);
  }

  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesAllocLocal(const OtherDomain &d, OutIter o, 
			       const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d, o, ctag);
  }


  // Find remote subdomains that touch on a given domain, and insert
  // the intersection of these subdomains into the given output
  // iterator.  Return the number of touching elements. This version
  // of touches can build either pointers or objects. 
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesRemote(const OtherDomain &, OutIter, 
		           const ConstructTag &) const 
  {
    return 0;
  }
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesAllocRemote(const OtherDomain &, OutIter, 
			        const ConstructTag &) const 
  {
    return 0;
  }

  // Find all/local/remote subdomains that touch on a given domain,
  // and insert the intersection of these subdomains into the given
  // output iterator.  Return the number of touching elements.  These
  // versions of touches can build only objects objects.
  
  template <class OtherDomain, class OutIter>
  inline int touches(const OtherDomain &d, OutIter o) const 
  {
    return touches(d, o, TouchesConstructNodeObj());
  }

  // same as above, but returns the touches on the allocated domains.

  template <class OtherDomain, class OutIter>
  inline int touchesAlloc(const OtherDomain &d, OutIter o) const 
  {
    return touchesAlloc(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesLocal(const OtherDomain &d, OutIter o) const 
  {
    return touchesLocal(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesAllocLocal(const OtherDomain &d, OutIter o) const 
  {
    return touchesAllocLocal(d, o, TouchesConstructNodeObj());
  }

  template <class OtherDomain, class OutIter>
  inline int touchesRemote(const OtherDomain &d, OutIter o) const 
  {
    return touchesRemote(d, o, TouchesConstructNodeObj());
  }
  
  template <class OtherDomain, class OutIter>
  inline int touchesAllocRemote(const OtherDomain &d, OutIter o) const 
  {
    return touchesRemote(d, o, TouchesConstructNodeObj());
  }

  //============================================================
  // Observer methods
  //============================================================

  // Respond to events generated by the LayoutData_t.
  // These are just passed on to our observers.

  virtual void notify(LayoutData_t &d, const ObserverEvent &event)
  {
    // We should only get this message from our LayoutData_t object

    PAssert(&d == pdata_m.rawPointer());
    Observable_t::notify(event);
  }


  //============================================================
  // Output
  //============================================================

  // Print a DynamicLayout on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const 
  {
    ostr << "DynamicLayout " << ID() << " on global domain " 
      << domain() << ":" << '\n';
    ostr << "   Total subdomains: " << sizeGlobal() << '\n';
    ostr << "   Local subdomains: " << sizeLocal() << '\n';
    ostr << "  Remote subdomains: " << sizeRemote() << '\n';
    ostr << "        Grid blocks: " << blocks() << '\n';
    const_iterator a;
    for (a = beginGlobal(); a != endGlobal(); ++a)
      ostr << "  Global subdomain = " << *a << '\n';
    for (a = beginLocal(); a != endLocal(); ++a)
      ostr << "   Local subdomain = " << *a << '\n';
    for (a = beginRemote(); a != endRemote(); ++a)
      ostr << "  Remote subdomain = " << *a << '\n';
    pdata_m->print(ostr);
  }

private:

  //============================================================
  // Data
  //============================================================

  // DynamicLayout stores its data in a RefCounted class to
  // simplify memory management.
  
  RefCountedPtr<LayoutData_t> pdata_m;
};


/**
 * The data object held by a DynamicLayoutView object.
 */

class DynamicLayoutViewData : public RefCounted
{
public:

  typedef DynamicLayout                     Layout_t;
  typedef DynamicLayoutView                 ViewLayout_t;

  typedef Interval<1>                       Domain_t;
  typedef Range<1>                          BaseDomain_t;
  typedef int                               Context_t;
  typedef Unique::Value_t                   ID_t;

  typedef Layout_t::Domain_t                AllocatedDomain_t;

  typedef Node<Domain_t,AllocatedDomain_t>  Value_t;
  typedef std::vector<Value_t *>            List_t;        // for convenience

  typedef DynamicLayoutViewData             LayoutData_t;

  //============================================================
  // Constructors
  //============================================================

  // Default constructor. Creates an empty view of an empty layout.
  
  DynamicLayoutViewData() 
  : id_m(Unique::get()), subdomainsComputed_m(false)
  { }
  
  // Other constructors.
  
  template <class DT>
  inline DynamicLayoutViewData(const Layout_t &layout, 
                               const Domain<1, DT> &dom)
  : id_m(Unique::get()),
    layout_m(layout), 
    domain_m(dom.unwrap().length()),
    baseDomain_m(dom.unwrap()),
    subdomainsComputed_m(false)
  { 
    // The layout passed in must be initialized.  

    PAssert(layout_m.initialized());

    // The domain we're passing in must be contained in the base
    // layout.

    PAssert(contains(layout_m.domain(), dom.unwrap()));
  }

  template <class DT>
  inline
  DynamicLayoutViewData(const ViewLayout_t &layout, const Domain<1, DT> &dom);

  // Destructor

  ~DynamicLayoutViewData() 
  {
    List_t::iterator a;
    for (a = all_m.begin(); a != all_m.end(); ++a)
      delete (*a);
  }

  // Return the globalID of the node containing the given element.

  inline int globalID(const Loc<1> &loc, Loc<1> &oloc) const
  {
    oloc = baseDomain_m.first() + baseDomain_m.stride() * loc.first();
    return layout_m.globalID(oloc);
  }

  inline int globalID(int i0, Loc<1> &oloc) const
  {
    oloc = baseDomain_m.first() + baseDomain_m.stride() * i0;
    return layout_m.globalID(oloc);
  }

  // Touches calculation.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, 
	      const ConstructTag &ctag) const 
  {
    // Transform the local domain to base coordinates.

    BaseDomain_t bd = Pooma::NoInit();
    localToBase(d, bd);

    // Run the touches function for our underlying layout

    std::vector<Node<BaseDomain_t,AllocatedDomain_t> > tnodes;
    int count = layout_m.touches(bd, std::back_inserter(tnodes));

    // Convert the domains back to the local coordinates, construct
    // appropriate return values, and push them onto the output list.

    Interval<1> ld = Pooma::NoInit();

    for (int i = 0; i < count; i++)
      {
        baseToLocal(tnodes[i].domain(), ld);
	*o++ = touchesConstruct(ld,
			        tnodes[i].allocated(), // Don't convert this!
			        tnodes[i].affinity(),
			        tnodes[i].context(),
			        tnodes[i].globalID(), 
			        tnodes[i].localID(), 
			        ctag);
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

    // We need to find the nodes that intersect with our base domain.
    // To do this, run the touches function for our underlying layout.

    std::vector<Node<BaseDomain_t,AllocatedDomain_t> > tnodes;
    int count = layout_m.touches(baseDomain_m, std::back_inserter(tnodes));

    // Now, run through the nodes we've got and patch the domains up.

    Domain_t ld = Pooma::NoInit();

    const Context_t myContext = Pooma::context();
    
    for (int i = 0; i < count; i++)
      {
        baseToLocal(tnodes[i].domain(), ld);
	Value_t *pt = touchesConstruct(ld,
                                       tnodes[i].allocated(), // Don't convert
                                       tnodes[i].affinity(),
			               tnodes[i].context(),
			               tnodes[i].globalID(),
			               tnodes[i].localID(),
			               TouchesConstructNodePtr());
	all_m.push_back(pt);
	
	// Sort these into local and remote.
	
        if (pt->context() == myContext
	    || pt->context() == -1 )
          local_m.push_back(pt);
        else
          remote_m.push_back(pt);
      }

    // Set flag indicating we've computed these subdomains.

    subdomainsComputed_m = true;
  }

  template <class Domain>
  void localToBase(const Domain &d, BaseDomain_t &bd) const
  {
    const BaseDomain_t &b = baseDomain_m;
    bd = Range<1>(b.first() + d.first()  * b.stride(),
                  b.first() + d.last()   * b.stride(),
                              d.stride() * b.stride());
  }
  
  void baseToLocal(const BaseDomain_t &bd, Interval<1> &d) const
  {
    const BaseDomain_t &b = baseDomain_m;
    PAssert(b.stride() == bd.stride());
    d = Interval<1>((bd.first() - b.first()) / b.stride(),
                    (bd.last()  - b.first()) / b.stride());
  }
  
  //============================================================
  // Data
  //============================================================

  // Our unique ID number.

  ID_t id_m;

  // A copy of the ultimate Layout object that we are viewing.

  DynamicLayout layout_m;

  // Domain of this view
  
  Domain_t domain_m;
  
  // Copy of the base domain - the domain used to subset the original
  // domain to obtain this view.
  
  BaseDomain_t baseDomain_m;
  
  // The list of all, local, and remote subdomains. Declared mutable
  // since these are evaluated in a lazy fashion in order to make
  // the taking of views inexpensive.

  // NOTE: These are subsets of the underlying Layout's lists. They
  // are NOT maps from global/local ID to nodes.
  
  mutable List_t all_m;
  mutable List_t local_m;
  mutable List_t remote_m;

  // Have we filled our subdomain lists?

  mutable bool subdomainsComputed_m;
};



/**
 * DynamicLayoutView is a Layout class that provides a view of an
 * existing DynamicLayout object. 
 *
 * To construct a DynamicLayoutView, you need an existing
 * DynamicLayout or a DynamicLayoutView and the subdomain that
 * is being viewed. This class does not have a useful default
 * constructor since it is based on an existing DynamicLayout.
 *
 * Once created, DynamicLayoutView has the same interface as
 * Layout (see Layout.h). It also provides this extra interface:
 *
 * int globalID(const Loc<1> &pos) : return the globalID
 *    of the node that contains the point.
 */

class DynamicLayoutView
{
public:
  
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumerations.
  
  enum { dimensions = 1 };

  // General public typedefs.

  typedef DynamicLayoutViewData                    LayoutData_t; 

  typedef LayoutData_t::Domain_t                   Domain_t;
  typedef LayoutData_t::BaseDomain_t               BaseDomain_t;
  typedef LayoutData_t::Context_t                  Context_t;
  typedef LayoutData_t::ID_t                       ID_t;
  
  typedef LayoutData_t::Layout_t                   Layout_t;
  typedef LayoutData_t::AllocatedDomain_t          AllocatedDomain_t;
  typedef LayoutData_t::Value_t                    Value_t;
  
  typedef LayoutData_t::List_t                     List_t;

  typedef DynamicLayoutView                        This_t;       // convenience
  typedef DynamicLayoutView                        ViewLayout_t; // convenience
  
  // Iterator through nodes. Basically the same as the vector iterator
  // except it dereferences automatically.
  
  typedef DerefIterator<Value_t>                   iterator;
  typedef ConstDerefIterator<Value_t>              const_iterator;


  //============================================================
  // Constructors
  //============================================================
  
  // Default constructor - creates an empty view of an empty layout.

  DynamicLayoutView()
  : pdata_m(new LayoutData_t())
  { }
  
  // Constructor building a DynamicLayoutView from a
  // DynamicLayout and a non-slice domain like an Interval<1> or
  // Range<1>.
  
  template <class DT>
  DynamicLayoutView(const Layout_t &layout, const Domain<1, DT> &dom)
  : pdata_m(new LayoutData_t(layout,dom))
  { }

  // Constructor building a DynamicLayoutView from another
  // DynamicLayoutView and a non-slice domain like an
  // Interval<1> or Range<1>.
  
  template <class DT>
  DynamicLayoutView(const ViewLayout_t &layout, const Domain<1, DT> &dom)
  : pdata_m(new LayoutData_t(layout,dom))
  { }

  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
  inline DynamicLayoutView(const This_t &model) 
    : pdata_m(model.pdata_m)
  { }
  
  inline This_t &operator=(const This_t &model) 
  {
    if (this != &model)
      {
        pdata_m = model.pdata_m;
      }
    return *this;
  }


  //============================================================
  // Destructor
  //============================================================

  // The actual data will be cleaned up by the DynamicLayoutData
  // destructor if all references to the data go away, so there is
  // nothing to do here.

  inline ~DynamicLayoutView() 
  { }


  //============================================================
  // Accessors
  //============================================================

  // Return ID value. Our ID is unique. However, use the one from the 
  // layout we're viewing as the baseID.
  
  inline ID_t ID() const { return pdata_m->id_m; }
  inline ID_t baseID() const { return pdata_m->layout_m.baseID(); }

  // Return that this layout has been initialized. There is no other
  // way to legally construct an object of this type.
  
  inline bool initialized() const { return true; }

  // Return the global domain in our coordinate system and in that of
  // our original domain.
  
  inline const Domain_t &domain() const 
  { 
    return pdata_m->domain_m; 
  }
  
  inline const BaseDomain_t &baseDomain() const 
  {
    return pdata_m->baseDomain_m; 
  }

  inline const Layout_t &baseLayout() const
  {
    return pdata_m->layout_m;
  }

  template <class DT>
  BaseDomain_t &localToBase(const Domain<1, DT> &dlocal, 
                            BaseDomain_t &base) const 
  {
    pdata_m->localToBase(dlocal, base);
    return base;
  }

  // Return the first index in the specified direction.
  // (Always zero since this is a zero-based engine.)
  
  inline int first(int) const { return 0; }

  // Compare to another Layout.  The layouts are the same if:
  //   1. They have the same base ID value.
  //   2. They have the same base domain.
  
  template <class Layout>
  inline bool operator==(const Layout &layout) const 
  {
    return (baseID() == layout.baseID() && 
            baseDomain() == layout.baseDomain());
  }

  // Compare for inequality.

  template <class Layout>
  inline bool operator!=(const Layout &layout) const
  {
    return !(*this == layout);
  }


  //============================================================
  // Patch lookup
  //============================================================

  // Return the globalID of the node containing the given element.

  int globalID(const Loc<1> &, Loc<1> &) const;
  int globalID(int, Loc<1> &) const;

  //============================================================
  // Iterators
  //============================================================

  // Return begin and end iterators for the list of all subdomains;
  // take care to compute subdomains first
  
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

  // Return begin and end iterators for the list of all subdomains;
  // take care to compute subdomains first
  
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

  // Return begin and end iterators for the list of all subdomains;
  // take care to compute subdomains first
  
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


  //============================================================
  // Touch methods
  //============================================================

  // Find all subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.
  // Return the number of touching elements. This version of touches
  // can build either pointers or objects.
  
  template <class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, 
              const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d,o,ctag);
  }

  // Find local subdomains that touch on a given domain, and insert
  // the intersection of these subdomains into the given output
  // iterator.  Return the number of touching elements. This version
  // of touches can build either pointers or objects. 

  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesLocal(const OtherDomain &d, OutIter o, 
                   const ConstructTag &ctag) const 
  {
    return pdata_m->touches(d, o, ctag);
  }

  // Find remote subdomains that touch on a given domain, and insert
  // the intersection of these subdomains into the given output
  // iterator.  Return the number of touching elements. This version
  // of touches can build either pointers or objects. 

  template <class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesRemote(const OtherDomain &, OutIter, 
                    const ConstructTag &) const 
  {
    return 0;
  }

  // Find all/local/remote subdomains that touch on a given domain,
  // and insert the intersection of these subdomains into the given
  // output iterator.  Return the number of touching elements.  These
  // versions of touches can build only objects objects.

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

  //============================================================
  // Output
  //============================================================
    
  // Print a DynamicLayoutView on an output stream

  template <class Ostream>
  void print(Ostream &ostr) const 
  {
    ostr << "DynamicLayoutView " << ID() << " on global domain " 
      << domain() << ":" << '\n';
    ostr << "   Base ID:          " << baseID() << '\n';
    ostr << "   Base domain:      " << baseDomain() << '\n';
    ostr << "   Total subdomains: " << sizeGlobal() << '\n';
    ostr << "   Local subdomains: " << sizeLocal() << '\n';
    ostr << "  Remote subdomains: " << sizeRemote() << '\n';
    const_iterator a;
    for (a = beginGlobal(); a != endGlobal(); ++a)
      ostr << "  Global subdomain = " << *a << '\n';
    for (a = beginLocal(); a != endLocal(); ++a)
      ostr << "   Local subdomain = " << *a << '\n';
    for (a = beginRemote(); a != endRemote(); ++a)
      ostr << "  Remote subdomain = " << *a << '\n';
  }

private:

  //============================================================
  // Private utility functions
  //============================================================

  // Fill our subdomain lists.

  void computeSubdomains() const { pdata_m->computeSubdomains(); }
   
  //============================================================
  // Data
  //============================================================

  // The data is stored in a RefCounted class to simplify memory
  // management.  This is probably not as important for ViewLayout
  // classes as for Layout classes, but we do it for
  // consistency. Currently ViewLayouts are not observed directly by
  // anyone. Of course, the Layout that we have a copy of is observed.

  RefCountedPtr<LayoutData_t> pdata_m;
  
};


//============================================================
// DynamicLayoutViewData inline method definitions
//============================================================

template <class DT>
DynamicLayoutViewData::DynamicLayoutViewData(const ViewLayout_t &layout,
					     const Domain<1, DT> &dom)
  : id_m(Unique::get()),
    layout_m(layout.baseLayout()), 
    domain_m(dom.unwrap().length()),
    subdomainsComputed_m(false)
{
  // The layout passed in must be initialized. 

  PAssert(layout_m.initialized());

  // The domain we're passing in must be contained in the base layout.

  PAssert(contains(layout.domain(), dom.unwrap()));

  // need to compute our base domain from given view and domain

  baseDomain_m = layout.localToBase(dom.unwrap(),baseDomain_m);
}


//============================================================
// DynamicLayoutView inline method definitions
//============================================================

//-----------------------------------------------------------------------------
//
// Return the globalID to the given element, expressed either as
// a single Loc or as an int.
//
// These convert to the coordinate in the base domain and then 
// ask the underlying layout to perform the calculation.
//-----------------------------------------------------------------------------

inline int
DynamicLayoutView::
globalID(const Loc<1> &loc, Loc<1> &oloc) const
{
  return pdata_m->globalID(loc,oloc);
}

inline int
DynamicLayoutView::
globalID(int i0, Loc<1> &oloc) const
{
  return pdata_m->globalID(i0,oloc);
}


//=============================================================================
// NewDomain1 traits classes for DynamicLayout and DynamicLayoutView
//=============================================================================

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a DynamicLayout.
//
//-----------------------------------------------------------------------------

template <>
struct NewDomain1<DynamicLayout>
{
  typedef DynamicLayout &Type_t;

  inline static Type_t combine(const DynamicLayout &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// This is so an array can be initialized with a DynamicLayoutView.
//
//-----------------------------------------------------------------------------

template <>
struct NewDomain1<DynamicLayoutView>
{
  typedef DynamicLayoutView &Type_t;

  inline static Type_t combine(const DynamicLayoutView &a)
    {
      return const_cast<Type_t>(a);
    }
};

//-----------------------------------------------------------------------------
//
// ostream inserters for DynamicLayout and DynamicLayoutView:
//
//-----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &ostr, 
			 const DynamicLayout &layout);

std::ostream &operator<<(std::ostream &ostr, 
			 const DynamicLayoutView &layout);


// } // namespace POOMA

#endif // POOMA_LAYOUT_DYNAMICLAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DynamicLayout.h,v $   $Author: richard $
// $Revision: 1.39 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
