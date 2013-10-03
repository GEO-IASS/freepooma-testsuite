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
// Functions:
//   makeUniformGrid utility function (in anonymous namespace - not global)
//   DynamicLayout non-inline non-template member function definitions.
//-----------------------------------------------------------------------------

#include "Layout/DynamicLayout.h"

#include "Domain/Contains.h"
#include "Domain/Grid.h"
#include "Partition/UniformGridPartition.h"
#include "Partition/GridPartition.h"
#include "Utilities/PAssert.h"

// For multiple-context synchronization:

#include "PETE/PETE.h" // for OpSumAssign
#include "Tulip/ReduceOverContexts.h"
#include "Tulip/PatchSizeSyncer.h"

#include <iostream>

//============================================================
// Utility functions
//============================================================

//-----------------------------------------------------------------------------
//
// Grid<1> <anonymous>::makeUniformGrid(const Interval<1> gdom, int nblocks)
//
// Utility function that constructs a near-uniform patch decomposition
// of the domain gdom, having nblocks subdomains. The decomposition is
// returned as a Grid<1> object. This Grid object is not used as a 
// "domain" per se. Rather, we use its blockIterator to generate the 
// patches that tile the global domain. This difference is important 
// to note since the blockIterator patches are Interval<1>(p1,p2-1),
// where p1 and p2 are adjacent points in the domain. Thus, we must 
// construct a grid-object whose upper end point is one more than the 
// last point of the input domain. That is, if we are subdividing the
// interval [0,19] into two blocks of 10 elements, then we want the 
// Grid's patch iterator to return [0,9] and [10,19]. This requires a 
// Grid with the points [0,10,20].
//
//-----------------------------------------------------------------------------

namespace { // anonymous - local linkage only

Grid<1> makeUniformGrid(const Interval<1> &gdom, int nblocks)
{
  PAssert(!gdom.empty());
  
  // First calculate the approximate block size and the remainder. 
  
  long blocksize = gdom.size() / nblocks;
  long remainder = gdom.size() % nblocks;
  
  if (remainder == 0)
    {
      // If the remainder is zero, then we can construct a completely 
      // uniform grid, which can be done by constucting it with a Range 
      // object (having one extra point).

      Range<1> ret(gdom.first(), gdom.last() + 1, blocksize);
      return Grid<1>(ret);
    }
  else
    {      
      // If the remainder is non-zero, we make the last "remainder"
      // sub-blocks one element longer than the rest. 
      
      IndirectionList<int> vertexlist(nblocks + 1);
      int j;
      vertexlist(0) = gdom.first();
      for (j = 1; j < nblocks + 1 - remainder; ++j)
        {
          vertexlist(j) = vertexlist(j-1) + blocksize;
	}
      for (j = nblocks + 1 - remainder; j < nblocks + 1; ++j)
        {
          vertexlist(j) = vertexlist(j-1) + blocksize + 1;
        }
      PAssert(vertexlist(nblocks) == gdom.last() + 1);
      return Grid<1>(vertexlist);
    }
}

} // close anonymous namespace

//============================================================
// DynamicLayoutData non-inline non-template method definitions
//============================================================

//-----------------------------------------------------------------------------
//
// DynamicLayoutData::DynamicLayoutData()
//
// Default constructor for DynamicLayoutData. Initializes this object
// to look like an "empty" layout, with no patches and an empty domain.  
// The "initialize" method can be used to complete the initialization.
//
//-----------------------------------------------------------------------------

DynamicLayoutData::DynamicLayoutData()
  : Observable_t(*this), 
    ID_m(Unique::get()), 
    dirtyLayout_m(true)
{ }

//-----------------------------------------------------------------------------
//
// DynamicLayoutData::~DynamicLayoutData()
//
// Destructor for DynamicLayoutData. This removes all the nodes and
// then the Observable destructor will tell all users of this layout that
// it is going away.
//
//-----------------------------------------------------------------------------

DynamicLayoutData::~DynamicLayoutData()
{
  // Delete existing nodes and clear all the lists.
  
  for (int i = 0; i < all_m.size(); ++i) delete all_m[i];
  all_m.clear();
  local_m.clear();
  remote_m.clear();
}

//-----------------------------------------------------------------------------
//
// Create new elements by extending the current domain of the specified 
// local patch by the requested number of elements. 'local' means on 
// this same context.  The patch is referred to by a local index, 
// from 0 ... # local patches - 1. If patch=-1, create elements in the 
// last local patch.
//
// All observers are notified of the change, then we change our
// domain value.
//
//-----------------------------------------------------------------------------

void DynamicLayoutData::create(CreateSize_t num, PatchID_t patch)
{
  PAssert(num >= 0);
  if (num == 0) return;

  // If the patch number is < 0, change it to the last local patch.

  if (patch < 0) patch = local_m.size() - 1;
  PAssert(patch < local_m.size());

  // Let all users know of the create request.

  notify(CreateEvent(num, patch));

  // Modify the domain for this local patch.  When sync is called,
  // everything else will get updated.

  addElements(local_m[patch]->domain(), num);

  // Note that we will need to rebuild things.

  dirtyLayout_m = true;
}

//-----------------------------------------------------------------------------
//
// Perform a "multiple patch" copy, using a list of IndirectionList's
// for a set of source patches, and an IndirectionList giving the
// patch ID for the source patches.  Copy data into the destination
// patch.  The source and desination patches must be specified, this
// is only for "zero-based" index lists.  If the last argument is
// true, storage is created at the end, otherwise elements are
// just copied to the end of the existing storage.
//
//-----------------------------------------------------------------------------

void DynamicLayoutData::
copy(const IndirectionList<IndirectionList<int> > &lists,
     const IndirectionList<int> &fromlist,
     PatchID_t toPatch,
     bool docreate)
{
  // If the toPatch number is < 0, change it to the last local patch.
  // Is this really a useful default???
  
  if (toPatch < 0) toPatch = local_m.size() - 1;
  PAssert(toPatch < local_m.size());

  // Let all users know of the copy request.

  notify(CopyPatchEvent(lists, fromlist, toPatch, docreate));

  // Modify the domain for this local patch.  When sync is called,
  // everything else will get updated.
  
  // Why is there a bool??? Why not just figure out if new allocation
  // is needed and if so, do it, if not, don't. (JAC)

  if (docreate)
    {
      int np = lists.size();
      int created = 0;
      for (int i = 0; i < np; ++i)
	created += lists(i).size();

      addElements(local_m[toPatch]->domain(), created);

      // Note that we will need to rebuild things.
      
      dirtyLayout_m = true;
    }
}

//-----------------------------------------------------------------------------
//
// Sync up the layout with any other contexts, taking into account
// that other contexts may have performed create/destroy operations.
// This will reset all the local domains to be properly contiguous,
// and let all engine's using this layout reset their domains.
//
// The multiple-context stuff is done by syncGlobalDomains below.
//
//-----------------------------------------------------------------------------

void DynamicLayoutData::sync()
{
  int nContexts = Pooma::contexts();
  
  // First check if the layout is (globally) dirty or not.
  
  if (nContexts == 1) // no communication required
    {
      if (!initialized() || !dirty()) return;
    }
  else
    {
      typedef ReduceOverContexts<int, OpAddAssign> GlobalSum_t;

      // Do a global reduction on the initialized and dirty flags

      int globalInitialized;
      GlobalSum_t(int(initialized())).broadcast(globalInitialized);

      // They'd better either all be initialized or not.

      PAssert(globalInitialized == 0 || globalInitialized == nContexts);

      if (globalInitialized == 0) return;

      int globalDirty;
      GlobalSum_t(int(dirty())).broadcast(globalDirty);

      if (globalDirty == 0) return;
    }
      
  // Recalculate and renumber the domains. If we are multi-context,
  // this does the global calculations to fix the global decomposition
  // of the current total domain.

  calcDomains();

  // Recalculate the domain maps, if necessary.  We need to do this
  // now since we'll need to call globalID() routines from threads in
  // later operations ... the alternative is to make checking and rebuilding
  // of the domain maps a mutually-exclusive operation.

  calcMaps();

  // The domains & maps are up-to-date, so clear our dirty flag.

  dirtyLayout_m = false;

  // Notify all the users that they can sync up their patches.

  notify(SyncEvent());
}

void DynamicLayoutData::syncGlobalDomains()
{
  // First we build a Grid<1> object that represents the local
  // patches. The points in this Grid are the "first" points
  // for every subdomain, plus one past the last point of the
  // last subdomain (consecutive pairs of points can be considered 
  // as begin-end pairs defining half open "intervals" in the 
  // STL sense). We return the total number of elements.

  int nlocal = local_m.size();
  IndirectionList<int> lgdata(nlocal+1);

  // Since we're not guaranteed that there are *any* local elements,
  // we construct a Grid that is zero-based and let the PatchSizeSyncer
  // figure out the final domains based on context ordering.
  
  long pos = 0;
  for (int i = 0; i < nlocal; ++i)
    {
      lgdata(i) = pos;
      pos += local_m[i]->domain().size();
    }
  lgdata(nlocal) = pos;
      
  Grid<1> localGrid(lgdata);

  // Now initialize a PatchSizeSyncer object with the local data
  // and call the calcGlobalGrid method to do the communication
  // and return a global Grid object that represents the
  // global patch decomposition.

  Grid<1> globalGrid;
  Context_t myContext = Pooma::context();
  Pooma::PatchSizeSyncer(myContext,localGrid).calcGlobalGrid(globalGrid);

  // The number of patches in a dynamic layout is fixed, so the number
  // of points in this grid had better match the number of patches
  // (+1).

  PAssert(globalGrid.size() == all_m.size() + 1);
  PAssert(domain_m.empty() || globalGrid.first() == domain_m.first());

  // Finally, use the consecutive points in the Grid to reset the
  // domains for all of our Nodes. This loop relies on the global
  // patch ordering being such that their subdomains are contiguous.
  // Note the special handling for empty patches.

  // JCC: The assumption about Node ordering of the original code here
  // does not work in general.  Instead, I am assuming that the global
  // Grid information is ordered by context.  Thus we have domains for 
  // all of the Nodes on context 0, followed by domains for all the 
  // Nodes on context 1, etc.  This ordering is independent of the 
  // partitioning scheme or context mapper used.
  // Later we should add a map between the node ordering in the all_m
  // node list and the ordering based on context number for efficiency.

  Context_t numContexts = Pooma::contexts();
  PatchID_t numNodes = all_m.size();
  int c, i, j = 0;
  for (c = 0; c < numContexts; ++c)
    {
      for (i = 0; i < numNodes; ++i)
	{
	  if (all_m[i]->context() == c)
	    {
	      int begin = globalGrid(j);
	      int end   = globalGrid(j+1);
	      PAssert(begin <= end);

	      Domain_t dom = Pooma::NoInit();
	      if (begin < end) 
		dom = Domain_t(begin, end - 1); // [begin,end) domain
	      else
		dom = Domain_t();               // empty domain

	      all_m[i]->setDomain(dom);
	      all_m[i]->setAllocated(dom);
	      j++;
	    }
	}
    }

  int begin = globalGrid.first();
  int end   = globalGrid.last();
  if (begin < end)
    domain_m = Domain_t(begin, end - 1);
  else
    domain_m = Domain_t();
}

//-----------------------------------------------------------------------------
//
// void DynamicLayoutData::calcDomains()
//
// Calculates the total domain of each patch and this total layout, since
// this can change due to dynamic operations.
//
//-----------------------------------------------------------------------------

void DynamicLayoutData::calcDomains()
{
  // This does not check the dirty flag - that should be done prior
  // to calling this. Wasn't a big deal for single-context stuff, but
  // now this is a global reduction, so try to only do it once.

  // We scan through the local domains, and adjust their starting
  // offsets to be contiguous. We will start everyone off at
  // domain_m.first(). This way we can skip recalculating the domains
  // if we're on a single context. If we're on multiple contexts, it
  // doesn't matter what the set of local domains start with prior to
  // calling syncGlobalDomains().

  int first = domain_m.first();
  CreateSize_t pos = first, len = 0;
  for (int i = 0; i < local_m.size(); ++i)
    {
      Domain_t dom = local_m[i]->domain();
      len = dom.length();
      if (len > 0)
	{
	  dom = Domain_t(pos, pos + len - 1);
	  pos += len;
	}

      // Give this new domain to the Node ... it will have the same
      // size, but a possibly different initial offset.

      local_m[i]->setDomain(dom);
      local_m[i]->setAllocated(dom);
    }

  // Update the remote and total domains...
  // Add a check here for a replicated mapping of Nodes.
  // In this case, no global synchronization is needed.

  if (Pooma::contexts() > 1 && all_m[0]->context() != -1) 
    {
      syncGlobalDomains();
    }
  else
    {
      // Just update the total domain.

      if (pos == first)
	domain_m = Domain_t();
      else
	domain_m = Interval<1>(first, pos - 1);
    }

}


//-----------------------------------------------------------------------------
//
// void DynamicLayoutData::calcMaps()
//
// Calculates the DomainMaps's for this object, based on the current
// settings for the blocks, since this can change due to dynamic operations.
//
//-----------------------------------------------------------------------------

void DynamicLayoutData::calcMaps() const
{
  // Initialize the map...

  // Clear out any existing info

  map_m.zap();

  // If this is empty, there is nothing to do.

  if (!domain_m.empty())
    {
      // Initialize the map and then add each non-empty subdomain
      // to the map along with its global ID.
      
      map_m.initialize(domain_m);

      for (int j = 0; j < all_m.size(); ++j)
	{
	  const Interval<1> &blockDom = all_m[j]->domain();
          PAssert(j == all_m[j]->globalID());
	  if (!blockDom.empty())
	    {
	      typedef DomainMap<Interval<1>,int>::Value_t Val_t;
	      map_m.insert(Val_t(blockDom, j));
	    }
	}

      // Update the DomainMap

      map_m.update();
    }
}

//-----------------------------------------------------------------------------
//
// globalID takes a position within the domain of the layout, and returns
// the global ID for that node.
//
//-----------------------------------------------------------------------------

int DynamicLayoutData::globalID(const Loc<1> &loc) const
{
  PAssert(!dirtyLayout_m);
  
  // Make sure the point is in our domain.

  PAssert(contains(domain_m, loc));

  // Find the position of the point.

  typedef DomainMapTouchIterator<Interval<1>,int> MapIterator_t;
  
  MapIterator_t dmti = (map_m.touch(Interval<1>(loc))).first;
  PAssert(dmti != MapIterator_t()); // Default constructor produces end iterator

  // Return the offset (dmti dereferences to an int, which is the globalID): 

  return *dmti;
}

int DynamicLayoutData::globalID(int i0) const
{
  // Call the Loc version.
  return globalID(Loc<1>(i0));
}

//============================================================
// DynamicLayout non-inline non-template method definitions
//============================================================

//-----------------------------------------------------------------------------
//
// Constructors and initialize methods
// 
// See comments in the class definition (in DynamicLayout.h)
//
//-----------------------------------------------------------------------------

DynamicLayout::DynamicLayout()
  : Observable<This_t>(*this),
    pdata_m(new LayoutData_t()) 
{ 
  pdata_m->attach(*this);
}
  
DynamicLayout::DynamicLayout(const Domain_t &gdom)
  : Observable<This_t>(*this), 
    pdata_m( new LayoutData_t( gdom, GridPartition<1>(), UniformMapper() ) )
{
  pdata_m->attach(*this);
}

DynamicLayout::DynamicLayout(const Domain_t &gdom, 
			     int blocks)
  : Observable<This_t>(*this)
{
  UniformMapper cmap(blocks);
  if (!gdom.empty()) 
    {
      Grid<1> grid = makeUniformGrid(gdom,blocks);
      GridPartition<1> gpar(grid);
      pdata_m = new LayoutData_t(gdom, gpar, cmap);
    }
  else
    {
      Loc<1> decomp(blocks);
      GridPartition<1> gpar(decomp);
      pdata_m = new LayoutData_t(gdom, gpar, cmap);
    }
  pdata_m->attach(*this);
}

DynamicLayout::DynamicLayout(const Grid<1> &grid)
  : Observable<This_t>(*this)
{
  Domain_t gdom(grid.first(),grid.last()-1);
  GridPartition<1> gpar(grid);
  UniformMapper cmap(gpar);
  pdata_m = new LayoutData_t(gdom, gpar, cmap);  
  pdata_m->attach(*this);
}

DynamicLayout::DynamicLayout(const This_t &model) 
  : Observable<This_t>(*this),
    pdata_m(model.pdata_m)
{ 
   pdata_m->attach(*this);
}

//-----------------------------------------------------------------------------
//
// assignment operator for DynamicLayout
//
//-----------------------------------------------------------------------------

DynamicLayout &DynamicLayout::operator=(const This_t &model)
{
  if (this != &model)
    {
      pdata_m->detach(*this);
      pdata_m = model.pdata_m;
      pdata_m->attach(*this);
    }

  return *this;
}


//-----------------------------------------------------------------------------
//
// Initialize methods for DynamicLayout
//
//-----------------------------------------------------------------------------

void DynamicLayout::initialize(const Domain_t &gdom)
{
  pdata_m->initialize(gdom, GridPartition<1>(), UniformMapper());
}

void DynamicLayout::initialize(const Domain_t &gdom,
			       int blocks)
{
  UniformMapper cmap(blocks);
  if (!gdom.empty()) 
    {
      Grid<1> grid = makeUniformGrid(gdom,blocks);
      GridPartition<1> gpar(grid);
      pdata_m->initialize(gdom, gpar, cmap);
    }
  else
    {
      Loc<1> decomp(blocks);
      GridPartition<1> gpar(decomp);
      pdata_m->initialize(gdom, gpar, cmap);
    }
}

void DynamicLayout::initialize(const Domain_t &gdom,
			       const Grid<1> &grid)
{
  GridPartition<1> gpar(grid);
  UniformMapper cmap(gpar);
  pdata_m->initialize(gdom, gpar, cmap);
}

void DynamicLayout::initialize(const Grid<1> &grid)
{
  Domain_t gdom(grid.first(),grid.last()-1);
  GridPartition<1> gpar(grid);
  UniformMapper cmap(gpar);
  pdata_m->initialize(gdom, gpar, cmap);  
}

//-----------------------------------------------------------------------------
//
// Ostream inserter definitions...
//
//-----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &ostr, 
			 const DynamicLayout &layout)
{
  layout.print(ostr);
  return ostr;
}

std::ostream &operator<<(std::ostream &ostr, 
			 const DynamicLayoutView &layout)
{
  layout.print(ostr);
  return ostr;
}


// } // namespace POOMA

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DynamicLayout.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
