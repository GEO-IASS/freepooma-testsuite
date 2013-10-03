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
//   GridLayout<Dim> template definitions.
//-----------------------------------------------------------------------------

#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/Grid.h"
#include "Layout/GridLayout.h"
#include "Partition/GridPartition.h"
#include "Threads/PoomaSmarts.h"
#include "Utilities/PAssert.h"
#include <vector>
#include <iostream>

//============================================================
// GridLayoutData non-inline method definitions
//============================================================

template<int Dim>
GridLayoutData<Dim>::GridLayoutData()
 : LayoutBaseData<Dim>(false,false,
		   GuardLayers_t(0),GuardLayers_t(0),
		   Interval<Dim>(),Interval<Dim>()),
    Observable<GridLayoutData>(*this) 
{
  for (int d=0; d < Dim; ++d)
    this->firsti_m[d] = this->firste_m[d] = blockStrides_m[d] = 0;
}


//-----------------------------------------------------------------------------
//
// template<int Dim>
// template<class Partitioner>
// GridLayoutData<Dim>::GridLayoutData(Domain_t, Partitioner)
//
// Constructor for GridLayoutData that takes a global domain and a
// partitioner.  The global domain can be empty.  This will invoke the
// initialize method to set up the layout state and invoke the partitioner.
//
//-----------------------------------------------------------------------------
template<int Dim>
template<class Partitioner> 
GridLayoutData<Dim>::GridLayoutData(const Grid<Dim> &gdom, 
				    const Partitioner &gpar,
				    const ContextMapper<Dim> &cmap)
  : LayoutBaseData<Dim>(false,false,
			GuardLayers_t(0),GuardLayers_t(0),
			Interval<Dim>(), Interval<Dim>()),
    Observable<GridLayoutData>(*this)  
{
  initialize(gdom, gpar, cmap);
} 

template<int Dim>
template<class Partitioner>
GridLayoutData<Dim>::GridLayoutData(const Domain_t &gdom,
				    const Partitioner &gpar,
				    const ContextMapper<Dim> &cmap)
  : LayoutBaseData<Dim>(false,false,
			GuardLayers_t(0),GuardLayers_t(0),
			gdom,gdom),
    Observable<GridLayoutData>(*this)  
{
  for (int d=0; d < Dim; ++d)
    this->firsti_m[d] = this->firste_m[d] = blockStrides_m[d] = 0;
  initialize(gdom, gpar, cmap);
}

template<int Dim>
GridLayoutData<Dim>::~GridLayoutData()
{
  for (typename List_t::iterator a = this->all_m.begin(); a != this->all_m.end(); ++a)
    delete (*a);
} 

template<int Dim>
template <class Partitioner>
inline void GridLayoutData<Dim>::initialize(const Grid<Dim> &gdom, 
		       const Partitioner &gpar,
		       const ContextMapper<Dim> &cmap)
{
  Domain_t idom = Pooma::NoInit();
  for (int d=0; d < Dim; ++d)
    idom[d] = Interval<1>(gdom[d].first(), gdom[d].last() - 1);
  initialize(idom, gpar, cmap);
}


template<int Dim>
template<class Partitioner>
inline void GridLayoutData<Dim>::initialize(const Domain_t &gdom,
				     const Partitioner &gpar,
				     const ContextMapper<Dim> &cmap)
{
  int i;

  // This will work with grid (and simpler) partitioners.

  CTAssert(Partitioner::gridded);

  // delete existing nodes and clear all the lists

  if (this->all_m.size() > 0)
    {
      for (i=0; i < this->all_m.size(); ++i)
	delete this->all_m[i];
      this->all_m.clear();
      this->local_m.clear();
      this->remote_m.clear();
    }

  // After this, we will need to rebuild things.

  dirtyLayout_m = true;

  // Initially, our total and owned domains are the same.

  this->domain_m = gdom;
  this->innerdomain_m = gdom;

  // Examine the partitioner for info about guard cells.  Change our
  // domains if necessary, and save guard cell info for later.

  this->hasInternalGuards_m = (gpar.hasInternalGuards() && gpar.maxSize() > 1);
  if (this->hasInternalGuards_m)
    {
      this->internalGuards_m = gpar.internalGuards();
    }

  this->hasExternalGuards_m = (gpar.hasExternalGuards() && ! this->domain_m.empty());
  if (this->hasExternalGuards_m)
    {
      this->externalGuards_m = gpar.externalGuards();
      GuardLayers<Dim>::addGuardLayers(this->domain_m, this->externalGuards_m);
    }

  // Get the number of blocks in each dimension from the partitioner.

  this->blocks_m = gpar.blocks();

  // Determine the initial offsets for each dimension.  This stays fixed,
  // even though the total size can change from dynamic operations.  Also
  // calculate the strides for converting an i,j,k index to a serialized
  // block index.

  for (i=0; i < Dim; ++i)
    {
      if (!this->domain_m[i].empty())
	this->firsti_m[i] = this->domain_m[i].first();
      blockStrides_m[i] = (i == 0 ? 1 :
			   blockStrides_m[i-1] * this->blocks_m[i-1].first());
    }

  // Invoke the partitioner, which will pass back the all_m list. 
  
  gpar.partition(this->innerdomain_m, this->all_m, cmap);

  typename List_t::iterator start = this->all_m.begin();
  typename List_t::iterator end   = this->all_m.end();
  
  for ( ; start!=end ;++start )
    {
      if( (*start)->context() == Pooma::context() 
	  ||(*start)->context() == -1 )          // per SASMITH's request
	this->local_m.push_back(*start);
      else
	this->remote_m.push_back(*start);
    }


  // Initially we calculate the domain maps.

  calcMaps();
  calcAllocMaps();

  // Calculate what we need to do in a fill-guard-cell operation.

  calcGCFillList();
}

//-----------------------------------------------------------------------------
//
// template<int Dim>
// void GridLayoutData<Dim>::initialize
//
// Used by an I/O or data management entity to initialize the layout based
// on detailed state information previously stored. As in the case of the
// initializer with the partitioner argument, this method will call 'addDomain'
// to add in the new domains it creates and will initialize
// guard cell info, etc.
//
//-----------------------------------------------------------------------------

template<int Dim>
void GridLayoutData<Dim>::initialize(const Domain_t& idom,
				 const List_t& nodes,
				 const Loc<Dim>& blocks,
				 bool hasIG, bool hasEG,
				 const GuardLayers_t& ig,
				 const GuardLayers_t& eg)
{
  int i;

  // delete existing nodes and clear all the lists

  if (this->all_m.size() > 0)
    {
      for (i=0; i < this->all_m.size(); ++i)
	delete this->all_m[i];
      this->all_m.clear();
      this->local_m.clear();
      this->remote_m.clear();
    }

  // After this, we will need to rebuild things.

  dirtyLayout_m  = true;

  // Initially, our total and owned domains are the same.

  this->domain_m = idom;
  this->innerdomain_m = idom;

  // Examine the info about guard cells.  Change our domains if necessary,
  // and save guard cell info for later.

  this->hasInternalGuards_m = hasIG;
  if (this->hasInternalGuards_m)
    {
      this->internalGuards_m = ig;
    }

  this->hasExternalGuards_m = (hasEG && ! this->domain_m.empty());
  if (this->hasExternalGuards_m)
    {
      this->externalGuards_m = eg;
      GuardLayers<Dim>::addGuardLayers(this->domain_m, this->externalGuards_m);
    }

  // Get the number of blocks in each dimension.

  this->blocks_m = blocks;

  // Determine the initial offsets for each dimension.  This stays fixed,
  // even though the total size can change from dynamic operations.  Also
  // calculate the strides for converting an i,j,k index to a serialized
  // block index.

  for (i=0; i < Dim; ++i)
    {
      if (!this->domain_m[i].empty())
	this->firsti_m[i] = this->domain_m[i].first();
      blockStrides_m[i] = (i == 0 ? 1 :
			   blockStrides_m[i-1] * this->blocks_m[i-1].first());
    }

  // Assign the given list of nodes to the total list.
  this->all_m= nodes;

  // Iterate through the complete list of nodes provided and assign to the
  // appropriate subcategories.

  typename List_t::iterator start = this->all_m.begin();
  typename List_t::iterator end   = this->all_m.end();
  
  for ( ; start!=end ;++start )
    {
      if( (*start)->context() == Pooma::context() ||
	  (*start)->context() == -1 )
	this->local_m.push_back(*start);
      else
	this->remote_m.push_back(*start);
    }

  // Calculate the domain maps.

  calcMaps();
  calcAllocMaps();

  // Calculate what we need to do in a fill-guard-cell operation.

  calcGCFillList();
}

template<int Dim>
void GridLayoutData<Dim>::sync()
{
  PAssert(Dim==1);  

  // Nothing to do if our layout is not dirty.

  if (!this->initialized() || !dirty())
    return;

  // Recalculate and renumber the domains.

  calcDomains();

  // Recalculate the domain maps, if necessary.  We need to do htis
  // now since we'll need to call globalID() routines from threads in
  // later operations ... the alternative is to make checking and rebuilding
  // of the domain maps a mutually-exclusive operation.

  calcMaps();

  // Do not recalculate allocated maps or guard cell fill lists
  // here; these can get rebuilt as needed.

  // Notify all the users that they can sync up their patches.

  this->notify(SyncEvent());
}


//-----------------------------------------------------------------------------
//
// template<int Dim>
// void GridLayout<Dim>::calcGCFillList()
//
// Calculates the cached information needed by MultiPatch Engine to fill
// the guard cells. 
//
//-----------------------------------------------------------------------------

template<int Dim>
void GridLayoutData<Dim>::calcGCFillList()
{
  int d;

  // We can just return if we don't need this.

  if (!this->initialized() || !this->hasInternalGuards_m)
    return;

  // Clear our the existing list

  this->gcFillList_m.clear();

  // We want to create the list in such a manner that
  // all communication in a particular direction is done first,
  // allowing parallelism with the least amount of contention
  // for patches. Thus we have an outer loop over numPatches, doing
  // the upward copies first, then the downward copies.

  int numPatches = this->all_m.size();
  this->gcFillList_m.reserve(2*Dim*this->local_m.size()); 

  // Make sure we have the same number of patches as blocks in the grid
  // (this is just a sanity check).

  PAssert(numPatches == blockStrides_m[Dim-1] * this->blocks_m[Dim-1].first());

  // Create an Interval for use in iterating over the "grid" blocks.

  Interval<Dim> grid;
  for (d=0; d < Dim; ++d)
    grid[d] = Interval<1>(this->blocks_m[d].first() + 1);

  // When we go to multiple contexts, this algorithm will still
  // work if the local's always form a block that is also
  // stored in fortran storage order.
    
  for (d=0; d < Dim; ++d) 
    {
      if (this->internalGuards_m.lower(d) > 0)
	{	
	  typename Interval<Dim>::blockIterator start = grid.beginBlock();
	  typename Interval<Dim>::blockIterator end = grid.endBlock();

	  for (; start != end; ++start)
	    {
	      // Get the serialized index of the patch we're looking at

	      int sourceID = start.index();

	      // Find the indices of the patch below this one

	      Loc<Dim> tmp = start.point();
	      tmp[d] += 1;	// we add one here since we're looking down
		
	      // If this is an edge in this direction/dimension, continue.

	      if (!(tmp[d] >= this->blocks_m[d] || tmp[d].first() < 0))
		{
		  // It's in the system ... find the serialized index
		  // of the patch below this one (where tmp points).

		  int destID = blockIndex(tmp);

		  // Only add in the guard cell fill info if both patches
		  // have non-empty domains.

		  if (!(this->all_m[sourceID]->domain().empty() ||
			this->all_m[destID]->domain().empty()))
		    {
		      // check that we didn't screw up the indexing.

		      PAssert(destID >= 0 && destID < numPatches);

		      // Create a domain that just covers the guard region.

		      Domain_t gcdom(this->all_m[sourceID]->allocated());
		      int max = this->all_m[sourceID]->domain()[d].last();
		      int min = max - this->internalGuards_m.lower(d) + 1;
		      gcdom[d] = Interval<1>(min, max);

		      // Now, push IDs and source into cache...
		
		      this->gcFillList_m.push_back(GCFillInfo_t(gcdom, sourceID, destID, d*2));
		    }
		}
	    }	
	}

      // now go in the other direction, 

      if (this->internalGuards_m.upper(d) > 0)
	{	
	  typename Interval<Dim>::blockIterator start = grid.beginBlock();
	  typename Interval<Dim>::blockIterator end = grid.endBlock();

	  for (; start != end; ++start)
	    {  
	      // Get the serialized index of the patch we're looking at

	      int sourceID = start.index();

	      // Find the indices of the patch above this one

	      Loc<Dim> tmp = start.point();
	      tmp[d] -= 1;	// we subtract one here since we're looking up

	      // If this is an edge in this direction/dimension, continue.

	      if (!(tmp[d] >= this->blocks_m[d] || tmp[d].first() < 0))
		{
		  // It's in the system ... find the serialized index
		  // of the patch below this one (where tmp points).

		  int destID = blockIndex(tmp);

		  // Only add in the guard cell fill info if both patches
		  // have non-empty domains.

		  if (!(this->all_m[sourceID]->domain().empty() ||
			this->all_m[destID]->domain().empty()))
		    {
		      // check that we didn't screw up the indexing.

		      PAssert(destID < numPatches);

		      // Create a domain that just covers the guard region.

		      Domain_t gcdom(this->all_m[sourceID]->allocated());  
		      int min = this->all_m[sourceID]->domain()[d].first();
		      int max = min + this->internalGuards_m.upper(d) - 1;
		      gcdom[d] = Interval<1>(min, max);  
		
		      // Now, push IDs and source into cache...
		
		      this->gcFillList_m.push_back(GCFillInfo_t(gcdom, sourceID, destID, d*2+1));
		    }
		}
	    }
	}
    }
}


//-----------------------------------------------------------------------------
//
// template<int Dim>
// void GridLayoutData<Dim>::calcDomains()
//
// Calculates the total domain of each patch and this total layout, since
// this can change due to dynamic operations.
//
//-----------------------------------------------------------------------------

template<int Dim>
void GridLayoutData<Dim>::calcDomains()
{
  // Don't need to do this if we're not initialized or not dirty.

  if (!this->initialized() || !dirty())
    return;

  // This will only work at present if there are no guard cells.

  PAssert(!(this->hasInternalGuards_m || this->hasExternalGuards_m));

  // And it will only work for 1D.

  PAssert(Dim == 1);

  // We scan through the domains, and adjust their starting offsets to be
  // contiguous.

  CreateSize_t sizes = this->firsti_m[0];
  bool allempty = true;
  for (int i=0; i < this->all_m.size(); ++i)
    {
      Domain_t dom = this->all_m[i]->domain();
      if (!dom[0].empty())
	{
	  dom[0] = Interval<1>(sizes, sizes + dom[0].length() - 1);
	  sizes += dom[0].length();
	  allempty = false;
	}

      // Give this new domain to the Node ... it will have the same size,
      // but a possibly different initial offset.

      this->all_m[i]->setDomain(dom);
      this->all_m[i]->setAllocated(dom);
    }

  // Update the total domains ... we make them the same since we assert
  // that this works only if there are no guard cells.

  if (allempty)
    this->domain_m = Domain_t();
  else
    this->domain_m = Interval<1>(this->firsti_m[0], sizes - 1);

  this->innerdomain_m = this->domain_m;

  // And finally clear our dirty flag.

  dirtyLayout_m = false;
}


//-----------------------------------------------------------------------------
//
// template<int Dim>
// void GridLayoutData<Dim>::calcMaps()
//
// Calculates the DomainMaps's for this object, based on the current
// settings for the blocks, since this can change due to dynamic operations.
//
//-----------------------------------------------------------------------------

template<int Dim>
void GridLayoutData<Dim>::calcMaps() const
{
  int i, j;

  // Don't need to do this if we're not initialized or not dirty.

  if (!this->initialized() || !dirty())
    return;

  // Initialize the maps in each dimension.

  for (i=0; i < Dim; ++i)
    {
      // Clear out any existing info

      map_m[i].zap();

      // If this dimension is empty, we can just quit here.

      if (this->domain_m[i].empty())
	continue;

      // Create a Loc used to get a serialized block index

      Loc<Dim> blockLoc(0);

      // For the owned domain.  The external guard
      // layers are considered part of the owned domain.

      map_m[i].initialize(Interval<1>(this->domain_m[i].first() -
				      this->externalGuards_m.lower(i),
				      this->domain_m[i].last() +
				      this->externalGuards_m.upper(i)));

      int b = this->blocks_m[i].first();
      for (j=0; j < b; ++j)
	{
	  // Get the serialized index of the block with grid position
	  // {0,...j,...0}, that is, all zero except for dimension d.

	  blockLoc[i] = Loc<1>(j);
	  int k = blockIndex(blockLoc);
	  Interval<1> blockDom = this->all_m[k]->domain()[i];

	  // Only add this block in if it is not empty

	  if (!blockDom.empty())
	    {
	      // Check for special case of edges, in which case the external
	      // guard layer is considered part of the owned domain.

	      int lo = (j == 0       ? this->externalGuards_m.lower(i) : 0);
	      int hi = (j == (b - 1) ? this->externalGuards_m.upper(i) : 0);

	      Interval<1> mval(blockDom.first() - lo, blockDom.last() + hi);
	      typename DomainMap<Interval<1>,int>::Value_t val(mval, j);

	      map_m[i].insert(val);
	    }
	}

      // update the DomainMap's

      map_m[i].update();
    }
}


//-----------------------------------------------------------------------------
//
// template<int Dim>
// void GridLayout<Dim>::calcAllocMaps()
//
// Calculates the DomainMaps's for this object, based on the current
// settings for the blocks, since this can change due to dynamic operations.
// This version sets up the allocated domain maps, which are possibly
// different if we have guard cells.
//
//-----------------------------------------------------------------------------

template<int Dim>
void GridLayoutData<Dim>::calcAllocMaps() const
{
  int i, j;

  // Don't need to do this if we're not initialized or not dirty.

  if (!this->initialized() || !dirty())
    return;

  // Initialize the maps in each dimension.

  for (i=0; i < Dim; ++i)
    {
      // Clear out any existing info

      mapAloc_m[i].zap();

      // If this dimension is empty, we can just quit here.

      if (this->domain_m[i].empty())
	continue;

      // Create a Loc used to get a serialized block index

      Loc<Dim> blockLoc(0);

      // For the allocated domain.

      mapAloc_m[i].initialize(Interval<1>(this->domain_m[i].first() -
					  this->externalGuards_m.lower(i),
					  this->domain_m[i].last() +
					  this->externalGuards_m.upper(i)));

      int b = this->blocks_m[i].first();
      for (j=0; j < b; ++j)
	{
	  // Get the serialized index of the block with grid position
	  // {0,...j,...0}, that is, all zero except for dimension d.

	  blockLoc[i] = Loc<1>(j);
	  int k = blockIndex(blockLoc);
	  Interval<1> blockDom = this->all_m[k]->domain()[i];

	  // Only add this block in if it is not empty

	  if (!blockDom.empty())
	    {
	      // Check for the special case for edge patches

	      int lo = (j==0 ?
			this->externalGuards_m.lower(i) :
			this->internalGuards_m.lower(i));
	      int hi = (j == (b - 1) ?
			this->externalGuards_m.upper(i) : 
			this->internalGuards_m.upper(i));

	      Interval<1> ival(blockDom.first() - lo, blockDom.last() + hi);
	      typename DomainMap<Interval<1>,int>::Value_t valAl(ival, j);

	      mapAloc_m[i].insert(valAl);
	    }
	}

      // update the DomainMap's

      mapAloc_m[i].update();
    }
}

//-----------------------------------------------------------------------------
//
// globalID takes a position within the domain of the layout, and returns
// the global ID for that node.
//
//-----------------------------------------------------------------------------

template<int Dim>
int GridLayoutData<Dim>::globalID(const Loc<Dim> &loc) const
{
  // Make sure the point is in our domain.

  PAssert(contains(this->domain_m, loc));

  // Find the i,j,k position of the point.

  Loc<Dim> point;
  for (int i=0; i < Dim; ++i)
    {
      // Find the position for this dimensio

      DomainMapTouchIterator<Interval<1>,int> dmti = 
	(map_m[i].touch(Interval<1>(loc[i]))).first;
      DomainMapTouchIterator<Interval<1>,int> baditerator;
      PAssert(dmti != baditerator);

      // Save the index for this dimension

     point[i] = *dmti;
    }

  // Serialize this point to get the global ID

  return blockIndex(point);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 1);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  return globalID(loc);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0, int i1) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 2);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  loc[1] = i1;
  return globalID(loc);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0, int i1, int i2) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 3);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  loc[1] = i1;
  loc[2] = i2;
  return globalID(loc);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 4);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  loc[1] = i1;
  loc[2] = i2;
  loc[3] = i3;
  return globalID(loc);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
				  int i4) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 5);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  loc[1] = i1;
  loc[2] = i2;
  loc[3] = i3;
  loc[4] = i4;
  return globalID(loc);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
				  int i4, int i5) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 6);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  loc[1] = i1;
  loc[2] = i2;
  loc[3] = i3;
  loc[4] = i4;
  loc[5] = i5;
  return globalID(loc);
}

template <int Dim>
int GridLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
				  int i4, int i5, int i6) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 7);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  loc[1] = i1;
  loc[2] = i2;
  loc[3] = i3;
  loc[4] = i4;
  loc[5] = i5;
  loc[6] = i6;
  return globalID(loc);
}


//============================================================
// touches operations
//============================================================

// Find all subdomains that touch on a given domain, and insert the
// intersection of these subdomains into the given output iterator.
// Return the number of touching elements. This version of touches
// can build either pointers or objects.

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int GridLayoutData<Dim>::touches(const OtherDomain &fulld, OutIter o,
				 const ConstructTag &ctag) const
{
  int i;

  // Make sure we have a valid layout

  PAssert(this->initialized());

  // Figure the type of the domain resulting from the intersection

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;

  // We only need to do touches for the overlapping domain.  If there is
  // nothing left, we can just return.

  OutDomain_t d = intersect(this->domain_m, fulld);
  if (d.empty())
    return 0;

  // Create an object of this output domain type for use later.

  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;
  
  int hiAxisIndex[Dim];
  int loAxisIndex[Dim];
  Loc<Dim> curnode;

  // Find the min and max block indices of the touching domains in each
  // dimension.

  for (i=0;i < Dim; ++i)
    {
      loAxisIndex[i] =
	*((map_m[i].touch(Interval<1>(d[i].first(),d[i].first()))).first);
      PAssert(loAxisIndex[i] >= 0 && 
	      loAxisIndex[i] <= this->blocks_m[i].first());

      hiAxisIndex[i] =
	*((map_m[i].touch(Interval<1>(d[i].last(),d[i].last() ))).first);
      PAssert(hiAxisIndex[i] >= 0 && 
	      hiAxisIndex[i] <= this->blocks_m[i].first());

      // if this is a reversed range, swap indexes. 

      if(loAxisIndex[i]>hiAxisIndex[i])
	{
	  int tmp = hiAxisIndex[i];
	  hiAxisIndex[i] = loAxisIndex[i];
	  loAxisIndex[i] = tmp;
	}

      curnode[i] = loAxisIndex[i];

    }

  // Go through all the blocks and output the values.

  int count=0;
  while(1)
    { 
      // curnode is the current block ... convert it to a serial index

      int nodeListIndex = blockIndex(curnode);

      // we can skip this block if it is empty

      if (!this->all_m[nodeListIndex]->domain().empty())
	{
	  // Make sure that block is OK ... this is a sanity check.

	  outDomain = intersect(fulld,this->all_m[nodeListIndex]->domain());
	  PAssert(!outDomain.empty());

	  // Output this block.

	  *o = touchesConstruct(outDomain,
				this->all_m[nodeListIndex]->allocated(),
				this->all_m[nodeListIndex]->affinity(),
				this->all_m[nodeListIndex]->context(),
				this->all_m[nodeListIndex]->globalID(),
				this->all_m[nodeListIndex]->localID(),
				ctag);
	  ++count;
	}

      // Move on to the next block

      curnode[0] += 1;
      for (i=0; i < Dim; ++i)
	{
	  if (curnode[i] == (hiAxisIndex[i]+1))
	    {
	      if (i==(Dim-1))
		break;

	      curnode[i] = loAxisIndex[i];
	      curnode[i+1] += 1;
	    }
	}

      // Are we at the end?

      if (curnode[Dim-1] == (hiAxisIndex[Dim-1]+1))
	break;
    }

  return count;
}


//-----------------------------------------------------------------------------
//
// Find all subdomains that touch on a given domain, and insert the
// intersection of these subdomains into the given output iterator.
// Return the number of touching elements. This version of touches
// can build either pointers or objects.
//
//-----------------------------------------------------------------------------

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int GridLayoutData<Dim>::touchesAlloc(const OtherDomain &fulld, OutIter o,
				      const ConstructTag &ctag) const 
{
  int i;

  // Make sure we have a valid layout

  PAssert(this->initialized());

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;

  // We only need to do touches for the overlapping domain.  If there is
  // nothing left, we can just return.

  OutDomain_t d = intersect(this->domain_m, fulld);
  if (d.empty())
    return 0;

  // Create an object of this output domain type for use later.

  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;
  
  int hiAxisIndex[Dim];
  int loAxisIndex[Dim];
  Loc<Dim> curnode;

  // Find the min and max block indices of the touching domains in each
  // dimension.

  for (i=0; i < Dim; ++i)
    {
      loAxisIndex[i] =
	*((mapAloc_m[i].touch(Interval<1>(d[i].first(),d[i].first()))).first);
      PAssert(loAxisIndex[i] >= 0 && 
	      loAxisIndex[i] < this->blocks_m[i].first());

      hiAxisIndex[i] = 
	*((mapAloc_m[i].touch(Interval<1>(d[i].last(),d[i].last()))).first);
      PAssert(hiAxisIndex[i] >= 0 && 
	      hiAxisIndex[i] < this->blocks_m[i].first());

      // if this is a reversed range, swap indexes. 

      if(loAxisIndex[i]>hiAxisIndex[i])
	{
	  int tmp = hiAxisIndex[i];
	  hiAxisIndex[i] = loAxisIndex[i];
	  loAxisIndex[i] = tmp;
	}

      curnode[i] = loAxisIndex[i];
    }

  // Go through all the blocks and output the values.

  int count=0;
  while(1)
    {
      // curnode is the current block ... convert it to a serial index

      int nodeListIndex = blockIndex(curnode);

      // we can skip this block if it is empty

      if (!this->all_m[nodeListIndex]->domain().empty())
	{
	  // Make sure that block is OK ... this is a sanity check.

	  outDomain = intersect(fulld,this->all_m[nodeListIndex]->allocated());

	  PAssert(!outDomain.empty());

	  // Output this block.

	  *o = touchesConstruct(outDomain,
				this->all_m[nodeListIndex]->allocated(),
				this->all_m[nodeListIndex]->affinity(),
				this->all_m[nodeListIndex]->context(),
				this->all_m[nodeListIndex]->globalID(),
				this->all_m[nodeListIndex]->localID(),
				ctag);
	  ++count;
	}
	
      // Move on to the next block

      curnode[0] += 1;
      for (i=0; i < Dim; ++i)
	{
	  if (curnode[i] == (hiAxisIndex[i]+1))
	    {
	      if (i == (Dim-1))
		break;

	      curnode[i] = loAxisIndex[i];
	      curnode[i+1] += 1;
	    }
	}

      // Are we at the end?

      if (curnode[Dim-1] == (hiAxisIndex[Dim-1]+1))
	break;
    }

  return count;
}

template<int Dim>
template<class Out>
void GridLayoutData<Dim>::print(Out & ostr)
{
  int i;
  ostr << " hasInternalGuards_m, hasExternalGuards_m "
       << this->hasInternalGuards_m << ' ' << this->hasExternalGuards_m
       << "\n internalGuards_m ";
  for (i=0; i<Dim; ++i)
    ostr << this->internalGuards_m.upper(i) << '-'
	 << this->internalGuards_m.lower(i) << ' ';
  ostr << "\n externalGuards_m ";
  for (i=0; i<Dim; ++i) 
    ostr << this->externalGuards_m.upper(i) << '-'
	 << this->externalGuards_m.lower(i) << ' ';
  ostr << '\n';
  FillIterator_t gstart = this->gcFillList_m.begin();
  FillIterator_t gend = this->gcFillList_m.end();
  ostr << " this->gcFillList_m\n";
  for(; gstart!=gend; ++gstart)
    ostr << "       "
	 << gstart->domain_m << ' '
	 << gstart->ownedID_m << ' '
	 << gstart->guardID_m << '\n';
  ostr << std::flush;
}



//============================================================
// GridLayout non-inline method definitions
//============================================================

//-----------------------------------------------------------------------------
//
// Constructors and initialize methods
// 
// See comments in the class definition (in GridLayout.h)
//
//-----------------------------------------------------------------------------

template <int Dim>
GridLayout<Dim>::GridLayout()
  : LayoutBase<Dim,GridLayoutData<Dim> >(new LayoutData_t()),
    Observable<This_t>(*this)
{ 
  this->pdata_m->attach(*this);
}
  
template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom,const DistributedTag & )
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t(gdom,
		  GridPartition<Dim>(Loc<Dim>(1)),
		  DistributedMapper<Dim>(GridPartition<Dim>(Loc<Dim>(1))))),
    Observable<This_t>(*this)     
{
  this->pdata_m->attach(*this); 
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom,const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t(gdom,
		  GridPartition<Dim>(Loc<Dim>(1)),
		  LocalMapper<Dim>())),
    Observable<This_t>(*this)     
{
  this->pdata_m->attach(*this); 
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const GuardLayers_t &gcs,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t(gdom,
		  GridPartition<Dim>(Loc<Dim>(1),gcs),  // single block
		  DistributedMapper<Dim>(
				GridPartition<Dim>(Loc<Dim>(1),gcs)))),
    Observable<This_t>(*this)     
{
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const GuardLayers_t &gcs,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t(gdom,
		  GridPartition<Dim>(Loc<Dim>(1),gcs), // single block
		  LocalMapper<Dim>())),
    Observable<This_t>(*this)   
{
  this->pdata_m->attach(*this);
}


template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Loc<Dim> &blocks,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >((LayoutData_t*) 0),
    Observable<This_t>(*this)
{
  if (!gdom.empty()) 
    this->pdata_m =
      new LayoutData_t( gdom,
			GridPartition<Dim>(makeRGrid(gdom,blocks)),
       	DistributedMapper<Dim>(GridPartition<Dim>(makeRGrid(gdom,blocks)))
			);
  else
    this->pdata_m = new LayoutData_t( gdom, 
				GridPartition<Dim>(blocks),
	DistributedMapper<Dim>(GridPartition<Dim>(blocks))
				);
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Loc<Dim> &blocks,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >((LayoutData_t*) 0),
    Observable<This_t>(*this)
{
  if (!gdom.empty()) 
    this->pdata_m = new LayoutData_t( gdom,
				GridPartition<Dim>(makeRGrid(gdom,blocks)),
				LocalMapper<Dim>());
  else
    this->pdata_m = new LayoutData_t( gdom,
				GridPartition<Dim>(blocks),
				LocalMapper<Dim>());
  this->pdata_m->attach(*this);
}


template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Loc<Dim> &blocks, 
			    const GuardLayers_t &gcs,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >((LayoutData_t*) 0),
    Observable<This_t>(*this)
{ 
  if (!gdom.empty()) 
    this->pdata_m = new LayoutData_t( gdom,
				GridPartition<Dim>(makeRGrid(gdom,blocks),gcs),
       DistributedMapper<Dim>(GridPartition<Dim>(makeRGrid(gdom,blocks),gcs)));
  else
    this->pdata_m = new LayoutData_t( gdom,
				GridPartition<Dim>(blocks,gcs),
       DistributedMapper<Dim>(GridPartition<Dim>(blocks,gcs)));
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Loc<Dim> &blocks, 
			    const GuardLayers_t &gcs,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >((LayoutData_t*) 0),
    Observable<This_t>(*this)
{ 
  if (!gdom.empty()) 
    this->pdata_m = new LayoutData_t(gdom,
			       GridPartition<Dim>(makeRGrid(gdom,blocks),gcs),
			       LocalMapper<Dim>());
  else
    this->pdata_m = new LayoutData_t(gdom,
			       GridPartition<Dim>(blocks,gcs),
			       LocalMapper<Dim>());
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Loc<Dim> &blocks, 
			    const GuardLayers_t &igcs, 
			    const GuardLayers_t &egcs,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >((LayoutData_t*) 0),
    Observable<This_t>(*this)
{
  if (!gdom.empty()) 
    this->pdata_m = new LayoutData_t( gdom,
      GridPartition<Dim>(makeRGrid(gdom,blocks),igcs,egcs),
      DistributedMapper<Dim>(
       GridPartition<Dim>(makeRGrid(gdom,blocks),igcs,egcs)));
  else
    this->pdata_m = new LayoutData_t( gdom,
      GridPartition<Dim>(blocks,igcs,egcs),
      DistributedMapper<Dim>(GridPartition<Dim>(blocks,igcs,egcs)));

  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Loc<Dim> &blocks, 
			    const GuardLayers_t &igcs, 
			    const GuardLayers_t &egcs,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >((LayoutData_t*) 0),
    Observable<This_t>(*this)
{
  if (!gdom.empty()) 
    this->pdata_m = new LayoutData_t(
		   gdom,
		   GridPartition<Dim>(makeRGrid(gdom,blocks),igcs,egcs),
		   LocalMapper<Dim>() );
  else
    this->pdata_m = new LayoutData_t( gdom,
				GridPartition<Dim>(blocks,igcs,egcs),
				LocalMapper<Dim>() );

  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Grid<Dim> &grid,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >( 
      new LayoutData_t( grid, 
			GridPartition<Dim>(grid),
			DistributedMapper<Dim>(GridPartition<Dim>(grid)))),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Grid<Dim> &grid,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >( 
      new LayoutData_t( grid, 
			GridPartition<Dim>(grid),
			LocalMapper<Dim>())),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Grid<Dim> &grid, 
			    const GuardLayers_t &gcs,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
( new LayoutData_t( grid, 
		    GridPartition<Dim>(grid,gcs),
		    DistributedMapper<Dim>(GridPartition<Dim>(grid,gcs) )) ),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const Grid<Dim> &grid, 
			    const GuardLayers_t &gcs,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
( new LayoutData_t( grid, 
		    GridPartition<Dim>(grid,gcs),
		    LocalMapper<Dim>()) ),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}


template <int Dim>
GridLayout<Dim>::GridLayout(const Grid<Dim> &grid, 
			    const GuardLayers_t &igcs, 
			    const GuardLayers_t &egcs,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t(grid, 
		  GridPartition<Dim>(grid,igcs,egcs),
		  DistributedMapper<Dim>(GridPartition<Dim>(grid,igcs,egcs)))),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}


template <int Dim>
GridLayout<Dim>::GridLayout(const Grid<Dim> &grid, 
			    const GuardLayers_t &igcs, 
			    const GuardLayers_t &egcs,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t( grid, 
		   GridPartition<Dim>(grid,igcs,egcs),
		   LocalMapper<Dim>() ) ),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
template <class Partitioner>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Partitioner &gpar,
			    const DistributedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
(new LayoutData_t(gdom, 
		  gpar,
		  DistributedMapper<Dim>(gpar)) ),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
template <class Partitioner>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Partitioner &gpar,
			    const ReplicatedTag &)
  : LayoutBase<Dim,GridLayoutData<Dim> >
      (new LayoutData_t(gdom, 
       		  	gpar,
			LocalMapper<Dim>()) ),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
template <class Partitioner>
GridLayout<Dim>::GridLayout(const Domain_t &gdom, 
			    const Partitioner &gpar,
			    const ContextMapper<Dim> &cmap)
  : LayoutBase<Dim,GridLayoutData<Dim> >(new LayoutData_t(gdom, gpar, cmap) ),
    Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template <int Dim>
GridLayout<Dim>::GridLayout(const This_t &model) 
  : LayoutBase<Dim,GridLayoutData<Dim> >(model.pdata_m),
    Observable<This_t>(*this)
{ 
   this->pdata_m->attach(*this);
}
  

//-----------------------------------------------------------------------------
//
// assignment operator for GridLayout
//
//-----------------------------------------------------------------------------

template <int Dim>
GridLayout<Dim> &GridLayout<Dim>::operator=(const This_t &model)
{
  if (this != &model)
    {
      this->pdata_m->detach(*this);
      this->pdata_m = model.pdata_m;
      this->pdata_m->attach(*this);
    }

  return *this;
}


//-----------------------------------------------------------------------------
//
// Initialize methods for GridLayout
//
//-----------------------------------------------------------------------------

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const DistributedTag &)
{
  this->pdata_m->initialize( gdom, 
		       GridPartition<Dim>(),
		       DistributedMapper<Dim>(GridPartition<Dim>()) );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const GuardLayers_t &gcs,
				 const DistributedTag &)
{
  this->pdata_m->initialize( gdom, 
		       GridPartition<Dim>(gcs),
		       DistributedMapper<Dim>(GridPartition<Dim>(gcs))  );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Loc<Dim> &blocks,
				 const DistributedTag &)
{
  if (!gdom.empty())
    this->pdata_m->initialize(
	   gdom, 
	   GridPartition<Dim>(makeRGrid(gdom,blocks)),
	   DistributedMapper<Dim>(GridPartition<Dim>(makeRGrid(gdom,blocks))));
  else 
    this->pdata_m->initialize( gdom, 
			 GridPartition<Dim>(blocks),
			 DistributedMapper<Dim>(GridPartition<Dim>(blocks)) );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Loc<Dim> &blocks,
				 const GuardLayers_t &gcs,
				 const DistributedTag &)
{
  if (!gdom.empty())
    this->pdata_m->initialize( 
	      gdom,
	      GridPartition<Dim>(makeRGrid(gdom,blocks),gcs),
	      DistributedMapper<Dim>(
                GridPartition<Dim>(makeRGrid(gdom,blocks),gcs)));
  else
    this->pdata_m->initialize( 
			gdom, 
			GridPartition<Dim>(blocks,gcs),
			DistributedMapper<Dim>(
			  GridPartition<Dim>(blocks,gcs)) );
  
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Loc<Dim> &blocks,
				 const GuardLayers_t &igcs,
				 const GuardLayers_t &egcs,
				 const DistributedTag &)
{
  if (!gdom.empty())
    this->pdata_m->initialize( 
      gdom,
      GridPartition<Dim>(makeRGrid(gdom,blocks),igcs,egcs),
      DistributedMapper<Dim>(
	GridPartition<Dim>(makeRGrid(gdom,blocks),igcs,egcs)));
  else 
    this->pdata_m->initialize( 
			gdom, 
			GridPartition<Dim>(blocks,igcs,egcs),
			DistributedMapper<Dim>(
			  GridPartition<Dim>(blocks,igcs,egcs)));
}

template <int Dim>
void GridLayout<Dim>::initialize(const Grid<Dim> &grid,
				 const DistributedTag &)
{
  this->pdata_m->initialize( grid, 
		       GridPartition<Dim>(grid),
		       DistributedMapper<Dim>(GridPartition<Dim>(grid)));
}

template <int Dim>
void GridLayout<Dim>::initialize(const Grid<Dim> &grid,
				 const GuardLayers_t &gcs,
				 const DistributedTag &)
{
  this->pdata_m->initialize( grid, 
		       GridPartition<Dim>(grid,gcs),
		       DistributedMapper<Dim>(GridPartition<Dim>(grid,gcs)) );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Grid<Dim> &grid,
				 const GuardLayers_t &igcs,
				 const GuardLayers_t &egcs,
				 const DistributedTag &)
{
  this->pdata_m->initialize( 
		      grid, 
		      GridPartition<Dim>(grid,igcs,egcs),
		      DistributedMapper<Dim>(
			GridPartition<Dim>(grid,igcs,egcs)) );
}

template <int Dim>
template <class Partitioner>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Partitioner &gpar,
				 const DistributedTag &)
{
  this->pdata_m->initialize(gdom, 
		      gpar,
		      DistributedMapper<Dim>(gpar));
}

// ReplicatedTag

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const ReplicatedTag &)
{
  this->pdata_m->initialize( gdom, GridPartition<Dim>(),LocalMapper<Dim>() );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const GuardLayers_t &gcs,
				 const ReplicatedTag &)
{
  this->pdata_m->initialize( gdom, GridPartition<Dim>(gcs),LocalMapper<Dim>() );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Loc<Dim> &blocks,
				 const ReplicatedTag &)
{
  if (!gdom.empty())
    this->pdata_m->initialize( 
			gdom, 
			GridPartition<Dim>(makeRGrid(gdom,blocks)),
			LocalMapper<Dim>() );
  else 
    this->pdata_m->initialize( 
			gdom, 
			GridPartition<Dim>(blocks),
			LocalMapper<Dim>() );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Loc<Dim> &blocks,
				 const GuardLayers_t &gcs,
				 const ReplicatedTag &)
{
  if (!gdom.empty())
    this->pdata_m->initialize( gdom,
			 GridPartition<Dim>(makeRGrid(gdom,blocks),gcs),
			 LocalMapper<Dim>() );
  else
    this->pdata_m->initialize( gdom, 
			 GridPartition<Dim>(blocks,gcs),
			 LocalMapper<Dim>() );
  
}

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Loc<Dim> &blocks,
				 const GuardLayers_t &igcs,
				 const GuardLayers_t &egcs,
				 const ReplicatedTag &)
{
  if (!gdom.empty())
    this->pdata_m->initialize( 
	      gdom,
	      GridPartition<Dim>(makeRGrid(gdom,blocks),igcs,egcs),
	      LocalMapper<Dim>() );
  else 
    this->pdata_m->initialize( 
			gdom, 
			GridPartition<Dim>(blocks,igcs,egcs),
			LocalMapper<Dim>() );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Grid<Dim> &grid,
				 const ReplicatedTag &)
{
  this->pdata_m->initialize( grid, 
		       GridPartition<Dim>(grid),
		       LocalMapper<Dim>() );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Grid<Dim> &grid,
				 const GuardLayers_t &gcs,
				 const ReplicatedTag &)
{
  this->pdata_m->initialize( 
		      grid, 
		      GridPartition<Dim>(grid,gcs),
		      LocalMapper<Dim>() );
}

template <int Dim>
void GridLayout<Dim>::initialize(const Grid<Dim> &grid,
				 const GuardLayers_t &igcs,
				 const GuardLayers_t &egcs,
				 const ReplicatedTag &)
{
  this->pdata_m->initialize( 
		      grid, 
		      GridPartition<Dim>(grid,igcs,egcs),
		      LocalMapper<Dim>() );
}


template <int Dim>
template <class Partitioner>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Partitioner &gpar,
				 const ReplicatedTag &)
{
  this->pdata_m->initialize(gdom, 
		      gpar,
		      LocalMapper<Dim>());
}

// This initializer is intented to be used by the I/O system

template <int Dim>
void GridLayout<Dim>::initialize(const Domain_t& idom,
				 const List_t& nodes,
				 const Loc<Dim>& blocks,
				 bool hasIG, bool hasEG,
				 const GuardLayers_t& ig,
				 const GuardLayers_t& eg)
{
  this->pdata_m->initialize(idom,nodes,blocks,hasIG,hasEG,ig,eg);
}



template <int Dim>
template <class Partitioner>
void GridLayout<Dim>::initialize(const Domain_t &gdom,
				 const Partitioner &gpar,
				 const ContextMapper<Dim> &cmap)
{
  this->pdata_m->initialize(gdom, gpar, cmap);
}

template <int Dim>
template <class Ostream>
void GridLayout<Dim>::print(Ostream &ostr) const
{
  ostr << "GridLayout " << this->ID() << " on global domain " 
       << this->domain() << ":" << '\n';
  ostr << "   Total subdomains: " << this->sizeGlobal() << '\n';
  ostr << "   Local subdomains: " << this->sizeLocal() << '\n';
  ostr << "  Remote subdomains: " << this->sizeRemote() << '\n';
  ostr << "        Grid blocks: " << this->blocks() << '\n';
  typename GridLayout<Dim>::const_iterator a;
  for (a = this->beginGlobal(); a != this->endGlobal(); ++a)
    ostr << "  Global subdomain = " << *a << '\n';
  for (a = this->beginLocal(); a != this->endLocal(); ++a)
    ostr << "   Local subdomain = " << *a << '\n';
  for (a = this->beginRemote(); a != this->endRemote(); ++a)
    ostr << "  Remote subdomain = " << *a << '\n';
  this->pdata_m->print(ostr);
}

template<int Dim, int Dim2>
template <class Ostream>
void GridLayoutView<Dim, Dim2>::print(Ostream &ostr) const 
{
  ostr << "GridLayoutView " << this->ID() << " on global domain " 
       << this->domain() << ':' << '\n';
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


// } // namespace POOMA

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: GridLayout.cpp,v $   $Author: richi $
// $Revision: 1.93 $   $Date: 2004/11/10 22:02:10 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
