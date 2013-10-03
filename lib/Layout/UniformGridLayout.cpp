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
//   UniformGridLayoutData<Dim> template definitions.
//   UniformGridLayout<Dim> template definitions.
//-----------------------------------------------------------------------------

#include "Threads/PoomaSmarts.h"
#include "Layout/UniformGridLayout.h"
#include "Utilities/PAssert.h"

#include <vector>


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//============================================================
// UniformGridLayoutData non-inline method definitions
//============================================================

//-----------------------------------------------------------------------------
//
// template <int Dim>
// template <class Partitioner,class Mapper>
// void UniformGridLayout<Dim>::
// UniformGridLayoutData(const Domain_t &gdom,
//                       const Partitioner &gpar,
//                       const Mapper &cmap)
// Originally, we provided a slew of constructors, mirroring those
// in UniformGridLayout. However, once we eliminated storage of a
// partitioner object, it became simpler to have the Layout construct
// the partitioner and then initialize its data object by passing 
// that partitioner to this constructor.
//
//-----------------------------------------------------------------------------

template <int Dim>
inline UniformGridLayoutData<Dim>::
UniformGridLayoutData() 
  : Observable<UniformGridLayoutData>(*this) 
{ 
  for (int i = 0; i < Dim; ++i)
    blockstride_m[i] = blocksizes_m[i] = 0;
}

template <int Dim>
template <class Partitioner>
UniformGridLayoutData<Dim>::
UniformGridLayoutData(const Domain_t &gdom, 
		      const Partitioner &gpar,
		      const ContextMapper<Dim> & cmap )
  : LayoutBaseData<Dim>(false,
			false,
			GuardLayers_t(0),
			GuardLayers_t(0),
			gdom,
			gdom),
    Observable<UniformGridLayoutData>(*this)
{
  // Figure out if we have guards to worry about.
    
  if (gpar.hasInternalGuards() && gpar.maxSize() > 1)
    {
      this->hasInternalGuards_m = true;
      this->internalGuards_m = gpar.internalGuards();
    }
      
  if (gpar.hasExternalGuards())
    {
      this->hasExternalGuards_m = true;
      this->externalGuards_m = gpar.externalGuards();
      GuardLayers<Dim>::addGuardLayers(this->domain_m,this->externalGuards_m);
    }
    
  // Do the partitioning. 
  // This initializes allDomain_m, firsti_m, etc.
      
  partition(gpar,cmap);

}

template <int Dim>
template <class Partitioner>
void UniformGridLayoutData<Dim>::partition(const Partitioner &gpar,
					   const ContextMapper<Dim> &cmap)
{
  int i;

  // In spite of being templated, this only works with uniform-grid
  // partitioners.
    
  CTAssert(Partitioner::uniform);

  // We must have something to partition, and the domain lists must be
  // empty.
  
  PAssert(this->domain_m.size() > 0);
  PAssert(this->innerdomain_m.size() > 0);
  PAssert(this->all_m.size() == 0);
  PAssert(this->local_m.size() == 0);
  PAssert(this->remote_m.size() == 0);

  // Save the first and block size info from the current domain.

  this->blocks_m = gpar.blocks();

  // Note, for the purposes of partitioning, we pretend like we're
  // only working with the inner domain. The total domain includes the
  // external guards, and those do not affect the partitioning.  

  blockstride_m[0] = 1;
  int blocks[Dim];
  for (i = 0; i < Dim; ++i)
  {
    this->firsti_m[i] = this->innerdomain_m[i].first();
    this->firste_m[i] = this->domain_m[i].first();
    blocks[i] = gpar.blocks()[i].first();
    allDomain_m[i] = Interval<1>(blocks[i]);
    blocksizes_m[i] = this->innerdomain_m[i].length() / blocks[i];
    if (i > 0)
      blockstride_m[i] = blockstride_m[i-1] * blocks[i-1];
  }

  // Invoke the partitioner.
  
  gpar.partition(this->innerdomain_m, this->all_m, cmap);

  // fill local and remote lists

  typename List_t::const_iterator start = this->all_m.begin();
  typename List_t::const_iterator end   = this->all_m.end();
  
  for ( ; start!=end ; ++start)
    {
      if ( (*start)->context() == Pooma::context()
	   || (*start)->context() == -1 )
	{ 
	  (*start)->localID() = this->local_m.size();
	  this->local_m.push_back(*start);
	}
      else
	this->remote_m.push_back(*start);
    }

  if (this->hasInternalGuards_m) 
    {
      this->gcFillList_m.clear();
      calcGCFillList();
    }
}
//-----------------------------------------------------------------------------
//
// template<int Dim>
// void UniformGridLayoutData<Dim>::initialize
//
// Used by an I/O or data management entity to initialize the layout based
// on detailed state information previously stored. As in the case of the
// initializer with the partitioner argument, this method will call 'addDomain'
// to add in the new domains it creates and will initialize
// guard cell info, etc.
//
//-----------------------------------------------------------------------------

template<int Dim>
void UniformGridLayoutData<Dim>::initialize(const Domain_t& idom,
				 const List_t& nodes,
				 const Loc<Dim>& ublocks,
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

  // Save the first and block size info from the current domain.

  this->blocks_m = ublocks;

  // Note, for the purposes of partitioning, we pretend like we're
  // only working with the inner domain. The total domain includes the
  // external guards, and those do not affect the partitioning.
  
  blockstride_m[0] = 1;
  int blocks[Dim];
  for (i = 0; i < Dim; ++i)
  {
    this->firsti_m[i] = this->innerdomain_m[i].first();
    this->firste_m[i] = this->domain_m[i].first();
    blocks[i] = ublocks[i].first();
    allDomain_m[i] = Interval<1>(blocks[i]);
    blocksizes_m[i] = this->innerdomain_m[i].length() / blocks[i];
    if (i > 0)
      blockstride_m[i] = blockstride_m[i-1] * blocks[i-1];
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

  if (this->hasInternalGuards_m) 
    {
      this->gcFillList_m.clear();
      calcGCFillList();
    }
}
//-----------------------------------------------------------------------------
//
// template <int Dim>
// void UniformGridLayoutData<Dim>::calcGCFillList()
//
// Calculates the cached information needed by MultiPatch Engine to
// fill the guard cells.
//
//-----------------------------------------------------------------------------

template <int Dim>
void UniformGridLayoutData<Dim>::calcGCFillList()
  {
    int d, p;

    // We want to create the list in such a manner that all
    // communication in a particular direction is done first, allowing
    // parallelism with the least amount of contention for
    // patches. Thus we have an outer loop over Dim, doing the upward
    // copies first, then the downward copies.
     
    int numPatches = this->all_m.size();
     
    this->gcFillList_m.reserve(2*Dim*this->local_m.size());
    
    for (d = 0; d < Dim; ++d)
      {
        // First we "send" up in every direction, meaning that we fill
        // the "lower" internal guard cells for domains that have
        // them.
        
        if (this->internalGuards_m.lower(d) > 0)
          {
            // We use a DomainIterator to figure out if we're at edges
            // as we iterate through the patches.
            
            // NOTE!!! Implicit in this is that all of the domains are
            // stored in fortran storage order in the all_m array.  Of
            // course, this is also only valid for single context
            // stuff.  When we go to multiple contexts, this algorithm
            // will still work if the local's always form a block that
            // is also stored in fortran storage order.
            
            typename Interval<Dim>::iterator pos = allDomain_m.begin();
            
            for (p = 0; p < numPatches; ++p, ++pos)
              {
                // Edge detection. If this element is at the upper
                // edge in the direction that we're sending, skip it
                // and continue.
               
                if ( (*pos)[d].first() == allDomain_m[d].last() ) continue;
                  
                // The destination ID is one step "up" in the "d"
                // direction, which is at an offset in all_m of
                // blockstride_m[d]:
                
                int sourceID = p;
                int destID   = p + blockstride_m[d];
                
                // Check that our destination is in range.

                PAssert(destID < numPatches);

                // We should never get here if we're at the last cell.

                PAssert(pos != allDomain_m.end()); 
                                      
                // Calculate the domain of the overlapping cells that
                // need to be communicated. This is the total domain
                // in all directions but "d", where it is just the top
                // guard-cell width of the source domain.
                
                // (This causes copying of some uninitialized data,
                // since the first direction includes guards [which
                // haven't been filled] in the perpendicular directions,
                // but that data later gets overwritten by good data.
                // Could change this to use more conservative sets
                // of domains, but then the accumulation would have to
                // happen in reverse order [I think???].)
                
                Domain_t gcdom(this->all_m[p]->allocated());
                
                int max = this->all_m[p]->domain()[d].last();
                int min = max - this->internalGuards_m.lower(d) + 1;
                                
                gcdom[d] = Interval<1>(min,max);  
                 
                // Now, push IDs and source into cache...
 		if (
		    this->all_m[sourceID]->context() == -1 || 
		    this->all_m[sourceID]->context() == Pooma::context() || 
 		    this->all_m[destID]->context() == Pooma::context()
		    )
                this->gcFillList_m.push_back(GCFillInfo_t(gcdom,sourceID,destID,d*2));
              }
          }

        // Next we "send" down in every direction, meaning that we
        // fill the "upper" internal guard cells for domains that have
        // them.

        if (this->internalGuards_m.upper(d) > 0)
          {
            typename Interval<Dim>::iterator pos = allDomain_m.begin();

            for (p = 0; p < numPatches; ++p, ++pos)
              {
                // Edge detection. If this element is at the lower
                // edge in the direction that we're sending, skip it
                // and continue.
               
                if ( (*pos)[d].first() == allDomain_m[d].first() ) continue;
                  
                // The destination ID is one step "down" in the "d"
                // direction, which is at an offset in all_m of
                // blockstride_m[d]:
                
                int sourceID = p;
                int destID   = p - blockstride_m[d];
                 
                // Check that destination is in range.

                PAssert(destID >= 0);

                // Calculate the domain of the overlapping cells that
                // need to be communicated. See comments above.

                Domain_t gcdom(this->all_m[p]->allocated());
                
                int min = this->all_m[p]->domain()[d].first();
                int max = min + this->internalGuards_m.upper(d) - 1;
                
                gcdom[d] = Interval<1>(min,max);  
                 
                // Now, push IDs and source into cache...
 		if (
		    this->all_m[sourceID]->context() == -1 || 
		    this->all_m[sourceID]->context() == Pooma::context() || 
 		    this->all_m[destID]->context() == Pooma::context()
		    )
		  this->gcFillList_m.push_back(GCFillInfo_t(gcdom,sourceID,destID,d*2+1));
              }
          }
      }
  }



//-----------------------------------------------------------------------------
//
// template <int Dim>
// template <class Partitioner>
// void UniformGridLayout<Dim>::
// repartition(const Partitioner &)
//
// Repartition the layout using a new Partitioner scheme.  The initial
// domain lists are cleared out, the partitioner is invoked, and then
// all the observers are notified.  This can only be done with a
// GridPartition partitioner.
//
//-----------------------------------------------------------------------------

template <int Dim>
template <class Partitioner>
bool UniformGridLayoutData<Dim>::
repartition(const Partitioner &p,
	    const ContextMapper<Dim>& cmap)
{
  // We can only repartition if we have been initialized to some domain.

  PAssert(this->domain_m.size() > 0);

  // Delete existing nodes and clear all the lists.

  for (int i = 0; i < this->all_m.size(); ++i)
    delete this->all_m[i];
    
  this->all_m.clear();
  this->local_m.clear();
  this->remote_m.clear();

  // Do the new partitioning.

  partition(p,cmap);

  if (this->hasInternalGuards_m) 
    {
      this->gcFillList_m.clear();
      calcGCFillList();
    }

  // Notify all users.

  this->notify(repartitionEvent);

  return true;
}


// Find all subdomains that touch on a given domain, and insert the
// intersection of these subdomains into the given output iterator.
// Return the number of touching elements. This version of touches
// can build either pointers or objects.

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int UniformGridLayoutData<Dim>::touches(const OtherDomain &d, OutIter o,
					const ConstructTag &ctag) const 
{
  int i, count = 0;

  // Make sure we have a valid touching domain.

  PAssert(this->initialized());
  PAssert(contains(this->domain_m, d));

  // Find the starting and ending grid positions, and store as an
  // Interval.

  Interval<Dim> box = Pooma::NoInit();
  for (i = 0; i < Dim; ++i) 
    {
      int a, b;
      if (!this->hasExternalGuards_m)
	{
	  a = (d[i].min() - this->firsti_m[i]) / blocksizes_m[i];
	  b = (d[i].max() - this->firsti_m[i]) / blocksizes_m[i];
	}
      else
	{
	  // If we're in the lower guards, this will fall
	  // through to give zero.
            
	  a = b = 0;
	  int pos = d[i].min();
	  int last = this->innerdomain_m[i].last();
	  int del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      a = del / blocksizes_m[i];
	    else
	      a = allDomain_m[i].last();
            
	  pos = d[i].max();
	  del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      b = del / blocksizes_m[i];
	    else
	      b = allDomain_m[i].last();
	}
      box[i] = Interval<1>(a, b);
    }

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;
  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;

  // Iterate through the Interval grid positions.

  typename Interval<Dim>::const_iterator boxiter = box.begin();
    
  while (boxiter != box.end()) 
    {
      // Calculate the linear position of the current node.
        
      int indx = (*boxiter)[0].first();
      for (i = 1; i < Dim; ++i)
	indx += blockstride_m[i] * (*boxiter)[i].first();

      // Get that node, intersect the domain with the requested one,
      // and write the result out to the output iterator.
        
      PAssert(indx >= 0 && indx < this->all_m.size());
        
      outDomain = intersect(d, this->all_m[indx]->domain());

      PAssert(!outDomain.empty());
        
      *o = touchesConstruct(outDomain,
			    this->all_m[indx]->allocated(),
			    this->all_m[indx]->affinity(),
			    this->all_m[indx]->context(),
			    this->all_m[indx]->globalID(),
			    this->all_m[indx]->localID(),
			    ctag);

      // Increment output iterator, count, and grid iterator.

      ++o;
      ++count;
      ++boxiter;
    }

  // Return the number of non-empty domains we found.
      
  return count;
}

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int UniformGridLayoutData<Dim>::touchesLocal(const OtherDomain &d, 
					     OutIter o, 
					     const ConstructTag &ctag) const 
{
  int i, count = 0;

  // Make sure we have a valid touching domain.

  PAssert(this->initialized());
  PAssert(contains(this->domain_m, d));

  // Find the starting and ending grid positions, and store as an
  // Interval.

  Interval<Dim> box = Pooma::NoInit();
  for (i = 0; i < Dim; ++i) 
    {
      int a, b;
      if (!this->hasExternalGuards_m)
	{
	  a = (d[i].min() - this->firsti_m[i]) / blocksizes_m[i];
	  b = (d[i].max() - this->firsti_m[i]) / blocksizes_m[i];
	}
      else
	{
	  // If we're in the lower guards, this will fall
	  // through to give zero.
            
	  a = b = 0;
	  int pos = d[i].min();
	  int last = this->innerdomain_m[i].last();
	  int del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      a = del / blocksizes_m[i];
	    else
	      a = allDomain_m[i].last();
            
	  pos = d[i].max();
	  del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      b = del / blocksizes_m[i];
	    else
	      b = allDomain_m[i].last();
	}
      box[i] = Interval<1>(a, b);
    }

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;
  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;

  // Iterate through the Interval grid positions.

  typename Interval<Dim>::const_iterator boxiter = box.begin();
    
  while (boxiter != box.end()) 
    {
      // Calculate the linear position of the current node.
        
      int indx = (*boxiter)[0].first();
      for (i = 1; i < Dim; ++i)
	indx += blockstride_m[i] * (*boxiter)[i].first();

      // Get that node, intersect the domain with the requested one,
      // and write the result out to the output iterator.
        
      PAssert(indx >= 0 && indx < this->local_m.size());
        
      outDomain = intersect(d, this->local_m[indx]->domain());

      PAssert(!outDomain.empty());
        
      *o = touchesConstruct(outDomain,
			    this->local_m[indx]->allocated(),
			    this->local_m[indx]->affinity(),
			    this->local_m[indx]->context(),
			    this->local_m[indx]->globalID(),
			    this->local_m[indx]->localID(),
			    ctag);

      // Increment output iterator, count, and grid iterator.

      ++o;
      ++count;
      ++boxiter;
    }

  // Return the number of non-empty domains we found.
      
  return count;
}

// Find all Remote subdomains that touch on a given domain, and insert the
// intersection of these subdomains into the given output iterator.
// Return the number of touching elements. This version of touches
// can build either pointers or objects.

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int UniformGridLayoutData<Dim>::touchesRemote(const OtherDomain &d, 
					      OutIter o, 
					      const ConstructTag &ctag) const 
{
  int i, count = 0;

  // Make sure we have a valid touching domain.

  PAssert(this->initialized());
  PAssert(contains(this->domain_m, d));

  // Find the starting and ending grid positions, and store as an
  // Interval.

  Interval<Dim> box = Pooma::NoInit();
  for (i = 0; i < Dim; ++i) 
    {
      int a, b;
      if (!this->hasExternalGuards_m)
	{
	  a = (d[i].min() - this->firsti_m[i]) / blocksizes_m[i];
	  b = (d[i].max() - this->firsti_m[i]) / blocksizes_m[i];
	}
      else
	{
	  // If we're in the lower guards, this will fall
	  // through to give zero.
            
	  a = b = 0;
	  int pos = d[i].min();
	  int last = this->innerdomain_m[i].last();
	  int del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      a = del / blocksizes_m[i];
	    else
	      a = allDomain_m[i].last();
            
	  pos = d[i].max();
	  del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      b = del / blocksizes_m[i];
	    else
	      b = allDomain_m[i].last();
	}
      box[i] = Interval<1>(a, b);
    }

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;
  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;

  // Iterate through the Interval grid positions.

  typename Interval<Dim>::const_iterator boxiter = box.begin();
    
  while (boxiter != box.end()) 
    {
      // Calculate the linear position of the current node.
        
      int indx = (*boxiter)[0].first();
      for (i = 1; i < Dim; ++i)
	indx += blockstride_m[i] * (*boxiter)[i].first();

      // Get that node, intersect the domain with the requested one,
      // and write the result out to the output iterator.
        
      PAssert(indx >= 0 && indx < this->remote_m.size());
        
      outDomain = intersect(d, this->remote_m[indx]->domain());

      PAssert(!outDomain.empty());
        
      *o = touchesConstruct(outDomain,
			    this->remote_m[indx]->allocated(),
			    this->remote_m[indx]->affinity(),
			    this->remote_m[indx]->context(),
			    this->remote_m[indx]->globalID(),
			    this->remote_m[indx]->localID(),
			    ctag);

      // Increment output iterator, count, and grid iterator.

      ++o;
      ++count;
      ++boxiter;
    }

  // Return the number of non-empty domains we found.
      
  return count;
}

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int UniformGridLayoutData<Dim>::touchesAlloc(const OtherDomain &d, OutIter o, 
					     const ConstructTag &ctag) const 
{
  // If there are no internal guard cells, then this calculation is the
  // same as the normal touches calculation.
    
  if (!this->hasInternalGuards_m) return touches(d,o,ctag);
    
  int i, count = 0;

  // Make sure we have a valid touching domain.

  PAssert(this->initialized());
  PAssert(contains(this->domain_m, d));

  // Find the starting and ending grid positions, and store as an
  // Interval.

  Interval<Dim> box = Pooma::NoInit();
  for (i = 0; i < Dim; ++i) 
    {
      int a, b;
        
      // This part is unchanged from touches. What we're going to do
      // here is simply extend the range in each direction by one block
      // and then let the intersection calculation below sort out if
      // there is actually an intersection. Otherwise we'd have to do
      // comparisons on all blocks up here and then still do the
      // intersection on the remaining blocks below. On average, this
      // should be slightly faster I think.
        
      if (!this->hasExternalGuards_m)
	{
	  a = (d[i].min() - this->firsti_m[i]) / blocksizes_m[i];
	  b = (d[i].max() - this->firsti_m[i]) / blocksizes_m[i];
	}
      else
	{
	  // If we're in the lower guards, this will fall
	  // through to give zero.
            
	  a = b = 0;
	  int pos = d[i].min();
	  int last = this->innerdomain_m[i].last();
	  int del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      a = del / blocksizes_m[i];
	    else
	      a = allDomain_m[i].last();
            
	  pos = d[i].max();
	  del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      b = del / blocksizes_m[i];
	    else
	      b = allDomain_m[i].last();
	}
          
      // Now we check that we're not at the ends of the domain for the
      // brick of blocks, and extend the region accordingly.
        
      if (a > 0) --a;
      if (b < allDomain_m[i].last()) ++b;
        
      box[i] = Interval<1>(a, b);
    }

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;
  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;

  // Iterate through the Interval grid positions.

  typename Interval<Dim>::const_iterator boxiter = box.begin();
    
  while (boxiter != box.end()) 
    {
      // Calculate the linear position of the current node.
        
      int indx = (*boxiter)[0].first();
      for (i = 1; i < Dim; ++i)
	indx += blockstride_m[i] * (*boxiter)[i].first();

      // Get that node, intersect the *allocated* domain with the 
      // requested one, and write the result out to the output iterator.
        
      PAssert(indx >= 0 && indx < this->all_m.size());
        
      outDomain = intersect(d, this->all_m[indx]->allocated());

      // We can no longer assert that outDomain is not empty since
      // we extended the search box without checking. Thus we now 
      // have an if tests around the output.
        
      if (!outDomain.empty())
	{
	  *o = touchesConstruct(outDomain,
				this->all_m[indx]->allocated(),
				this->all_m[indx]->affinity(),
				this->all_m[indx]->context(),
				this->all_m[indx]->globalID(),
				this->all_m[indx]->localID(),
				ctag);
	}

      // Increment output iterator, count, and grid iterator.

      ++o;
      ++count;
      ++boxiter;
    }

  // Return the number of non-empty domains we found.
      
  return count;
}

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int UniformGridLayoutData<Dim>::touchesAllocLocal(const OtherDomain &d, OutIter o, 
					          const ConstructTag &ctag) const 
{
  // If there are no internal guard cells, then this calculation is the
  // same as the normal touches calculation.
    
  if (!this->hasInternalGuards_m) return touches(d,o,ctag);
    
  int i, count = 0;

  // Make sure we have a valid touching domain.

  PAssert(this->initialized());
  PAssert(contains(this->domain_m, d));

  // Find the starting and ending grid positions, and store as an
  // Interval.

  Interval<Dim> box = Pooma::NoInit();
  for (i = 0; i < Dim; ++i) 
    {
      int a, b;
        
      // This part is unchanged from touches. What we're going to do
      // here is simply extend the range in each direction by one block
      // and then let the intersection calculation below sort out if
      // there is actually an intersection. Otherwise we'd have to do
      // comparisons on all blocks up here and then still do the
      // intersection on the remaining blocks below. On average, this
      // should be slightly faster I think.
        
      if (!this->hasExternalGuards_m)
	{
	  a = (d[i].min() - this->firsti_m[i]) / blocksizes_m[i];
	  b = (d[i].max() - this->firsti_m[i]) / blocksizes_m[i];
	}
      else
	{
	  // If we're in the lower guards, this will fall
	  // through to give zero.
            
	  a = b = 0;
	  int pos = d[i].min();
	  int last = this->innerdomain_m[i].last();
	  int del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      a = del / blocksizes_m[i];
	    else
	      a = allDomain_m[i].last();
            
	  pos = d[i].max();
	  del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      b = del / blocksizes_m[i];
	    else
	      b = allDomain_m[i].last();
	}
          
      // Now we check that we're not at the ends of the domain for the
      // brick of blocks, and extend the region accordingly.
        
      if (a > 0) --a;
      if (b < allDomain_m[i].last()) ++b;
        
      box[i] = Interval<1>(a, b);
    }

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;
  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;

  // Iterate through the Interval grid positions.

  typename Interval<Dim>::const_iterator boxiter = box.begin();
    
  while (boxiter != box.end()) 
    {
      // Calculate the linear position of the current node.
        
      int indx = (*boxiter)[0].first();
      for (i = 1; i < Dim; ++i)
	indx += blockstride_m[i] * (*boxiter)[i].first();

      // Get that node, intersect the *allocated* domain with the 
      // requested one, and write the result out to the output iterator.
        
      PAssert(indx >= 0 && indx < this->local_m.size());
        
      outDomain = intersect(d, this->local_m[indx]->allocated());

      // We can no longer assert that outDomain is not empty since
      // we extended the search box without checking. Thus we now 
      // have an if tests around the output.
        
      if (!outDomain.empty())
	{
	  *o = touchesConstruct(outDomain,
				this->local_m[indx]->allocated(),
				this->local_m[indx]->affinity(),
				this->local_m[indx]->context(),
				this->local_m[indx]->globalID(),
				this->local_m[indx]->localID(),
				ctag);
	}

      // Increment output iterator, count, and grid iterator.

      ++o;
      ++count;
      ++boxiter;
    }

  // Return the number of non-empty domains we found.
      
  return count;
}

template <int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int UniformGridLayoutData<Dim>::touchesAllocRemote(const OtherDomain &d, OutIter o, 
					           const ConstructTag &ctag) const 
{
  // If there are no internal guard cells, then this calculation is the
  // same as the normal touches calculation.
    
  if (!this->hasInternalGuards_m) return touches(d,o,ctag);
    
  int i, count = 0;

  // Make sure we have a valid touching domain.

  PAssert(this->initialized());
  PAssert(contains(this->domain_m, d));

  // Find the starting and ending grid positions, and store as an
  // Interval.

  Interval<Dim> box = Pooma::NoInit();
  for (i = 0; i < Dim; ++i) 
    {
      int a, b;
        
      // This part is unchanged from touches. What we're going to do
      // here is simply extend the range in each direction by one block
      // and then let the intersection calculation below sort out if
      // there is actually an intersection. Otherwise we'd have to do
      // comparisons on all blocks up here and then still do the
      // intersection on the remaining blocks below. On average, this
      // should be slightly faster I think.
        
      if (!this->hasExternalGuards_m)
	{
	  a = (d[i].min() - this->firsti_m[i]) / blocksizes_m[i];
	  b = (d[i].max() - this->firsti_m[i]) / blocksizes_m[i];
	}
      else
	{
	  // If we're in the lower guards, this will fall
	  // through to give zero.
            
	  a = b = 0;
	  int pos = d[i].min();
	  int last = this->innerdomain_m[i].last();
	  int del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      a = del / blocksizes_m[i];
	    else
	      a = allDomain_m[i].last();
            
	  pos = d[i].max();
	  del = pos - this->firsti_m[i];
            
	  if (del >= 0)
	    if (pos <= last)
	      b = del / blocksizes_m[i];
	    else
	      b = allDomain_m[i].last();
	}
          
      // Now we check that we're not at the ends of the domain for the
      // brick of blocks, and extend the region accordingly.
        
      if (a > 0) --a;
      if (b < allDomain_m[i].last()) ++b;
        
      box[i] = Interval<1>(a, b);
    }

  // Figure the type of the domain resulting from the intersection.

  typedef typename 
    IntersectReturnType<Domain_t,OtherDomain>::Type_t OutDomain_t;
  OutDomain_t outDomain = Pooma::NoInit();

  // Generate the type of the node pushed on the output iterator.

  typedef Node<OutDomain_t,Domain_t> OutNode_t;

  // Iterate through the Interval grid positions.

  typename Interval<Dim>::const_iterator boxiter = box.begin();
    
  while (boxiter != box.end()) 
    {
      // Calculate the linear position of the current node.
        
      int indx = (*boxiter)[0].first();
      for (i = 1; i < Dim; ++i)
	indx += blockstride_m[i] * (*boxiter)[i].first();

      // Get that node, intersect the *allocated* domain with the 
      // requested one, and write the result out to the output iterator.
        
      PAssert(indx >= 0 && indx < this->remote_m.size());
        
      outDomain = intersect(d, this->remote_m[indx]->allocated());

      // We can no longer assert that outDomain is not empty since
      // we extended the search box without checking. Thus we now 
      // have an if tests around the output.
        
      if (!outDomain.empty())
	{
	  *o = touchesConstruct(outDomain,
				this->remote_m[indx]->allocated(),
				this->remote_m[indx]->affinity(),
				this->remote_m[indx]->context(),
				this->remote_m[indx]->globalID(),
				this->remote_m[indx]->localID(),
				ctag);
	}

      // Increment output iterator, count, and grid iterator.

      ++o;
      ++count;
      ++boxiter;
    }

  // Return the number of non-empty domains we found.
      
  return count;
}


// } // namespace POOMA


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformGridLayout.cpp,v $   $Author: richard $
// $Revision: 1.43 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo

