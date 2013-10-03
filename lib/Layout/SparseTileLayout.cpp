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
//   SparseTileLayout<Dim> template definitions.
//-----------------------------------------------------------------------------

#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/DomainRemoveOverlap.h"
#include "Layout/SparseTileLayout.h"
#include "Partition/GridPartition.h"
#include "Threads/PoomaSmarts.h"
#include "Utilities/PAssert.h"
#include <vector>
#include <map>
#include <iostream>

//============================================================
// SparseTileLayoutData non-inline method definitions
//============================================================

//-----------------------------------------------------------------------------
//
// template<int Dim>
// SparseTileLayoutData<Dim>::SparseTileLayoutData()
//
// Default constructor for SparseTileLayoutData, this sets up this object
// to look like an "empty" layout, with no patches and an empty domain
// with no guard cells.  The "initialize" method can be used to change
// it to a different state.
//
//-----------------------------------------------------------------------------

template<int Dim>
SparseTileLayoutData<Dim>::SparseTileLayoutData()
  : Observable<SparseTileLayoutData>(*this) 
{
}

template<int Dim>
SparseTileLayoutData<Dim>::SparseTileLayoutData(const Domain_t &boundingbox,
						const PatchList_t & PatchList,
						const ContextMapper<Dim> &cmap)
 : LayoutBaseData<Dim>(false,false,
		       GuardLayers_t(0),GuardLayers_t(0),
		       boundingbox,boundingbox),
   Observable<SparseTileLayoutData>(*this)
{
  this->blocks_m = Loc<Dim>();
  initialize(boundingbox,PatchList,cmap);
}


template<int Dim>
SparseTileLayoutData<Dim>::SparseTileLayoutData(const Domain_t &boundingbox,
						const GuardLayers_t & globalGL,
						const PatchList_t & PatchList,
						const ContextMapper<Dim> &cmap)
  : LayoutBaseData<Dim>(true,true,
			globalGL,globalGL,
			boundingbox,boundingbox),
  Observable<SparseTileLayoutData>(*this)
{
  this->blocks_m = Loc<Dim>();
  initialize(boundingbox,globalGL,PatchList,cmap);
}

template<int Dim>
SparseTileLayoutData<Dim>::SparseTileLayoutData(const Domain_t &boundingbox,
						const GuardLayers_t & internalGL,
						const GuardLayers_t & externalGL,
						const PatchList_t & PatchList,
						const ContextMapper<Dim> &cmap)
 : LayoutBaseData<Dim>(true,true,
			internalGL,externalGL,
		       boundingbox,boundingbox),
  Observable<SparseTileLayoutData>(*this)
{
  this->blocks_m = Loc<Dim>();
  initialize(boundingbox,internalGL,externalGL,PatchList,cmap );
}

template<int Dim>
SparseTileLayoutData<Dim>::SparseTileLayoutData(const Domain_t &boundingbox)
  : LayoutBaseData<Dim>(false, false,
			GuardLayers_t(),GuardLayers_t(),
			boundingbox,boundingbox),
  Observable<SparseTileLayoutData>(*this)
{
  this->blocks_m = Loc<Dim>();
  initialize(boundingbox);
}

template<int Dim>
SparseTileLayoutData<Dim>::SparseTileLayoutData(const Domain_t &boundingbox,
						const GuardLayers_t & internalGL,
						const GuardLayers_t & externalGL)
  :  LayoutBaseData<Dim>(true, true,
			internalGL,externalGL,
			boundingbox,boundingbox),
  Observable<SparseTileLayoutData>(*this)
{
  this->blocks_m = Loc<Dim>();
  initialize(boundingbox,internalGL,externalGL);
}

template<int Dim>
template<class Partitioner>
SparseTileLayoutData<Dim>::SparseTileLayoutData(const Domain_t &bbox,
						const Partitioner &gpar,
						const ContextMapper<Dim> &cmap)
   :  LayoutBaseData<Dim>(false,false,
			  GuardLayers_t(0),GuardLayers_t(0),
			  bbox,bbox),
    Observable<SparseTileLayoutData>(*this)
{
  this->blocks_m = Loc<Dim>();
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
    
  initialize(bbox, gpar,cmap);
}


template<int Dim>
SparseTileLayoutData<Dim>::~SparseTileLayoutData()
{
  for (typename List_t::iterator a = this->all_m.begin(); a != this->all_m.end(); ++a)
    delete (*a);
}

template<int Dim>
void SparseTileLayoutData<Dim>::syncPatch()
{
  typename List_t::iterator start = this->all_m.begin();
  typename List_t::iterator end = this->all_m.end();
  for ( ; start != end ; ++start)
    if ( (*start)->context() == Pooma::context()
	 ||(*start)->context() == -1 )
      this->local_m.push_back(*start);
    else
      this->remote_m.push_back(*start);
  
   // calculate the lookup maps
   calcMaps();

   calcAllocMaps();
   
   // Calculate what we need to do in a fill-guard-cell operation.
   
   calcGCFillList();
}

template<int Dim>
void SparseTileLayoutData<Dim>::calcMaps()
{

  // Don't need to do this if we're not initialized

  if (!this->initialized())
    return;

  map_m.zap();

  map_m.initialize(this->domain_m);

  typename List_t::const_iterator start = this->all_m.begin();
  typename List_t::const_iterator end   = this->all_m.end();
  int i=0;

  for ( ; start != end ; ++start, ++i )
    {
      pidx_t tmp((*start)->globalID(),i);
      // for Node, domain() returns the owned domain.
      typename DomainMap<Domain_t,pidx_t>::Value_t val((*start)->domain(),tmp);
      
      map_m.insert(val);
    }
  map_m.update();
}

template<int Dim>
void SparseTileLayoutData<Dim>::calcAllocMaps()
{
 // Don't need to do this if we're not initialized 

  if (!this->initialized())
    return;

  mapAloc_m.zap();

  mapAloc_m.initialize(this->domain_m);
 
  typename List_t::const_iterator start = this->all_m.begin();
  typename List_t::const_iterator end   = this->all_m.end();

  int i=0;

  for ( ; start!=end ; ++start , ++i)
    {   
      pidx_t tmp((*start)->globalID(),i);
      typename DomainMap<Domain_t,pidx_t>::Value_t val((*start)->allocated(),tmp);
						   
      mapAloc_m.insert(val);
    }
  mapAloc_m.update();

}



template<int Dim>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &gdom)
{
  this->blocks_m = Loc<Dim>();
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

  for(i=0;i<Dim;++i)
    this->firste_m[i] = this->firsti_m[i] = this->domain_m[i].first();

  // Initially, our total and owned domains are the same.

  this->domain_m = gdom;
  this->innerdomain_m = gdom;

  // Examine the partitioner for info about guard cells.  Change our
  // domains if necessary, and save guard cell info for later.

  this->hasInternalGuards_m = this->hasExternalGuards_m = false;
  this->internalGuards_m = this->externalGuards_m = GuardLayers_t();

}

template<int Dim>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &gdom,
					   const GuardLayers_t & globalGL)
{
  this->blocks_m = Loc<Dim>();
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

  this->domain_m = gdom;
  this->innerdomain_m = gdom;

  for(i=0;i<Dim;++i)
    this->firsti_m[i] = this->domain_m[i].first();

  this->hasInternalGuards_m = this->hasExternalGuards_m=true;
  this->internalGuards_m = this->externalGuards_m = globalGL;

  GuardLayers<Dim>::addGuardLayers(this->domain_m, this->externalGuards_m);
  
  for(i=0;i<Dim;++i)
    this->firste_m[i]=this->domain_m[i].first();
}


template<int Dim>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &gdom,
				     const GuardLayers_t & internalGL,
				     const GuardLayers_t & externalGL)
{
  this->blocks_m = Loc<Dim>();
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

  this->domain_m = gdom;
  this->innerdomain_m = gdom;

  for(i=0;i<Dim;++i)
    this->firsti_m[i]=this->domain_m[i].first();

  this->hasInternalGuards_m = this->hasExternalGuards_m = true;
  this->internalGuards_m = internalGL;
  this->externalGuards_m = externalGL;

  GuardLayers<Dim>::addGuardLayers(this->domain_m, this->externalGuards_m);
  
  for(i=0;i<Dim;++i)
    this->firste_m[i] = this->domain_m[i].first();
}

template<int Dim>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &bbox,
					   const PatchList_t &plist,
					   const ContextMapper<Dim> &cmap)
{
  this->blocks_m = Loc<Dim>();
  initialize(bbox);
  TilePartition<Dim> gpar(plist);
  gpar.partition(bbox,this->all_m,cmap);
  syncPatch();
}


template<int Dim>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &bbox,
					   const GuardLayers_t & globalGL,
					   const PatchList_t &plist,
					   const ContextMapper<Dim> &cmap)
{
  this->blocks_m = Loc<Dim>();
  initialize(bbox,globalGL);
  TilePartition<Dim> gpar(bbox,plist,globalGL);
  gpar.partition(bbox,this->all_m,cmap);
  syncPatch();
}

template<int Dim>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &bbox,
					   const GuardLayers_t & internalGL,
					   const GuardLayers_t & externalGL,
					   const PatchList_t &plist,
					   const ContextMapper<Dim> &cmap)
{
  this->blocks_m = Loc<Dim>();
  initialize(bbox,internalGL,externalGL);
  TilePartition<Dim> gpar(bbox,plist,internalGL,externalGL);
  gpar.partition(bbox,this->all_m,cmap);

  syncPatch();

}


template<int Dim>
template<class Partitioner>
void SparseTileLayoutData<Dim>::initialize(const Domain_t &bbox,
					   const Partitioner &gpar,
					   const ContextMapper<Dim> &cmap)
{
  this->blocks_m = Loc<Dim>();
  initialize(bbox,gpar.internalGuards(),gpar.externalGuards());
  gpar.partition(bbox,this->all_m,cmap);
  syncPatch();
}

template<int Dim>
void SparseTileLayoutData<Dim>::calcGCFillList()
{
  if(!this->initialized() || !this->hasInternalGuards_m)
    return;
  
  this->gcFillList_m.clear();
  gcBorderFillList_m.clear();

  typedef Node<Domain_t,AllocatedDomain_t> NNode_t;
  typedef std::vector<NNode_t>             TouchList_t;
  TouchList_t tlist;

  // first we do the internal overlap regions
  typename List_t::iterator start = this->all_m.begin();
  typename List_t::iterator end = this->all_m.end();

  for ( ; start!=end; ++start)
    {
      touches((*start)->allocated(),
	      std::back_inserter(tlist),
	      TouchesConstructNodeObj());

      // now pack the tlist into the GCFillInfo object
      // The if test is to remove the self-touch entry
      typename TouchList_t::iterator GCLstart = tlist.begin();
      typename TouchList_t::iterator GCLend = tlist.end();
      
      for( ; GCLstart != GCLend ;++GCLstart)
	{
	  if(GCLstart->globalID() == (*start)->globalID()) 
	    {
	      tlist.erase(GCLstart);
	      break;
	    }
	}

      GCLstart = tlist.begin();
      GCLend = tlist.end();
      for ( ; GCLstart!=GCLend ; ++GCLstart )
	{
	  //  removed the external guard layer area. 
	  this->gcFillList_m.push_back(GCFillInfo_t((*GCLstart).domain(),
						    (*GCLstart).globalID(),
						    (*start)->globalID()));
	}
      tlist.clear();
    }

  // Next, generate the list of all the internalGuardLayer
  // regions, and 'subtract' the above list from it to produce
  // a list of internal guard layers that will be filled externally.
 
  std::vector<GCBorderFillInfo> bfv;

  start = this->all_m.begin();
  for ( ; start!=this->all_m.end(); ++start)
    {
      for(int d=0;d<Dim;++d)
	{
	  Domain_t gcdom( (*start)->allocated());
	  
	  int max = (*start)->allocated()[d].last();
	  int min = max - this->internalGuards_m.upper(d) + 1;
	  gcdom[d] = Interval<1>(min, max);
	  gcdom = intersect(this->innerdomain_m,gcdom);
	  if (gcdom.size()>0) 
	    {
	      bfv.push_back(GCBorderFillInfo(gcdom, (*start)->globalID() ));
	    }
	  // now do the other side of this dimension. 
	  gcdom = (*start)->allocated();
	  min = (*start)->allocated()[d].first();
	  max = min + this->internalGuards_m.lower(d) -1;
 	  gcdom[d] = Interval<1>(min, max);
	  gcdom = intersect(this->innerdomain_m,gcdom);
	  
	  if (gcdom.size() > 0)
	  {
	    bfv.push_back(GCBorderFillInfo(gcdom,(*start)->globalID())); 
	  }
	} 
    }

  // remove overlap of GCFillInfo on GCBorderFillInfo
  std::vector<Domain_t> temp2,temp3,temp4;

  std::vector<GCBorderFillInfo> temp;

  BorderFillIterator_t bst = bfv.begin();
  BorderFillIterator_t ben = bfv.end();

  for ( ; bst != ben ; ++bst)
    {
      FillIterator_t gst = this->beginFillList();
      FillIterator_t gen = this->endFillList(); 
      temp2.clear();
      temp2.push_back(bst->domain());
      
      for ( ; gst!=gen ; ++gst )
	{
	  typename std::vector<Domain_t>::iterator ts = temp2.begin();
	  for ( ; ts != temp2.end() ; ++ts )
	    {
	      temp3 = DomainRemoveOverlap(*ts,gst->domain_m);
	      temp4.insert(temp4.end(),temp3.begin(),temp3.end());
	    }
	  temp2 = temp4;
	  temp4.clear();
	}
      
      typename std::vector<Domain_t>::iterator ts = temp2.begin();
      for( ; ts != temp2.end(); ++ts)
	temp.push_back(GCBorderFillInfo(*ts, bst->patchID() ));
    }
  
  //gcBorderFillList_m.clear();
  gcBorderFillList_m = temp;
}
//-----------------------------------------------------------------------------
//
// globalID takes a position within the domain of the layout, and returns
// the global ID for that node.
//
//-----------------------------------------------------------------------------

template<int Dim>
int SparseTileLayoutData<Dim>::globalID(const Loc<Dim> &loc) const
{
  // Make sure the point is in our domain.

  PAssert(contains(this->domain_m, loc));

  DomainMapTouchIterator<Interval<Dim>,pidx_t> dmti = 
	(map_m.touch(Interval<Dim>(loc))).first;
  DomainMapTouchIterator<Interval<Dim>,pidx_t> baditerator;
  
  PInsist(dmti!=baditerator,"Bad location requested in SparseTileLayout");
  return (*dmti).first;
}

template <int Dim>
int SparseTileLayoutData<Dim>::globalID(int i0) const
{
  // Make sure the point is in our domain.

  PAssert(Dim == 1);

  // Call the Loc version.

  Loc<Dim> loc;
  loc[0] = i0;
  return globalID(loc);
}

template <int Dim>
int SparseTileLayoutData<Dim>::globalID(int i0, int i1) const
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
int SparseTileLayoutData<Dim>::globalID(int i0, int i1, int i2) const
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
int SparseTileLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3) const
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
int SparseTileLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
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
int SparseTileLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
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
int SparseTileLayoutData<Dim>::globalID(int i0, int i1, int i2, int i3,
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

template<int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int SparseTileLayoutData<Dim>::touches(const OtherDomain &fulld, 
				       OutIter o, 
				       const ConstructTag &ctag) const
{
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
  // Ack!!! OutDomain_t is a Range, but DomainMap::touches requires an Interval.
  typename DomainMap<Interval<Dim>,pidx_t>::Touch_t dmti =
    map_m.touch(Interval<Dim>(d));
  typename DomainMap<Interval<Dim>,pidx_t>::touch_iterator a;

  int count = 0;

  for( a = dmti.first ;a != dmti.second; ++a)
    {
      int nodeListIndex = (*a).second;

      outDomain =  intersect(a.domain(), fulld);
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
  return count;
}

template<int Dim>
template <class OtherDomain, class OutIter, class ConstructTag>
int SparseTileLayoutData<Dim>::touchesAlloc(const OtherDomain &fulld, 
					    OutIter o, 
					    const ConstructTag &ctag) const
{
  int i;

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

  typename DomainMap<Interval<Dim>,pidx_t>::Touch_t dmti = map_m.touch(d);
  typename DomainMap<Interval<Dim>,pidx_t>::touch_iterator a;

  int count = 0;

  for( a = dmti.first ;a != dmti.second; ++a)
    {
      int nodeListIndex = (*a).second;

     
      outDomain =  intersect(a.domain(), fulld);
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
  return count;
}

// print the SparseTileLayoutData object
template<int Dim>
template<class Out>
void SparseTileLayoutData<Dim>::print(Out & o) const
{
  int i;
  o<< " SparseTileLayoutData<"<<Dim<<">: "<<std::endl;
  o<< " ID_m " << this->ID_m << std::endl;
  o<< " domain_m " << this->domain_m <<std::endl;
  o<< " innerdomain_m " << this->innerdomain_m <<std::endl;
  o<< " all_m : " << std::endl;
  typename List_t::const_iterator start = this->all_m.begin();
  typename List_t::const_iterator end   = this->all_m.end();
  for ( ; start!=end ; ++start)
    o<< (*start)->globalID()<<" "<<
      (*start)->domain()<<" "<< 
      (*start)->allocated()<<" "
     <<std::endl; 
  o<< " local_m : " << std::endl;
  start = this->local_m.begin();
  end = this->local_m.end();
  for ( ; start!=end ; ++start)
    o<< (*start)->globalID()<<" "<<
      (*start)->localID()<<" " <<
      (*start)->domain()<<" "<< 
      (*start)->allocated()<<" "
     <<std::endl;
  o<< " firste_m[Dim] " ;
  for ( i=0;i<Dim;++i) o<< this->firste_m[i]<<" ";
  o<< std::endl;
  o<< " firsti_m[Dim] " ;
  for ( i=0;i<Dim;++i) o<< this->firsti_m[i]<<" ";
  o<< std::endl;
  o<< " hasInternalGuards_m, hasExternalGuards_m " <<
    this->hasInternalGuards_m <<" " << this->hasExternalGuards_m <<std::endl;
  o<< " internalGuards_m " ;
   for ( i=0;i<Dim;++i) 
     o<< this->internalGuards_m.upper(i)<<"-"<<this->internalGuards_m.lower(i)<<" ";
   o<<std::endl;
  o<< " externalGuards_m " ;
   for ( i=0;i<Dim;++i) 
     o<< this->externalGuards_m.upper(i)<<"-"<<this->externalGuards_m.lower(i)<<" ";
   o<<std::endl;
   
   FillIterator_t gstart = this->gcFillList_m.begin();
   FillIterator_t gend = this->gcFillList_m.end();
   o<< " gcFillList_m " <<std::endl;
   for( ; gstart!=gend ; ++gstart)
     o<<"       "
      <<gstart->domain_m<<" "
      <<gstart->ownedID_m<<" "
      <<gstart->guardID_m<<std::endl;

   BorderFillIterator_t bgstart = gcBorderFillList_m.begin();
   BorderFillIterator_t bgend =  gcBorderFillList_m.end();
   o<< " gcBorderFillList_m " <<std::endl;
   for( ; bgstart!=bgend ; ++bgstart)
     o<<"       "
      <<bgstart->domain()<<" "
      <<bgstart->patchID()<<std::endl;
}



//-----------------------------------------------------------------------
// SparseTileLayout member function
//-----------------------------------------------------------------------


template<int Dim>
void SparseTileLayout<Dim>::syncPatch()
{
  this->pdata_m->syncPatch();
}

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout()
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >(new LayoutData_t()),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(Domain_t & boundingbox,
					const PatchList_t &patchlist,
					const ReplicatedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(boundingbox,patchlist,LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(Domain_t & boundingbox,
					const PatchList_t &patchlist,
					const DistributedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(boundingbox,patchlist,DistributedMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}
  template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &boundingbox,
					const GuardLayers_t & globalGL,
					const PatchList_t  & patchlist,
					const DistributedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
  (new LayoutData_t(boundingbox,globalGL,patchlist,DistributedMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}
template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &boundingbox,
					const GuardLayers_t & globalGL,
					const PatchList_t  & patchlist,
					const ReplicatedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
  (new LayoutData_t(boundingbox,globalGL,patchlist,LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & internalGL,
		   const GuardLayers_t & externalGL,
		   const PatchList_t  & patchlist,
		   const DistributedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
              (new LayoutData_t(boundingbox,
				internalGL,
				externalGL,
				patchlist,
				DistributedMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & internalGL,
		   const GuardLayers_t & externalGL,
		   const PatchList_t  & patchlist,
		   const ReplicatedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
              (new LayoutData_t(boundingbox,
				internalGL,
				externalGL,
				patchlist,
				LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &boundingbox)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(boundingbox)),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}
  
 template<int Dim> 
SparseTileLayout<Dim>:: SparseTileLayout(const Domain_t &boundingbox,
		   const GuardLayers_t & globalGL)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(boundingbox,globalGL)),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}  

template<int Dim> 
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &boundingbox,
					const GuardLayers_t & internalGL,
					const GuardLayers_t & externalGL)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(boundingbox,internalGL,externalGL)),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
} 


template <int Dim>
template <class Partitioner>
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &bbox, 
					const Partitioner &gpar,
					const DistributedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(bbox, gpar,DistributedMapper<Dim>(gpar))),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}
template <int Dim>
template <class Partitioner>
SparseTileLayout<Dim>::SparseTileLayout(const Domain_t &bbox, 
					const Partitioner &gpar,
					const ReplicatedTag &)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >
    (new LayoutData_t(bbox, gpar,LocalMapper<Dim>())),
  Observable<This_t>(*this)
{
  this->pdata_m->attach(*this);
}
  // Copy constructor & assignment operator
  // Shallow copies with reference counting.
  
template<int Dim>
SparseTileLayout<Dim>::SparseTileLayout(const This_t &model)
  : LayoutBase<Dim,SparseTileLayoutData<Dim> >(model.pdata_m),
  Observable<This_t>(*this)
{ 
   this->pdata_m->attach(*this);
}
