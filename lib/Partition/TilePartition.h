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

#ifndef POOMA_PARTITION_TILE_PARTITION_H
#define POOMA_PARTITION_TILE_PARTITION_H

//-----------------------------------------------------------------------------
// Classes:
// TilePartition<Dim>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Partition
 * @brief
 * A layout partitioner that will break a given global
 * domain into blocks specified by a domain Tile. 
 *
 * The user must provide consistent information; if the subdomain bounds
 * are greater than the global domain, an insist failure will result.
 * 
 * Secondarily, if Global and Internal guard cells are specified, those
 * internal guard cell regions may not span more than the adjacent patch;
 * violations will result in an insist failure. 
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Layout/Node.h"
#include "Partition/ContextMapper.h"
#include "Layout/GuardLayers.h"
#include "Utilities/PAssert.h"
#include <iosfwd>


///////////////////////////////////////////////////////////////////////////
// namespace POOMA {

template<int Dim>
class DefaultTPmapper
  : public ContextMapper<Dim>
{ 
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef Interval<Dim>                       Domain_t;
  typedef Node<Domain_t>                      Value_t;
  typedef std::vector<Value_t *>              List_t;

 

  template<class Partitioner>
  DefaultTPmapper(const Partitioner & gp)
  {
  }

  void map(const List_t & templist) const;
     
  // member data
private:

};

template<int Dim>
void DefaultTPmapper<Dim>::map(const List_t & templist) const
{
  int ncontexts = Pooma::contexts();
  int npc = templist.size()/ncontexts;
  if (templist.size()%ncontexts!=0) ++npc;

  typename List_t::const_iterator start = templist.begin();
  typename List_t::const_iterator end = templist.end();
  int c = 0;
  int p = 0;
  for ( ; start!=end ; ++start)
    {
      (*start)->context() = p;
      if(p == Pooma::context())
	(*start)->localID() = c;

      ++c;
      if (c > npc)
	{
	  ++p;
	  c = 0;
	}
    }

  int affinityMax = Smarts::concurrency();
  int idMax = 0;

  start = templist.begin();
  for ( ; start != end ; ++start)
    if((*start)->context()==Pooma::context())
      {
	(*start)->localID()=idMax;
	++idMax;
      }
  start = templist.begin();
  for ( ; start != end ; ++start)
    { 
      if((*start)->context()==Pooma::context())
	(*start)->affinity() = static_cast<int>( affinityMax * 
						 ( (*start)->localID() /
			       static_cast<double>(idMax) ) );
    }
}



/**
 * TilePartition is a layout partitioner; it is created with the info
 * needed to partition a global domain into a (possibly sparse) 
 * list of Dim dimensional non-overlapping patches. 
 *
 *
 * A TilePartition object is constructed with the following information:
 *   TilePartition())
 *   Creates an empty partition, with no guard cells
 * 
 *  TilePartition(const PatchList_t &plist)
 * Creates patches from plist
 *
 *  TilePartition(const PatchList_t &plist, 
 *                       const GuardLayers<Dim> &gcs)
 * Same as above, with internal and external guard cell sizes set to gcs.
 * 
 *                       
 * TilePartition(const PatchList_t &plist, 
 *                       const GuardLayers<Dim> &igcs,
 *                       const GuardLayers<Dim> &egcs)
 * Same as above, with internal and external guard cell sizes specified 
 * independently.
 */

template<int Dim>
class TilePartition
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef LocalMapper<Dim>                  DefaultMapper_t;
  typedef Interval<Dim>                     Domain_t;
  typedef std::vector<Domain_t>             PatchList_t;
  typedef Node<Domain_t>                    Value_t;
  typedef std::vector<Value_t *>            List_t;

  // These traits are used by the layouts to determine if
  // a particular partitioner is appropriate for that layout.

  enum { uniform  = false };
  enum { gridded  = false };
  enum { tile     = true };
  enum { general  = true };  // <<<<< Is this correct? ? ? 

  enum { dimensions = Dim };

  //============================================================
  // Constructors
  //============================================================

  // The default constructor
  TilePartition() 
    : hasInternalGuards_m(false), 
      hasExternalGuards_m(false),
      num_m(0),
      internalGuards_m(0),
      externalGuards_m(0)
  {
  }

  TilePartition(const PatchList_t &pList)
    : hasInternalGuards_m(false), 
      hasExternalGuards_m(false),
      internalGuards_m(0),
      externalGuards_m(0),
      tile_m(pList)
  {
    num_m = pList.size();
  }
 
  TilePartition(const PatchList_t &pList,
		const GuardLayers<Dim> &gcs)
    : hasInternalGuards_m(true), 
      hasExternalGuards_m(false), 
      internalGuards_m(gcs), 
      externalGuards_m(gcs),
      tile_m(pList)
      
  { 
    num_m = tile_m.size();
  }

  TilePartition(const PatchList_t &pList,
		const GuardLayers<Dim> &igcs, 
		const GuardLayers<Dim> &egcs)
    : hasInternalGuards_m(true), 
      hasExternalGuards_m(true), 
      internalGuards_m(igcs), 
      externalGuards_m(egcs),
      tile_m(pList)
  {
    num_m = tile_m.size();
  }

  // copy constructor

  TilePartition(const TilePartition<Dim> & b)
    : hasInternalGuards_m(b.hasInternalGuards_m),
      hasExternalGuards_m(b.hasExternalGuards_m),
      internalGuards_m(b.internalGuards_m),
      externalGuards_m(b.externalGuards_m),
      tile_m(b.tile_m),
      num_m(b.num_m)
  {
  }

  //============================================================
  // Destructor
  //============================================================

  // The destructor for this partitioner does not have to do anything.
  
  ~TilePartition() { }

  //============================================================
  // Operators
  //============================================================

  // Assignment operator.

  TilePartition<Dim> &operator=(const TilePartition<Dim> &g)
  {
    if (this != &g)
      {
	hasInternalGuards_m = g.hasInternalGuards_m;
	hasExternalGuards_m = g.hasExternalGuards_m;
	internalGuards_m    = g.internalGuards_m;
	externalGuards_m    = g.externalGuards_m;
	tile_m              = const_cast<TilePartition<Dim>&>(g).tileList();
	num_m               = g.maxSize();
	if (!hasInternalGuards_m)
	  internalGuards_m = GuardLayers<Dim>(0);
	if (!hasExternalGuards_m)
	  externalGuards_m = GuardLayers<Dim>(0);
      }
    return *this;
  }

  //============================================================
  // Accessors
  //============================================================

  int maxSize() const { return num_m; }

  PatchList_t tileList() { return tile_m; }
  
  // Guard layer info:

  bool hasGuards() const { return hasInternalGuards_m||hasExternalGuards_m;  }

  bool hasCustomEdgeGuards() const 
  {
    if (hasInternalGuards_m&&!hasExternalGuards_m) return true;
    if (!hasInternalGuards_m&&hasExternalGuards_m) return true;
    if (hasInternalGuards_m&&hasExternalGuards_m
       &&(internalGuards_m!=externalGuards_m)) return true;
    return false;
  }

  bool hasInternalGuards() const { return hasInternalGuards_m; }

  bool hasExternalGuards() const { return hasExternalGuards_m;}

  const GuardLayers<Dim> &internalGuards() const 
  { 
    return internalGuards_m; 
  }
  
  const GuardLayers<Dim> &externalGuards() const 
  { 
    return externalGuards_m; 
  }
 
  //============================================================
  // Partition methods
  //============================================================

  template<class D>
  int partition(const D &bbox,List_t &all,const ContextMapper<Dim> &cmapper) const;
  
  template<class D>
  int partition(const D &bbox,List_t &all) const
  {
    return partition(bbox,all,LocalMapper<Dim>(*this));
  }
 
  template<class Out>
  void print(Out &o) const;

private:

  // Do we have guard layers?
  
  bool hasInternalGuards_m;
  
  // Are the external guards different from the internal?
  
  bool hasExternalGuards_m;
  
  // Specification of internal guard layers.
  
  GuardLayers<Dim> internalGuards_m;
  
  // Specification of external guard layers.
  
  GuardLayers<Dim> externalGuards_m;

  // The total number of blocks to create.
  
  int num_m;
  
  // The list of tiles or patches. 

  PatchList_t tile_m;

};

template<int Dim>
template<class D>
int TilePartition<Dim>::partition(const D &bbox, List_t &all,
				  const ContextMapper<Dim> &cmapper) const
{
  typedef typename DomainTraits<Domain_t>::Element_t Element_t;

  // Make sure we have the right dimensionality.

  CTAssert(Dim == DomainTraits<Domain_t>::dimensions);

  typename PatchList_t::const_iterator start = tile_m.begin();
  typename PatchList_t::const_iterator end = tile_m.end();
 
  while (start!=end)
    {
      Domain_t o = Pooma::NoInit();	
      o = * start;
      Domain_t oo = o;
      Domain_t a = Pooma::NoInit();	
      a = * start;

      if (hasInternalGuards()||hasExternalGuards())
	{
	  for (int i=0;i<Dim;i++)
	    {
	      if (oo[i].first() == bbox[i].first())
		{
		  if (hasExternalGuards())
		    {
		      o[i]=Interval<1>(o[i].first()-externalGuards().lower(i),
				       o[i].last());
		      a[i]=Interval<1>(a[i].first()-externalGuards().lower(i),
				       a[i].last());
		    }
		  if (hasInternalGuards() && oo[i].last() != bbox[i].last() )
		    a[i]=Interval<1>(a[i].first(),
				     a[i].last()+internalGuards().upper(i));
		}
	      if (oo[i].last()== bbox[i].last())
		{
		  if (hasExternalGuards())
		    {
		      o[i]=Interval<1>(o[i].first(),
				       o[i].last()+externalGuards().upper(i));
		      a[i]=Interval<1>(a[i].first(),
				       a[i].last()+externalGuards().upper(i));
		    }
		  if (hasInternalGuards()&&(oo[i].first() != bbox[i].first()))
		    a[i]=Interval<1>(a[i].first()-internalGuards().lower(i),
				     a[i].last());
		}
	      if (oo[i].first()!=bbox[i].first() &&
		  oo[i].last() != bbox[i].last() &&
		  hasInternalGuards()) // it's a fully internal patch
		a[i]=Interval<1>(o[i].first()-internalGuards().lower(i),
				 o[i].last()+internalGuards().upper(i));

	    }
	}

      // Add the domain to the layout.

      Value_t * node = new Value_t(o,a,-1,all.size(),-1);
      all.push_back(node);

      // Increment our counters and iterators.

      ++start;
    }

  cmapper.map(all);

  return all.size();
}

template<int Dim>
template<class Out>
void TilePartition<Dim>::print(Out &o) const
{
  int i;
  o << "TilePartition<" << Dim << ">:" << std::endl;
  o << "  hasInternalGuards_m  hasExternalGuards_m = ";
  o << hasInternalGuards_m<< " "<<hasExternalGuards_m<< std::endl;
  o << "  internalGuards_m:" << std::endl;
  o << "      upper       ";
  for (i=0; i < Dim; ++i)
    o << internalGuards_m.upper(i) << " ";
  o << std::endl;
  o << "      lower       ";
  for (i=0; i < Dim; ++i)
    o << internalGuards_m.lower(i) << " ";
  o << std::endl;
  o << "  externalGuards_m:" << std::endl;
  o << "      upper       ";
  for (i=0; i < Dim; ++i)
    o << externalGuards_m.upper(i) << " ";
  o << std::endl;
  o << "      lower       ";
  for (i=0; i < Dim; ++i)
    o << externalGuards_m.lower(i) << " ";
  o << std::endl;
  o << "  num_m = " << num_m << std::endl;
}


/// A specialization of the Inform traits used to say that node has
/// a print method.

template <int Dim>
std::ostream &operator<<(std::ostream &o, const TilePartition<Dim> &gp)
{
  gp.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTITION_GRID_PARTITION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TilePartition.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:17:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
