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

#ifndef POOMA_PARTITION_GRID_PARTITION_H
#define POOMA_PARTITION_GRID_PARTITION_H

//-----------------------------------------------------------------------------
// Classes:
// GridPartition<Dim>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Partition
 * @brief
 * A layout partitioner that will break a given global
 * domain into blocks specified by a domain Grid. 
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
#include "Domain/Grid.h"
#include "Layout/GuardLayers.h"
#include "Utilities/PAssert.h"
#include "Partition/ContextMapper.h"
#include "Partition/DistributedMapper.h"
#include <iosfwd>
#include <vector>
#include <list>

///////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/// Utility function to return a Grid from a global domain and a loc blocks

template <int Dim>
Grid<Dim> makeRGrid(const Interval<Dim> & gdom, const Loc<Dim> & blocks);

// forward declarations

template<int Dim>
class UniformGridPartition;



/**
 * GridPartition is a layout partitioner; it is created with the info
 * needed to partition a global domain into subdomains using N 1-dim upper
 * sub-domain specifications along each axis, or any of the specifiers used
 * for the UniformGridPartition. 
 *
 * GridPartition inherits from UniformGridPartition.
 *
 * A GridPartition object is constructed with the following information:
 *   GridPartition()
 *   Creates one partition, with no guard cells
 *
 *  GridPartition(const Loc<Dim> &n, int p=-1)
 * Creates n[i] blocks along each i'th dimension
 *
 *  GridPartition(const Loc<Dim> &n, 
 *                       const GuardLayers<Dim> &gcs)
 * Same as above, with internal and external guard cell sizes set to gcs.
 * 
 *                       
 * GridPartition(const Loc<Dim> &n, 
 *                       const GuardLayers<Dim> &igcs,
 *                       const GuardLayers<Dim> &egcs)
 * Same as above, with internal and external guard cell sizes specified 
 * independently.
 *
 *  GridPartition(const Grid<Dim> &g)
 * Partitions according to the Grid object.
 *
 *  GridPartition(const Grid<Dim> &g, 
 *                       const GuardLayers<Dim> &gcs)
 * Same as above, with internal and external guard cell sizes set to gcs.
 * 
 *                       
 * GridPartition(const Grid<Dim> &g, 
 *                       const GuardLayers<Dim> &igcs,
 *                       const GuardLayers<Dim> &egcs)
 * Same as above, with internal and external guard cell sizes specified 
 * independently.
 */

template<int Dim>
class GridPartition
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================
  // **** This DefaultMapper_t typedef needs to be changed...!
  typedef LocalMapper<Dim>        DefaultMapper_t;
  typedef Interval<Dim>           Domain_t;
  typedef Node<Domain_t>          Value_t;
  typedef std::vector<Value_t *>  List_t;  

  // These traits are used by the layouts to determine if
  // a particular partitioner is appropriate for that layout.

  enum { uniform  = false };
  enum { gridded  = true  };
  enum { tile     = false };
  enum { general  = false };

  enum { dimensions = Dim };

  //============================================================
  // Constructors
  //============================================================


  GridPartition(const Grid<Dim> &g)
    : hasInternalGuards_m(false), 
      hasExternalGuards_m(false),
      internalGuards_m(0),
      externalGuards_m(0),
      grid_m(g)
  {
    num_m=1;
    for (int i=0;i<Dim;i++)
      {
	blocks_m[i] = Loc<1>(grid_m[i].size()-1);
	num_m*=blocks_m[i].first();
      }
    
  }
  
  GridPartition(const Grid<Dim> &g, 
		const GuardLayers<Dim> &gcs)
    : hasInternalGuards_m(true),
      hasExternalGuards_m(true),
      internalGuards_m(gcs),
      externalGuards_m(gcs),
      grid_m(g)
  {
    num_m=1;
    for (int i=0;i<Dim;i++)
      {
	blocks_m[i] = Loc<1>(grid_m[i].size()-1);
	num_m*=blocks_m[i].first();
      }
  }
                   
  GridPartition(const Grid<Dim> &g, 
		const GuardLayers<Dim> &igcs,
		const GuardLayers<Dim> &egcs)
    : hasInternalGuards_m(true), 
      hasExternalGuards_m(true), 
      internalGuards_m(igcs), 
      externalGuards_m(egcs),
      grid_m(g)
  {
    num_m=1;
    for (int i=0;i<Dim;i++)
      {
	blocks_m[i] = Loc<1>(grid_m[i].size()-1);
	num_m*=(grid_m[i].size()-1);
      }
  }

  // The default constructor
  GridPartition() 
    : hasInternalGuards_m(false), 
      hasExternalGuards_m(false),
      internalGuards_m(0),
      externalGuards_m(0),
      num_m(1)
  {  
    for (int i=0;i<Dim;++i)
      blocks_m[i]=Loc<1>(1);
  }

  GridPartition(const Loc<Dim> &a)
    : blocks_m(a),
      hasInternalGuards_m(false),
      hasExternalGuards_m(false),
      internalGuards_m(0),
      externalGuards_m(0)
  {
    num_m = blocks_m[0].first();
    for (int d=1; d < Dim; ++d)
      num_m *= blocks_m[d].first();
  }
 
  GridPartition(const Loc<Dim> &a,
		const GuardLayers<Dim> &gcs)
    : blocks_m(a),
      hasInternalGuards_m(true), 
      hasExternalGuards_m(true), 
      internalGuards_m(gcs), 
      externalGuards_m(gcs)
  { 
    num_m = blocks_m[0].first();
    for (int d=1; d < Dim; ++d) 
      num_m *= blocks_m[d].first();
  }

  GridPartition(const Loc<Dim> &a,
		const GuardLayers<Dim> &igcs, 
		const GuardLayers<Dim> &egcs)
    : blocks_m(a),
      hasInternalGuards_m(true), 
      hasExternalGuards_m(true), 
      internalGuards_m(igcs), 
      externalGuards_m(egcs)
  {
    num_m = blocks_m[0].first();
    for (int d=1; d < Dim; ++d) 
      num_m *= blocks_m[d].first();
  }

  // copy constructor

  GridPartition(const GridPartition<Dim> & b)
    : blocks_m(b.blocks_m),
      hasInternalGuards_m(b.hasInternalGuards_m),
      hasExternalGuards_m(b.hasExternalGuards_m),
      internalGuards_m(b.internalGuards_m),
      externalGuards_m(b.externalGuards_m),
      num_m(b.num_m),
      grid_m(b.grid_m)
  {
  }

  // copy from a GP
  
 GridPartition(const UniformGridPartition<Dim> & b)
   : blocks_m(b.blocks_m),
     hasInternalGuards_m(b.hasInternalGuards_m),
     hasExternalGuards_m(b.hasExternalGuards_m),
     internalGuards_m(b.internalGuards_m),
     externalGuards_m(b.externalGuards_m),
     num_m(b.num_m)
  {}

      
  //============================================================
  // Destructor
  //============================================================

  // The destructor for this partitioner does not have to do anything.
  
  ~GridPartition() { }

  //============================================================
  // Operators
  //============================================================

  // Assignment operator.

  GridPartition<Dim> &operator=(const GridPartition<Dim> &g)
  {
    if (this != &g)
      {
	hasInternalGuards_m = g.hasInternalGuards_m;
	hasExternalGuards_m = g.hasExternalGuards_m;
	internalGuards_m    = g.internalGuards_m;
	externalGuards_m    = g.externalGuards_m;
	blocks_m            = g.blocks();
	num_m               = g.maxSize();
	grid_m              = g.grid();
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

  // Return all the blocks for the dimensions, as a Loc.
  
  const Loc<Dim> &blocks() const { return blocks_m; }
  
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

  const Grid<Dim> &grid() const { return grid_m;}

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
  int partition(const D &domain,
		List_t & all,
		const ContextMapper<Dim> &cmapper ) const;
  
  template<class D>
  int partition(const D &domain,
		List_t & all) const
  {
    return partition(domain,all,DefaultMapper_t(*this));
  }

  template<class Out>
  void print(Out &o) const;

private:

  // The number of blocks along each dimension.
  
  Loc<Dim> blocks_m;
  
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
  
  // Grid object for  Grid Layouts
  
  Grid<Dim> grid_m;
};


template<int Dim>
template<class D>
int GridPartition<Dim>::partition(const D &domain,
				  List_t & all,
				  const ContextMapper<Dim> &cmapper ) const
{
  typedef typename DomainTraits<Domain_t>::Element_t Element_t;

  // Make sure we have the right dimensionality.

  CTAssert(Dim == DomainTraits<D>::dimensions);
  CTAssert(Dim == DomainTraits<Domain_t>::dimensions);

  Grid<Dim> tgrid = grid();

  // if an empty domain was passed, we make empty patches
  if (domain.empty())
    {
      int np = 1;
      for (int i=0;i<Dim;++i)
	np*=blocks()[i].first();

      int start=0;
      Domain_t o;
      Domain_t a;
      while (start<np)
	{    
	  int gid = all.size();
	  int lid = -1;
	  Value_t *node = new Value_t(o, a, -1, gid, lid);
	  all.push_back(node);
	  ++start;
	}
      cmapper.map(all);
      return maxSize(); 
    }


  if (tgrid.empty()&&!domain.empty() )
    {
      tgrid = makeRGrid(domain,blocks_m);
    }

  typename Grid<Dim>::blockIterator start = tgrid.beginBlock();
  typename Grid<Dim>::blockIterator end = tgrid.endBlock();
 
  while (start!=end)
    {
      Loc<Dim> idx = start.point();
      Domain_t o = Pooma::NoInit();	
      o = * start;
      Domain_t a = Pooma::NoInit();	
      a = * start;

      // Calculate the guard cell specification for this domain.

      if (hasInternalGuards()||hasExternalGuards())
	{
	  for (int i=0;i<Dim;i++)
	    {
	      if (idx[i]==0)
		{
		  if (hasExternalGuards())
		    {
		      o[i]=Interval<1>(o[i].first()-externalGuards().lower(i),
				       o[i].last());
		      a[i]=Interval<1>(a[i].first()-externalGuards().lower(i),
				       a[i].last());
		    }
		  if (hasInternalGuards() && idx[i]!=(blocks()[i].first()-1))
		    a[i]=Interval<1>(a[i].first(),
				     a[i].last()+internalGuards().upper(i));
		}
	      if (idx[i]==blocks()[i].first()-1)
		{
		  if (hasExternalGuards())
		    {
		      o[i]=Interval<1>(o[i].first(),
				       o[i].last()+externalGuards().upper(i));
		      a[i]=Interval<1>(a[i].first(),
				       a[i].last()+externalGuards().upper(i));
		    }
		  if (hasInternalGuards()&&(idx[i]!=0))
		    a[i]=Interval<1>(a[i].first()-internalGuards().lower(i),
				     a[i].last());
		}
	      if (idx[i]!=0&&
		  idx[i]!=(blocks()[i].first()-1)&&
		  hasInternalGuards()) // it's a fully internal patch
		a[i]=Interval<1>(o[i].first()-internalGuards().lower(i),
				 o[i].last()+internalGuards().upper(i));

	    }
	}

      // Add the domain to the layout.
      int gid = all.size();
      int lid = -1;
      Value_t *node = new Value_t(o, a, -1, gid, lid);
      all.push_back(node);

      // Increment our counters and iterators.

      ++start;
    }

  cmapper.map(all);
  return maxSize();
}

template<int Dim>
template<class Out>
void GridPartition<Dim>::print(Out &o) const
{
  int i;
  o << "GridPartition<" << Dim << ">:" << std::endl;
  o << "  blocks_m = " << blocks_m << std::endl;
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
  o << "  grid_m = ";
  if (grid_m.empty() )
    o << "(empty)" << std::endl;
  else
    o << grid_m << std::endl;
}


/// A specialization of the Inform traits used to say that node has
/// a print method.

template <int Dim>
std::ostream &operator<<(std::ostream &o, const GridPartition<Dim> &gp)
{
  gp.print(o);
  return o;
}

/// Utility function to create a near-uniform partitioning of a domain. 
/// The result is represented by a Grid object. 

template <int Dim>
Grid<Dim> makeRGrid(const Interval<Dim> &gdom, const Loc<Dim> &blocks)
{
  Grid<Dim> ret;
  for (int i=0;i<Dim;++i)
    {

      if (gdom[i].size()%blocks[i].first()==0&& !gdom[i].empty())
	{
	  ret[i]=Grid<1>(Range<1>(gdom[i].first(),
				  gdom[i].last()+1, 
				  (int) gdom[i].size()/blocks[i].first()));
	}
      else if (!gdom[i].empty() )
	{
	  IndirectionList<int> rret(blocks[i].first()+1);
	  rret(0) = gdom[i].first();
	  for (int j = 1;j<blocks[i].first()+1;++j)
	    {
	      int temp = (int) gdom[i].size()/blocks[i].first();
	      rret(j)= rret(j-1) +  temp ;
	      rret(j) += (j > (blocks[i].first() - 
			       (gdom[i].size()- temp*blocks[i].first()) ) );
	    }
	  ret[i]=Grid<1>(rret);
	}
      else
	{
	  // do nothing.
	}
    }
  return ret;
}

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTITION_GRID_PARTITION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: GridPartition.h,v $   $Author: richard $
// $Revision: 1.32 $   $Date: 2004/11/01 18:17:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
