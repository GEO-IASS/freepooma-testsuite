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

#ifndef POOMA_PARTITION_UNIFORMGRIDPARTITION_H
#define POOMA_PARTITION_UNIFORMGRIDPARTITION_H

//-----------------------------------------------------------------------------
// Classes:
// UniformGridPartition<Dim>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Partition
 * @brief
 * A layout partitioner that will break a given
 * global domain into N equally-sized blocks, where the user specifies
 * how many subdivisions to make along each dimension (S_i).
 *
 * Thus, N = Prod(S_i).  The user must provide consistent information; if the
 * global domain given to this partitioner does not have the proper
 * size in each dimension to allow it to be divided by S_i evenly, it
 * will produce an assertion failure.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Layout/GuardLayers.h"
#include "Partition/ContextMapper.h"
#include "Partition/DistributedMapper.h"
#include "Utilities/PAssert.h"

#include <vector>
#include <list>

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * UniformGridPartition is a layout partitioner; it is created with
 * the information needed to partition a global domain into subdomains
 * using a grid layout.  All the blocks in a UniformGridLayout will
 * have the same size, and the domain that is given to a
 * UniformGridPartition object must be of the proper size to be
 * divided evenly into these sized blocks.
 *
 * A UniformGridPartition object is constructed with the following
 * information:
 *  - UniformGridPartition()
 *       Default constructor, makes a single block. 
 *  - UniformGridPartition(const GuardLayers<Dim> &guards)
 *       this is like the previous constructor, but takes a guard
 *       layer specification for external guards around the block.
 *  - UniformGridPartition(const Loc<Dim> &s) 
 *       the first argument is a Loc which contains the number of
 *       blocks to make along each dimension.  The total number of
 *       blocks is then the product of the Loc elements.
 *  - UniformGridPartition(const Loc<Dim> &s, 
 *                        const GuardLayers<Dim> &guards) 
 *  - UniformGridPartition(const Loc<Dim> &s, 
 *                        const GuardLayers<Dim> &internal_guards,
 *                        const GuardLayers<Dim> &external_guards) 
 *       These are like the second one, but specify guard layers
 *       for the partitioning. If only a single guard layer object is
 *       specified, then this specification is used for both the
 *       internal and external guard layers. These can be specified
 *       separately with the last constructor.
 *
 * UniformGridPartition objects, once created, are generally given to
 * Layout objects to initialize or repartition them.  The Layout, in
 * turn, will sets things up and call
 *
 *      UniformGridPartition::partition(domain, layoutData,mapper)
 *
 * where
 *  - domain     == global domain to partition
 *  - layoutData == the layout's data object that is being
 *                  partitioned.
 *  - mapper     == the ContextMapper to be used for this partitioning. 
 *
 * The partition method should divide up the global domain, and then
 * for each subdomain call 
 *
 *      layoutData.addDomain(domain[n], internal_guards, 
 *                           external_guards)
 * 
 * where domain[n] is the Nth subdomain, internal_guards and 
 * external_guards are the internal and external guard layer
 * specification for that particular domain (the partitioner must
 * determine if this domain is on an edge and assign external
 * vs. internal guards accordingly. Note that the partitioner partitions the
 * internal domain, not the total domain including the external
 * guards; i.e. the set of domain objects added to the layout is
 * independent of the external or internal guard specification. The
 * partitioner must know about the guards, however, since it needs to
 * check that the partitioning is consistent with the guards, and
 * since it is up to the partitioner to calculate the guards object
 * that is passed back to the layout so that IT can build Node objects
 * with the appropriate owned and allocated domains.
 *
 * If UniformGridPartition is asked to partition a global domain that
 * is empty, it will generate a set of N empty subdomains, with no
 * guard cells.  This is a special case, useful for the Particles
 * UniformLayout.
 *
 * UniformGridPartition also has accessor method to query how it will
 * divide the space: maxSize(), blocks(), along with a set
 * of accessors for examining guard layer data.
 */

template<int Dim>
class UniformGridPartition
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef LocalMapper<Dim>           DefaultMapper_t;
  typedef Interval<Dim>           Domain_t;
  typedef Node<Domain_t>          Value_t;
  typedef std::vector<Value_t *>  List_t;  
  // These traits are used by the layouts to determine if
  // a particular partitioner is appropriate for that layout.
  
  enum { uniform  = true  };
  enum { gridded  = true  };
  enum { tile     = false };
  enum { general  = false };

  enum { dimensions = Dim };

  //============================================================
  // Constructors
  //============================================================

  // The default constructor partitions the domain by just making one
  // single block.
  
  UniformGridPartition();

  // Ditto, but specifying external guards via a GuardLayers object:
  
  UniformGridPartition(const GuardLayers<Dim> &gcs);

  // All UniformGridPartition constructors take the following
  // arguments:
  //   a Loc<Dim> used to specify the number of blocks in each dimension.
  
  UniformGridPartition(const Loc<Dim> &a);

  // Ditto, but specifying guards (internal and external) via a 
  // GuardLayers object:
  
  UniformGridPartition(const Loc<Dim> &a, 
                       const GuardLayers<Dim> &gcs );

  // Ditto, with separate specifications of internal and external guards.
  
  UniformGridPartition(const Loc<Dim> &a, 
                       const GuardLayers<Dim> &igcs,
                       const GuardLayers<Dim> &egcs);

  // The copy constructor
  
  UniformGridPartition(const UniformGridPartition<Dim> &b);

  //============================================================
  // Destructor
  //============================================================

  // The destructor for this partitioner does not have to do anything.
  
  ~UniformGridPartition() { }


  //============================================================
  // Operators
  //============================================================

  // Assignment operator.
  
  UniformGridPartition<Dim> &
  operator=(const UniformGridPartition<Dim> &g) 
  {
    if (this != &g)
      {
	blocks_m              = g.blocks();
	hasGuards_m           = g.hasGuards_m;
	hasCustomEdgeGuards_m = g.hasCustomEdgeGuards_m;
	internalGuards_m      = g.internalGuards_m;
	externalGuards_m      = g.externalGuards_m;
	num_m                 = g.maxSize();
      }
    return *this;
  }
  

  //============================================================
  // Accessors
  //============================================================

  // Return the maximum number of subdomains this will generate.  The
  // actual number may end up being less, if the domain to be
  // partitioned is too small.
  
  int maxSize() const { return num_m; }

  // Return all the blocks for the dimensions, as a Loc.
  
  const Loc<Dim> &blocks() const { return blocks_m; }

  // Guard layer info:
  
  bool hasGuards() const 
  { 
    PAssert(hasGuards_m == (hasInternalGuards() || hasExternalGuards()));
    return hasGuards_m; 
  }

  bool hasInternalGuards() const 
  { 
    return hasGuards_m && internalGuards_m != 0;
  }

  bool hasExternalGuards() const 
  { 
    return hasGuards_m && externalGuards_m != 0; 
  }

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

  // For the given global domain, partition it into subdomains and put
  // the results in the provided layout object by calling
  // 'layoutData.addDomainList(List_t &templist)'.  Return the
  // total number of subdomains added.
  
  template<class D>
  int partition(const D &domain, 
		List_t & all,
		const ContextMapper<Dim>& cmapper) const;

  template<class D>
  int partition(const D &domain, List_t & list) const 
  {
    return partition(domain,list,DefaultMapper_t(*this));
  }

protected:

  // The number of blocks along each dimension.
  
  Loc<Dim> blocks_m;
  
  // Do we have guard layers?
  
  bool hasGuards_m;
  
  // Are the external guards different from the internal?
  
  bool hasCustomEdgeGuards_m;
  
  // Specification of internal guard layers.
  
  GuardLayers<Dim> internalGuards_m;
  
  // Specification of external guard layers.
  
  GuardLayers<Dim> externalGuards_m;

  // The total number of blocks to create.
  
  int num_m;
  
  // Calculate num_m from blocks_m:
  
  void calcNum()
  {
    num_m = blocks_m[0].first();
    for (int d = 1; d < Dim; ++d) 
      {
	num_m *= blocks_m[d].first();
      }
  }
    
};


//============================================================
// UniformGridPartition inline method definitions
//============================================================

template<int Dim>
template<class D>
int UniformGridPartition<Dim>::partition(const D &domain, 
					 List_t & all,
					 const ContextMapper<Dim>& cmapper) const 
{
  // The type info for domain we should be creating for the layout.

  typedef typename DomainTraits<Domain_t>::Element_t Element_t;

  // Make sure we have the right dimensionality.

  CTAssert(Dim == DomainTraits<D>::dimensions);
  CTAssert(Dim == DomainTraits<Domain_t>::dimensions);

  // This will only work with UnitStride domains

  CTAssert(DomainTraits<D>::unitStride == 1);
  CTAssert(DomainTraits<Domain_t>::unitStride == 1);

  // make sure the list is empty

  PAssert(all.size() == 0);

  // Cache the origin of the domain and make sure the domain is
  // properly sized. Also, build a domain corresponding to the
  // number of blocks in each direction for iterating over below.

  Element_t origin[Dim];
  Element_t sizes[Dim];
  Interval<Dim> bdomain = Pooma::NoInit(); // dummy initializer

  int i;

  for (i = 0; i < Dim; ++i) 
    {
      if (!domain.empty())
	{
	  int gcwidth = 
	    (internalGuards_m.lower(i) > internalGuards_m.upper(i)) ?
	    internalGuards_m.lower(i) : internalGuards_m.upper(i);

	  PInsist((domain[i].length() % blocks()[i].first()) == 0,
		  "All the blocks in a grid must be the same size.");

	  origin[i]  = domain[i].first();
	  sizes[i]   = domain[i].length() / blocks()[i].first();

	  PInsist(sizes[i] >= gcwidth,
		  "Block sizes too small for guard layer specification.");
	}
      bdomain[i] = Interval<1>(blocks()[i].first());
    }

  // Loop over all the blocks, creating new domains. 

  typename Interval<Dim>::const_iterator it = bdomain.begin();
  while (it != bdomain.end()) 
    {
      // Start with an initially empty domain and empty guard cells.

      Domain_t owned;
      GuardLayers<Dim> iguards(0);
      GuardLayers<Dim> eguards(0);

      // Calculate the subdomain, if the global domain is not empty.
      // If it is, we just use the empty domain.

      if (!domain.empty())
	{
	  Loc<Dim> pos = *it;
	  for (i = 0; i < Dim; ++i) 
	    {
	      int position = pos[i].first();
	      Element_t a = origin[i] + sizes[i]*position;
	      Element_t b = a + sizes[i] - 1;
	      typedef typename 
		DomainTraits<Domain_t>::OneDomain_t OneDomain_t;
	      owned[i] = OneDomain_t(a, b);
	    }

	  // Calculate the internal and external guard layer specifications
	  // for this domain.

	  if (hasGuards_m)
	    {
	      iguards = internalGuards_m;

	      // Check if we're at an edge, and, if so, use the
	      // external specification for that edge.

	      for (int d = 0; d < Dim; ++d)
		{
		  int position = pos[d].first();
		  if ( position == bdomain[d].first() ) 
		    {
		      eguards.lower(d) = externalGuards_m.lower(d);
		      iguards.lower(d) = 0;
		    }
		  if ( position == bdomain[d].last() ) 
		    {
		      eguards.upper(d) = externalGuards_m.upper(d);
		      iguards.upper(d) = 0;
		    }
		}
	    }
	}
      typename Value_t::ID_t gid = all.size();
      typename Value_t::ID_t lid = (-1);

      // Create a new Node object to store the subdomain data.

      GuardLayers<Dim>::addGuardLayers(owned,eguards);

      Domain_t allocated = owned;

      GuardLayers<Dim>::addGuardLayers(allocated,iguards);

      Value_t *node = new Value_t(owned, allocated, -1, gid, lid);

      all.push_back(node);

      // Increment our counters and iterators.

      ++it;
    }

  cmapper.map(all);

  // At the end, return # of domains created.

  return num_m;
}


//-----------------------------------------------------------------------------
//
// Constructors & assignment operator.
//
//-----------------------------------------------------------------------------

template <int Dim>
inline UniformGridPartition<Dim>::
UniformGridPartition()
: hasGuards_m(false),
  hasCustomEdgeGuards_m(false),
  num_m(1)
{
  blocks_m = 1;
}

template <int Dim>
inline UniformGridPartition<Dim>::
UniformGridPartition(const GuardLayers<Dim> &gcs)
: hasGuards_m(gcs != 0), 
  hasCustomEdgeGuards_m(gcs != 0), 
  externalGuards_m(gcs),
  num_m(1) 
{
  blocks_m = 1;
}

template <int Dim>
inline UniformGridPartition<Dim>::
UniformGridPartition(const Loc<Dim> &a)
: blocks_m(a),
  hasGuards_m(false), 
  hasCustomEdgeGuards_m(false)
{
  calcNum();
}
  
template <int Dim>
inline UniformGridPartition<Dim>::
UniformGridPartition(const Loc<Dim> &a, 
                     const GuardLayers<Dim> &gcs)
: blocks_m(a),
  hasGuards_m(gcs != 0), 
  hasCustomEdgeGuards_m(false), 
  internalGuards_m(gcs), 
  externalGuards_m(gcs)
{
  calcNum();
}

template <int Dim>
inline UniformGridPartition<Dim>::
UniformGridPartition(const Loc<Dim> &a, 
                     const GuardLayers<Dim> &igcs,
                     const GuardLayers<Dim> &egcs)
: blocks_m(a),
  hasGuards_m(igcs != 0 || egcs != 0), 
  hasCustomEdgeGuards_m(igcs != egcs), 
  internalGuards_m(igcs), 
  externalGuards_m(egcs)
{
  calcNum();
}

template <int Dim>
inline UniformGridPartition<Dim>::
UniformGridPartition(const UniformGridPartition<Dim> &b)
: blocks_m(b.blocks_m), 
  hasGuards_m(b.hasGuards_m), 
  hasCustomEdgeGuards_m(b.hasCustomEdgeGuards_m), 
  internalGuards_m(b.internalGuards_m), 
  externalGuards_m(b.externalGuards_m),
  num_m(b.num_m) 
{ }

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTITION_UNIFORMGRIDPARTITION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformGridPartition.h,v $   $Author: richard $
// $Revision: 1.31 $   $Date: 2004/11/01 18:17:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
