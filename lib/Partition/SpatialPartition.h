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

#ifndef POOMA_PARTITION_SPATIAL_PARTITION_H
#define POOMA_PARTITION_SPATIAL_PARTITION_H

//-----------------------------------------------------------------------------
// Classes:
// SpatialPartition<Layout>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Partition
 * @brief
 * A layout partitioner that will generate a set of
 * initially empty domains for insertion into a layout, based on the
 * information from another layout.
 *
 * The generated domains will all be
 * Dim-dimensional, regardless of the dimensionality of the given layout,
 * where Dim is the dimensionality of a layout object for which this object
 * is requested to generate new patches.
 * The same number of patches will be generated on each context as there
 * are in the reference layout, with the same memory affinity.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include "Layout/Node.h"
#include "Layout/GuardLayers.h"
#include "Utilities/PAssert.h"
#include <iosfwd>



///////////////////////////////////////////////////////////////////////////
// namespace POOMA {
class DefaultSPmapper
  : public ContextMapper<1>
{ 
  // this mapper class is a 'no op', as the partitioner that uses it
  // SpatialPartition is so specialized that it can and does do both
  // the partitioning and mapping in one shot. 

public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef Interval<1>                         Domain_t;
  typedef Node<Domain_t>                      Value_t;
  typedef std::vector<Value_t *>              List_t;

  template<class ReferenceLayout>
  DefaultSPmapper(const ReferenceLayout & L)
  {
  }

  DefaultSPmapper()
  {
  }
 
  inline void map(const List_t &) const
  {
  }

private:
};


/**
 * SpatialPartition is a layout partitioner; it is created with the info
 * needed to partition a global domain into subdomains based on information
 * in a second "reference" layout.
 *
 * This partitioner will generate initially empty patches, and so is only
 * useful if it is asked to partition an empty global domain.  It is primarily
 * used in situations where you need a second layout to have the same number
 * of patches as a reference domain, with the same memory affinity for each
 * patch.  The size of the new patches, and their dimension, can be different
 * than for the reference layout.
 *
 * The template parameter of SpatialLayout is the type of reference
 * layout; an instance of an object of this type must be provided in the
 * constructor, for later use in the "partition" method.  This can be any
 * layout-like type that follows the normal layout interface rules.
 *
 * The main interface to SpatialPartition is via its "partition" method,
 * where it is provided the layout object that will get new domains after
 * they are created and assigned to contexts, plus the global domain to
 * partition.  The global domain to partition in this case should be an
 * empty domain, if it is not it is an error.  For this partitioner,
 * the generated domains will have the same
 * dimensionality as the provided layout (but not necessarily the same
 * dimensionality as the reference layout).  "partition" will go through
 * the patches of the reference layout, and generate the same number of
 * patches for the layout being partitioned.  The generated patches will
 * have empty domains and no guard cells.  They will have the same memory
 * affinity as the corresponding patches in the reference domain.
 */

template<class ReferenceLayout>
class SpatialPartition
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef DefaultSPmapper           DefaultMapper_t;
  typedef ReferenceLayout           Layout_t;
  typedef Interval<1>               Domain_t;
  typedef Node<Domain_t>            Value_t;
  typedef std::vector<Value_t *>    List_t;

  // A useful typedef to refer to the type of this object, and other
  // typedefs.

  typedef SpatialPartition<ReferenceLayout>         This_t;
  typedef GuardLayers<1>                            GuardLayers_t;

  // A typedef for the type of reference layout

  typedef ReferenceLayout                           RefLayout_t;

  // Enumerated characteristics of this partitioner.

  enum { uniform = false };
  enum { gridded = true  };
  enum { tile    = false };
  enum { general = false };

  enum { dimensions = 1 };


  //============================================================
  // Constructors
  //============================================================

  // The main way to construct this partitioner is with a reference
  // layout object; this is used to generate patches in the partition
  // method.

  SpatialPartition(const RefLayout_t &layout)
    : reference_m(layout)
    {
    }

  // Copy constructor

  SpatialPartition(const This_t &sp)
    : reference_m(sp.reference_m)
    {
    }

  // Assignment operator.

  This_t &operator=(const This_t &sp)
    {
      if (this != &sp)
	reference_m = sp.reference_m;

      return *this;
    }


  //============================================================
  // Destructor
  //============================================================

  // The destructor for this partitioner does not have to do anything.
  
  ~SpatialPartition() { }


  //============================================================
  // Accessors
  //============================================================

  // Return the maximum number of subdomains this will generate.  This
  // is the same as the global number of patches in the reference layout.

  int maxSize() const
    {
      return reference().sizeGlobal();
    }

  // Return the number of blocks that will be generated in each dimension.
  // This partitioner will generate a single list of blocks, so this
  // will return a 1D Loc set equal to the number of blocks in the
  // reference layout.

  Loc<1> blocks() const
    {
      return Loc<1>(maxSize());
    }

  // Return which context these subdomains will be assigned to, or a
  // number < 0 if they will be assigned to all contexts.
  
  int context() const
    {
      return (-1);
    }

  // Guard layer info

  bool hasGuards() const
    {
      return false;
    }

  bool hasCustomEdgeGuards() const 
    {
      return false;
    }

  bool hasInternalGuards() const
    {
      return false;
    }

  bool hasExternalGuards() const
    {
      return false;
    }

  inline GuardLayers_t internalGuards() const
    {
      return GuardLayers_t();
    }

  inline GuardLayers_t externalGuards() const
    {
      return GuardLayers_t();
    }

  // Our reference layout.

  const RefLayout_t &reference() const
    {
      return reference_m;
    }


  //============================================================
  // Partition methods
  //============================================================

  template<class D,class Dom>
  int partition(const D &domain, 
		std::vector< Node<Dom> *> & all, 
		const ContextMapper<1> & cmapper ) const;

  template<class D,class Dom> 
  int partition(const D & domain, 
	        std::vector< Node<Dom> * > & all ) const
    {
      return partition(domain,all,DefaultSPmapper(reference_m));
    }

  //============================================================
  // I/O/
  //============================================================

  template<class Out>
  void print(Out &o) const;

private:
  // The reference layout

  RefLayout_t reference_m;

  // Make default constructor private, since we don't need to use it.

  SpatialPartition();
};

template<class ReferenceLayout>
template<class D,class Dom>
int SpatialPartition<ReferenceLayout>::partition(const D &domain,
						 std::vector< Node<Dom> *> & all, 
						 const ContextMapper<1> & cmapper) const
{
  // Make sure we have the right dimensionality between the provided
  // domain and the layout's domain.  We do NOT need to have the
  // dimensionality match between the provided layout and reference
  // layout, however.

  CTAssert(DomainTraits<D>::dimensions ==
	   DomainTraits<Dom>::dimensions);

  // For now, this will only work with a 1D provided layout.

  CTAssert(DomainTraits<D>::dimensions == 1);

  // The provided domain must actually be empty for this to work,
  // since we generate empty domains.

  PAssert(domain.empty());

  // Loop through the patches in the reference domain now.  For each one,
  // just add in an empty domain, assigned to the same context
  // as the current node that we're iterating over.

  typename RefLayout_t::const_iterator refpatch, 
    endref = reference().endGlobal();

  for (refpatch=reference().beginGlobal(); refpatch != endref; ++refpatch)
    { 
      Node<Dom> *node = new Node<Dom>(Dom(), 
				      Dom(),
				      refpatch->context(),
				      refpatch->globalID(),
				      refpatch->localID());
      all.push_back(node);
    }

  // This call to the mapper for the case of DefaultSPMapper 
  // is a no-op, and could be commented out. The above loop
  // does both the partitioning and the mapping. The structure
  // of the DefaultSPmapper is maintained for this special case
  // even though it's redundant.

  cmapper.map(all);

  // Return the number of domains we added in.

  return maxSize();
}

template<class ReferenceLayout>
template<class Out>
void SpatialPartition<ReferenceLayout>::print(Out &o) const
{
  o << "SpatialPartitioner:\n";
  o << "  reference layout = " << reference() << "\n";
  o << "  maximum patches = " << maxSize() << "\n";
}


/// A specialization of the Inform traits used to say that node has
/// a print method.

template <class L>
std::ostream &operator<<(std::ostream &o, const SpatialPartition<L> &sp)
{
  sp.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTITION_SPATIAL_PARTITION_H
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SpatialPartition.h,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:17:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo

