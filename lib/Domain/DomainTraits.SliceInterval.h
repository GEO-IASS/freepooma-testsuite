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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_SLICE_INTERVAL_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_SLICE_INTERVAL_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<SliceInterval<Dim,SliceDim>>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<SliceInterval<Dim,SliceDim>> is a specialization of the
 * general DomainTraits class, for the case of SliceInterval domain objects.
 *
 * It defines the general behavior of SliceInterval, including its typedef
 * and static data characteristics, how to store data for a SliceInterval, etc.
 * It is used by the SliceDomain base class of SliceInterval to implement most
 * of the public interface.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Domain/Interval.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int TotalDim, int SliceDim> class SliceInterval;


/**
 * DomainTraits<SliceInterval<Dim,SliceDim>> stores the characteristics and
 * much of the implementation details for SliceInterval domain objects.
 * A SliceInterval represents a set of two domain objects, one a "total"
 * domain, and the other a "slice" domain which is a subset of the total.
 * SliceInterval stores these two domains as Interval<> objects.
 *
 * The DomainTraits for slice domains is quite a bit simpler than the
 * DomainTraits for regular domains.  This is because SliceDomains have
 * a much simpler interface than regular domains, and are not intended for
 * direct user manipulation.  The DomainTraits for SliceDomain subclasses
 * like SliceInterval includes the following interface:
 *
 * static const int dimensions = # of total dimensions
 * static const int sliceDimensions = # of slice dimensions
 *
 * typedef SliceInterval<TotalDim,SliceDim> Domain_t;
 * typedef Interval<SliceDim>               SliceDomain_t;
 * typedef Interval<TotalDim>               TotalDomain_t;
 * typedef Interval<1>                      OneDomain_t;
 *
 * static OneDomain_t &getDomain(Domain_t &d, int n);
 * static OneDomain_t &getSliceDomain(Domain_t &d, int n);
 */

template<int TotalDim, int SliceDim>
struct DomainTraits< SliceInterval<TotalDim,SliceDim> >
{
  // necessary static data
  enum { domain          = true };
  enum { dimensions      = TotalDim,
	 sliceDimensions = SliceDim };
  enum { unitStride      = true };
  enum { singleValued    = false };
  enum { wildcard        = false };

  // necessary typedefs
  typedef SliceInterval<TotalDim,SliceDim> Domain_t;
  typedef SliceInterval<TotalDim,SliceDim> NewDomain1_t;
  typedef Interval<SliceDim>               SliceDomain_t;
  typedef Interval<TotalDim>               TotalDomain_t;
  typedef Interval<1>                      OneDomain_t;
  typedef Interval<1>                      PointDomain_t;

  // get the Nth element of the total domain, and return a OneDomain_t
  // object with it.
  static OneDomain_t &getDomain(Domain_t &d, int n) {
    return d.totalDomain()[n];
  }
  static const OneDomain_t &getDomain(const Domain_t &d,int n) {
    return d.totalDomain()[n];
  }

  // get the Nth element of the sliced domain, and return a OneDomain_t
  // object with it
  static OneDomain_t &getSliceDomain(Domain_t &d, int n) {
    return d.sliceDomain()[n];
  }
  static const OneDomain_t &getSliceDomain(const Domain_t &d, int n) {
    return d.sliceDomain()[n];
  }
  
  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a OneDomain_t,
  // since this is not a single-valued domain.
  static PointDomain_t &getPointDomain(Domain_t &d, int n) {
    return getDomain(d, n);
  }
  static const PointDomain_t &getPointDomain(const Domain_t &d, int n) {
    return getDomain(d, n);
  }

  // set the given dimension as ignorable
  static void cantIgnoreDomain(Domain_t &d, int n) {
    d.cantIgnoreDomain(n);
  }

  // get the ignore status of the given dimension
  static bool getIgnorable(const Domain_t &d, int n) {
    return d.ignorable(n);
  }

  // set the ignore status of the given dimension
  static void setIgnorable(Domain_t &d, int n, bool i) {
    d.ignorable(n) = i;
  }
  
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_SLICE_INTERVAL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.SliceInterval.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
