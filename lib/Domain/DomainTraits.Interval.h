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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_INTERVAL_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_INTERVAL_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<Interval<N>>
// DomainChangeDim<Interval<N>,Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<Interval<N>> is a specialization of the general DomainTraits
 * class, for the case of Interval domain objects.
 *
 * It defines the general
 * behavior of Interval, including its typedef and static data
 * characteristics, how to store data for a Interval, etc.  It is used by the
 * Domain base class of Interval to implement most of the public interface.
 *
 * DomainTraits<Interval<Dim>> stores the characteristics and much of the
 * implementation details for Interval domain objects.  An Interval represents
 * a sequence of numbers [a, a+1, ... b], with a hard-coded stride of +1.
 * Thus, it is unit-strided, but not necessarily single-valued.
 *
 * A general version of DomainTraits<Interval<Dim>> is defined here, which
 * only includes the basic information to make Interval<Dim> look like an
 * array of Interval<1> objects.  DomainTraits<Interval<1>> is a more specific
 * specialization which provides most of the necessary interface information
 * for items which need to know about Interval.  Since most of the interface
 * for a domain object is only available for 1D versions of that domain
 * object, the Interval<1> specialization defines more interface functions than
 * the Interval<Dim> case.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Utilities/UninitializedVector.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template <int Dim> class Loc;
template <> class Loc<1>;
template <int Dim> class Interval;
template <> class Interval<1>;
template <int Dim> class Range;
template <> class Range<1>;


/**
 * DomainTraits<Interval<Dim>>:
 * The specialization of DomainTraits for Interval, for dimensions greater than
 * one.
 */

template<int Dim>
struct DomainTraits< Interval<Dim> >
  : public DomainTraitsDomain<Interval<Dim>, int, Dim>
{
  // convenience typedef 
  typedef DomainTraitsDomain<Interval<Dim>, int, Dim>     Base_t;

  // necessary typedefs
  typedef typename Base_t::Element_t                      Element_t;
  typedef typename Base_t::Domain_t                       Domain_t;
  typedef typename Base_t::NewDomain1_t                   NewDomain1_t;
  typedef Interval<1>   OneDomain_t;
  typedef Interval<1>   PointDomain_t;
  typedef Interval<Dim> BlockDomain_t;
  typedef Loc<Dim>      AskDomain_t;
  typedef Interval<Dim> AddResult_t;
  typedef Range<Dim>    MultResult_t;

  // type for storage of this domain's data
  typedef UninitializedVector<OneDomain_t,Dim,Element_t> Storage_t;

  // necessary static data
  enum { domain = Base_t::domain };
  enum { dimensions = Base_t::dimensions,
	 sliceDimensions = Dim };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = true };
  enum { wildcard = false };

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  inline
  static OneDomain_t &getDomain(Domain_t &d, int n) { return d[n]; }
  inline
  static const OneDomain_t &getDomain(const Domain_t &d,int n) { return d[n]; }

  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a OneDomain_t,
  // since this is not a single-valued domain.
  inline
  static PointDomain_t &getPointDomain(Domain_t &d, int n) {
    return getDomain(d, n);
  }
  inline
  static const PointDomain_t &getPointDomain(const Domain_t &d, int n) {
    return getDomain(d, n);
  }
  
  // Domains get the chance to do special initialization.
  inline
  static void initializeStorage(Storage_t &dom) { dom.initialize(); }
};


/**
 * DomainTraits<Interval<1>>:
 * The specialization of DomainTraits for Interval, for dimension == 1.
 */

template<>
struct DomainTraits< Interval<1> >
  : public DomainTraitsDomain<Interval<1>, int, 1>
{
  // necessary typedefs
  typedef Interval<1> OneDomain_t;
  typedef Interval<1> PointDomain_t;
  typedef Interval<1> BlockDomain_t;
  typedef Loc<1>      AskDomain_t;
  typedef Interval<1> AddResult_t;
  typedef Range<1>    MultResult_t;

  // 1D necessary typedefs.  Interval requires two pieces of
  // information, the begin point and the length.  If length==0, this is
  // empty.  If the object is not empty, the stride is always 1,
  // and d[0] <= d[1].
  typedef Element_t Storage_t[2];

  // necessary static data
  enum { dimensions = 1,
         sliceDimensions = 1 };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = true };
  enum { wildcard = false };

  // return size, endpoint, stride, and loop information
  inline
  static Element_t first(const Storage_t &d)  { return d[0]; }
  inline
  static Element_t last(const Storage_t &d)   { return d[0] + d[1] - 1; }
  inline
  static Element_t stride(const Storage_t &)  { return 1; }
  inline
  static Element_t length(const Storage_t &d) { return d[1]; }
  inline
  static Element_t min(const Storage_t &d)    { return d[0]; }
  inline
  static Element_t max(const Storage_t &d)    { return d[0] + d[1] - 1; }
  inline
  static bool      empty(const Storage_t &d)  { return (d[1] < 1); }
  inline
  static int       loop(const Storage_t &)    { return 0; }

  // get the Nth value of the domain, where value # 0 is first(), etc.
  inline
  static Element_t elem(const Storage_t &d, int n) { return d[0] + n; }

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  inline
  static OneDomain_t &getDomain(Domain_t &d, int) { return d; }
  inline
  static const OneDomain_t &getDomain(const Domain_t &d, int) { return d; }

  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a OneDomain_t,
  // since this is not a single-valued domain.
  inline
  static PointDomain_t &getPointDomain(Domain_t &d, int n) {
    return getDomain(d, n);
  }
  inline
  static const PointDomain_t &getPointDomain(const Domain_t &d, int n) {
    return getDomain(d, n);
  }
  
  // Domains get the chance to do special initialization.
  // Domains get the chance to do special initialization.  Interval's are
  // initialized to have length 0 and, just to avoid having a random value,
  // to start at 0 (although, for a length 0 domain, the endpoints are
  // actually undefined).
  inline
  static void initializeStorage(Storage_t &dom) {
    dom[0] = 0;			// first
    dom[1] = 0;			// length
  }

  // change this domain object to the given one.  If things do not
  // match properly, assert a compile-time or run-time error.
  // For Interval, we must have:
  // 1) the same dimensions==1
  // 2) stride of newdom == 1
  template<class T>
  inline
  static void setDomain(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    PAssert(DomainTraits<T>::getStride(newdom) == 1);
    dom[0] = DomainTraits<T>::getFirst(newdom);
    dom[1] = DomainTraits<T>::getLength(newdom);
  }

  // a specialized version of setDomain which accepts begin & end values.
  // For Interval, we must have generally have begval <= endval, 
  // since the stride is hardcoded as +1. However, it seems overly restrictive
  // to disable the creation of zero-length Intervals. Hence, the slightly more
  // complicated PAssert.
  template<class T1, class T2>
  inline
  static void setDomain(Storage_t &dom, const T1 &begval, const T2 &endval) {
    CTAssert(DomainTraits<T1>::dimensions == 1);
    CTAssert(DomainTraits<T2>::dimensions == 1);
    CTAssert(DomainTraits<T1>::singleValued);
    CTAssert(DomainTraits<T2>::singleValued);
    dom[0] = begval;
    dom[1] = (endval - begval + 1);
    PAssert(begval <= endval || dom[1] == 0);
  }

  // change the loop variable for this object.  For Interval, this is a no-op.
  inline
  static void setLoop(Storage_t &, int) { }

  // change the value of this 1D domain given a user-supplied reference
  // domain and a wildcard.
  template<class UT, class T>
  inline
  static void setWildcardDomain(Storage_t &dom, const UT &u, const T &newdom) {
    CTAssert(DomainTraits<T>::wildcard);
    CTAssert(DomainTraits<T>::dimensions == 1);
    CTAssert(DomainTraits<UT>::dimensions == 1);
    dom[0] = newdom.first(u);	// uses wildcard version of first
    dom[1] = newdom.length(u);	// uses wildcard version of length
  }

  //
  // compare this domain type to the given domain.  For the comparisons
  // to be meaningful for Interval, we must have:
  // 1) the same dimensions==1
  // 2) both must not be empty
  //

  // 'isLessThan' returns true if dom < newdom
  template<class T>
  static bool isLessThan(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    PAssert(!(dom[1] < 1 || DomainTraits<T>::getEmpty(newdom)));
    return (dom[1] < DomainTraits<T>::getLength(newdom) ||
	    (dom[1] == DomainTraits<T>::getLength(newdom) &&
	     (dom[0] < DomainTraits<T>::getFirst(newdom) ||
	      (dom[0] == DomainTraits<T>::getFirst(newdom) &&
	       DomainTraits<T>::getStride(newdom) > 1))));
  }

  // 'isEqualTo' returns true if dom == newdom
  template<class T>
  static bool isEqualTo(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    return ((dom[1] == 0 && DomainTraits<T>::getLength(newdom) == 0) ||
	    (dom[0] == DomainTraits<T>::getFirst(newdom) &&
	     dom[1] == DomainTraits<T>::getLength(newdom) &&
	     DomainTraits<T>::getStride(newdom) == 1));
  }

  //
  // arithmetic accumulation operators.  These only work with
  // other domain objects with the following characteristics:
  // 1) they are singleValue'd
  // 2) they have dimensions == 1
  //
  // Note that for Intervals, we do NOT allow *= or /=.  You
  // must convert an Interval to a Range before doing multiplicative ops.
  //

  // addAccum means dom[0] += newdom
  template<class T>
  inline
  static void addAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom[0] += DomainTraits<T>::getFirst(newdom);
  }

  // subtractAccum means dom[0] -= newdom
  template<class T>
  inline
  static void subtractAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom[0] -= DomainTraits<T>::getFirst(newdom);
  }
};


/**
 * DomainChangeDim<T, int> is used to convert from a domain of one dimension
 * to another dimension (the second template parameter).
 * For Interval<Dim1>, it changes from Dim1 to Dim2.
 */

template<int Dim1, int Dim2>
struct DomainChangeDim<Interval<Dim1>, Dim2>
{
  // the type of the old and new domain
  typedef Interval<Dim1> OldType_t;
  typedef Interval<Dim2> NewType_t;

  // static data for old and new dimensions
  enum { oldDim = Dim1,
	 newDim = Dim2 };
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_INTERVAL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.Interval.h,v $   $Author: richard $
// $Revision: 1.30 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
