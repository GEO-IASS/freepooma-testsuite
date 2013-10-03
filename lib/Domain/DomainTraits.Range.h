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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_RANGE_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_RANGE_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<Range<N>>
// DomainChangeDim<Range<N>,Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<Range<N>> is a specialization of the general DomainTraits
 * class, for the case of Range domain objects.
 *
 * It defines the general
 * behavior of Range, including its typedef and static data
 * characteristics, how to store data for a Range, etc.  It is used by the
 * Domain base class of Range to implement most of the public interface.
 *
 * DomainTraits<Range<Dim>> stores the characteristics and much of the
 * implementation details for Range domain objects.  A Range represents
 * a sequence of numbers [a, a+s, a+2s, ... b], with a run-time stride s.
 *
 * A general version of DomainTraits<Range<Dim>> is defined here, which
 * only includes the basic information to make Range<Dim> look like an
 * array of Range<1> objects.  DomainTraits<Range<1>> is a more specific
 * specialization which provides most of the necessary interface information
 * for items which need to know about Range.  Since most of the interface
 * for a domain object is only available for 1D versions of that domain
 * object, the Range<1> specialization defines more interface functions than
 * the Range<Dim> case.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Utilities/UninitializedVector.h"


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
 * DomainTraits<Range<Dim>>:
 * The specialization of DomainTraits for Range, for dimensions greater than
 * one.
 */

template<int Dim>
struct DomainTraits< Range<Dim> >
  : public DomainTraitsDomain<Range<Dim>, int, Dim>
{
  // convenience typedef 
  typedef DomainTraitsDomain<Range<Dim>, int, Dim>     Base_t;

  // necessary typedefs
  typedef typename Base_t::Element_t                   Element_t;
  typedef typename Base_t::Domain_t                    Domain_t;
  typedef typename Base_t::NewDomain1_t                NewDomain1_t;
  typedef Range<1>      OneDomain_t;
  typedef Range<1>      PointDomain_t;
  typedef Interval<Dim> BlockDomain_t;
  typedef Loc<Dim>      AskDomain_t;
  typedef Range<Dim>    AddResult_t;
  typedef Range<Dim>    MultResult_t;

  // type for storage of this domain's data
  typedef UninitializedVector<OneDomain_t,Dim,Element_t> Storage_t;

  // necessary static data
  enum { domain = Base_t::domain };
  enum { dimensions = Base_t::dimensions,
	 sliceDimensions = Dim };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = false };
  enum { wildcard = false };

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  static OneDomain_t &getDomain(Domain_t &d, int n) { return d[n]; }
  static const OneDomain_t &getDomain(const Domain_t &d,int n) { return d[n]; }

  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a OneDomain_t,
  // since this is not a single-valued domain.
  static PointDomain_t &getPointDomain(Domain_t &d, int n) {
    return getDomain(d, n);
  }
  static const PointDomain_t &getPointDomain(const Domain_t &d, int n) {
    return getDomain(d, n);
  }
  
  // Domains get the chance to do special initialization.
  static void initializeStorage(Storage_t &dom) { dom.initialize(); }
};


/**
 * DomainTraits<Range<1>>:
 * The specialization of DomainTraits for Range, for dimension == 1.
 */

template<>
struct DomainTraits< Range<1> >
  : public DomainTraitsDomain<Range<1>, int, 1>
{
  // necessary typedefs
  typedef Range<1>    OneDomain_t;
  typedef Range<1>    PointDomain_t;
  typedef Interval<1> BlockDomain_t;
  typedef Loc<1>      AskDomain_t;
  typedef Range<1>    AddResult_t;
  typedef Range<1>    MultResult_t;

  // 1D necessary typedefs.  Range requires three pieces of
  // information, the begin point, the length, and the stride.
  // If length==0, this is empty.
  typedef Element_t Storage_t[3];

  // necessary static data
  enum { dimensions = 1,
         sliceDimensions = 1 };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = false };
  enum { wildcard = false };

  // return size, endpoint, stride, and loop information
  static Element_t first(const Storage_t &d)  { return d[0]; }
  static Element_t last(const Storage_t &d)   { return d[0] + (d[1]-1)*d[2]; }
  static Element_t stride(const Storage_t &d) { return d[2]; }
  static Element_t length(const Storage_t &d) { return d[1]; }
  static Element_t min(const Storage_t &d)    {
    return (d[2] > 0 ? d[0] : d[0] + (d[1]-1)*d[2]);
  }
  static Element_t max(const Storage_t &d)    {
    return (d[2] < 0 ? d[0] : d[0] + (d[1]-1)*d[2]);
  }
  static bool      empty(const Storage_t &d)  { return (d[1] < 1); }
  static int       loop(const Storage_t &)    { return 0; }

  // get the Nth value of the domain, where value # 0 is first(), etc.
  static Element_t elem(const Storage_t &d, int n) { return d[0] + n*d[2]; }

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  static OneDomain_t &getDomain(Domain_t &d, int) { return d; }
  static const OneDomain_t &getDomain(const Domain_t &d, int) { return d; }

  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a OneDomain_t,
  // since this is not a single-valued domain.
  static PointDomain_t &getPointDomain(Domain_t &d, int n) {
    return getDomain(d, n);
  }
  static const PointDomain_t &getPointDomain(const Domain_t &d, int n) {
    return getDomain(d, n);
  }
  
  // Domains get the chance to do special initialization.  Range's are
  // initialized to have length 0 and, just to avoid having a random value,
  // to start at 0 (although, for a length 0 domain, the endpoints are
  // actually undefined) and have stride = 1.
  static void initializeStorage(Storage_t &dom) {
    dom[0] = 0;			// first
    dom[1] = 0;			// length
    dom[2] = 1;			// stride
  }

  // change this domain object to the given one.  If things do not
  // match properly, assert a compile-time or run-time error.
  // For Range, we must have:
  // 1) the same dimensions==1
  template<class T>
  static void setDomain(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    dom[0] = DomainTraits<T>::getFirst(newdom);
    dom[1] = DomainTraits<T>::getLength(newdom);
    dom[2] = DomainTraits<T>::getStride(newdom);
  }

  // a specialized version of setDomain which accepts begin & end values.
  // The stride is set to + or - 1.
  template<class T1, class T2>
  static void setDomain(Storage_t &dom, const T1 &begval, const T2 &endval) {
    CTAssert(DomainTraits<T1>::dimensions == 1);
    CTAssert(DomainTraits<T2>::dimensions == 1);
    CTAssert(DomainTraits<T1>::singleValued);
    CTAssert(DomainTraits<T2>::singleValued);
    Element_t strideval = (endval < begval ? -1 : 1);
    dom[0] = begval;
    dom[1] = (endval - begval)/strideval + 1;
    dom[2] = strideval;
  }

  // a specialized version of setDomain which accepts begin & end values and
  // a stride.  For Range, we must have (endval - begval) % stride == 0, so
  // that the endpoints are consistent with the stride.
  // NOTE: The endpoint restriction has been removed; if the endpoint is
  // not consistent, it will be truncated.  To put this back to the
  // original method, uncomment the PAssert below.
  template<class T1, class T2, class T3>
  static void setDomain(Storage_t &dom, const T1 &begval, const T2 &endval,
			const T3 &strideval) {
    CTAssert(DomainTraits<T1>::dimensions == 1);
    CTAssert(DomainTraits<T2>::dimensions == 1);
    CTAssert(DomainTraits<T3>::dimensions == 1);
    CTAssert(DomainTraits<T1>::singleValued);
    CTAssert(DomainTraits<T2>::singleValued);
    CTAssert(DomainTraits<T3>::singleValued);
    //PAssert(strideval != 0 && ((endval - begval) % strideval)==0);
    dom[0] = begval;
    dom[1] = (endval - begval)/strideval + 1;
    dom[2] = strideval;
  }

  // change the loop variable for this object.  For Range, this is a no-op.
  static void setLoop(Storage_t &, int) { }

  // change the value of this 1D domain given a user-supplied reference
  // domain and a wildcard.
  template<class UT, class T>
  static void setWildcardDomain(Storage_t &dom, const UT &u, const T &newdom) {
    CTAssert(DomainTraits<T>::wildcard);
    CTAssert(DomainTraits<T>::dimensions == 1);
    CTAssert(DomainTraits<UT>::dimensions == 1);
    dom[0] = newdom.first(u);	// wildcard first method
    dom[1] = newdom.length(u);	// wildcard length method
    dom[2] = newdom.stride(u);	// wildcard stride method
  }

  //
  // compare this domain type to the given domain.  For the comparisons
  // to be meaningful for Range, we must have:
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
	       dom[2] < DomainTraits<T>::getStride(newdom)))));
  }

  // 'isEqualTo' returns true if dom == newdom
  template<class T>
  static bool isEqualTo(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    return ((dom[1] == 0 && DomainTraits<T>::getLength(newdom) == 0) ||
	    (dom[0] == DomainTraits<T>::getFirst(newdom) &&
	     dom[1] == DomainTraits<T>::getLength(newdom) &&
	     dom[2] == DomainTraits<T>::getStride(newdom)));
  }

  //
  // arithmetic accumulation operators.  These only work with
  // other domain objects with the following characteristics:
  // 1) they are singleValue'd
  // 2) they have dimensions == 1
  //

  // addAccum means dom[0] += newdom
  template<class T>
  static void addAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom[0] += DomainTraits<T>::getFirst(newdom);
  }

  // subtractAccum means dom[0] -= newdom
  template<class T>
  static void subtractAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom[0] -= DomainTraits<T>::getFirst(newdom);
  }

  // multiplyAccum means dom[0] *= newdom and dom[2] *= newdom
  template<class T>
  static void multiplyAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom[0] *= DomainTraits<T>::getFirst(newdom);
    dom[2] *= DomainTraits<T>::getFirst(newdom);
  }

  // divideAccum means dom[0] /= newdom and dom[2] /= newdom
  template<class T>
  static void divideAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom[0] /= DomainTraits<T>::getFirst(newdom);
    dom[2] /= DomainTraits<T>::getFirst(newdom);
  }
};


/**
 * DomainChangeDim<T, int> is used to convert from a domain of one dimension
 * to another dimension (the second template parameter).
 * For Range<Dim1>, it changes from Dim1 to Dim2.
 */

template<int Dim1, int Dim2>
struct DomainChangeDim<Range<Dim1>, Dim2>
{
  // the type of the old and new domain
  typedef Range<Dim1> OldType_t;
  typedef Range<Dim2> NewType_t;

  // static data for old and new dimensions
  enum { oldDim = Dim1,
	 newDim = Dim2 };
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_RANGE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.Range.h,v $   $Author: richard $
// $Revision: 1.27 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
