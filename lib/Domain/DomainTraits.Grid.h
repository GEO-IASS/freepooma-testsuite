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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_GRID_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_GRID_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<Grid<N>>
// DomainChangeDim<Grid<N>,Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<Grid<N>> is a specialization of the general DomainTraits
 * class, for the case of Grid domain objects.
 *
 * It defines the general
 * behavior of Grid, including its typedef and static data
 * characteristics, how to store data for a Grid, etc.  It is used by the
 * Domain base class of Grid to implement most of the public interface.
 *
 *
 * DomainTraits<Grid<Dim>> stores the characteristics and much of the
 * implementation details for Grid domain objects.  A Grid represents
 * a sequence of numbers [a0, a1, ... aN] for each dimension; the numbers
 * a0 ... aN can be any list, as long as they are sorted in ascending
 * or descending order.  Data is stored internally for each dimension
 * using an IndirectionList<int>, and the total domain is the tensor
 * product of the 1D lists.
 *
 * A general version of DomainTraits<Grid<Dim>> is defined here, which
 * only includes the basic information to make Grid<Dim> look like an
 * array of Grid<1> objects.  DomainTraits<Grid<1>> is a more specific
 * specialization which provides most of the necessary interface information
 * for items which need to know about Grid.  Since most of the interface
 * for a domain object is only available for 1D versions of that domain
 * object, the Grid<1> specialization defines more interface functions than
 * the Grid<Dim> case.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Domain/IndirectionList.h"
#include "Utilities/UninitializedVector.h"



//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


template <int Dim> class Loc;
template <> class Loc<1>;
template <int Dim> class Interval;
template <> class Interval<1>;
template <int Dim> class Grid;
template <> class Grid<1>;


/**
 * DomainTraits<Grid<Dim>>:
 * The specialization of DomainTraits for Grid, for dimensions greater than
 * one.
 */

template<int Dim>
struct DomainTraits< Grid<Dim> >
  : public DomainTraitsDomain<Grid<Dim>, int, Dim>
{
  // convenience typedef 
  typedef DomainTraitsDomain<Grid<Dim>, int, Dim>     Base_t;

  // necessary typedefs
  typedef typename Base_t::Element_t                  Element_t;
  typedef typename Base_t::Domain_t                   Domain_t;
  typedef typename Base_t::NewDomain1_t               NewDomain1_t;
  typedef Grid<1>       OneDomain_t;
  typedef Grid<1>       PointDomain_t;
  typedef Interval<Dim> BlockDomain_t;
  typedef Loc<Dim>      AskDomain_t;
  typedef Grid<Dim>     AddResult_t;
  typedef Grid<Dim>     MultResult_t;

  // type for storage of this domain's data
  typedef void *        ElemPtr_t;
  typedef UninitializedVector<OneDomain_t,Dim,ElemPtr_t> Storage_t;

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
 * DomainTraits<Grid<1>>:
 * The specialization of DomainTraits for Grid, for dimension == 1.
 */

template<>
struct DomainTraits< Grid<1> >
  : public DomainTraitsDomain<Grid<1>, int, 1>
{
  // necessary typedefs
  typedef Grid<1>     OneDomain_t;
  typedef Grid<1>     PointDomain_t;
  typedef Interval<1> BlockDomain_t;
  typedef Loc<1>      AskDomain_t;
  typedef Grid<1>     AddResult_t;
  typedef Grid<1>     MultResult_t;

  // 1D necessary typedefs.  Grid stores data in an IndirectionList<Element_t>.
  typedef IndirectionList<Element_t> Storage_t;

  // necessary static data
  enum { dimensions = 1,
         sliceDimensions = 1 };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = false };
  enum { wildcard = false };

  // return size, endpoint, stride, and loop information
  static Element_t first(const Storage_t &d)  { return d.first(); }
  static Element_t last(const Storage_t &d)   { return d.last(); }
  static Element_t stride(const Storage_t &d) { return d.stride(); }
  static Element_t length(const Storage_t &d) { return d.length(); }
  static Element_t min(const Storage_t &d)    { return d.min(); }
  static Element_t max(const Storage_t &d)    { return d.max(); }
  static bool      empty(const Storage_t &d)  { return d.empty(); }
  static int       loop(const Storage_t &)    { return 0; }

  // get the Nth value of the domain, where value # 0 is first(), etc.
  static Element_t elem(const Storage_t &d, int n) { return d(n); }

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
  
  // Domains get the chance to do special initialization.  Grid's
  // start out with an empty domain already, though, and don't need
  // to do anything extra.
  static void initializeStorage(Storage_t &) { }

  // change this domain object to the given one.  If things do not
  // match properly, assert a compile-time or run-time error.
  // For Grid, we must have:
  // 1) the same dimensions==1
  template<class T>
  static void setDomain(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    dom = Storage_t(DomainTraits<T>::getFirst(newdom),
		    DomainTraits<T>::getStride(newdom),
		    DomainTraits<T>::getLength(newdom));
  }

  // change this domain object to the given Grid<>.  This is
  // a special version, since we dont' want to use first/length/stride
  // queries, we want to copy over the IL contents directly.
  template<int Dim>
  static void setDomain(Storage_t &dom, const Grid<Dim> &newdom) {
    CTAssert(Dim == 1);
    dom = newdom.storage();
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
    dom = Storage_t(begval, strideval, (endval - begval)/strideval + 1);
  }

  // a specialized version of setDomain which accepts begin & end values and
  // a stride.  For Grid, we must have (endval - begval) % stride == 0, so
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
    dom = Storage_t(begval, strideval, (endval - begval)/strideval + 1);
  }

  // change this domain object to the given IndirectionList.  This is
  // a special version, since we dont' want to use first/length/stride
  // queries, we want to copy over the IL contents directly.
  static void setDomain(Storage_t &dom, const Storage_t &newdom) {
    dom = newdom;
  }

  // change the loop variable for this object.  For Grid, this is a no-op.
  static void setLoop(Storage_t &, int) { }

  // change the value of this 1D domain given a user-supplied reference
  // domain and a wildcard.
  template<class UT, class T>
  static void setWildcardDomain(Storage_t &dom, const UT &u, const T &newdom) {
    CTAssert(DomainTraits<T>::wildcard);
    CTAssert(DomainTraits<T>::dimensions == 1);
    CTAssert(DomainTraits<UT>::dimensions == 1);
    dom = Storage_t(newdom.first(u), newdom.stride(u), newdom.length(u));
  }

  //
  // compare this domain type to the given domain.  For the comparisons
  // to be meaningful for Grid, we must have:
  // 1) the same dimensions==1
  // 2) both must not be empty
  //

  // 'isLessThan' returns true if dom < newdom
  template<class T>
  static bool isLessThan(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    PAssert(!(dom.empty() || DomainTraits<T>::getEmpty(newdom)));
    return (dom.first() < DomainTraits<T>::getFirst(newdom) ||
	    (dom.first() == DomainTraits<T>::getFirst(newdom) &&
	     (dom.last() < DomainTraits<T>::getLast(newdom) ||
	      (dom.last() == DomainTraits<T>::getLast(newdom) &&
	       dom.length() < DomainTraits<T>::getLength(newdom)))));
  }

  // 'isEqualTo' returns true if dom == newdom
  template<class T>
  static bool isEqualTo(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    return ((dom.empty() && DomainTraits<T>::getEmpty(newdom)) ||
	    (dom.first() == DomainTraits<T>::getFirst(newdom) &&
	     dom.last() == DomainTraits<T>::getLast(newdom) &&
	     dom.length() == DomainTraits<T>::getLength(newdom)));
  }

  //
  // arithmetic accumulation operators.  These only work with
  // other domain objects with the following characteristics:
  // 1) they are singleValue'd
  // 2) they have dimensions == 1
  //

  // addAccum means add newdom to all elements
  template<class T>
  static void addAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom += DomainTraits<T>::getFirst(newdom);
  }

  // subtractAccum means subtract newdom from all elements
  template<class T>
  static void subtractAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom -= DomainTraits<T>::getFirst(newdom);
  }

  // multiplyAccum means multiply all elements by newdom
  template<class T>
  static void multiplyAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom *= DomainTraits<T>::getFirst(newdom);
  }

  // divideAccum means divide all elements by newdom
  template<class T>
  static void divideAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom /= DomainTraits<T>::getFirst(newdom);
  }
};


/**
 * DomainChangeDim<T, int> is used to convert from a domain of one dimension
 * to another dimension (the second template parameter).
 * For Grid<Dim1>, it changes from Dim1 to Dim2.
 */

template<int Dim1, int Dim2>
struct DomainChangeDim<Grid<Dim1>, Dim2>
{
  // the type of the old and new domain
  typedef Grid<Dim1> OldType_t;
  typedef Grid<Dim2> NewType_t;

  // static data for old and new dimensions
  enum { oldDim = Dim1,
	 newDim = Dim2 };
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_GRID_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.Grid.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
