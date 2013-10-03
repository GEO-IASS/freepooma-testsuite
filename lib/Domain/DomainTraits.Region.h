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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_REGION_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_REGION_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<Region<N,T>> : public DomainTraitsRegion<N,T>
// DomainChangeDim<Region<N,T>,Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<Region<N,T>> is a specialization of the general DomainTraits
 * class, for the case of Region domain objects.
 *
 * It defines the general
 * behavior of Region, including its typedef and static data
 * characteristics, how to store data for a Region, etc.  It is used by the
 * Domain base class of Region to implement most of the public interface.
 *
 * DomainTraits<Region<Dim,T>> stores the characteristics and much of the
 * implementation details for Region domain objects.  A Region represents
 * a continuous region of values in an N-dimensional space, by storing
 * the endpoints in each dimension which define an N-dimensional rectangle.
 * There is no stride associated with a Region; when asked, it reports a
 * stride equal to the width of the 1D area.
 *
 * A general version of DomainTraitsRegion<Dim,T> is defined which
 * only includes the basic information to make Region<Dim> look like an
 * array of Region<1> objects.  DomainTraitsRegion<1,T> is a more specific
 * specialization which provides most of the necessary interface information
 * for items which need to know about Region.  Since most of the interface
 * for a domain object is only available for 1D versions of that domain
 * object, the Region<1,T> specialization defines more interface functions than
 * the Region<Dim,T> case.
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

template <int Dim, class T> class Region;
template <class T> class Region<1,T>;
template <> class Region<1,POOMA_DEFAULT_POSITION_TYPE>;


/**
 * DomainTraits<Region<Dim,T>>:
 * The traits for an N-dimensional Region domain.
 */

template<int Dim, class T>
struct DomainTraits< Region<Dim,T> >
  : public DomainTraitsDomain<Region<Dim,T>, T, Dim>
{
  // convenience typedef 
  typedef DomainTraitsDomain<Region<Dim,T>, T, Dim>     Base_t;

  // necessary typedefs
  typedef typename Base_t::Element_t                    Element_t;
  typedef typename Base_t::Domain_t                     Domain_t;
  typedef typename Base_t::NewDomain1_t                 NewDomain1_t;
  typedef Region<1,T>   OneDomain_t;
  typedef Region<1,T>   PointDomain_t;
  typedef Region<Dim,T> BlockDomain_t;
  typedef Region<Dim,T> AskDomain_t;
  typedef Region<Dim,T> AddResult_t;
  typedef Region<Dim,T> MultResult_t;

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
 * DomainTraits<Region<Dim,T>>:
 * The traits for an 1-dimensional Region domain.
 */

template<class T>
struct DomainTraits< Region<1,T> >
  : public DomainTraitsDomain<Region<1,T>, T, 1>
{
  // convenience typedef 
  typedef DomainTraitsDomain<Region<1,T>, T, 1>         Base_t;

  // necessary typedefs
  typedef typename Base_t::Element_t                    Element_t;
  typedef typename Base_t::Domain_t                     Domain_t;
  // necessary typedefs
  typedef Region<1,T> OneDomain_t;
  typedef Region<1,T> BlockDomain_t;
  typedef Region<1,T> AskDomain_t;
  typedef Region<1,T> AddResult_t;
  typedef Region<1,T> MultResult_t;

  // 1D necessary typedefs.  Region requires two pieces of
  // information, the begin point and the length.  If length==0, this is
  // just a point.  For the iterator, we need to know the current position
  // and the stride (which will be the width of the domain).
  typedef Element_t Storage_t[2];
  typedef Element_t IteratorStorage_t[2];

  // necessary static data
  enum { domain = Base_t::domain };
  enum { dimensions = Base_t::dimensions,
	 sliceDimensions = 1 };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = true };
  enum { wildcard = false };

  // return size, endpoint, stride, and loop information given the storage
  static Element_t first(const Storage_t &d)  { return d[0]; }
  static Element_t last(const Storage_t &d)   { return d[0] + d[1]; }
  static Element_t stride(const Storage_t &d) { return d[1]; }
  static Element_t length(const Storage_t &d) { return d[1]; }
  static Element_t min(const Storage_t &d)    {
    return (length(d) >= 0 ? first(d) : last(d));
  }
  static Element_t max(const Storage_t &d)    {
    return (length(d) >= 0 ? last(d) : first(d));
  }
  static bool      empty(const Storage_t &d)  { return false; }
  static int       loop(const Storage_t &)    { return 0; }

  // get the Nth value of the domain, where value # 0 is first(), etc.
  // for Region, this is only useful for n=0 and n=1.
  static Element_t elem(const Storage_t &d, int n) { return d[0] + n*d[1]; }

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  static OneDomain_t &getDomain(Domain_t &d, int) { return d; }
  static const OneDomain_t &getDomain(const Domain_t &d, int) { return d; }

  // Domains get the chance to do special initialization.  Region's are
  // initialized to have length 0 and, just to avoid having a random value,
  // to start at 0 (although, for a length 0 domain, the endpoints are
  // actually undefined).
  static void initializeStorage(Storage_t &dom) {
    dom[0] = 0;			// first
    dom[1] = 0;			// length
  }

  // change this domain object to the given one.  If things do not
  // match properly, assert a compile-time or run-time error.
  // For Region, we must have:
  // 1) the same dimensions==1
  // 2) stride of newdom == 1
  template<class DT>
  static void setDomain(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::dimensions == 1);
    dom[0] = DomainTraits<DT>::getFirst(newdom);
    dom[1] = DomainTraits<DT>::getLast(newdom) - dom[0];
  }

  // a specialized version of setDomain which accepts begin & end values.
  // For Region, we must have begval <= endval, since the stride is
  // hardcoded as +1.
  static void setDomain(Storage_t &dom, Element_t begval, Element_t endval) {
    dom[0] = begval;
    dom[1] = (endval - begval);
  }

  // change the loop variable for this object.  For Region, this is a no-op.
  static void setLoop(Storage_t &, int) { }

  // change the value of this 1D domain given a user-supplied reference
  // domain and a wildcard.
  template<class UT, class DT>
  static void setWildcardDomain(Storage_t &dom, const UT &u, const DT &newdom)
  {
    CTAssert(DomainTraits<DT>::wildcard);
    CTAssert(DomainTraits<DT>::dimensions == 1);
    CTAssert(DomainTraits<UT>::dimensions == 1);
    dom[0] = newdom.first(u); 	           // wildcard min
    dom[1] = newdom.last(u) - dom[0];      // wildcard max
  }

  //
  // compare this domain type to the given domain.  For the comparisons
  // to be meaningful for Region, we must have:
  // 1) the same dimensions==1
  // 2) both must not be empty
  //

  // 'isLessThan' returns true if dom < newdom
  template<class DT>
  static bool isLessThan(const Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::dimensions == 1);
    return (dom[1] < DomainTraits<DT>::getLength(newdom) ||
	    (dom[1] == DomainTraits<DT>::getLength(newdom) &&
	     dom[0]  < DomainTraits<DT>::getFirst(newdom)));
  }

  // 'isEqualTo' returns true if dom == newdom
  template<class DT>
  static bool isEqualTo(const Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::dimensions == 1);
    return (dom[0] == DomainTraits<DT>::getFirst(newdom) &&
	    dom[1] == DomainTraits<DT>::getLength(newdom));
  }

  //
  // arithmetic accumulation operators.  These only work with
  // other domain objects with the following characteristics:
  //   1) they are singleValue'd
  //   2) they have dimensions == 1
  //

  // addAccum means dom[0] += newdom
  template<class DT>
  static void addAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued);
    CTAssert(DomainTraits<DT>::dimensions == 1);
    dom[0] += DomainTraits<DT>::getFirst(newdom);
  }

  // subtractAccum means dom[0] -= newdom
  template<class DT>
  static void subtractAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued);
    CTAssert(DomainTraits<DT>::dimensions == 1);
    dom[0] -= DomainTraits<DT>::getFirst(newdom);
  }

  // multiplyAccum means dom[0] *= newdom and dom[1] *= newdom
  template<class DT>
  static void multiplyAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued &&
	     DomainTraits<DT>::dimensions == 1);
    dom[0] *= DomainTraits<DT>::getFirst(newdom);
    dom[1] *= DomainTraits<DT>::getFirst(newdom);
  }

  // divideAccum means dom[0] /= newdom and dom[1] /= newdom
  template<class DT>
  static void divideAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued &&
	     DomainTraits<DT>::dimensions == 1);
    dom[0] /= DomainTraits<DT>::getFirst(newdom);
    dom[1] /= DomainTraits<DT>::getFirst(newdom);
  }

  //
  // Iterator operations.  These functions work with the iterator
  // storage type to perform initialization and increment/decrement of the
  // iterator.  By putting this here, we can specialize for the cases
  // where we know the stride is a fixed number.
  //

  // initialize the iterator storage to the values from a domain.
  // For Region, this just means setting the current value of the iterator
  // to the beginning point of the interval we're referring to.
  static void initializeIterator(const Storage_t &d, IteratorStorage_t &i) {
    i[0] = d[0];
    i[1] = d[1];
  }

  // initialize the iterator storage to the values from a domain.
  // This version of initialize sets the resulting storage i to point to
  // d + 2*length(d2), generally used to set up an end iterator.
  // For Region, we know the length is d2[1].
  static void initializeIterator(const Storage_t &d1, const Storage_t &d2,
				 IteratorStorage_t &i) {
    i[0] = d1[0] + d2[1] + d2[1];
    i[1] = d2[1];
  }

  // copy the values from the first iterator storage into another
  static void copyIterator(IteratorStorage_t d, IteratorStorage_t &i) {
    i[0] = d[0];
    i[1] = d[1];
  }

  // return the current value of an iterator from the given iterator storage
  static Element_t currentIterator(IteratorStorage_t i) { return i[0]; }

  // compare for equality the two iterators from their storage
  static bool compareIterator(IteratorStorage_t a, IteratorStorage_t b) {
    return (a[0] == b[0] && a[1] == b[1]);
  }

  // increment or decrement the given iterator's storage.  For Region, we
  // know that there is nothing to do, since it has only one point.  But we
  // must still increment/decrement by one, since our end iterator is one past
  // the end of the domain.
  static void incrementIterator(IteratorStorage_t &i) { i[0] += i[1]; }
  static void decrementIterator(IteratorStorage_t &i) { i[0] -= i[1]; }
};


/**
 * DomainTraits<Region<1,POOMA_DEFAULT_POSITION_TYPE>>:
 * The traits for an 1-dimensional Region of type POOMA_DEFAULT_POSITION_TYPE.
 */

template<>
struct DomainTraits< Region<1,POOMA_DEFAULT_POSITION_TYPE> >
  : public DomainTraitsDomain<Region<1,POOMA_DEFAULT_POSITION_TYPE>, 
                              POOMA_DEFAULT_POSITION_TYPE, 1>
{
  // convenience typedef 
  typedef POOMA_DEFAULT_POSITION_TYPE T;

  // necessary typedefs
  typedef Region<1,T> OneDomain_t;
  typedef Region<1,T> BlockDomain_t;
  typedef Region<1,T> AskDomain_t;
  typedef Region<1,T> AddResult_t;
  typedef Region<1,T> MultResult_t;

  // 1D necessary typedefs.  Region requires two pieces of
  // information, the begin point and the length.  If length==0, this is
  // just a point.  For the iterator, we need to know the current position
  // and the stride (which will be the width of the domain).
  typedef Element_t Storage_t[2];
  typedef Element_t IteratorStorage_t[2];

  // necessary static data
  enum { dimensions = 1,
	 sliceDimensions = 1 };
  enum { loopAware = false };
  enum { singleValued = false };
  enum { unitStride = true };
  enum { wildcard = false };

  // return size, endpoint, stride, and loop information given the storage
  static Element_t first(const Storage_t &d)  { return d[0]; }
  static Element_t last(const Storage_t &d)   { return d[0] + d[1]; }
  static Element_t stride(const Storage_t &d) { return d[1]; }
  static Element_t length(const Storage_t &d) { return d[1]; }
  static Element_t min(const Storage_t &d)    {
    return (length(d) >= 0 ? first(d) : last(d));
  }
  static Element_t max(const Storage_t &d)    {
    return (length(d) >= 0 ? last(d) : first(d));
  }
  static bool      empty(const Storage_t &)  { return false; }
  static int       loop(const Storage_t &)    { return 0; }

  // get the Nth value of the domain, where value # 0 is first(), etc.
  // for Region, this is only useful for n=0 and n=1.
  static Element_t elem(const Storage_t &d, int n) { return d[0] + n*d[1]; }

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  static OneDomain_t &getDomain(Domain_t &d, int) { return d; }
  static const OneDomain_t &getDomain(const Domain_t &d, int) { return d; }

  // Domains get the chance to do special initialization.  Region's are
  // initialized to have length 0 and, just to avoid having a random value,
  // to start at 0 (although, for a length 0 domain, the endpoints are
  // actually undefined).
  static void initializeStorage(Storage_t &dom) {
    dom[0] = 0;			// first
    dom[1] = 0;			// length
  }

  // change this domain object to the given one.  If things do not
  // match properly, assert a compile-time or run-time error.
  // For Region, we must have:
  // 1) the same dimensions==1
  // 2) stride of newdom == 1
  template<class DT>
  static void setDomain(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::dimensions == 1);
    dom[0] = DomainTraits<DT>::getFirst(newdom);
    dom[1] = DomainTraits<DT>::getLast(newdom) - dom[0];
  }

  // a specialized version of setDomain which accepts begin & end values.
  // For Region, we must have begval <= endval, since the stride is
  // hardcoded as +1.
  static void setDomain(Storage_t &dom, Element_t begval, Element_t endval) {
    dom[0] = begval;
    dom[1] = (endval - begval);
  }

  // change the loop variable for this object.  For Region, this is a no-op.
  static void setLoop(Storage_t &, int) { }

  // change the value of this 1D domain given a user-supplied reference
  // domain and a wildcard.
  template<class UT, class DT>
  static void setWildcardDomain(Storage_t &dom, const UT &u, const DT &newdom)
  {
    CTAssert(DomainTraits<DT>::wildcard);
    CTAssert(DomainTraits<DT>::dimensions == 1);
    CTAssert(DomainTraits<UT>::dimensions == 1);
    dom[0] = newdom.first(u); 	           // wildcard min
    dom[1] = newdom.last(u) - dom[0];      // wildcard max
  }

  //
  // compare this domain type to the given domain.  For the comparisons
  // to be meaningful for Region, we must have:
  // 1) the same dimensions==1
  // 2) both must not be empty
  //

  // 'isLessThan' returns true if dom < newdom
  template<class DT>
  static bool isLessThan(const Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::dimensions == 1);
    return (dom[1] < DomainTraits<DT>::getLength(newdom) ||
	    (dom[1] == DomainTraits<DT>::getLength(newdom) &&
	     dom[0]  < DomainTraits<DT>::getFirst(newdom)));
  }

  // 'isEqualTo' returns true if dom == newdom
  template<class DT>
  static bool isEqualTo(const Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::dimensions == 1);
    return (dom[0] == DomainTraits<DT>::getFirst(newdom) &&
	    dom[1] == DomainTraits<DT>::getLength(newdom));
  }

  //
  // arithmetic accumulation operators.  These only work with
  // other domain objects with the following characteristics:
  //   1) they are singleValue'd
  //   2) they have dimensions == 1
  //

  // addAccum means dom[0] += newdom
  template<class DT>
  static void addAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued);
    CTAssert(DomainTraits<DT>::dimensions == 1);
    dom[0] += DomainTraits<DT>::getFirst(newdom);
  }

  // subtractAccum means dom[0] -= newdom
  template<class DT>
  static void subtractAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued);
    CTAssert(DomainTraits<DT>::dimensions == 1);
    dom[0] -= DomainTraits<DT>::getFirst(newdom);
  }

  // multiplyAccum means dom[0] *= newdom and dom[1] *= newdom
  template<class DT>
  static void multiplyAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued &&
	     DomainTraits<DT>::dimensions == 1);
    dom[0] *= DomainTraits<DT>::getFirst(newdom);
    dom[1] *= DomainTraits<DT>::getFirst(newdom);
  }

  // divideAccum means dom[0] /= newdom and dom[1] /= newdom
  template<class DT>
  static void divideAccum(Storage_t &dom, const DT &newdom) {
    CTAssert(DomainTraits<DT>::singleValued &&
	     DomainTraits<DT>::dimensions == 1);
    dom[0] /= DomainTraits<DT>::getFirst(newdom);
    dom[1] /= DomainTraits<DT>::getFirst(newdom);
  }

  //
  // Iterator operations.  These functions work with the iterator
  // storage type to perform initialization and increment/decrement of the
  // iterator.  By putting this here, we can specialize for the cases
  // where we know the stride is a fixed number.
  //

  // initialize the iterator storage to the values from a domain.
  // For Region, this just means setting the current value of the iterator
  // to the beginning point of the interval we're referring to.
  static void initializeIterator(const Storage_t &d, IteratorStorage_t &i) {
    i[0] = d[0];
    i[1] = d[1];
  }

  // initialize the iterator storage to the values from a domain.
  // This version of initialize sets the resulting storage i to point to
  // d + 2*length(d2), generally used to set up an end iterator.
  // For Region, we know the length is d2[1].
  static void initializeIterator(const Storage_t &d1, const Storage_t &d2,
				 IteratorStorage_t &i) {
    i[0] = d1[0] + d2[1] + d2[1];
    i[1] = d2[1];
  }

  // copy the values from the first iterator storage into another
  static void copyIterator(IteratorStorage_t d, IteratorStorage_t &i) {
    i[0] = d[0];
    i[1] = d[1];
  }

  // return the current value of an iterator from the given iterator storage
  static Element_t currentIterator(IteratorStorage_t i) { return i[0]; }

  // compare for equality the two iterators from their storage
  static bool compareIterator(IteratorStorage_t a, IteratorStorage_t b) {
    return (a[0] == b[0] && a[1] == b[1]);
  }

  // increment or decrement the given iterator's storage.  For Region, we
  // know that there is nothing to do, since it has only one point.  But we
  // must still increment/decrement by one, since our end iterator is one past
  // the end of the domain.
  static void incrementIterator(IteratorStorage_t &i) { i[0] += i[1]; }
  static void decrementIterator(IteratorStorage_t &i) { i[0] -= i[1]; }
};


/**
 * DomainChangeDim<T, int> is used to convert from a domain of one dimension
 * to another dimension (the second template parameter).
 * For Region<Dim1>, it changes from Dim1 to Dim2.
 */

template<int Dim1, int Dim2, class T>
struct DomainChangeDim<Region<Dim1,T>, Dim2>
{
  // the type of the old and new domain
  typedef Region<Dim1,T> OldType_t;
  typedef Region<Dim2,T> NewType_t;

  // static data for old and new dimensions
  enum { oldDim = Dim1,
	 newDim = Dim2 };
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_REGION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.Region.h,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
