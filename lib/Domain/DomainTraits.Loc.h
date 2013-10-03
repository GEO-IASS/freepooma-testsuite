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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_LOC_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_LOC_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<Loc<N>>
// DomainChangeDim<Loc<N>,Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<Loc<N>> is a specialization of the general DomainTraits
 * class, for the case of Loc domain objects.
 *
 * It defines the general
 * behavior of Loc, including its typedef and static data characteristics,
 * how to store data for a Loc, etc.  It is used by the Domain base class
 * of Loc to implement most of the public interface.
 *
 * DomainTraits<Loc<Dim>> stores the characteristics and much of the
 * implementation details for Loc domain objects.  A Loc acts like a
 * single integer point in N-dimensional space, so it is a single-valued,
 * unit-stride domain.
 *
 * A general version of DomainTraits<Loc<Dim>> is defined here, which
 * only includes the basic information to make Loc<Dim> look like an
 * array of Loc<1> objects.  DomainTraits<Loc<1>> is a more specific
 * specialization which provides most of the necessary interface information
 * for items which need to know about Loc.  Since most of the interface
 * for a domain object is only available for 1D versions of that domain
 * object, the Loc<1> specialization defines more interface functions than
 * the Loc<Dim> case.
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


/**
 * DomainTraits<Loc<Dim>>:
 * The specialization of DomainTraits for Loc, for dimensions greater than
 * one.
 */

template<int Dim>
struct DomainTraits< Loc<Dim> >
  : public DomainTraitsDomain<Loc<Dim>, int, Dim>
{
  // convenience typedef 
  typedef DomainTraitsDomain<Loc<Dim>, int, Dim>     Base_t;

  // necessary typedefs
  typedef typename Base_t::Element_t                 Element_t;
  typedef typename Base_t::Domain_t                  Domain_t;
  typedef typename Base_t::NewDomain1_t              NewDomain1_t;
  typedef Loc<1>        OneDomain_t;
  typedef Loc<1>        PointDomain_t;
  typedef Interval<Dim> BlockDomain_t;
  typedef Loc<Dim>      AskDomain_t;
  typedef Loc<Dim>      AddResult_t;
  typedef Loc<Dim>      MultResult_t;

  // type for storage of this domain's data
  typedef UninitializedVector<OneDomain_t,Dim,Element_t> Storage_t;

  // necessary static data
  enum { domain = Base_t::domain };
  enum { dimensions = Base_t::dimensions,
	 sliceDimensions = 0 };
  enum { loopAware = false };
  enum { singleValued = true };
  enum { unitStride = true };
  enum { wildcard = false };

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  inline
  static OneDomain_t &getDomain(Domain_t &d, int n) {
    return d[n];
  }
  inline
  static const OneDomain_t &getDomain(const Domain_t &d, int n) {
    return d[n];
  }

  // convert from the Nth element of the domain to a single point, and
  // return a PointDomain_t.
  inline
  static PointDomain_t &getPointDomain(Domain_t &d, int n) {
    return d[n];
  }
  inline
  static const PointDomain_t &getPointDomain(const Domain_t &d, int n) {
    return d[n];
  }

  // Domains get the chance to do special initialization.
  static void initializeStorage(Storage_t &dom) { dom.initialize(); }


  // addAccum means dom += newdom
  template<class T>
  inline
  static void addAccum(Storage_t &dom, const T &newdom) 
  {
    CTAssert(DomainTraits<T>::singleValued &&
	     (DomainTraits<T>::dimensions == 1 || 
	      DomainTraits<T>::dimensions == dimensions ) );

    if (DomainTraits<T>::dimensions > 1)
      for (int i = 0;i< DomainTraits<T>::dimensions ; ++i)
	dom[i] += DomainTraits<T>::getFirst(newdom[i]);
    else
      for (int i = 0;i< dimensions ; ++i)
	dom[i] += DomainTraits<T>::getFirst(newdom[0]);
  }

  // subtractAccum means dom -= newdom
  template<class T>
  inline
  static void subtractAccum(Storage_t &dom, const T &newdom) 
  {
    CTAssert(DomainTraits<T>::singleValued &&
	     (DomainTraits<T>::dimensions == 1 || 
	      DomainTraits<T>::dimensions == dimensions ) );
    if (DomainTraits<T>::dimensions > 1)
      for (int i = 0;i< DomainTraits<T>::dimensions ; ++i)
	dom[i] -= DomainTraits<T>::getFirst(newdom[i]);
    else
      for (int i = 0;i< dimensions ; ++i)
	dom[i] -= DomainTraits<T>::getFirst(newdom);
  }

  template<class T>
  static void multiplyAccum(Storage_t &dom, const T &newdom) 
  {
    CTAssert(DomainTraits<T>::singleValued &&
	     (DomainTraits<T>::dimensions == 1 || 
	      DomainTraits<T>::dimensions == dimensions ) );
    if (DomainTraits<T>::dimensions > 1)
      for (int i = 0;i< DomainTraits<T>::dimensions ; ++i)
	dom[i] *= DomainTraits<T>::getFirst(newdom[i]);
    else
      for (int i = 0;i< dimensions ; ++i)
	dom[i] *= DomainTraits<T>::getFirst(newdom);

  }
  
  template<class T>
  static void divideAccum(Storage_t &dom, const T &newdom) 
  {
    CTAssert(DomainTraits<T>::singleValued &&
	     (DomainTraits<T>::dimensions == 1 || 
	      DomainTraits<T>::dimensions == dimensions ) );
    if (DomainTraits<T>::dimensions > 1)
      for (int i = 0;i< DomainTraits<T>::dimensions ; ++i)
	dom[i] /= DomainTraits<T>::getFirst(newdom[i]);
    else
      for (int i = 0;i< dimensions ; ++i)
	dom[i] /= DomainTraits<T>::getFirst(newdom);

  }

};


/**
 * DomainTraits<Loc<1>>:
 * The specialization of DomainTraits for Loc, for dimension == 1.
 */

template<>
struct DomainTraits< Loc<1> >
  : public DomainTraitsDomain<Loc<1>, int, 1>
{
  // necessary typedefs
  typedef Loc<1>      OneDomain_t;
  typedef Loc<1>      PointDomain_t;
  typedef Interval<1> BlockDomain_t;
  typedef Loc<1>      AskDomain_t;
  typedef Loc<1>      AddResult_t;
  typedef Loc<1>      MultResult_t;

  // 1D necessary typedefs.  Loc's store just a single integer, which is the
  // point.  They cannot represent empty domains, and always have length == 1,
  // stride == 1.
  typedef Element_t   Storage_t;

  // necessary static data
  enum { dimensions = 1,
         sliceDimensions = 0 };
  enum { loopAware = false };
  enum { singleValued = true };
  enum { unitStride = true };
  enum { wildcard = false };

  // return size, endpoint, stride, and loop information from the storage
  // data for this domain
  inline
  static Element_t first(Storage_t d)    { return d; }
  inline
  static Element_t last(Storage_t d)     { return d; }
  inline
  static Element_t stride(Storage_t)     { return 1; }
  inline
  static Element_t length(Storage_t)     { return 1; }
  inline
  static Element_t min(Storage_t d)      { return d; }
  inline
  static Element_t max(Storage_t d)      { return d; }
  inline
  static bool      empty(Storage_t)      { return false; }
  inline
  static int       loop(Storage_t)       { return 0; }

  // get the Nth value of the domain, where value # 0 is first(), etc.
  inline
  static Element_t elem(Storage_t d, int) { return d; }

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  inline
  static OneDomain_t &getDomain(Domain_t &d, int) {
    return d;
  }
  inline
  static const OneDomain_t &getDomain(const Domain_t &d, int) {
    return d;
  }

  // convert from the Nth element of the domain to a single point, and
  // return a PointDomain_t.
  inline
  static PointDomain_t &getPointDomain(Domain_t &d, int) {
    return d;
  }
  inline
  static const PointDomain_t &getPointDomain(const Domain_t &d, int) {
    return d;
  }

  // Domains get the chance to do special initialization.
  // 1D Loc's are initialized to zero.
  inline
  static void initializeStorage(Storage_t &dom) {
    dom = 0;
  }

  // change this domain object to the given one.  If things do not
  // match properly, assert a compile-time or run-time error.
  // For Loc, we must have:
  // 1) the same dimensions==1
  // 2) length() == 1 for the new domain
  // NOTE: this could be changed to a more rigorous compile-time check
  // that the copied domain is singleValued, if that is considered important
  // for performance.
  template<class T>
  inline
  static void setDomain(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    PAssert(DomainTraits<T>::getLength(newdom) == 1);
    dom = DomainTraits<T>::getFirst(newdom);
  }

  // change the loop variable for this object.  For Loc, this is a no-op.
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
    dom = newdom.first(u);	// uses wildcard version of first()
  }

  //
  // compare this domain type to the given domain.  For the comparisons
  // to be meaningful for Loc, we must have:
  // 1) the same dimensions==1
  // 2) length() == 1 for the new domain
  // NOTE: this could be changed to a more rigorous compile-time check
  // that the copied domain is singleValued, if that is considered important
  // for performance.
  //

  // 'isLessThan' returns true if dom < newdom
  template<class T>
  static bool isLessThan(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    PAssert(DomainTraits<T>::getLength(newdom) == 1);
    return (dom < DomainTraits<T>::getFirst(newdom));
  }

  // 'isEqualTo' returns true if dom == newdom
  template<class T>
  static bool isEqualTo(const Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    PAssert(DomainTraits<T>::getLength(newdom) == 1);
    return (dom == DomainTraits<T>::getFirst(newdom));
  }

  //
  // arithmetic accumulation operators.  These only work with
  // other domain objects with the following characteristics:
  // 1) they are singleValue'd
  // 2) they have dimensions == 1 -or- the have dimension 
  // equal to the dimension of *this
  //
  // Note that for Locs, we do NOT allow *= or /=.  You
  // must convert a Loc to a Range before doing multiplicative operations.
  //

  // addAccum means dom += newdom
  template<class T>
  inline
  static void addAccum(Storage_t &dom, const T &newdom) {
   
    CTAssert(DomainTraits<T>::singleValued &&
    	     DomainTraits<T>::dimensions == 1);
   dom += DomainTraits<T>::getFirst(newdom);

  }

  // subtractAccum means dom -= newdom
  template<class T>
  inline
  static void subtractAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom -= DomainTraits<T>::getFirst(newdom);
  }

  template<class T>
  static void multiplyAccum(Storage_t &dom, const T &newdom) {
    CTAssert(DomainTraits<T>::singleValued &&
	     DomainTraits<T>::dimensions == 1);
    dom *= DomainTraits<T>::getFirst(newdom);
  }
  
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
 * For Loc<Dim1>, it changes from Dim1 to Dim2.
 */

template<int Dim1, int Dim2>
struct DomainChangeDim<Loc<Dim1>, Dim2>
{
  // the type of the old and new domain
  typedef Loc<Dim1> OldType_t;
  typedef Loc<Dim2> NewType_t;

  // static data for old and new dimensions
  enum { oldDim = Dim1,
	 newDim = Dim2 };
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_LOC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.Loc.h,v $   $Author: richard $
// $Revision: 1.34 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
