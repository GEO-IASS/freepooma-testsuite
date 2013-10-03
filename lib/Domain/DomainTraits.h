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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_H

//-----------------------------------------------------------------------------
// Class:
// DomainTraits<T>
// DomainTraitsDomain<DomT, T, Dim>
// DOmainTraitsScalar<DomT, T>
// DomainChangeDim<T, Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits is a traits class for domain objects. Includes DomainChangeDim.
 *
 * This traits class
 * is used to specialize the Domain and DomainBase base classes for all
 * domain objects to do the proper action for that domain type, and to store
 * the proper information.
 *
 * DomainChangeDim is a simple struct used to convert a domain of type T
 * with a certain number of dimensions Dim1, to the same general type of
 * domain but with a different number of dimensions Dim2.  It defines two
 * typedefs 'OldType_t' and 'NewType_t' with the type of domain with Dim1
 * and Dim2, respectively.
 *
 * DomainTraits<T> is a traits class which provides all the specific
 * information and functionality to specialize how a given domain object
 * such as Loc<N> or Range<1> should store its data and interact with other
 * domain objects.  The general DomainTraits class is empty; when a new
 * domain object is added, a new, partially specialized DomainTraits class
 * must be written to provide the information on how that domain class
 * should act.  This domain class, along with a new class derived from
 * Domain<DomainTraits<NewDomainType<N>>, N> and the information in
 * NewDomain on how to combine this domain with others, is all that is
 * needed to add a new Domain type into POOMA.
 *
 * DomainTraits must be able to provide information about how the
 * domain object can be combined with others.  There are different ways
 * in which domain objects can be combined (see NewDomain.h for a list).
 *
 * Each DomainTraits class that is specialized for a specific domain object
 * should do the following:
 * -# DomainTraits should be templated on the main type, e.g., 'Loc<N>'.
 * -# DomainTraits should define the following typedefs and enums:
 *    - typedef Element_t : the type of data for the domain (e.g., int)
 *    - typedef Size_t : the type of data for the domain's size (e.g., long)
 *    - typedef Domain_t  : the domain type (what DomainTraits is templated on
 *    - typedef OneDomain_t : type of 1D domain which should be used to
 *                            combine this domain with others
 *    - typedef PointDomain_t : type of 1D domain that is returned if we try
 *                              to convert the domain to a single point.  Only
 *                              singleValue'd domains can really be converted
 *                              to a point; others just return OneDomain_t
 *    - typedef AskDomain_t : the type of domain returned when asking questions
 *                            about the domain, like "what are the first vals?"
 *    - typedef Storage_t : how the data for a domain should be stored.  An
 *                          instance of type Storage_t is kept in DomainBase.
 *    - typedef NewDomain1_t : the type of NewDomain1<Dom>::SliceType_t; 
 *                             the same as Domain_t, except for Scalar Domains.
 *    - typedef AddResult_t : the type of domain which should result if you
 *                            add or subtract this domain to a Loc or int.
 *    - typedef MultResult_t : the type of domain which should result if you
 *                             multiply or divide this domain by a Loc or int.
 *    - enum { domain } : if true, the type derived from Domain, 
 *                        but if false, the type is not really a 
 *                        domain (such as int)
 *    - enum { dimensions = # of dimensions }
 *    - enum { sliceDimensions = # of slice dimensions }
 *    - enum { loopAware } : whether domain stores info about which
 *                           loop variable it refers to
 *    - enum { singleValued } : whether domain refers to single point
 *    - enum { unitStride } : whether domain is unit-stride sequence
 * -# DomainTraits should include the following static methods, for
 *             multidimensional domains (that is, where Dim > 1)
 *    - static OneDomain_t &getDomain(Domain_t& d, int dim)
 *    - static const OneDomain_t &getDomain(const Domain_t& d, int dim)
 *        ==> return d[dim], for use in combining elements of this type
 *            of domain with others.
 *    - static PointDomain_t &getPointDomain(Domain_t& d, int dim)
 *    - static const PointDomain_t &getPointDomain(const Domain_t& d, int dim)
 *        ==> for single-valued domains, return the Nth 1D domain converted
 *            to a single point (such as an int or Loc, referring to a point).
 *            for a non-single-valued domain, just return the Nth 1D domain.
 *            this is used to distinguish between times when you need the
 *            domain as a single point, or converted to a larger domain (which
 *            can happen with 'getDomain').
 *    - static void initializeStorage(Storage_t &s)
 *        ==> initialize the given storage object to an initial, default state.
 * -# DomainTraits for a 1D domain type should include the static methods:
 *    - static OneDomain_t &getDomain(Domain_t& d, int dim)
 *    - static const OneDomain_t &getDomain(const Domain_t& d, int)
 *        ==> return d[0], for use in combining elements of this type
 *            of domain with others.  For the 1D case, this just returns
 *            the domain object back again, regardless of dim, so that 1D
 *            objects can be combined with all N 1D domains from another
 *            multidimensional domain.
 *    - static PointDomain_t &getPointDomain(Domain_t& d, int dim)
 *    - static const PointDomain_t &getPointDomain(const Domain_t& d, int dim)
 *        ==> same as for the N-dim case
 *    - static void initializeStorage(Storage_t &s)
 *        ==> initialize the given storage object to an initial, default state.
 *    - static Element_t first(const Storage_t &d)
 *        ==> return the first endpoint of the domain, using the info in d
 *    - static Element_t last(const Storage_t &d)
 *        ==> return the last endpoint of the domain, using the info in d
 *    - static Size_t length(const Storage_t &d)
 *        ==> return the length of the domain, using the info in d
 *    - static Element_t stride(const Storage_t &d)
 *        ==> return the stride of the domain, using the info in d
 *    - static Element_t min(const Storage_t &d)
 *        ==> return the minimum endpoint of the domain, using the info in d
 *    - static Element_t max(const Storage_t &d)
 *        ==> return the maximum endpoint of the domain, using the info in d
 *    - static bool empty(const Storage_t &d)
 *        ==> return true if the length of the domain is 0
 *    - static int loop(const Storage_t &d)
 *        ==> return which loop variable the domain refers to
 *    - static Element_t elem(const Storage_t &d, int n)
 *        ==> return the Nth element in the 1D domain; the 0th elem is
 *            first(), the length()-1'th elem is last(), etc.
 *    - template<class T> static void setDomain(Storage_t &, const T &newdom)
 *    - static void setDomain(Storage_t &dom, int newdom) {
 *        ==> change the domain to the values in newdom, which should be
 *            a domain object with the proper interface, or an integer.
 *    - static void setDomain(Storage_t &dom, ** necessary args **)
 *        ==> change the domain to the values from a specialized list of
 *            information, which is domain-type specific.  For example, for
 *            Range, this takes three extra args: first, last, stride
 *    - static void setLoop(Storage_t &, int loop)
 *        ==> change which loop variable is referred to
 *    - template<class T> static bool isLessThan(const Storage_t &, const T &)
 *    - static bool isLessThan(const Storage_t &dom, int newdom)
 *        ==> return if the domain object from the given storage is less than
 *            the second argument, either templated or specialized to int.
 *    - template<class T> static bool isEqualTo(const Storage_t &, const T &)
 *    - static bool isEqualTo(const Storage_t &dom, int newdom)
 *        ==> return if the domain object from the given storage is equal to
 *            the second argument, either templated or specialized to int.
 *    - template<class T> static void addAccum(Storage_t &, const T &);
 *    - static void addAccum(Storage_t &dom, int newdom)
 *        ==> add the value from the second argument to this domain.
 *    - template<class T> static void subtractAccum(Storage_t &, const T &);
 *    - static void subtractAccum(Storage_t &dom, int newdom)
 *        ==> subtract the value of the second argument from this domain.
 *    - template<class T> static void multiplyAccum(Storage_t &, const T &);
 *    - static void multiplyAccum(Storage_t &dom, int newdom)
 *        ==> multiply the value of this domain by the second argument
 *    - template<class T> static void divideAccum(Storage_t &, const T &);
 *    - static void divideAccum(Storage_t &dom, int newdom)
 *        ==> divide the value of this domain by the second argument
 * -# DomainTraits for all 1D domains and any other scalars need to also
 *    provide the following static methods, which all take as an argument
 *    the domain object (or scalar) itself:
 *    - static Element_t getFirst(const Domain_t &)
 *    - static Element_t getLast(const Domain_t &)
 *    - static Element_t getStride(const Domain_t &)
 *    - static Size_t    getLength(const Domain_t &)
 *    - static Size_t    getSize(const Domain_t &)
 *    - static Element_t getMin(const Domain_t &)
 *    - static Element_t getMax(const Domain_t &)
 *    - static bool      getEmpty(const Domain_t &)
 *    - static int       getLoop(const Domain_t &)
 *    .
 *    All of these methods, for domain objects, should just return the
 *    value of the relevant function (e.g., getFirst(d) returns d.first()).
 *    For non-domain objects such as scalars, they just return the scalar
 *    value itself or the relevant value (such as 'false' for getEmpty())
 *    The base classes 'DomainTraitsDomain' and 'DomainTraitsScalar' are
 *    provided for users to inherit from which provide implementations of
 *    these static methods already, along with a couple other domain trait
 *    settings (described next).
 *
 * There are two general classes of DomainTraits:
 *   -# DomainTraits for Domain classes, which provide traits for actual Domain
 *      subclasses.  When defining traits of of this form, you should have
 *      the partially-specialized class derive from
 *      DomainTraitsDomain<DomT, T, Dim>, where DomT is the type for Domain_t,
 *      T is the type for Element_t, and Dim is the dimension.  In this case,
 *      the base class DomainTraitsDomain provided part of the DomainTraits
 *      interface, including:
 *        - For N-D object:  Element_t, Domain_t, domain, and dimensions
 *        - For 1-D objects: Element_t, Domain_t, domain, and dimensions
 *                         get* methods (see list above), which just call
 *                         the corresponding method in the provide Domain
 *                         argument.
 *   -# Non-Domain DomainTraits, for things such as scalars.  For all types
 *      which do not have specific DomainTraits specializations, the default
 *      behavior is to treat them as scalars, and so the default version of
 *      DomainTraits<T> is derived from DomainTraitsScalar<T, T>.  The
 *      class DomainTraitsScalar is templated on the type for Domain_t and
 *      Element_t.  This is a separate class (and not just
 *      the default implementation of DomainTraits) because for the special
 *      case of integral types, there is a separate DomainTraits used to
 *      convert integral scalars to Locs or Intervals when needed.
 *      If you need to provide special traits for
 *      any other scalars, then have the specialization of DomainTraits
 *      inherit from DomainTraitsScalar.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

// Helper struct to compute proper type for Size_t typedef in DomainTraits

template <class T>
struct SizeTypePromotion
{
  typedef T Type_t;
};

template <>
struct SizeTypePromotion<int>
{
  typedef long Type_t;
};

template <>
struct SizeTypePromotion<float>
{ 
  typedef double Type_t;
};


/**
 * DomainTraitsDomain<DomT, T, Dim> can act as a base class for the partially-
 * specialized versions of DomainTraits for domain-like classes, that is,
 * classes which are derived from Domain.  It it templated on the types to
 * use for Domain_t and Element_t, respectively, and the dimension of the
 * domain.  The N-D version just defines Domain_t, Element_t, and static data
 * domain and dimensions, while the 1-D version also provided implementations
 * of the static get* methods (such as getFirst, getMin, etc).  Just inherit
 * from DomainTraitsDomain if you're defining traits for a new Domain subclass.
 *
 * N-dimensional version of DomainTraitsDomain
 */
template<class DomT, class T, int Dim>
struct DomainTraitsDomain
{
  // necessary typedefs
  typedef typename SizeTypePromotion<T>::Type_t Size_t;
  typedef    T Element_t;
  typedef DomT Domain_t;
  typedef DomT NewDomain1_t;

  // necessary static data
  enum { domain = true };
  enum { dimensions = Dim };

  // Put this here since no standard domain will ever override.
  static bool      getIgnorable(const Domain_t &, int)   { return false; }
};

/**
 * 1-dimensional specialized version of DomainTraitsDomain
 */
template<class DomT, class T>
struct DomainTraitsDomain<DomT, T, 1>
{
  // necessary typedefs
  typedef typename SizeTypePromotion<T>::Type_t Size_t;
  typedef    T Element_t;
  typedef DomT Domain_t;
  typedef DomT NewDomain1_t;

  // necessary static data
  enum { domain = true };
  enum { dimensions = 1 };

  // for the domain type, return the first, last, stride, and other
  // characteristics.  Since this base class is for domains, we know
  // we can call the same method for the provided argument
  inline
  static Element_t getFirst(const Domain_t &d)  { return d.first(); }
  inline
  static Element_t getLast(const Domain_t &d)   { return d.last(); }
  inline
  static Element_t getStride(const Domain_t &d) { return d.stride(); }
  inline
  static Size_t    getLength(const Domain_t &d) { return d.length(); }
  inline
  static Size_t    getSize(const Domain_t &d)   { return d.size(); }
  inline
  static Element_t getMin(const Domain_t &d)    { return d.min(); }
  inline
  static Element_t getMax(const Domain_t &d)    { return d.max(); }
  inline
  static bool      getEmpty(const Domain_t &d)  { return d.empty(); }
  inline
  static int       getLoop(const Domain_t &d)   { return d.loop(); }
  inline
  static Element_t getElem(const Domain_t &d, int n)   { return d.elem(n); }
  inline
  static bool      getIgnorable(const Domain_t &, int)   { return false; }
};


/**
 * DomainTraitsScalar<DomT, T, NewDom1T> can act as a base class for partially
 * specialized versions of DomainTraits for non-domain classes and types,
 * such as the basic scalar types. It it templated on the types to use for
 * Domain_t, Element_t, and NewDomain1_t.  It is used
 * as the base class for the DomainTraits<int> specialization, and for
 * DomainTraits<T> in general for all types T that do not have any other
 * specialized traits defined.  It provides definitions for most of the
 * standard traits settings, and implementations of the static get* methods.
 * For scalars, get* functions mostly just return back the same scalar,
 * except that:
 *   - for a scalar, the stride and length are always 1 and integers
 *   - for a scalar, getEmpty() is always false
 *   - for a scalar, getLoop() always returns 0 as an integer
 */

template<class DomT, class T, class NewDom1T>
struct DomainTraitsScalar
{
  // necessary typedefs
  typedef DomT      Domain_t;
  typedef DomT      OneDomain_t;
  typedef NewDom1T  NewDomain1_t;
  typedef T         PointDomain_t;
  typedef T         Element_t;
  typedef int       Size_t;

  // necessary static data
  enum { domain = false };
  enum { dimensions = 1,
	 sliceDimensions = 0 };
  enum { loopAware = false };
  enum { singleValued = true };
  enum { unitStride = true };
  enum { wildcard = false };

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).
  inline
  static OneDomain_t getDomain(T d, int) { return OneDomain_t(d); }

  // convert the scalar to a single point domain.  Since this is for
  // scalars, this means just return back the scalar.
  inline
  static PointDomain_t getPointDomain(T d, int) { return d; }

  // for the scalar type, return the first, last, stride, and other
  // characteristics
  inline
  static Element_t getFirst(const Element_t &d) { return d; }
  inline
  static Element_t getLast(const Element_t &d)  { return d; }
  inline
  static int       getStride(const Element_t &) { return 1; }
  inline
  static Size_t    getLength(const Element_t &) { return 1; }
  inline
  static Size_t    getSize(const Element_t &d)  { return 1; }
  inline
  static Element_t getMin(const Element_t &d)   { return d; }
  inline
  static Element_t getMax(const Element_t &d)   { return d; }
  inline
  static bool      getEmpty(const Element_t &)  { return false; }
  inline
  static int       getLoop(const Element_t &)   { return 0; }
  inline
  static Element_t getElem(const Element_t &d, int)   { return d; }
};


/**
 * So now, finally, we can define the default version of DomainTraits<T>
 * which just inherits from DomainTraitsScalar<T, T, T>
 */

template<class T>
struct DomainTraits : public DomainTraitsScalar<T, T, T> 
{ 
  // convenience typedef
  typedef DomainTraitsScalar<T, T, T>       Base_t;

  // necessary typedefs
  typedef typename Base_t::Domain_t         Domain_t;
  typedef typename Base_t::OneDomain_t      OneDomain_t;
  typedef typename Base_t::NewDomain1_t     NewDomain1_t;
  typedef typename Base_t::PointDomain_t    PointDomain_t;
  typedef typename Base_t::Element_t        Element_t;
  typedef typename Base_t::Size_t           Size_t;

  // necessary static data
  enum { domain = Base_t::domain };
  enum { dimensions = Base_t::dimensions,
	 sliceDimensions = Base_t::sliceDimensions };
  enum { loopAware = Base_t::loopAware };
  enum { singleValued = Base_t::singleValued };
  enum { unitStride = Base_t::unitStride };
  enum { wildcard = Base_t::wildcard };
};


/**
 * DomainChangeDim is a struct which is templated on a domain of type T,
 * which has some original number of dimensions oldDim, and a new number
 * of dimensions Dim.  It is used to determine a new domain type which is
 * of the same basic type but with a different number of dimensions (basically,
 * it changes the dimension of T from oldDim --> Dim).  It should define
 * the following typedefs and static data:
 *  - typedef ... OldType_t;               // the previous domain type, T<oldDim>
 *  - typedef ... NewType_t;               // the new domain type, T<Dim>
 *  - static const int oldDim = ... ;      // the original number of dimensions
 *  - static const int newDim = Dim;       // the new number of dimensions
 *
 * We define here the basic form of the struct, but empty.  Specialized
 * domain types such define a partially specialized version of this struct
 * with the above typedefs and static data relevant to that domain type.
 */

template<class T, int Dim>
struct DomainChangeDim {
  // the type of the old and new domain
  typedef T OldType_t;
  typedef T NewType_t;

  // static data for old and new dimensions
  enum { oldDim = Dim,
	 newDim = Dim };
};


/// global function templates for invoking DomainTraits::setDomain()

template <class Dom, class Storage, class T1, class T2>
inline void
setDomain(const Dom&, Storage& data, const T1& beg, const T2& end) {
  DomainTraits<Dom>::setDomain(data, beg, end);
}

// include specializations of DomainTraits for built-in integral types
#include "Domain/DomainTraits.int.h"

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.h,v $   $Author: richard $
// $Revision: 1.25 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
