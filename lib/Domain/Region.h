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
//-----------------------------------------------------------------------------
// Class:
// Region<int, T = double>
//-----------------------------------------------------------------------------

#ifndef POOMA_DOMAIN_REGION_H
#define POOMA_DOMAIN_REGION_H

/** @file
 * @ingroup Domain
 * @brief
 * Region is a general type of continuous domain, which refers to all points
 * between two endpoints a and b.
 *
 * It is basically an array of Region<1> objects.
 * It is templated on the number of dimensions, and the data type used to
 * store the values (generally double or float, but possibly any other type).
 * The macro POOMA_DEFAULT_POSITION_TYPE defines the type for a default
 * parameter value for the floating point type; if this macro is not defined,
 * double is used.  So you can construct a Region<N>, and there will be a
 * default type T = double used.  The user can override what the default
 * type should be by defining POOMA_DEFAULT_POSITION_TYPE when their
 * application is built.
 *
 * Region defers most of its implementation to the Domain<DomainTraits<Region>>
 * base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/DomainTraits.Region.h"
#include "Domain/NewDomain.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Region<N,T> is a domain representing a set of N continuous 1D regions, one
 * for each dimension N.  The regions have endpoints [a,b], and Region refers
 * to the rectangular region defined by these endpoints.
 *
 * You can construct a Region object using other domain objects.
 * The constructors accept up to 7 domain objects of various types.
 * Domain types are, for example, Loc, Region, Interval. An int, double, or
 * float may also be used
 * in the constructor for Region; it acts like a Loc<1> object
 * in that context.  The domain arguments for the Region
 * constructors are combined together to form a single domain object with
 * a dimension equal to the sum of the arguments dimensions; for example,
 * if you try to create a Region<3,T> from a Loc<2> and an Interval<1>, e.g.
 *   Region<3,T> a(Loc<2>(1,2), Interval<1>(3,5));
 * the Loc<2> and Interval arguments are combined into a (2+1=3) dimensional
 * domain object, used to initialize the Region<3>.  The number of dimensions
 * for the arguments must be <= number of dimensions of the newly constructed
 * Region.
 *
 * Since Region has possibly floating-point values, when constructed from
 * integer domains such as Interval it just takes the first and last values of
 * their domain as the endpoints.  When using a Region to construct an integer
 * domain, the standard first(), length(), stride() methods will be called,
 * and for Region this will just lead to integer conversion (truncation) to
 * get the integer values.
 *
 * For Region<1,T>, the list of constructors is limited to just the following:
 *  - Region<1,T> a() - default constructor, which creates a Region for just the
 *                     origin, e.g. 0.
 *  - Region<1,T> a(n) - sets the Region to the sequence [0 ... n]
 *  - Region<1,T> a(m,n) - sets the Region to the sequence [m ... n]
 *  - Region<1,T> a(m,n,s) - sets the Region to the sequence [m ... n];
 *                          s is ignored
 *
 * The default Region<1,T> constructor initializes the Region to the origin.
 *
 * In addition to the constructors, Region has the following public
 * interface, similar to all domain objects.  There are two classes of
 * interface methods, one class which includes methods which any Region<N,T>
 * object has, regardless of dimensions, the other class which includes extra
 * interface methods that are available for just Region<1,T> objects.
 *
 * Region<N,T> interface:
 *  - T size() - return the 'volume' of the domain, which is the product
 *      of the lengths of the N 1D Regions
 *  - bool empty() - always false here
 *  - Region<1,T> operator[](int N) - return the Nth Region<1,T> in a
 *      multidimensional Region<M,T>.  For Region<1,T> objects, this just
 *      returns the object back.
 *  - comparison operators: <, >, !=, ==, <=, >= : compare a Region<N,T> to
 *      another domain object.  The compared domains must have the same
 *      number of dimensions.
 *  - arithmetic accumulation operators +=, -=, *=, /= : add or subtract in a
 *      given domain.  The added domain must have the same number of
 *      dimensions, or a dimension of 1 (in which case, the same value
 *      is used for all dimensions), and be known to be single-valued (which
 *      is true for Loc and int's).  Note that for Region, *= and /= ARE
 *      allowed, since Region can have its stride changed at run time.  *=
 *      and /= result in scaling of the endpoints and stride, which leaves
 *      the length (and size) the same.  += and -= shift the beginning
 *      endpoints by the given values, also leaving the length and size the
 *      same.  Negation of a Region negates the endpoints and stride.
 *  - binary arithmetic operators +, -, *, / : for + and -, adding a Region
 *      to another Loc or int returns a new Region.  For * and /, scaling
 *      by a Loc or int also returns a Region object, since the stride may
 *      change.
 *  - increment/decrement operator ++, -- : only prefix versions of ++ and --
 *      are provided; they act just like += length() and -= length() ops.
 *
 * Region<1,T> interface:
 * all the methods for Region<N,T> are also available for Region<1,T>. Plus:
 *  - T length() - number of elements (including endpoints) of the domain.
 *     Really, this should be either 1 or infinity, but it defined here
 *     somewhat differently, as just the distance between the endpoints.  So
 *     a length of zero really means this just refers to one point, not that
 *     this is empty.  A Region cannot be empty, it must refer to SOME point(s).
 *  - T first() - the beginning endpoint.
 *  - T last() - the ending endpoint.
 *  - T stride() - here, the same as the length()
 *  - T min(), T max() - min or max of the endpoints.
 *  - bool empty() - always false for a Region
 *  - Region<1,T>::iterator begin() and end() - return iterators for the 1D
 *      domain.  These act like (at least) forward iterators.
 *
 * For the special case of Region<1,T>, there is a specialization given
 * after the general case that has different constructors.
 *
 * Region inherits much of its activity from Domain<DomainTraits<Region>>
 */

template<int Dim, class T = POOMA_DEFAULT_POSITION_TYPE>
class Region : public Domain<Dim, DomainTraits<Region<Dim,T> > >
{
  // convenience typedefs
  typedef DomainTraits< Region<Dim,T> >         DT_t;
  typedef Domain<Dim, DT_t>                     Base_t;

public:
  // Typedefs obtained from the base class and DomainTraits

  typedef typename Base_t::iterator             iterator;
  typedef typename Base_t::const_iterator       const_iterator;
  typedef typename Base_t::blockIterator        blockIterator;
  typedef typename Base_t::const_blockIterator  const_blockIterator;
  
  typedef typename DT_t::Element_t            Element_t;
  typedef typename DT_t::Domain_t             Domain_t;
  typedef typename DT_t::OneDomain_t          OneDomain_t;
  typedef typename DT_t::BlockDomain_t        BlockDomain_t;
  typedef typename DT_t::AskDomain_t          AskDomain_t;
  typedef typename DT_t::AddResult_t          AddResult_t;
  typedef typename DT_t::MultResult_t         MultResult_t;
  typedef typename DT_t::Storage_t            Storage_t;

  // duplicate static data from traits class

  enum { domain          = DT_t::domain };
  enum { dimensions      = DT_t::dimensions,
	 sliceDimensions = DT_t::sliceDimensions };
  enum { loopAware       = DT_t::loopAware };
  enum { singleValued    = DT_t::singleValued };
  enum { unitStride      = DT_t::unitStride };
  enum { wildcard        = DT_t::wildcard };

  //
  // Constructors.
  //

  // default constructor : initialize to refer to the origin
  Region() { }

  // NoInit constructor
  Region(const Pooma::NoInit &e) 
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(e) {
  }

  // copy constructor
  Region(const Region<Dim,T> &a)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain1<Region<Dim,T> >::fill(*this, a);
  }

  // templated constructors, taking from 1 to 7 different domain objects
  // (or integers).
  template<class T1>
  explicit Region(const T1 &a)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  template<class T1, class T2>
  Region(const T1 &a, const T2 &b)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain2<T1,T2>::fill(*this, a, b);
  }

  template<class T1, class T2, class T3>
  Region(const T1 &a, const T2 &b, const T3 &c)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain3<T1,T2,T3>::fill(*this, a, b, c);
  }

  template<class T1, class T2, class T3, class T4>
  Region(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain4<T1,T2,T3,T4>::fill(*this, a, b, c, d);
  }

  template<class T1, class T2, class T3, class T4, class T5>
  Region(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain5<T1,T2,T3,T4,T5>::fill(*this, a, b, c, d, e);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6>
  Region(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
	const T6 &f)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain6<T1,T2,T3,T4,T5,T6>::fill(*this, a, b, c, d, e, f);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6, class T7>
  Region(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
	const T6 &f, const T7 &g)
    : Domain<Dim, DomainTraits<Region<Dim,T> > >(Pooma::NoInit()) {
    NewDomain7<T1,T2,T3,T4,T5,T6,T7>::fill(*this, a, b, c, d, e, f, g);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Region() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T1>
  Region<Dim,T> &operator=(const T1 &newdom) {
    return NewDomain1<T1>::fill(*this, newdom);
  }

  Region<Dim,T> &operator=(const Region<Dim,T> &newdom) {
    return NewDomain1<Region<Dim,T> >::fill(*this, newdom);
  }

protected:

private:

};


/**
 * Region<1> is a 1D specialization of Region<N>; for the 1D case,
 * there are only a restricted set of constructors available.
 * For the special case of Region<1>, the following constructors
 * are defined:
 *  - Region<1> a() - default constructor, which creates an EMPTY Region
 *  - Region<1> a(n) - sets the Region to the sequence [0 ... (n-1)], stride 1
 *  - Region<1> a(m,n) - sets the Region to be [m ... n], stride 1 or -1
 *  - Region<1> a(m,n,s) - sets the Region to the sequence [m ... n], stride s
 *  - Region<1> a(Domain d) : a Region copied from d, which must be a
 *     1D domain object.
 *
 * Since partially-specialized classes cannot have default template
 * parameters (WHY???? Because The Standards Committee Is Evil), there
 * is a version of the 1D specialization of Region for general type T, and
 * a further specialization to 1D and type POOMA_DEFAULT_POSITION_TYPE.
 */

// specialization for 1D Region of type T
template<class T>
class Region<1,T> : public Domain<1, DomainTraits<Region<1,T> > >
{
  // convenience typedefs
  typedef DomainTraits< Region<1,T> >           DT_t;
  typedef Domain<1, DT_t>                       Base_t;

public:
  // Typedefs obtained from the base class and DomainTraits

  typedef typename Base_t::iterator             iterator;
  typedef typename Base_t::const_iterator       const_iterator;
  typedef typename Base_t::blockIterator        blockIterator;
  typedef typename Base_t::const_blockIterator  const_blockIterator;
  
  typedef typename DT_t::Element_t            Element_t;
  typedef typename DT_t::Domain_t             Domain_t;
  typedef typename DT_t::OneDomain_t          OneDomain_t;
  typedef typename DT_t::BlockDomain_t        BlockDomain_t;
  typedef typename DT_t::AskDomain_t          AskDomain_t;
  typedef typename DT_t::AddResult_t          AddResult_t;
  typedef typename DT_t::MultResult_t         MultResult_t;
  typedef typename DT_t::Storage_t            Storage_t;

  // duplicate static data from traits class

  enum { domain          = DT_t::domain };
  enum { dimensions      = DT_t::dimensions,
	 sliceDimensions = DT_t::sliceDimensions };
  enum { loopAware       = DT_t::loopAware };
  enum { singleValued    = DT_t::singleValued };
  enum { unitStride      = DT_t::unitStride };
  enum { wildcard        = DT_t::wildcard };

  //
  // Constructors.
  //

  // default constructor
  Region() { }

  // NoInit constructor
  Region(const Pooma::NoInit &e) 
    : Domain<1, DomainTraits<Region<1,T> > >(e) {
  }

  // copy constructor
  Region(const Region<1,T> &a)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    NewDomain1<Region<1,T> >::fill(*this, a);
  }

  // general argument constructor, to copy from a different domain type
  template<class T1>
  explicit Region(const T1 &a)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  // initialize from a single value: sets endpoints to [0..n].
  // domain_m is the domain information storage kept in the base class.
  Region(Element_t n)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    DomainTraits<Region<1,T> >::setDomain(this->domain_m, 0, n);
  }

  // initialize from a set of endpoints: sets endpoints to [m..n].
  // domain_m is the domain information storage kept in the base class.
  Region(Element_t m, Element_t n)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    DomainTraits<Region<1,T> >::setDomain(this->domain_m, m, n);
  }

  // initialize from a set of endpoints and with a given stride.
  // for Region, the stride is ignored.
  // domain_m is the domain information storage kept in the base class.
  Region(Element_t m, Element_t n, Element_t)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    DomainTraits<Region<1,T> >::setDomain(this->domain_m, m, n);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Region() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T1>
  Region<1,T> &operator=(const T1 &newdom) {
    return NewDomain1<T1>::fill(*this, newdom);
  }

  Region<1,T> &operator=(const Region<1,T> &newdom) {
    return NewDomain1<Region<1,T> >::fill(*this, newdom);
  }

};


// specialization for 1D Region of type POOMA_DEFAULT_POSITION_TYPE
template<>
class Region<1,POOMA_DEFAULT_POSITION_TYPE>
  : public Domain<1, DomainTraits<Region<1,POOMA_DEFAULT_POSITION_TYPE> > >
{
  // typedef for the element type, to make things easier here
  typedef POOMA_DEFAULT_POSITION_TYPE T;
  // convenience typedef
  typedef DomainTraits< Region<1,T> >  DT_t;

public:
  // typedefs from DomainTraits

  typedef DT_t::Element_t            Element_t;
  typedef DT_t::Domain_t             Domain_t;
  typedef DT_t::OneDomain_t          OneDomain_t;
  typedef DT_t::BlockDomain_t        BlockDomain_t;
  typedef DT_t::AskDomain_t          AskDomain_t;
  typedef DT_t::AddResult_t          AddResult_t;
  typedef DT_t::MultResult_t         MultResult_t;
  typedef DT_t::Storage_t            Storage_t;

  // duplicate static data from traits class

  enum { domain          = DT_t::domain };
  enum { dimensions      = DT_t::dimensions,
	 sliceDimensions = DT_t::sliceDimensions };
  enum { loopAware       = DT_t::loopAware };
  enum { singleValued    = DT_t::singleValued };
  enum { unitStride      = DT_t::unitStride };
  enum { wildcard        = DT_t::wildcard };

  //
  // Constructors.
  //

  // default constructor
  Region() { }

  // NoInit constructor
  Region(const Pooma::NoInit &e) 
    : Domain<1, DomainTraits<Region<1,T> > >(e) {
  }

  // copy constructor
  Region(const Region<1,T> &a)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    NewDomain1<Region<1,T> >::fill(*this, a);
  }

  // general argument constructor, to copy from a different domain type
  template<class T1>
  explicit Region(const T1 &a)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  // initialize from a single value: sets range to [0..n-1].  Must
  // have n >= 0. domain_m is the domain information storage kept in
  // the base class.
  Region(Region<1,T>::Element_t n)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    DomainTraits<Region<1,T> >::setDomain(this->domain_m, 0, n);
  }

  // initialize from a set of endpoints: sets range to [m ..n].
  // domain_m is the domain information storage kept in the base class.
  Region(Region<1,T>::Element_t m,
	Region<1,T>::Element_t n)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    DomainTraits<Region<1,T> >::setDomain(this->domain_m, m, n);
  }

  // initialize from a set of endpoints and with a given stride.
  // domain_m is the domain information storage kept in the base class.
  Region(Region<1,T>::Element_t m,
	Region<1,T>::Element_t n,
	Region<1,T>::Element_t)
    : Domain<1, DomainTraits<Region<1,T> > >(Pooma::NoInit()) {
    DomainTraits<Region<1,T> >::setDomain(this->domain_m, m, n);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Region() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T1>
  Region<1,T> &operator=(const T1 &newdom) {
    return NewDomain1<T1>::fill(*this, newdom);
  }

  Region<1,T> &operator=(const Region<1,T> &newdom) {
    return NewDomain1<Region<1,T> >::fill(*this, newdom);
  }

};

#endif     // POOMA_DOMAIN_REGION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Region.h,v $   $Author: richard $
// $Revision: 1.28 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
