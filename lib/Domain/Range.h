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

#ifndef POOMA_DOMAIN_RANGE_H
#define POOMA_DOMAIN_RANGE_H

//-----------------------------------------------------------------------------
// Class:
// Range<int>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Range is a general type of integer domain, which refers to a set of points
 * a, a+s, a+2s, ..., b.
 *
 * It has a run-time specified stride value s.  It
 * is basically an array of Range<1> objects.
 * Range defers most of its implementation to the Domain<DomainTraits<Range>>
 * base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/DomainTraits.Range.h"
#include "Domain/NewDomain.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template <int Dim> class Range;

template<int Dim>
inline
void fillRangeScalar(Range<Dim> &r, const int &a);

/**
 * Range<N> is a domain representing a set of N numeric sequences, one
 * for each dimension N.  The sequences have endpoints [a,b], with a stride
 * s.
 *
 * You can construct a Range object using other domain objects.
 * The constructors accept up to 7 domain objects of various types.
 * Domain types are, for example, Loc, Range, Interval. int may also be used
 * in the constructor for Range; it acts like a Loc<1> object
 * in that context.  The domain arguments for the Range
 * constructors are combined together to form a single domain object with
 * a dimension equal to the sum of the arguments dimensions; for example,
 * if you try to create a Range<3> from a Loc<2> and an Interval<1>, e.g.
 *   Range<3> a(Loc<2>(1,2), Interval<1>(3,5));
 * the Loc<2> and Interval arguments are combined into a (2+1=3) dimensional
 * domain object, used to initialize the Range<3>.  The number of dimensions
 * for the arguments must be <= number of dimensions of the newly constructed
 * Range.
 *
 * For Range<1>, the list of constructors is limited to just the following:
 *  - Range<1> a() - default constructor, which creates an EMPTY Range
 *  - Range<1> a(n) - sets the Range to the sequence [0 ... (n-1)], stride 1
 *  - Range<1> a(m,n) - sets the Range to the sequence [m ... n], stride 1 or -1
 *  - Range<1> a(m,n,s) - sets the Range to the sequence [m ... n], stride s
 *
 * The default Range<1> constructor initializes the Range to be empty,
 * that is, to have a length() of zero.  In that case, the endpoints are
 * undefined, as is any operation involving the Range.
 *
 * In addition to the constructors, Range has the following public
 * interface, similar to all domain objects.  There are two classes of
 * interface methods, one class which includes methods which any Range<N>
 * object has, regardless of dimensions, the other class which includes extra
 * interface methods that are available for just Range<1> objects.
 *
 * Range<N> interface:
 *  - long size() - return the 'volume' of the domain, which is the product
 *      of the lengths of the N 1D Ranges
 *  - bool empty() - return if any of the Range<1> objects have length == 0
 *  - Range<1> operator[](int N) - return the Nth Range<1> in a
 *      multidimensional Range<M>.  For Range<1> objects, this just
 *      returns the object back.
 *  - comparison operators: <, >, !=, ==, <=, >= : compare a Range<N> to
 *      another domain object.  The compared domains must have the same
 *      number of dimensions.
 *  - arithmetic accumulation operators +=, -=, *=, /= : add or subtract in a
 *      given domain.  The added domain must have the same number of
 *      dimensions, or a dimension of 1 (in which case, the same value
 *      is used for all dimensions), and be known to be single-valued (which
 *      is true for Loc and int's).  Note that for Range, *= and /= ARE
 *      allowed, since Range can have its stride changed at run time.  *=
 *      and /= result in scaling of the endpoints and stride, which leaves
 *      the length (and size) the same.  += and -= shift the beginning
 *      endpoints by the given values, also leaving the length and size the
 *      same.  Negation of a Range negates the endpoints and stride.
 *  - binary arithmetic operators +, -, *, / : for + and -, adding a Range
 *      to another Loc or int returns a new Range.  For * and /, scaling
 *      by a Loc or int also returns a Range object, since the stride may
 *      change.
 *  - increment/decrement operator ++, -- : only prefix versions of ++ and --
 *      are provided; they act just like += 1 and -= 1 operations.
 *
 * Range<1> interface:
 * all of the methods for Range<N> are also available for Range<1>. Plus:
 *  - long length() - number of elements (including endpoints) of the domain.
 *     for a non-unit-stride Range, the length refers to the number of
 *     strided points (including the endpoints), NOT the difference between
 *     the first and last endpoint.  That is, length = (end-beg)/stride + 1,
 *     NOT (end-beg) + 1.
 *  - int first() - the beginning endpoint.
 *  - int last() - the ending endpoint.
 *  - int min(), int max() - min or max of the endpoints.
 *  - Range<1>::iterator begin() and end() - return iterators for the 1D
 *      domain.  These act like (at least) forward iterators.
 *
 * For the special case of Range<1>, there is a specialization given
 * after the general case that has different constructors.
 *
 * Range inherits much of its activity from Domain<DomainTraits<Range>>
 */

template<int Dim>
class Range : public Domain<Dim, DomainTraits<Range<Dim> > >
{
  // convenience typedefs
  typedef DomainTraits< Range<Dim> >            DT_t;
  typedef Domain<Dim, DT_t>                     Base_t;

public:
  // Typedefs from parent class and DomainTraits

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

  // default constructor : initialize to an empty domain
  Range() { }

  // copy constructor
  Range(const Range<Dim> &a)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain1<Range<Dim> >::fill(*this, a);
  }

  // Uninitialized constructor
  Range(const Pooma::NoInit &a)
    : Domain<Dim, DomainTraits<Range<Dim> > >(a) 
  { }

  // templated constructors, taking from 1 to 7 different domain objects
  // (or integers).
  template<class T1>
  explicit Range(const T1 &a)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  template<class T1, class T2>
  Range(const T1 &a, const T2 &b)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain2<T1,T2>::fill(*this, a, b);
  }

  template<class T1, class T2, class T3>
  Range(const T1 &a, const T2 &b, const T3 &c)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain3<T1,T2,T3>::fill(*this, a, b, c);
  }

  template<class T1, class T2, class T3, class T4>
  Range(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain4<T1,T2,T3,T4>::fill(*this, a, b, c, d);
  }

  template<class T1, class T2, class T3, class T4, class T5>
  Range(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain5<T1,T2,T3,T4,T5>::fill(*this, a, b, c, d, e);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6>
  Range(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
	const T6 &f)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain6<T1,T2,T3,T4,T5,T6>::fill(*this, a, b, c, d, e, f);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6, class T7>
  Range(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
	const T6 &f, const T7 &g)
    : Domain<Dim, DomainTraits<Range<Dim> > >(Pooma::NoInit()) {
    NewDomain7<T1,T2,T3,T4,T5,T6,T7>::fill(*this, a, b, c, d, e, f, g);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Range() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  Range<Dim> &operator=(const T &newdom) {
    return NewDomain1<T>::fill(*this, newdom);
  }

  Range<Dim> &operator=(const Range<Dim> &newdom) {
    return NewDomain1<Range<Dim> >::fill(*this, newdom);
  }

  Range<Dim> &operator=(const int a) {
    fillRangeScalar(*this,a);
    return *this;
  }

protected:

private:

};


/**
 * Range<1> is a 1D specialization of Range<N>; for the 1D case,
 * there are only a restricted set of constructors available.
 * For the special case of Range<1>, the following constructors
 * are defined:
 *  - Range<1> a() - default constructor, which creates an EMPTY Range
 *  - Range<1> a(n) - sets the Range to the sequence [0 ... (n-1)], stride 1
 *  - Range<1> a(m,n) - sets the Range to the sequence [m ... n], stride 1 or -1
 *  - Range<1> a(m,n,s) - sets the Range to the sequence [m ... n], stride s
 *  - Range<1> a(Domain d) : a Range copied from d, which must be a
 *     1D domain object.
 */

template<>
class Range<1> : public Domain<1, DomainTraits<Range<1> > >
{
  // convenience typedef
  typedef DomainTraits< Range<1> >  DT_t;

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
  Range() { }

  // copy constructor
  Range(const Range<1> &a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    NewDomain1<Range<1> >::fill(*this, a);
  }

  // Uninitialized constructor
  Range(const Pooma::NoInit &a)
    : Domain<1, DomainTraits<Range<1> > >(a) 
  { }

  // general argument constructor, to copy from a different domain type
  template<class T1>
  explicit Range(const T1 &a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  // initialize from a single scalar.  must specialize to all scalars.
  Range(char a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - 1);
  }
  Range(unsigned char a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - 1);
  }
  Range(short a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    short s = (a < 0 ? -1 : 1);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - s);
  }
  Range(unsigned short a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - 1);
  }
  Range(int a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    int s = (a < 0 ? -1 : 1);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - s);
  }
  Range(unsigned int a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - 1);
  }
  Range(long a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    long s = (a < 0 ? -1 : 1);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - s);
  }
  Range(unsigned long a)
    : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
    PAssert(a != 0);
    DomainTraits<Range<1> >::setDomain(domain_m, 0, a - 1);
  }

  // initialize from a set of endpoints: sets range to [m ..n].
  // domain_m is the domain information storage kept in the base class.
  template<class T1, class T2>
  Range(const T1 &m, const T2 &n);

  // initialize from a set of endpoints and with a given stride.
  // domain_m is the domain information storage kept in the base class.
  template<class T1, class T2, class T3>
  Range(const T1 &m, const T2 &n, const T3 &s);

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Range() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  Range<1> &operator=(const T &newdom) {
    return NewDomain1<T>::fill(*this, newdom);
  }

  Range<1> &operator=(const Range<1> &newdom) {
    return NewDomain1<Range<1> >::fill(*this, newdom);
  }

};

// initialize from a set of endpoints: sets range to [m ..n].
// domain_m is the domain information storage kept in the base class.
// Not sure why this isn't required everywhere.
#if defined(__MWERKS__)
template<>
#endif
template <class T1, class T2>
inline
Range<1>::Range(const T1 &m, const T2 &n)
  : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
  DomainTraits<Range<1> >::setDomain(domain_m, m, n);
}

// initialize from a set of endpoints and with a given stride.
// domain_m is the domain information storage kept in the base class.
// Not sure why this isn't required everywhere.
#if defined(__MWERKS__)
template<>
#endif
template <class T1, class T2, class T3>
inline
Range<1>::Range(const T1 &m, const T2 &n, const T3 &s)
  : Domain<1, DomainTraits<Range<1> > >(Pooma::NoInit()) {
  DomainTraits<Range<1> >::setDomain(domain_m, m, n, s);
}


template<int Dim>
inline
void fillRangeScalar(Range<Dim> &r, const int &a)
{
  for (int i=0; i < Dim; ++i)
    r[i]=Range<1>(a);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_RANGE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Range.h,v $   $Author: richard $
// $Revision: 1.24 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
