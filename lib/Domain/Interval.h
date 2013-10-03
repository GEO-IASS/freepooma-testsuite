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

#ifndef POOMA_DOMAIN_INTERVAL_H
#define POOMA_DOMAIN_INTERVAL_H

//-----------------------------------------------------------------------------
// Class:
// Interval<int>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Interval is a very simple type of domain, which refers to a set of points
 * a, a+1, a+2, ..., b.
 *
 * It has a hard-coded stride of 1.  Interval<N>
 * is basically an array of Interval<1> objects.
 * Interval defers most of its implementation to the
 * Domain<DomainTraits<Interval>> base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/DomainTraits.Interval.h"
#include "Domain/NewDomain.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Interval<N> is a domain representing a set of N numeric sequences, one
 * for each dimension N.  The sequences have endpoints [a,b], and a
 * hard-coded stride of +1.
 *
 * You can construct an Interval object using other domain objects.
 * The constructors accept up to 7 domain objects of various types.
 * Domain types are, for example, Loc, Range, Interval. int may also be used
 * in the constructor for Interval; it acts like a Loc<1> object
 * in that context.  The domain arguments for the Interval
 * constructors are combined together to form a single domain object with
 * a dimension equal to the sum of the arguments' dimensions; for example,
 * if you try to create an Interval<3> from a Loc<2> and an Interval<1>, e.g.
 *   Interval<3> a(Loc<2>(1,2), Interval<1>(3,5));
 * the Loc<2> and Interval arguments are combined into a (2+1=3) dimensional
 * domain object, used to initialize the Interval<3>.  The number of dimensions
 * for the arguments must be <= number of dimensions of the newly constructed
 * Interval.
 *
 * For Interval<1>, the list of constructors is limited to just the following:
 *  - Interval<1> a() - default constructor, which creates an EMPTY Interval
 *  - Interval<1> a(n) - sets the Interval to the sequence [0 ... (n-1)]
 *  - Interval<1> a(m,n) - sets the Interval to the sequence [m ... n]
 *
 * The default Interval<1> constructor initializes the Interval to be empty,
 * that is, to have a length() of zero.  In that case, the endpoints are
 * undefined, as is any operation involving the Interval.
 *
 * In addition to the constructors, Interval has the following public
 * interface, similar to all domain objects.  There are two classes of
 * interface methods, one class which includes methods which any Interval<N>
 * object has, regardless of dimensions, the other class which includes extra
 * interface methods that are available for just Interval<1> objects.
 *
 * Interval<N> interface:
 *  - long size() - return the 'volume' of the domain, which is the product
 *      of the lengths of the N 1D Intervals
 *  - bool empty() - return if any of the Interval<1> objects have length == 0
 *  - Interval<1> operator[](int N) - return the Nth Interval<1> in a
 *      multidimensional Interval<M>.  For Interval<1> objects, this just
 *      returns the object back.
 *  - comparison operators: <, >, !=, ==, <=, >= : compare an Interval <N> to
 *      another domain object.  The compared domains must have the same
 *      number of dimensions.
 *  - arithmetic accumulation operators +=, -= : add or subtract in a
 *      given domain.  The added domain must have the same number of
 *      dimensions, or a dimension of 1 (in which case, the same value
 *      is used for all dimensions), and be known to be single-valued (which
 *      is true for Loc and int's).  Note that for Interval, *= and /= are NOT
 *      allowed, since these operators would change the stride and that is
 *      not allowed for Interval (it has a hard-coded stride of 1).
 *      The negation operator (operator-) is also NOT allowed for Interval.
 *  - binary arithmetic operators +, -, *, / : for + and -, adding an Interval
 *      to another Loc or int returns a new Interval.  For * and /, scaling
 *      by a Loc or int returns a Range object, since the stride may
 *      change.
 *  - increment/decrement operator ++, -- : only prefix versions of ++ and --
 *      are provided; they act just like += 1 and -= 1 operations.
 *
 * Interval<1> interface:
 * all of the methods for Interval<N> are also available for Interval<1>.
 *  - int length() - number of elements (including endpoints) of the domain.
 *  - int first() - the beginning endpoint.
 *  - int last() - the ending endpoint.
 *  - int min(), int max() - min or max of the endpoints.
 *  - Interval<1>::iterator begin() and end() - return iterators for the 1D
 *      domain.  These act like (at least) forward iterators.
 *
 * Interval inherits much of its activity from Domain<DomainTraits<Interval>>.
 *
 * For the special case of Interval<1>, there is a specialization given
 * after the general case that has different constructors (listed above).
 */

template<int Dim>
class Interval : public Domain<Dim, DomainTraits<Interval<Dim> > >
{
  // convenience typedefs
  typedef DomainTraits< Interval<Dim> >         DT_t;
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
  Interval() { }

  // copy constructor
  Interval(const Interval<Dim> &a)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain1<Interval<Dim> >::fill(*this, a);
  }

  // Uninitialized constructor
  Interval(const Pooma::NoInit &a)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(a) 
  { }
  
  // templated constructors, taking from 1 to 7 different domain objects
  // (or integers).
  template<class T1>
  explicit Interval(const T1 &a)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  template<class T1, class T2>
  Interval(const T1 &a, const T2 &b)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain2<T1,T2>::fill(*this, a, b);
  }

  template<class T1, class T2, class T3>
  Interval(const T1 &a, const T2 &b, const T3 &c)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain3<T1,T2,T3>::fill(*this, a, b, c);
  }

  template<class T1, class T2, class T3, class T4>
  Interval(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain4<T1,T2,T3,T4>::fill(*this, a, b, c, d);
  }

  template<class T1, class T2, class T3, class T4, class T5>
  Interval(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain5<T1,T2,T3,T4,T5>::fill(*this, a, b, c, d, e);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6>
  Interval(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
      const T6 &f)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain6<T1,T2,T3,T4,T5,T6>::fill(*this, a, b, c, d, e, f);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6, class T7>
  Interval(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
      const T6 &f, const T7 &g)
    : Domain<Dim, DomainTraits<Interval<Dim> > >(Pooma::NoInit()) {
    NewDomain7<T1,T2,T3,T4,T5,T6,T7>::fill(*this, a, b, c, d, e, f, g);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Interval() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  Interval<Dim> &operator=(const T &newdom) {
    return NewDomain1<T>::fill(*this, newdom);
  }

  Interval<Dim> &operator=(const Interval<Dim> &newdom) {
    return NewDomain1<Interval<Dim> >::fill(*this, newdom);
  }

protected:

private:

};


/**
 * Interval<1> is a 1D specialization of Interval<N>; for the 1D case,
 * there are only a restricted set of constructors available.
 * For the special case of Interval<1>, the following constructors
 * are defined:
 *  - Interval<1> a() - default constructor, which creates an EMPTY Interval
 *  - Interval<1> a(n) - sets the Interval to the sequence [0 ... (n-1)]
 *  - Interval<1> a(m,n) - sets the Interval to the sequence [m ... n]
 *  - Interval<1> a(Domain d) - an Interval copied from d, which must be a
 *     1D domain object.
 */

template<>
class Interval<1> : public Domain<1, DomainTraits<Interval<1> > >
{
  // convenience typedef
  typedef DomainTraits< Interval<1> >  DT_t;

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
  Interval() { }

  // copy constructor
  Interval(const Interval<1> &a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    NewDomain1<Interval<1> >::fill(*this, a);
  }

  // Uninitialized constructor
  Interval(const Pooma::NoInit &a)
    : Domain<1, DomainTraits<Interval<1> > >(a) 
  { }

  // general argument constructor, to copy from a different domain type
  template<class T1>
  explicit Interval(const T1 &a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    NewDomain1<T1>::fill(*this, a);
  }

  // initialize from a single scalar.  must specialize to all scalars.
  Interval(char a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(unsigned char a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(short a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(unsigned short a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(int a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(unsigned int a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(long a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }
  Interval(unsigned long a)
    : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
    DomainTraits<Interval<1> >::setDomain(domain_m, 0, a - 1);
  }

  // initialize from a set of endpoints: sets interval to [m ..n].  Must
  // have m <= n.
  template<class T1, class T2>
  Interval(const T1 &m, const T2 &n);

  // initialize from three integers: if the stride is not 1,
  // it is an error.
  template<class T1, class T2, class T3>
  Interval(const T1 &m, const T2 &n, const T3 &s);

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Interval() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  Interval<1> &operator=(const T &newdom) {
    return NewDomain1<T>::fill(*this, newdom);
  }

  Interval<1> &operator=(const Interval<1> &newdom) {
    return NewDomain1<Interval<1> >::fill(*this, newdom);
  }

};

// initialize from a set of endpoints: sets interval to [m ..n].  Must
// have m <= n.
// Not sure why this isn't required everywhere.
#if defined(__MWERKS__)
template<>
#endif
template <class T1, class T2>
inline
Interval<1>::Interval(const T1 &m, const T2 &n)
  : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
  DomainTraits<Interval<1> >::setDomain(domain_m, m, n);
}

// initialize from three integers: if the stride is not 1,
// it is an error.
// Not sure why this isn't required everywhere.
#if defined(__MWERKS__)
template<>
#endif
template <class T1, class T2, class T3>
inline
Interval<1>::Interval(const T1 &m, const T2 &n, const T3 &s)
  : Domain<1, DomainTraits<Interval<1> > >(Pooma::NoInit()) {
  PAssert(s == 1);
  DomainTraits<Interval<1> >::setDomain(domain_m, m, n);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_INTERVAL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Interval.h,v $   $Author: richard $
// $Revision: 1.24 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
