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

#ifndef POOMA_DOMAIN_LOC_H
#define POOMA_DOMAIN_LOC_H

//-----------------------------------------------------------------------------
// Class:
// Loc<int>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Loc<N> is a very simple type of domain, which refers to just one point.
 *
 * It acts very much like an N-dimensional vector of integers.  It can be
 * used to refer to a single point along a sequence of points in a domain.
 * Loc defers most of its implementation to the Domain<DomainTraits<Loc>>
 * base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/DomainTraits.Loc.h"
#include "Domain/NewDomain.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * fillLocStorage descriptions:
 *   fillLocStorage is a function (actually a set of overloaded functions)
 * that copies data from a given domain into a Loc.  It will modify the
 * second argument (the loc) starting at the index given by the first
 * argument.  The first argument will be incremented by the number of
 * dimensions filled in to the Loc.  This is defines for 1 ... 7 arguments.
 */

template<int Dim, class T1>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a)
{
  PAssert(currIndex >= 0 && (currIndex + DomainTraits<T1>::dimensions) <= Dim);
  for (int i=0; i < DomainTraits<T1>::dimensions; ++i)
    loc[currIndex++].setDomain(DomainTraits<T1>::getPointDomain(a, i));
  return currIndex;
}

template<int Dim, class T1, class T2>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a, const T2 &b)
{
  currIndex = fillLocStorage(currIndex, loc, a);
  currIndex = fillLocStorage(currIndex, loc, b);
  return currIndex;
}

template<int Dim, class T1, class T2, class T3>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a, const T2 &b,
		   const T3 &c)
{
  currIndex = fillLocStorage(currIndex, loc, a, b);
  currIndex = fillLocStorage(currIndex, loc, c);
  return currIndex;
}

template<int Dim, class T1, class T2, class T3, class T4>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a, const T2 &b,
		   const T3 &c, const T4 &d)
{
  currIndex = fillLocStorage(currIndex, loc, a, b, c);
  currIndex = fillLocStorage(currIndex, loc, d);
  return currIndex;
}

template<int Dim, class T1, class T2, class T3, class T4, class T5>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a, const T2 &b,
		   const T3 &c, const T4 &d, const T5 &e)
{
  currIndex = fillLocStorage(currIndex, loc, a, b, c, d);
  currIndex = fillLocStorage(currIndex, loc, e);
  return currIndex;
}

template<int Dim, class T1, class T2, class T3, class T4, class T5, class T6>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a, const T2 &b,
		   const T3 &c, const T4 &d, const T5 &e, const T6 &f)
{
  currIndex = fillLocStorage(currIndex, loc, a, b, c, d, e);
  currIndex = fillLocStorage(currIndex, loc, f);
  return currIndex;
}

template<int Dim, class T1, class T2, class T3, class T4, class T5, class T6,
         class T7>
inline
int fillLocStorage(int currIndex, Loc<Dim> &loc, const T1 &a, const T2 &b,
		   const T3 &c, const T4 &d, const T5 &e, const T6 &f,
		   const T7 &g)
{
  currIndex = fillLocStorage(currIndex, loc, a, b, c, d, e, f);
  currIndex = fillLocStorage(currIndex, loc, g);
  return currIndex;
}


/**
 * Full description of CopyLocStorage:
 *   CopyLocStorage is a simple struct with one static member 'copy' that
 * copies data out of a given single domain into the given Loc.  If the
 * data is greater than 1D, only the dimensions of the given domain are copied.
 * If the data is 1D, then the value of the domain is copied into all 'Dim'
 * dimensions of the Loc.  So there is a 1D specialization of CopyLocStorage.
 */

template<int Dim, class T, int DimT, bool wildcard>
struct CopyLocStorageImpl
{
  // for wildcard initializers, nothing to do, just leave uninitialized
  inline
  static void copy(Loc<Dim> &, const T &) { }
};

template<int Dim, class T, int DimT>
struct CopyLocStorageImpl<Dim, T, DimT, false>
{
  inline
  static void copy(Loc<Dim> &loc, const T &a) {
    CTAssert(DomainTraits<T>::dimensions == DimT);
    CTAssert(DomainTraits<T>::dimensions <= Dim);
    fillLocStorage(0, loc, a);
  }
};

template<int Dim, class T>
struct CopyLocStorageImpl<Dim, T, 1, false>
{
  inline
  static void copy(Loc<Dim> &loc, const T &a) {
    CTAssert(DomainTraits<T>::dimensions == 1);
    for (int i=0; i < Dim; ++i)
      loc[i].setDomain(DomainTraits<T>::getPointDomain(a, 0));
  }
};

template<int Dim, class T>
struct CopyLocStorage
{
  inline
  static void copy(Loc<Dim> &loc, const T &a) {
    CopyLocStorageImpl<Dim, T, DomainTraits<T>::dimensions,
                       DomainTraits<T>::wildcard>::copy(loc, a);
  }
};



/**
 * Loc<N> is a domain representing a single N-dimensional point.  It has
 * a stride of one, endpoints which are the same, and a single element.
 * Otherwise, it acts very much like all other domain objects such as Range
 * or Interval.
 *
 * You can construct a Loc object using other Loc objects or integers.
 * The constructor for Loc accepts up to 7 domain objects of various types.
 * Domain types are, for example, Loc, Range, Interval. int may also be used
 * in the constructor for Loc; it acts like a Loc<1> object in that context.
 * It is illegal to try to construct a Loc from other domain objects which
 * have lengths greater than one.  The domain arguments for the Loc
 * constructors are combined together to form a single domain object with
 * a dimension equal to the sum of the arguments dimensions; for example,
 * if you try to create a Loc<3> from a Loc<2> and an int, as in
 *   Loc<3> a(Loc<2>(1,2), 3);
 * the Loc<2> and int arguments are combined into a (2+1=3) dimensional
 * domain object, used to initialize the Loc<3>.  The number of dimensions
 * for the arguments must be <= number of dimensions of the newly constructed
 * Loc.
 *
 * The default Loc<1> constructor initializes the point to zero; the default
 * Loc<N> constructor initializes a set of Loc<1> objects all to zero.
 *
 * In addition to the constructors, Loc has the following public interface,
 * similar to all domain objects.  There are two classes of interface
 * methods, one class which includes methods which any Loc<N> object has,
 * regardless of dimensions, the other class which includes extra interface
 * methods that are available for just Loc<1> objects.
 *
 * Loc<N> interface:
 *  - long size() - return the 'volume' of the domain; for Loc, this is 1.
 *  - bool empty() - return if any of the Loc<1> objects have length == 0
 *  - Loc<1> operator[](int N) - return the Nth Loc<1> in a multidimensional
 *      Loc<M>.  For Loc<1> objects, this just returns the object back.
 *  - comparison operators: <, >, !=, ==, <=, >= : compare a Loc<N> to
 *      another domain object.  The compared domains must have the same
 *      number of dimensions.
 *  - arithmetic accumulation operators +=, -= : add or subtract in a
 *      given domain.  The added domain must have the same number of
 *      dimensions, or a dimension of 1 (in which case, the same value
 *      is used for all dimensions), and be known to be single-valued (which
 *      is true for Loc and int's). 
 *  - binary arithmetic operators +, -, *, / : for + and -, adding a Loc
 *      to another Loc or int returns a new Loc object.  For * and /, scaling
 *      by a Loc or int returns a Range object, since the stride may
 *      change.
 *  - increment/decrement operator ++, -- : only prefix versions of ++ and --
 *      are provided; they act just like += 1 and -= 1 operations.
 *
 * Loc<1> interface:
 * all of the methods for Loc<N> are also available for Loc<1>.  Plus:
 *  - int length() - number of elements (including endpoints) of the domain.
 *      For Loc<1>, this is always 1.
 *  - int first() - the beginning endpoint, for Loc<1> just the point itself.
 *  - int last() - the ending endpoint, for Loc<1> just the point itself.
 *  - int min(), int max() - min or max of the endpoints.
 *  - Loc<1>::iterator begin() and end() - return iterators for the 1D domain.
 *      These act like (at least) forward iterators.
 *
 * Loc inherits much of its activity from Domain<DomainTraits<Loc>,Dim>.
 * Domain is a base class that uses the template argument as a traits class
 * for Loc to specialize it's behavior.  The definition of Loc<N> itself
 * only contains Loc<N> constructors, destructor, and operator= methods.
 * All other interface methods are in Domain, or in DomainBase (from
 * which Domain inherits).
 */

template<int Dim>
class Loc : public Domain<Dim, DomainTraits<Loc<Dim> > >
{
  // convenience typedefs
  typedef DomainTraits< Loc<Dim> >              DT_t;
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

  // default constructor: initialize the points to zero.
  inline
  Loc() { }

  // copy constructor
  inline
  Loc(const Loc<Dim> &a)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    fillLocStorage(0, *this, a);
  }

  // Uninitialized constructor
  Loc(const Pooma::NoInit &a)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(a) 
  { }

  // templated constructors, taking from 1 to 7 different domain objects
  // (or integers).
  template<class T1>
  inline
  explicit Loc(const T1 &a)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CopyLocStorage<Dim, T1>::copy(*this, a);
  }

  template<class T1, class T2>
  inline
  Loc(const T1 &a, const T2 &b)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CTAssert(Dim >= (DomainTraits<T1>::dimensions +
		     DomainTraits<T2>::dimensions));
    fillLocStorage(0, *this, a, b);
  }

  template<class T1, class T2, class T3>
  inline
  Loc(const T1 &a, const T2 &b, const T3 &c)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CTAssert(Dim >= (DomainTraits<T1>::dimensions +
		     DomainTraits<T2>::dimensions +
		     DomainTraits<T3>::dimensions));
    fillLocStorage(0, *this, a, b, c);
  }

  template<class T1, class T2, class T3, class T4>
  inline
  Loc(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CTAssert(Dim >= (DomainTraits<T1>::dimensions +
		     DomainTraits<T2>::dimensions +
		     DomainTraits<T3>::dimensions +
		     DomainTraits<T4>::dimensions));
    fillLocStorage(0, *this, a, b, c, d);
  }

  template<class T1, class T2, class T3, class T4, class T5>
  inline
  Loc(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CTAssert(Dim >= (DomainTraits<T1>::dimensions +
		     DomainTraits<T2>::dimensions +
		     DomainTraits<T3>::dimensions +
		     DomainTraits<T4>::dimensions +
		     DomainTraits<T5>::dimensions));
    fillLocStorage(0, *this, a, b, c, d, e);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6>
  inline
  Loc(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
      const T6 &f)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CTAssert(Dim >= (DomainTraits<T1>::dimensions +
		     DomainTraits<T2>::dimensions +
		     DomainTraits<T3>::dimensions +
		     DomainTraits<T4>::dimensions +
		     DomainTraits<T5>::dimensions +
		     DomainTraits<T6>::dimensions));
    fillLocStorage(0, *this, a, b, c, d, e, f);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6, class T7>
  inline
  Loc(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
      const T6 &f, const T7 &g)
    : Domain<Dim, DomainTraits<Loc<Dim> > >(Pooma::NoInit()) {
    CTAssert(Dim >= (DomainTraits<T1>::dimensions +
		     DomainTraits<T2>::dimensions +
		     DomainTraits<T3>::dimensions +
		     DomainTraits<T4>::dimensions +
		     DomainTraits<T5>::dimensions +
		     DomainTraits<T6>::dimensions +
		     DomainTraits<T7>::dimensions));
    fillLocStorage(0, *this, a, b, c, d, e, f, g);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  inline
  ~Loc() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  inline
  Loc<Dim> &operator=(const T &newdom) {
    CopyLocStorage<Dim, T>::copy(*this, newdom);
    return *this;
  }

  inline
  Loc<Dim> &operator=(const Loc<Dim> &newdom) {
    fillLocStorage(0, *this, newdom);
    return *this;
  }

  // Override the standard DomainBase print() function with something
  // appropriate for Loc<>, which doesn't have separate first/last values, and
  // doesn't have a stride: print a domain to a stream, in the format "["
  // first, first, ... first "]"

  template<class Out>
  void print(Out &o) const {
    const Domain_t &d = this->unwrap();
    o << "[";
    for (int i=0; i < Dim; ++i) {
      o << d[i].first();
      if (i < (Dim-1))
        o << ",";
    }
    o << "]";
  }

protected:

private:

};


/**
 * Loc<1> is a 1D specialization of Loc<N>; for the 1D case,
 * there are only a restricted set of constructors available.
 * For the special case of Loc<1>, the following constructors
 * are defined:
 *  - Loc<1> a() - default constructor, which creates a Loc<1> set to zero.
 *  - Loc<1> a(n) - sets the Loc<1> to the point n
 */

template<>
class Loc<1> : public Domain<1, DomainTraits<Loc<1> > >
{
  // convenience typedef
  typedef DomainTraits< Loc<1> >  DT_t;

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
  inline
  Loc() { }

  // copy constructor
  inline
  Loc(const Loc<1> &a)
    : Domain<1, DomainTraits<Loc<1> > >(Pooma::NoInit()) {
    fillLocStorage(0, *this, a);
  }

  // Uninitialized constructor
  Loc(const Pooma::NoInit &a)
    : Domain<1, DomainTraits<Loc<1> > >(a) 
  { }

  // general argument constructor, to copy from a different domain type.
  // This sets the Loc to the given point.
  template<class T1>
  inline
  explicit Loc(const T1 &a)
    : Domain<1, DomainTraits<Loc<1> > >(Pooma::NoInit()) {
    CopyLocStorage<1, T1>::copy(*this, a);
  }

  // initialize from two integers: if they are not equal, it is an error.
  template<class T1, class T2>
  inline
  Loc(const T1 &a, const T2 &b)
    : Domain<1, DomainTraits<Loc<1> > >(Pooma::NoInit()) {
    CTAssert(DomainTraits<T1>::dimensions == 1);
    CTAssert(DomainTraits<T2>::dimensions == 1);
    CTAssert(DomainTraits<T1>::singleValued);
    CTAssert(DomainTraits<T2>::singleValued);
    PAssert(a == b);
    fillLocStorage(0, *this, a);
  }

  // initialize from three integers: if the first two are not equal,
  // it is an error.
  template<class T1, class T2, class T3>
  inline
  Loc(const T1 &a, const T2 &b, const T3 &c)
    : Domain<1, DomainTraits<Loc<1> > >(Pooma::NoInit()) {
    CTAssert(DomainTraits<T1>::dimensions == 1);
    CTAssert(DomainTraits<T2>::dimensions == 1);
    CTAssert(DomainTraits<T3>::dimensions == 1);
    CTAssert(DomainTraits<T1>::singleValued);
    CTAssert(DomainTraits<T2>::singleValued);
    CTAssert(DomainTraits<T3>::singleValued);
    PAssert(a == b);
    fillLocStorage(0, *this, a);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  inline
  ~Loc() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  inline
  Loc<1> &operator=(const T &newdom) {
    CopyLocStorage<1, T>::copy(*this, newdom);
    return *this;
  }

  inline
  Loc<1> &operator=(const Loc<1> &newdom) {
    fillLocStorage(0, *this, newdom);
    return *this;
  }

  template<class Out>
  void print(Out &o) const {
    const Domain_t &d = this->unwrap();
    o << "[";
    o << d[0].first();
    o << "]";
  }
};


/// print a Loc to a stream, in the format
///   "[" first, first, ... first "]"
/// This overrides the more general function in DomainBase.h

template<int Dim>
std::ostream& operator<<(std::ostream &o, const Loc<Dim> &loc)
{
  loc.print(o);
  return o;
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_LOC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Loc.h,v $   $Author: richard $
// $Revision: 1.30 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
