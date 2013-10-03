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

#ifndef POOMA_DOMAIN_GRID_H
#define POOMA_DOMAIN_GRID_H

//-----------------------------------------------------------------------------
// Class:
// Grid<int Dim>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Grid is a general type of integer domain, which refers to a set of points
 * a0, a1, ... aN for each dimension.
 *
 * The points can be any ascending or
 * descending sequence, there is no fixed stride.  This is basically a set
 * of Dim IndirectionList<int>'s, one for each dimension; the total domain
 * is the tensor product of these lists.  Grid<Dim> is basically an array
 * of Grid<1> objects.
 *
 * Grid defers most of its implementation to the Domain<DomainTraits<Grid>>
 * base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/DomainTraits.Grid.h"
#include "Domain/NewDomain.h"
#include "Domain/Loc.h"        // needed for use of operator<<
#include <iosfwd>


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * Grid<N> is a domain representing a set of N numeric sequences, one
 * for each dimension N.  The sequences are lists of ascending or descending
 * numbers, but without any fixed stride - the list for each dimension is
 * an IndirectionList<int>.
 *
 * You can construct a Grid object using other domain objects.
 * The constructors accept up to 7 domain objects of various types.
 * Domain types are, for example, Loc, Grid, Interval. int may also be used
 * in the constructor for Grid; it acts like a Loc<1> object
 * in that context.  The domain arguments for the Grid
 * constructors are combined together to form a single domain object with
 * a dimension equal to the sum of the argument's dimensions; for example,
 * if you try to create a Grid<3> from a Loc<2> and an Interval<1>, e.g.
 *   Grid<3> a(Loc<2>(1,2), Interval<1>(3,5));
 * the Loc<2> and Interval arguments are combined into a (2+1=3) dimensional
 * domain object, used to initialize the Grid<3>.  The number of dimensions
 * for the arguments must be <= number of dimensions of the newly constructed
 * Grid.
 *
 * Grid, unlike other domain objects, can also be constructed with
 * IndirectionList objects.  IndirectionList's look like 1D domain objects
 * to the constructor of Grid, so multiple lists can be used.  Grid's can
 * also be used in this same way to construct other Grid's.
 *
 * For Grid<1>, the list of constructors includes the following:
 *  - Grid<1> a() - default constructor, which creates an EMPTY Grid
 *  - Grid<1> a(n) - sets the Grid to the sequence [0 ... (n-1)], stride 1
 *  - Grid<1> a(m,n) - sets the Grid to the sequence [m ... n], stride 1 or -1
 *  - Grid<1> a(m,n,s) - sets the Grid to the sequence [m ... n], stride s
 *
 * The default Grid<1> constructor initializes the Grid to be empty,
 * that is, to have a length() of zero.  In that case, the endpoints are
 * undefined, as is any operation involving the Grid.
 *
 * In addition to the constructors, Grid has the following public
 * interface, similar to all domain objects.  There are two classes of
 * interface methods, one class which includes methods which any Grid<N>
 * object has, regardless of dimensions, the other class which includes extra
 * interface methods that are available for just Grid<1> objects.
 *
 * Grid<N> interface:
 *  - long size() - return the 'volume' of the domain, which is the product
 *      of the lengths of the N 1D Grids
 *  - bool empty() - return if any of the Grid<1> objects have length == 0
 *  - Grid<1> operator[](int N) - return the Nth Grid<1> in a
 *      multidimensional Grid<M>.  For Grid<1> objects, this just
 *      returns the object back.
 *  - comparison operators: <, >, !=, ==, <=, >= : compare a Grid<N> to
 *      another domain object.  The compared domains must have the same
 *      number of dimensions.
 *  - arithmetic accumulation operators +=, -=, *=, /= : add or subtract in a
 *      given domain.  The added domain must have the same number of
 *      dimensions, or a dimension of 1 (in which case, the same value
 *      is used for all dimensions), and be known to be single-valued (which
 *      is true for Loc and int's).  Note that for Grid, *= and /= ARE
 *      allowed, since Grid can have its stride changed at run time.  *=
 *      and /= result in scaling of the endpoints and stride, which leaves
 *      the length (and size) the same.  += and -= shift the beginning
 *      endpoints by the given values, also leaving the length and size the
 *      same.  Negation of a Grid negates the endpoints and stride.
 *  - binary arithmetic operators +, -, *, / : for + and -, adding a Grid
 *      to another Loc or int returns a new Grid.  For * and /, scaling
 *      by a Loc or int also returns a Grid object, since the stride may
 *      change.
 *  - increment/decrement operator ++, -- : only prefix versions of ++ and --
 *      are provided; they act just like += 1 and -= 1 operations.
 *
 * Grid<1> interface:
 * all of the methods for Grid<N> are also available for Grid<1>. Plus:
 *  - int length() - number of elements (including endpoints) of the domain.
 *     for a non-unit-stride Grid, the length refers to the number of
 *     strided points (including the endpoints), NOT the difference between
 *     the first and last endpoint.  That is, length = (end-beg)/stride + 1,
 *     NOT (end-beg) + 1.
 *  - int first() - the beginning endpoint.
 *  - int last() - the ending endpoint.
 *  - int min(), int max() - min or max of the endpoints.
 *  - Grid<1>::iterator begin() and end() - return iterators for the 1D
 *      domain.  These act like (at least) forward iterators.
 *
 * For the special case of Grid<1>, there is a specialization given
 * after the general case that has different constructors.
 *
 * Grid inherits much of its activity from Domain<DomainTraits<Grid>>
 */

template<int Dim>
class Grid : public Domain<Dim, DomainTraits<Grid<Dim> > >
{
  // convenience typedefs
  typedef DomainTraits< Grid<Dim> >             DT_t;
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
  Grid() { }

  // copy constructor
  Grid(const Grid<Dim> &a) {
    NewDomain1<Grid<Dim> >::fill(*this, a);
  }

  // templated constructors, taking from 1 to 7 different domain objects
  // (or integers).
  template<class T1>
  explicit Grid(const T1 &a) {
    NewDomain1<T1>::fill(*this, a);
  }

  template<class T1, class T2>
  Grid(const T1 &a, const T2 &b) {
    NewDomain2<T1,T2>::fill(*this, a, b);
  }

  template<class T1, class T2, class T3>
  Grid(const T1 &a, const T2 &b, const T3 &c) {
    NewDomain3<T1,T2,T3>::fill(*this, a, b, c);
  }

  template<class T1, class T2, class T3, class T4>
  Grid(const T1 &a, const T2 &b, const T3 &c, const T4 &d) {
    NewDomain4<T1,T2,T3,T4>::fill(*this, a, b, c, d);
  }

  template<class T1, class T2, class T3, class T4, class T5>
  Grid(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e) {
    NewDomain5<T1,T2,T3,T4,T5>::fill(*this, a, b, c, d, e);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6>
  Grid(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
       const T6 &f) {
    NewDomain6<T1,T2,T3,T4,T5,T6>::fill(*this, a, b, c, d, e, f);
  }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6, class T7>
  Grid(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
       const T6 &f, const T7 &g) {
    NewDomain7<T1,T2,T3,T4,T5,T6,T7>::fill(*this, a, b, c, d, e, f, g);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Grid() { 
    for (int i=0;i<Dim;++i)
      this->operator[](i).~OneDomain_t();
  }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  Grid<Dim> &operator=(const T &newdom) {
    return NewDomain1<T>::fill(*this, newdom);
  }

  Grid<Dim> &operator=(const Grid<Dim> &newdom) {
    return NewDomain1<Grid<Dim> >::fill(*this, newdom);
  }

  //
  // I/O/
  //

  // print a Grid<N> to a stream, in the format
  //   "[" value1,value2,...,valueN "]"

  template<class Out>
  void print(Out &o) const {
    iterator p    = this->begin();
    iterator pend = this->end();
    o << "[";				
    while (p != pend)			
      {					
	o << *p;				
	++p;				
	if (p != pend)			
	  o << ",";			
      }					
    o << "]";				
  }

protected:

private:

};


/**
 * Grid<1> is a 1D specialization of Grid<N>; for the 1D case,
 * there are only a restricted set of constructors available.
 * For the special case of Grid<1>, the following constructors
 * are defined:
 *  - Grid<1> a() - default constructor, which creates an EMPTY Grid
 *  - Grid<1> a(n) - sets the Grid to the sequence [0 ... (n-1)], stride 1
 *  - Grid<1> a(m,n) - sets the Grid to the sequence [m ... n], stride 1 or -1
 *  - Grid<1> a(m,n,s) - sets the Grid to the sequence [m ... n], stride s
 *  - Grid<1> a(Domain d) : a Grid copied from d, which must be a
 *     1D domain object.
 */

template<>
class Grid<1> : public Domain<1, DomainTraits<Grid<1> > >
{
  // convenience typedef
  typedef DomainTraits< Grid<1> >  DT_t;

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
  Grid() { }

  // copy constructor
  Grid(const Grid<1> &a) {
    NewDomain1<Grid<1> >::fill(*this, a);
  }

  // general argument constructor, to copy from a different domain type
  template<class T1>
  explicit Grid(const T1 &a) {
    NewDomain1<T1>::fill(*this, a);
  }

  // initialize from a single scalar.  must specialize to all scalars.
  Grid(char a) {
    PAssert(a != 0);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - 1);
  }
  Grid(unsigned char a) {
    PAssert(a != 0);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - 1);
  }
  Grid(short a) {
    PAssert(a != 0);
    short s = (a < 0 ? -1 : 1);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - s);
  }
  Grid(unsigned short a) {
    PAssert(a != 0);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - 1);
  }
  Grid(int a) {
    PAssert(a != 0);
    int s = (a < 0 ? -1 : 1);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - s);
  }
  Grid(unsigned int a) {
    PAssert(a != 0);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - 1);
  }
  Grid(long a) {
    PAssert(a != 0);
    long s = (a < 0 ? -1 : 1);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - s);
  }
  Grid(unsigned long a) {
    PAssert(a != 0);
    DomainTraits<Grid<1> >::setDomain(domain_m, 0, a - 1);
  }

  // initialize from a set of endpoints: sets range to [m ..n].
  // domain_m is the domain information storage kept in the base class.
  template<class T1, class T2>
  Grid(const T1 &m, const T2 &n);

  // initialize from a set of endpoints and with a given stride.
  // domain_m is the domain information storage kept in the base class.
  template<class T1, class T2, class T3>
  Grid(const T1 &m, const T2 &n, const T3 &s);

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~Grid() { }

  //
  // operator=, templated to allow assignment from other domain objects
  // this uses the same mechanism as the constructors to fill in to this
  // object the data from the given object.  If the new object has too
  // few dimensions, this will only change the first M dimensions of this
  // object, where M is the number of dimensions for newdom
  //

  template<class T>
  Grid<1> &operator=(const T &newdom) {
    return NewDomain1<T>::fill(*this, newdom);
  }

  Grid<1> &operator=(const Grid<1> &newdom) {
    return NewDomain1<Grid<1> >::fill(*this, newdom);
  }

  //
  // A special function used to initialize one Grid from another.  For
  // this, we need non-modifiable access to the storage.
  //

  const Storage_t &storage() const { return domain_m; }

  //
  // I/O/
  //

  // print a Grid<N> to a stream, in the format
  //   "[" value1,value2,...,valueN "]"

  template<class Out>
  void print(Out &o) const {
    iterator p    = begin();
    iterator pend = end();
    o << "[";				
    while (p != pend)			
      {					
	o << *p;				
	++p;				
	if (p != pend)			
	  o << ",";			
      }					
    o << "]";				
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
Grid<1>::Grid(const T1 &m, const T2 &n) {
  DomainTraits<Grid<1> >::setDomain(domain_m, m, n);
}

// initialize from a set of endpoints and with a given stride.
// domain_m is the domain information storage kept in the base class.
// Not sure why this isn't required everywhere.
#if defined(__MWERKS__)
template<>
#endif
template <class T1, class T2, class T3>
inline
Grid<1>::Grid(const T1 &m, const T2 &n, const T3 &s) {
  DomainTraits<Grid<1> >::setDomain(domain_m, m, n, s);
}


/// print a Grid<N> to a stream, in the format
///   "[" value1,value2,...,valueN "]"

template<int Dim>
std::ostream& operator<<(std::ostream &o, const Grid<Dim> &grid)
{
  grid.print(o);
  return o;
}

//////////////////////////////////////////////////////////////////////
//
// Specialization of the CHEETAH Serialize class for Grid<1>.
//
//////////////////////////////////////////////////////////////////////

#if POOMA_MESSAGING

#include "Tulip/Messaging.h"

namespace Cheetah {

template<>
class Serialize<CHEETAH, Grid<1> >
{
public:

  typedef Grid<1>                    Grid_t;
  typedef Grid_t::Element_t          Element_t;
  typedef IndirectionList<Element_t> List_t;

  static inline long
  size(const Grid_t &a)
  {
    return sizeof(int) + a.length() * sizeof(Element_t);
  }

  static inline int
  pack(const Grid_t &a, char *buffer)
  {
    *reinterpret_cast<int *>(buffer) = a.length();

    long length = a.length() * sizeof(Element_t);

    memcpy(buffer + sizeof(int), &a.storage()(0), length);
    
    return sizeof(int) + length;
  }

  static inline int
  unpack(Grid_t* &a, char *buffer)
  {
    int length = *reinterpret_cast<int *>(buffer);

    List_t list(length);

    length *= sizeof(Element_t);

    memcpy(&list(0), buffer + sizeof(int), length);

    a = new Grid_t(list);

    return sizeof(int) + length;
  }

  static inline void
  cleanup(Grid_t* a)
  {
    delete a;
  }
};

} // namespace Cheetah

#endif     // POOMA_MESSAGING

#endif     // POOMA_DOMAIN_GRID_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Grid.h,v $   $Author: richard $
// $Revision: 1.19 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
