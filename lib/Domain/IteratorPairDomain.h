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
// Classes
//   IteratorPairDomain<Iter>
//-----------------------------------------------------------------------------

#ifndef POOMA_DOMAIN_ITERATOR_PAIR_DOMAIN_H
#define POOMA_DOMAIN_ITERATOR_PAIR_DOMAIN_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"

#include <iterator>
#include <iosfwd>

#if POOMA_NONSTANDARD_ITERATOR
//-----------------------------------------------------------------------------
// Hack alert!!!
// Microsoft doesn't provide these specializations (not surprising since VC++ 
// doesn't do partial specializations) and Intel just uses the Microsoft 
// headers. We need this in the code below since we use int* as an iterator
// in some of the code that uses IteratorPairDomain. 
//-----------------------------------------------------------------------------

namespace std {
  // Copied from the SGI STL
  template <class _Tp>
  struct iterator_traits<_Tp*> {
    typedef random_access_iterator_tag iterator_category;
    typedef _Tp                         value_type;
    typedef ptrdiff_t                   difference_type;
    typedef _Tp*                        pointer;
    typedef _Tp&                        reference;
  };

  template <class _Tp>
  struct iterator_traits<const _Tp*> {
    typedef random_access_iterator_tag iterator_category;
    typedef _Tp                         value_type;
    typedef ptrdiff_t                   difference_type;
    typedef const _Tp*                  pointer;
    typedef const _Tp&                  reference;
  };
} // namespace std

#endif

namespace Pooma {

/** @file
 * @ingroup Domain
 * @brief
 * IteratorPairDomain<Iter> wraps a pair of iterators and provides a
 * subset of services provided by other 1D domains, and in particular by
 * IndirectionList.
 *
 * It is meant to be used in places where
 * IndirectionList might be used, but it should be cheaper since it is
 * just storing a pair of iterators. 
 */

/**
 * Like IndirectionList<T>, IteratorPairDomain<T_Iter> provides access
 * to an arbitrary sequence of T values. Rather than storing these in
 * a DataBlockPtr, IteratorPairDomain simply stores a [begin,end) pair
 * of iterators from some container. This is obviously cheaper, more
 * flexible, and less coupled to other types (IndirectionList is
 * implemented with a DataBlockPtr, which is not a particularly nice
 * dependency). When using asynchronous threads for evaluation, one
 * may need to take care with the type of iterators used. For
 * instance, since IndirectionList uses a DataBlockPtr, its data is
 * reference counted and will not go away if the original copy goes
 * out of scope. This would not be true for, say,
 * IteratorPairDomain<int*>. On the other hand, one can accomplish
 * basically the same result wiht IteratorPairDomain<DataBlockPtr>,
 * since DataBlockPtr is a valid iterator type.
 *
 * Unlike IndirectionList, IteratorPairDomain does not currently
 * support various arithmetic operations or various domain-domain
 * operations. These may be added at a future date, though I think
 * domain modifying operations are very dangerous when combined with
 * shallow copy semantics. (e.g., IndirectionList does support such
 * operations. However, things like
 *
 *     I1 = I2; 
 *     I2 += 1;
 *
 * don't do the same thing for IndirectionLists as they do for
 * Interval and Range since the second statement modifies I2. This is
 * A Bad Thing (TM,IMHO,JAC)).
 *
 * Note that many uses of IteratorPairDomains (and IndirectionLists)
 * assume that the points are ordered. However, there is nothing about
 * the abstraction that dictates this - it is simply a usage pattern.
 *
 * IteratorPairDomain<Iter> interface:
 *
 * Types:
 *  - This_t             - IteratorPairDomain<Iter>
 *  - IterTraits_t       - std::iterator_traits<Iter>
 *  - Element_t          - type returned by *Iter if return by value
 *  - ElementRef_t       - type returned by *Iter if return by reference
 *  - Size_t             - integral type used for sizes
 *  - iterator           - Iterator type (Iter)
 *
 * Constants:
 *  - dimensions         - number of dimensions (1)
 *  - loopAware          - is it loop aware? (false [n/a])
 *  - singleValued       - is it a single value (false)
 *  - unitStride         - is it unit stride (false)
 *
 * Constructors:
 *  - IteratorPairDomain()     
 *  - IteratorPairDomain(Iter begin, Iter end)
 *
 * Member functions:
 *  - This_t operator[](int) 
 *                        returns *this - for compatibilitiy with 
 *                        multidimensional domains.
 *  - ElementRef_t operator()(Size_t i) 
 *  - Element_t    operator()(Size_t i) const
 *                        returns the i'th element of the underlying
 *                        container. 
 *  - Size_t length()    - number of elements of the domain. 
 *  - Size_t size()      - same as above. 
 *  - bool empty()       - true if size() == 0
 *  - bool initialized() - true if size() > 0
 *  - Element_t first()  - the first point in the domain.
 *  - Element_t last()   - the last point in the domain.
 *  - Iter begin()       - iterator to start of sequence.
 *  - Iter end()         - iterator to the end of the sequence.
 *  - Element_t max()    - returns the largest element (assumes operator<())
 *  - Element_t min()    - returns the smallest element (assumes operator<())
 */

template <class Iter>
class IteratorPairDomain
{
public:
  //---------------------------------------------------------------------------
  // Typedefs and class constants.
  //---------------------------------------------------------------------------

  typedef IteratorPairDomain<Iter>            This_t;
  typedef Iter                                iterator;
  typedef std::iterator_traits<Iter>          IterTraits_t;
  typedef typename IterTraits_t::value_type   Element_t;
  typedef typename IterTraits_t::reference    ElementRef_t;
  typedef long                                Size_t;

  enum { dimensions   = 1 };
  enum { loopAware    = false };
  enum { singleValued = false };
  enum { unitStride   = false };

  //---------------------------------------------------------------------------
  // Constructors, destructor, etc.
  //---------------------------------------------------------------------------

  // Default constructor : initialize to an empty domain

  IteratorPairDomain() : size_m(0) { }

  // Primary constructor taking two iterators:

  IteratorPairDomain(Iter begin, Iter end)
    : begin_m(begin), end_m(end)
  {
    size_m = std::distance(begin_m,end_m);
  }

  // Copy constructors

  IteratorPairDomain(const This_t &a)
    : begin_m(a.begin_m), end_m(a.end_m), size_m(a.size_m)
  { }

  template <class OtherIter>
  IteratorPairDomain(const IteratorPairDomain<OtherIter> &a)
    : begin_m(a.begin()), end_m(a.end()), size_m(a.size())
  { }

  // Copy assignment.

  This_t &operator=(const This_t &model) 
  {
    begin_m = model.begin_m;
    end_m   = model.end_m;
    size_m  = model.size_m;
    return *this;
  }

  //---------------------------------------------------------------------------
  // Basic domain operations.
  //---------------------------------------------------------------------------

  This_t &operator[](int i)
  {
    PAssert(i == 0);
    return *this;
  }

  const This_t &operator[](int i) const
  {
    PAssert(i == 0);
    return *this;
  }

  // I tried delegating this to a templated member function that
  // avoided the temporary (i.e. that just did "return *(begin_m+i);"
  // for random-access iterators), but CodeWarrior didn't like the
  // specialization syntax, and Intel just plain ignored the
  // specialization.

  Element_t operator()(Size_t i) const 
  {
    PAssert(i >= 0 && i < size_m);
    Iter pos = begin_m;
    std::advance(pos,i);
    return *pos;
  }

  ElementRef_t operator()(Size_t i) 
  { 
    PAssert(i >= 0 && i < size_m);
    Iter pos = begin_m;
    std::advance(pos,i);
    return *pos;
  }

  Size_t length() const    { return size_m; }
  Size_t size() const      { return size_m; }
  bool empty() const       { return size_m == 0; }
  bool initialized() const { return size_m != 0; }

  Element_t first() const  { return *begin_m; }

  Element_t last() const
  {
    Iter lpos = begin_m;
    std::advance(lpos,size_m-1);
    return *lpos;
  }

  Element_t min() const
  {
    Iter pos = begin_m;
    Element_t result = *pos++;
    for (Size_t i = 1; i < size_m; ++i)
      {
        if (*pos < result) result = *pos;
        ++pos;
      }                               
    return result;
  }

  Element_t max() const
  {
    Iter pos = begin_m;
    Element_t result = *pos++;
    for (Size_t i = 1; i < size_m; ++i)
      {
        if (result < *pos) result = *pos;
        ++pos;
      }
    return result;
  }

  //---------------------------------------------------------------------------
  // Iterator production
  //---------------------------------------------------------------------------

  Iter begin() const { return begin_m; }
  Iter end()   const { return end_m; }

  //---------------------------------------------------------------------------
  // I/O
  //---------------------------------------------------------------------------

  // Print an IteratorPairDomain to a stream, in the format
  //   "[" val1, val2, ... , valN "]"

  template <class Out>
  void print(Out &o) const
  {
    o << "[";

    Iter pos = begin_m;
    o << *pos++;
    for (Size_t i = 1; i < size_m; ++i)
      {
        o << "," << *pos++;
      }
    o << "]";
  }

private:
  //---------------------------------------------------------------------------
  // Data
  //---------------------------------------------------------------------------

  // We don't manage the storage - rather we just manage a pair
  // of iterators into someone elses storage:

  Iter begin_m;
  Iter end_m;

  // Number of elements (for convenience)

  Size_t size_m;
};

//-----------------------------------------------------------------------------
// operator<<
//   - ostream inserter - simply calls the print method defined above.
//-----------------------------------------------------------------------------

template <class Iter>
std::ostream& 
operator<<(std::ostream &o, const IteratorPairDomain<Iter> &list)
{
  list.print(o);
  return o;
}

} // namespace Pooma

#endif  // POOMA_DOMAIN_ITERATOR_PAIR_DOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IteratorPairDomain.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
