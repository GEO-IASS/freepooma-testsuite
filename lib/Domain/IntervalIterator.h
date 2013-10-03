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

#ifndef POOMA_DOMAIN_INTERVAL_ITERATOR_H
#define POOMA_DOMAIN_INTERVAL_ITERATOR_H

//-----------------------------------------------------------------------------
// Classes: 
//   IntervalIterator
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Domain
 * @brief
 *   IntervalIterator - Iterates through Interval<1> points.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Utilities/PAssert.h"

#include <iterator> // std::random_access_iterator_tag
#include <stddef.h> // ptrdiff_t

//-----------------------------------------------------------------------------
// Open POOMA namespace:
//-----------------------------------------------------------------------------

// namespace POOMA {

/**
 * A random access iterator class that iterates through all of the points
 * of an Interval<1>, returning ints when dereferenced.
 */

class IntervalIterator
{
public:

  //============================================================
  // Typedefs
  //============================================================

  typedef IntervalIterator                    This_t;
  typedef Interval<1>                         Domain_t;
  typedef ptrdiff_t                           Value_t;

  typedef std::random_access_iterator_tag     iterator_category;
  typedef ptrdiff_t                           value_type;
  typedef ptrdiff_t                           difference_type;
  typedef const ptrdiff_t*                    pointer;
  typedef const ptrdiff_t&                    reference;
  
  //============================================================
  // Constructors
  //============================================================

  // The main IntervalIterator constructor stores the given domain
  // and initializes the iterator to the beginning.

  IntervalIterator(const Domain_t &d, int initial_pos = 0)
  : domain_m(d), val_m(d.first() + initial_pos)
  { }
    
  // Copy constructor.

  IntervalIterator(const This_t &it)
  : domain_m(it.domain_m), val_m(it.val_m)
  { }

  // The default constructor constructs an end iterator for an empty
  // domain.

  IntervalIterator() : val_m(1) { };

  //============================================================
  // Accessors
  //============================================================

  // Dereference operator. Returns const ref to internal Loc.

  inline const Value_t &operator*() const { PAssert(!done()); return val_m; }

  // Member selection operator. Not really useful (ints have no 
  // members to invoke), but part of the required interface. 

  inline const Value_t *operator->() const { PAssert(!done()); return &val_m; }

  // Equality tests.
  // This only tests that the iterators have the same value.
  // It does not test whether the underlying domains are the same.
  
  inline bool operator==(const This_t &i) const { return val_m == i.val_m; }
  inline bool operator!=(const This_t &i) const { return val_m != i.val_m; }
  inline bool operator< (const This_t &i) const { return val_m <  i.val_m; }
  inline bool operator<=(const This_t &i) const { return val_m <= i.val_m; }
  inline bool operator> (const This_t &i) const { return val_m >  i.val_m; }
  inline bool operator>=(const This_t &i) const { return val_m >= i.val_m; }
  
  inline Value_t operator[](int n) const { return *(*this + n); }
  
  This_t operator+(int n) const
  {
    IntervalIterator ret(*this);
    ret += n;
    return ret;
  }

  This_t operator-(int n) const
  {
    IntervalIterator ret(*this);
    ret -= n;
    return ret;
  }

  ptrdiff_t
  operator-(const This_t &it) const
  {
    PAssert(domain_m == it.domain_m);
    return val_m - it.val_m;
  }
    
  // At-end (false) test.
  // Returns true if this iterator is at-end.

  bool done() const { return (val_m > domain_m.last()); }

  //============================================================
  // Mutators
  //============================================================

  // Assignment operator.

  This_t &operator=(const This_t &it)
  {
    if (&it != this)
      {
        domain_m = it.domain_m;
        val_m    = it.val_m;
      }
    return *this;
  }

  // Pre-increment operator. 

  This_t &operator++() { increment();   return *this; }
  This_t &operator--() { increment(-1); return *this; }

  // Post-increment operator.
  // This has to make a copy, so prefer the above if possible.

  This_t operator++(int)
  {
    IntervalIterator save(*this);
    increment();
    return save;
  }

  This_t operator--(int)
  {
    IntervalIterator save(*this);
    increment(-1);
    return save;
  }

  inline This_t operator+=(int n) { increment(n);  return *this; }
  inline This_t operator-=(int n) { increment(-n); return *this; }

private:

  //============================================================
  // Data members.
  //============================================================

  // The domain we're iterating over.

  Domain_t domain_m;

  // Our current value.
  
  Value_t val_m;
    
  //============================================================
  // Implementation functions
  //============================================================

  // Increment iterator.

  inline void increment()
  {
    PAssert(!done());
    ++val_m;
  }
  
  inline void increment(int n)
  {
//    PAssert(!done()); // fails if we try to decrement from end, which is legal
    val_m += n;
  }
};

inline IntervalIterator operator+(int n, const IntervalIterator &it)
{
  IntervalIterator ret(it);
  ret += n;
  return ret;
}

// } // namespace POOMA

#endif // POOMA_DOMAIN_INTERVAL_ITERATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IntervalIterator.h,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
