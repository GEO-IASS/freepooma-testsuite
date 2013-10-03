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

#ifndef POOMA_DOMAIN_INDIRECTIONLIST_ITERATOR_H
#define POOMA_DOMAIN_INDIRECTIONLIST_ITERATOR_H

//-----------------------------------------------------------------------------
// Classes: 
//   IndirectionListIterator
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Domain
 * @brief
 * IndirectionListIterator - Iterates through IndirectionList<T> elements.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/IndirectionList.h"
#include "Utilities/PAssert.h"

#include <iterator> // std::random_access_iterator_tag
#include <stddef.h> // ptrdiff_t

//-----------------------------------------------------------------------------
// Open POOMA namespace:
//-----------------------------------------------------------------------------

// namespace POOMA {

/**
 * A random access iterator class that iterates through all of the elements
 * of an IndirectionList<T>, returning T's when dereferenced.
 */

template <class T>
class IndirectionListIterator
{
public:

  //============================================================
  // Typedefs
  //============================================================

  typedef IndirectionListIterator<T>          This_t;
  typedef IndirectionList<T>                  Domain_t;
  typedef T                                   Value_t;

  typedef std::random_access_iterator_tag     iterator_category;
  typedef T                                   value_type;
  typedef ptrdiff_t                           difference_type;
  typedef const T*                            pointer;
  typedef const T&                            reference;

  //============================================================
  // Constructors
  //============================================================

  // The main IndirectionListIterator constructor stores the given domain
  // and initializes the iterator to the beginning.

  IndirectionListIterator(const Domain_t &d, int initial_pos = 0)
  : domain_m(d), pos_m(initial_pos)
  {
    if (!done())
      val_m = domain_m(pos_m);
  }
    
  // Copy constructor.

  IndirectionListIterator(const This_t &it)
  : domain_m(it.domain_m), pos_m(it.pos_m), val_m(it.val_m)
  { }

  // The default constructor constructs an end iterator for an empty domain.

  IndirectionListIterator() { }

  //============================================================
  // Accessors
  //============================================================

  // Dereference operator. Returns const ref to T.

  inline const Value_t & operator*() const { PAssert(!done()); return val_m; }

  // Member selection operator. 

  inline const Value_t * operator->() const { PAssert(!done()); return &val_m; }

  // Equality tests.
  // This only tests that the iterators have the same position.
  // It does not test whether the underlying domains are the same.
  
  inline bool operator==(const This_t &i) const { return pos_m == i.pos_m; }
  inline bool operator!=(const This_t &i) const { return pos_m != i.pos_m; }
  inline bool operator< (const This_t &i) const { return pos_m <  i.pos_m; }
  inline bool operator<=(const This_t &i) const { return pos_m <= i.pos_m; }
  inline bool operator> (const This_t &i) const { return pos_m >  i.pos_m; }
  inline bool operator>=(const This_t &i) const { return pos_m >= i.pos_m; }
  
  inline Value_t operator[](int n) const
  {
    return domain_m(pos_m + n);
  }
  
  This_t operator+(int n) const
  {
    This_t ret(*this);
    ret += n;
    return ret;
  }

  This_t operator-(int n) const
  {
    This_t ret(*this);
    ret -= n;
    return ret;
  }

  ptrdiff_t
  operator-(const This_t &it) const
  {
    //    PAssert(domain_m == it.domain_m);
    //    Cannot check this currently for IndirectionLists!
    return pos_m - it.pos_m;
  }
    
  // At-end (false) test.
  // Returns true if this iterator is at-end.

  bool done() const { return (pos_m > domain_m.size()-1); }

  //============================================================
  // Mutators
  //============================================================

  // Assignment operator.

  This_t &operator=(const This_t &it)
  {
    if (&it != this)
      {
        domain_m = it.domain_m;
        pos_m    = it.pos_m;
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
    This_t save(*this);
    increment();
    return save;
  }

  This_t operator--(int)
  {
    This_t save(*this);
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

  // Our current position in the domain.
  
  int pos_m;
    
  // Our current value.
  
  Value_t val_m;
    
  //============================================================
  // Implementation functions
  //============================================================

  // Increment iterator.

  inline void increment()
  {
    PAssert(!done());
    ++pos_m;
    if (!done())
      val_m = domain_m(pos_m);
    else
      val_m = Value_t();
  }
  
  inline void increment(int n)
  {
    //    PAssert(!done()); 
    //    fails if we try to decrement from end, which is legal
    pos_m += n;
    if (!done())
      val_m = domain_m(pos_m);
    else
      val_m = Value_t();
  }
};

template <class T>
inline
IndirectionListIterator<T> operator+(int n,
  const IndirectionListIterator<T> &it)
{
  IndirectionListIterator<T> ret(it);
  ret += n;
  return ret;
}

// } // namespace POOMA

#endif // POOMA_DOMAIN_INDIRECTIONLIST_ITERATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IndirectionListIterator.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
