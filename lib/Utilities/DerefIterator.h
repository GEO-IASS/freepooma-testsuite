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

#ifndef POOMA_UTILITIES_DEREFITERATOR_H
#define POOMA_UTILITIES_DEREFITERATOR_H

//-----------------------------------------------------------------------------
// Classes: 
//   DerefIterator
//   ConstDerefIterator
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Utilities
 * @brief
 * STL style iterators for lists of pointers.
 *
 * Unlike vector<T*>::iterator,
 * these automatically dereference themselves and, in the process,
 * maintain const correctness.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include <vector>


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

template <class T>
class ConstDerefIterator;

/**
 * DerefIterator<T> and ConstDerefIterator<T> are STL-style iterators that
 * are used to properly handle iterating through lists of pointers. Not only
 * is this a convenience, as these iterators automatically dereference
 * themselves, it also solves a problem with const correctness. If one has
 * vector<T*>::const_iterator, this only keeps the user from modifying the
 * pointer, not from modifying the object that is pointed to. What one really
 * wants is vector<const T*>::const_iterator, but that is not something
 * one can get from vector<T*>. 
 */

template <class T>
class DerefIterator
{
public:

  // Convenience typedefs
  
  typedef T                                   Value_t;
  typedef std::vector<T*>                     List_t;
  typedef DerefIterator<Value_t>              iterator;
  typedef ConstDerefIterator<Value_t>         const_iterator;
  
  // Required iterator typedefs
  
  typedef std::random_access_iterator_tag     iterator_category;
  typedef Value_t                             value_type;
  typedef ptrdiff_t                           difference_type;
  typedef Value_t*                            pointer;
  typedef Value_t&                            reference;

  friend class ConstDerefIterator<T>;

  friend iterator operator+ (const difference_type n, const iterator &iter)
  {
    return iterator(iter.p_m + n);
  }

protected:

  typename List_t::iterator p_m;

public:

  inline DerefIterator() { }
  inline DerefIterator(typename List_t::iterator x) : p_m(x) { }

  inline reference operator*() const {
    return **p_m; 
  }
  inline pointer operator->() const {
    return *p_m; 
  }
  inline iterator &operator++() {
    ++p_m;
    return *this;
  }
  inline iterator operator++(int) {
    iterator tmp = *this;
    ++p_m;
    return tmp;
  }
  inline iterator &operator--() {
    --p_m;
    return *this;
  }
  inline iterator operator--(int) {
    iterator tmp = *this;
    --p_m;
    return tmp;
  }
  inline iterator &operator+=(const difference_type i) {
    p_m += i;
    return *this;
  }
  inline iterator &operator-=(const difference_type i) {
    p_m -= i;
    return *this;
  }
  inline iterator operator+(const difference_type i) const {
    return iterator(p_m + i);
  }
  inline iterator operator-(const difference_type i) const {
    return iterator(p_m - i);
  }
  inline difference_type operator-(const iterator &x) const {
    return p_m - x.p_m;
  }
  inline difference_type operator-(const const_iterator &x) const {
    return p_m - x.p_m;
  }
  inline reference operator[](const difference_type i) const {
    return **(p_m + i); 
  }

  inline bool operator==(const iterator &x) const {
    return p_m == x.p_m;
  }
  inline bool operator==(const const_iterator &x) const {
    return p_m == x.p_m;
  }
  inline bool operator<(const iterator &x) const {
    return p_m < x.p_m;
  }
  inline bool operator<(const const_iterator &x) const {
    return p_m < x.p_m;
  }
  inline bool operator!= (const iterator &y) const
  { return ! (*this == y); }
  inline bool operator>  (const iterator &y) const
  { return (y < *this); }
  inline bool operator<= (const iterator &y) const
  { return  ! (y < *this); }
  inline bool operator>= (const iterator &y) const
  { return  ! (*this < y); }
  inline bool operator!= (const const_iterator &y) const
  { return  ! (*this == y); }
  inline bool operator>  (const const_iterator &y) const
  { return  (y < *this); }
  inline bool operator<= (const const_iterator &y) const
  { return  ! (y < *this); }
  inline bool operator>= (const const_iterator &y) const
  { return  ! (*this < y); }
};

template <class T>
class ConstDerefIterator
{
public:

  // Convience typedefs
  
  typedef T                                   Value_t;
  typedef std::vector<T*>                     List_t;
  typedef DerefIterator<Value_t>              iterator;
  typedef ConstDerefIterator<Value_t>         const_iterator;
  
  // Required typedefs
  
  typedef std::random_access_iterator_tag     iterator_category;
  typedef Value_t                             value_type;
  typedef ptrdiff_t                           difference_type;
  typedef const Value_t*                      pointer;
  typedef const Value_t&                      reference;

  friend class DerefIterator<T>;

  friend const_iterator operator+ (const difference_type n, 
    const const_iterator &iter) 
  {
    return const_iterator(iter.p_m + n);
  }

protected:

  typename List_t::const_iterator p_m;

public:

  inline ConstDerefIterator() { }
  inline ConstDerefIterator(const typename List_t::const_iterator &x) : p_m(x) { }
  inline ConstDerefIterator(const typename List_t::iterator &x) : p_m(x) { }
  inline ConstDerefIterator(const iterator &x) : p_m(x.p_m) { }
  inline ConstDerefIterator(const const_iterator &x) : p_m(x.p_m) { }

  inline reference operator*() const {
    return **p_m; 
  }
  inline pointer operator->() const {
    return *p_m; 
  }
  inline const_iterator &operator++() {
    ++p_m;
    return *this;
  }
  inline const_iterator operator++(int) {
    const_iterator tmp = *this;
    ++p_m;
    return tmp;
  }
  inline const_iterator &operator--() {
    --p_m;
    return *this;
  }
  inline const_iterator operator--(int) {
    const_iterator tmp = *this;
    --p_m;
    return tmp;
  }
  inline const_iterator &operator+=(const difference_type i) {
    p_m += i;
    return *this;
  }
  inline const_iterator &operator-=(const difference_type i) {
    p_m -= i;
    return *this;
  }
  inline const_iterator operator+(const difference_type i) const {
    return const_iterator(p_m + i);
  }
  inline const_iterator operator-(const difference_type i) const {
    return const_iterator(p_m - i);
  }
  inline difference_type operator-(const const_iterator &x) const {
    return p_m - x.p_m;
  }
  inline difference_type operator-(const iterator &x)
    const {
    return p_m - x.p_m;
  }
  inline reference operator[](const difference_type i) const {
    return **(p_m + i); 
  }

  inline bool operator==(const const_iterator &x) const {
    return p_m == x.p_m;
  }
  inline bool operator==(const iterator &x) const {
    return p_m == x.p_m;
  }
  inline bool operator<(const const_iterator &x) const {
    return p_m < x.p_m;
  }
  inline bool operator<(const iterator &x) const {
    return p_m < x.p_m;
  }
  inline bool operator!= (const const_iterator &y) const
  { return ! (*this == y); }
  inline bool operator>  (const const_iterator &y) const
  { return (y < *this); }
  inline bool operator<= (const const_iterator &y) const
  { return  ! (y < *this); }
  inline bool operator>= (const const_iterator &y) const
  { return  ! (*this < y); }
  inline bool operator!= (const iterator &y) const
  { return  ! (*this == y); }
  inline bool operator>  (const iterator &y) const
  { return  (y < *this); }
  inline bool operator<= (const iterator &y) const
  { return  ! (y < *this); }
  inline bool operator>= (const iterator &y) const
  { return  ! (*this < y); }
};

// } // namespace POOMA

#endif // POOMA_UTILITIES_DEREFITERATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DerefIterator.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
