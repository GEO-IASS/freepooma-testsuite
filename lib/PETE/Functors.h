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

#ifndef PETE_PETE_FUNCTORS_H
#define PETE_PETE_FUNCTORS_H

///////////////////////////////////////////////////////////////////////////////
//
// WARNING: THIS FILE IS FOR INTERNAL PETE USE. DON'T INCLUDE IT YOURSELF
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Include files.
//-----------------------------------------------------------------------------

#include <iterator>

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   LeafFunctor<LeafType, LeafTag>
//
// DESCRIPTION
//   LeafType is the type of something at the leaf of the expression tree.
//   LeafTag specifies the operation being applied.
//
//   LeafFunctors are used by ForEach to apply operations to the leaves of the
//   expression tree. Typical functors are evaluators, counters, etc.
//   Users define functors for use with ForEach by specializing
//   the struct LeafFunctor<LeafType, LeafTag> for the user defined Functor and
//   any Leaf types that are necessary.
//
//   This isn't a functor in the conventional sense since it isn't invoked
//   operator() on a LeafFunctor object. Instead, the static apply member
//   function is called without an object. In a lot of ways, this is better
//   and more flexible than a regular functor.
//
//   LeafFunctor specializations must define the following:
//
//      typedef ... Type_t; 
//         - the return type of the functor.
//      static Type_t apply(const LeafType &l, const LeafTag &f) {} 
//         - evaluates the functor on leaf l.
//
//-----------------------------------------------------------------------------

template<class LeafType, class LeafTag>
struct LeafFunctor
{ };

template<class LeafType, class LeafTag>
inline
typename LeafFunctor<LeafType, LeafTag>::Type_t
leafFunctor(const LeafType &leaf, const LeafTag &tag)
{
  return LeafFunctor<LeafType, LeafTag>::apply(leaf, tag);
}

//-----------------------------------------------------------------------------
//
// CLASS NAMES
//   EvalLeaf1-7, LeafFunctor<Scalar<T>, EvalLeaf1-7 >
//
// DESCRIPTION
//   EvalLeaf1 through EvalLeaf7 are used to evaluate leaves using 1 to 7
//   integer indices. We supply the tags and scalar versions here. Users must
//   supply specializations for their own containers.
//
//-----------------------------------------------------------------------------

// 1D

struct EvalLeaf1
{
  int i1_m;
  inline EvalLeaf1(int i1) : i1_m(i1) { }
  inline int val1() const { return i1_m; }
};
  
template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf1>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf1 &) 
  {
    return s.value();
  }
};

// 2D

struct EvalLeaf2
{
  int i1_m, i2_m;
  inline EvalLeaf2(int i1, int i2) : i1_m(i1), i2_m(i2) { }
  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
};
  
template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf2>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf2 &) 
  {
    return s.value();
  }
};

// 3D

struct EvalLeaf3
{
  int i1_m, i2_m, i3_m;
  inline EvalLeaf3(int i1, int i2, int i3) 
    : i1_m(i1), i2_m(i2), i3_m(i3) { }
  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
};
  
template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf3>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf3 &) 
  {
    return s.value();
  }
};

// 4D

struct EvalLeaf4
{
  int i1_m, i2_m, i3_m, i4_m;
  inline EvalLeaf4(int i1, int i2, int i3, int i4) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4) { }
  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
};
  
template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf4>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf4 &) 
  {
    return s.value();
  }
};

// 5D

struct EvalLeaf5
{
  int i1_m, i2_m, i3_m, i4_m, i5_m;
  inline EvalLeaf5(int i1, int i2, int i3, int i4, int i5) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5) { }
  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  inline int val5() const { return i5_m; }
};

template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf5>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf5 &)
  {
    return s.value();
  }
};

// 6D

struct EvalLeaf6
{
  int i1_m, i2_m, i3_m, i4_m, i5_m, i6_m;
  inline EvalLeaf6(int i1, int i2, int i3, int i4, int i5, int i6) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5), i6_m(i6) { }
  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  inline int val5() const { return i5_m; }
  inline int val6() const { return i6_m; }
};
  
template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf6>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf6 &) 
  {
    return s.value();
  }
};

// 7D

struct EvalLeaf7
{
  int i1_m, i2_m, i3_m, i4_m, i5_m, i6_m, i7_m;
  inline EvalLeaf7(int i1, int i2, int i3, int i4, int i5, int i6,
    int i7) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5), i6_m(i6), i7_m(i7) { }
  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  inline int val5() const { return i5_m; }
  inline int val6() const { return i6_m; }
  inline int val7() const { return i7_m; }
};
  
template<class T>
struct LeafFunctor<Scalar<T>, EvalLeaf7>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const EvalLeaf7 &) 
  {
    return s.value();
  }
};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//   IncrementLeaf, LeafFunctor<{T, Scalar<T>}, IncrementLeaf >
//
// DESCRIPTION
//   A leaf-tag and functor used to increment an iterator for scalars and
//   leaves made up of STL iterators.
//
//-----------------------------------------------------------------------------

struct IncrementLeaf
{ };
  
template<class T>
struct LeafFunctor<T, IncrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const T &cl, const IncrementLeaf &) 
  {
    T &l = const_cast<T &>(cl);
    ++l;
    return 0;
  }
};

#if defined(__MWERKS__)

// Workaround for screwy CWPro 4.1 bug.

template <class T>
struct LeafFunctor<const T*, IncrementLeaf> 
{
  typedef int Type_t;
  inline static
  Type_t apply(const T* & const ci, const IncrementLeaf &)
  {
    T* &i = const_cast<T* &>(ci);
    ++i;
    return 0;
  }
};

#endif
  
template<class T>
struct LeafFunctor<Scalar<T>, IncrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const Scalar<T> &, const IncrementLeaf &) 
  {
    return 0;
  }
};


//-----------------------------------------------------------------------------
//
// CLASS NAME
//   DecrementLeaf, LeafFunctor<{T, Scalar<T>}, DecrementLeaf >
//
// DESCRIPTION
//   A leaf-tag and functor used to decrement an iterator for scalars and
//   leaves made up of STL iterators.
//
//-----------------------------------------------------------------------------

struct DecrementLeaf
{ };
  
template<class T>
struct LeafFunctor<T, DecrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const T &cl, const DecrementLeaf &) 
  {
    T &l = const_cast<T &>(cl);
    --l;
    return 0;
  }
};

#if defined(__MWERKS__)
// Workaround for screwy CWPro 4.1 bug.
template <class T>
struct LeafFunctor<const T*, DecrementLeaf> 
{
  typedef int Type_t;
  inline static
  Type_t apply(const T* & const ci, const IncrementLeaf &)
  {
    T* &i = const_cast<T* &>(ci);
    --i;
    return 0;
  }
};
#endif
  
template<class T>
struct LeafFunctor<Scalar<T>, DecrementLeaf>
{
  typedef int Type_t;
  inline static
  Type_t apply(const Scalar<T> &, const DecrementLeaf &) 
  {
    return 0;
  }
};

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   DereferenceLeaf, LeafFunctor<{T, Scalar<T>}, DereferenceLeaf >
//
// DESCRIPTION
//   A leaf-tag and functor used to dereference an iterator for scalars and
//   leaves made up of STL iterators.
//
//-----------------------------------------------------------------------------

struct DereferenceLeaf
{ };

template<class ForwardIterator>
struct LeafFunctor<ForwardIterator, DereferenceLeaf>
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Type_t;
  inline static
  Type_t apply(const ForwardIterator &i, const DereferenceLeaf &)
  {
    return *i;
  }
};

#if defined(__MWERKS__)
// Workaround for screwy CWPro 4.1 bug.
template <class T>
struct LeafFunctor<const T*, DereferenceLeaf> 
{
  typedef T Type_t;
  inline static
  Type_t apply(const T *i, const DereferenceLeaf &)
  {
    return *i;
  }
};
#endif
  
template<class T>
struct LeafFunctor<Scalar<T>, DereferenceLeaf>
{
  typedef T Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const DereferenceLeaf &) 
  {
    return s.value();
  }
};


#endif // PETE_PETE_FUNCTORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Functors.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:56 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
