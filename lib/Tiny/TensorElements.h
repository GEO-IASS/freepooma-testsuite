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

#ifndef POOMA_TINY_TENSOR_ELEMENTS_H
#define POOMA_TINY_TENSOR_ELEMENTS_H

//-----------------------------------------------------------------------------
// Class: 
//    TensorElem
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// Trait classes for getting elements of Tiny objects at compile time.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class E> class Tensor;
template<int D, class T, class E> class TensorEngine;
class Antisymmetric;
class Symmetric;
class Diagonal;



/**
 * Trivial class for returning something that can't assigned into:
 */

struct Unwritable
{
public:
  // Assignment is not allowed; make it a no-op:
  // 1) Assign from any other type:
  template<class T>
  void operator=(const T&) { }
  // 2) Assign from another Unwritable:
  void operator=(const Unwritable&) { }
};


/**
 * Returns true or false (compile-time value, really an enum = 0 or 1) for
 * whether element (I,J) of a Tensor type is writable. In general, for example
 * for Tensors using Full for their EngineTag parameter, this is always
 * true. For a Diagonal Tensor, for example, it's false except when I=J.
 */

// Generic (EngineTag) case, all elements writable:
template<int D, class E, int I, int J>
class Writable
{
public:
  enum { value = 1 };
};

// Antisymmetric case, only elements (i,j) with i>j writable:
template<int D, int I, int J>
class Writable<D, Antisymmetric, I, J>
{
public:
  enum { value = (I > J) };
};

// Symmetric case, only elements (i,j) with i>=j writable:
template<int D, int I, int J>
class Writable<D, Symmetric, I, J>
{
public:
  enum { value = (I >= J) };
};

// Diagonal case, only elements (i,j) with i=j writable:
template<int D, int I, int J>
class Writable<D, Diagonal, I, J>
{
public:
  enum { value = (I == J) };
};



/**
 * The "boolean" (really an int) parameter "B" in TensorEngineElem<> flags
 * whether the I,J element is writable, according to the engine type. For Full
 * engines, all are writable, for things like Antisymmetric, only some are
 * writable. The external class Writable<> answers the question whether B is
 * true or false, down in the TensorElem<> implementation below.
 */

template<int D, class T, class E, int I, int J, int B=1>
struct TensorEngineElem;

// Partial specialization for B=true (B=1); allows getting writable references:
template<int D, class T, class E, int I, int J>
struct TensorEngineElem<D, T, E, I, J, 1>
{
  typedef TensorEngine<D,T,E> V;
  typedef typename V::Element_t         Element_t;
  typedef typename V::CTConstElementRef_t ConstElementRef_t;
  typedef typename V::CTElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return x.template getIJ<I,J>(); }
  static ElementRef_t      get(      V& x) { return x.template getIJ<I,J>(); }
};

// Partial specialization for B=false (B=0); returns dummy unwritable
// references (instances of the trivial Unwritable class, for which
// operator=(T) does nothing:
template<int D, class T, class E, int I, int J>
struct TensorEngineElem<D, T, E, I, J, 0>
{
  typedef TensorEngine<D,T,E> V;
  typedef typename V::Element_t         Element_t;
  typedef typename V::CTConstElementRef_t ConstElementRef_t;
  typedef Unwritable ElementRef_t;
  static ConstElementRef_t get(const V& x) { return x.template getIJ<I,J>(); }
  static ElementRef_t      get(      V& x) { return Unwritable(); }
};



/**
 * The general templates for the class TensorElem.
 * TensorElem should be specialized for Tensor-like classes.
 *
 * The general definition is for scalars which cannot be subscripted.  We also
 * have specializations for tensors with arbitrary engines which just use
 * operator() with both integers.  This is the fallback if a given engine type
 * doesn't specify anything else.
 */

template<class V, int I, int J>
struct TensorElem
{
  typedef       V  Element_t;
  typedef const V& ConstElementRef_t;
  typedef       V& ElementRef_t;
  static ConstElementRef_t get(const V& x) { return x; }
  static      ElementRef_t get(      V& x) { return x; }
};

template<int D, class T, class E, int I, int J>
struct TensorElem< Tensor<D,T,E> , I , J>
{
  typedef Tensor<D,T,E> V;
  typedef TensorEngineElem<D,T,E,I,J,Writable<D,E,I,J>::value> TE;
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x.engine()); }
  static      ElementRef_t get(V& x)       { return TE::get(x.engine()); }
};

template<int D, class T, class E, int I, int J>
struct TensorElem< TensorEngine<D,T,E> , I , J>
{
  // empty default specialization for TensorEngine - these need to specialize
  // TensorElem for not falling back to the scalar default.
};



//-----------------------------------------------------------------------------
//
// TensorAssign
//
// Template metaprogram for copying out of one tensor and into another.
// Input:
//   The tensor we're writing into.
//   Something to copy out of.
//
// Evaluate by recursing on the quadrants of the tensor.
//
//-----------------------------------------------------------------------------

//
// The general case of copying divides the tensor into quadrants
// and calls copy on each quadrant.
// This will be applied if each axis of the tensor is larger than 1.
//

template<class T1, class T2, class Op, int B1, int L1, int B2, int L2>
struct TensorAssign
{
  enum { B11=B1 , L11=L1/2 , B12=B1+L1/2 , L12 = L1-L1/2 };
  enum { B21=B2 , L21=L2/2 , B22=B2+L2/2 , L22 = L2-L2/2 };
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      TensorAssign<T1,T2,Op,B11,L11,B21,L21>::apply(x,y,op);
      TensorAssign<T1,T2,Op,B12,L12,B21,L21>::apply(x,y,op);
      TensorAssign<T1,T2,Op,B11,L11,B22,L22>::apply(x,y,op);
      TensorAssign<T1,T2,Op,B12,L12,B22,L22>::apply(x,y,op);
    }
};

//
// The case for a column.
// Divide the column in two, and recurse.
//

template<class T1, class T2, class Op, int B1, int L1, int B2>
struct TensorAssign<T1,T2,Op,B1,L1,B2,1>
{
  enum { B11=B1 , L11=L1/2 , B12=B1+L1/2 , L12 = L1-L1/2 };
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      TensorAssign<T1,T2,Op,B11,L11,B2,1>::apply(x,y,op);
      TensorAssign<T1,T2,Op,B12,L12,B2,1>::apply(x,y,op);
    }
};

//
// The case for a row.
// Divide the row in two, and recurse.
//

template<class T1, class T2, class Op, int B1, int B2, int L2>
struct TensorAssign<T1,T2,Op,B1,1,B2,L2>
{
  enum { B21=B2 , L21=L2/2 , B22=B2+L2/2 , L22 = L2-L2/2 };
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      TensorAssign<T1,T2,Op,B1,1,B21,L21>::apply(x,y,op);
      TensorAssign<T1,T2,Op,B1,1,B22,L22>::apply(x,y,op);
    }
};

//
// The case for a single element.
// Just do it.
//

template<class T1, class T2, class Op, int B1, int B2>
struct TensorAssign<T1,T2,Op,B1,1,B2,1>
{
  static void apply(T1& x, const T2& y,Op op=Op())
    {
      op(TensorElem<T1,B1,B2>::get(x), TensorElem<T2,B1,B2>::get(y));
    }
};

//
// The case for a two by two block.
// Just do it.
//

template<class T1, class T2, class Op, int B1, int B2>
struct TensorAssign<T1,T2,Op,B1,2,B2,2>
{
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      op(TensorElem<T1,B1  ,B2  >::get(x), TensorElem<T2,B1  ,B2  >::get(y));
      op(TensorElem<T1,B1+1,B2  >::get(x), TensorElem<T2,B1+1,B2  >::get(y));
      op(TensorElem<T1,B1  ,B2+1>::get(x), TensorElem<T2,B1  ,B2+1>::get(y));
      op(TensorElem<T1,B1+1,B2+1>::get(x), TensorElem<T2,B1+1,B2+1>::get(y));
    }
};

//
// The case for a three by three block.
// Just do it.
//

template<class T1, class T2, class Op, int B1, int B2>
struct TensorAssign<T1,T2,Op,B1,3,B2,3>
{
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      op(TensorElem<T1,B1  ,B2  >::get(x), TensorElem<T2,B1  ,B2  >::get(y));
      op(TensorElem<T1,B1+1,B2  >::get(x), TensorElem<T2,B1+1,B2  >::get(y));
      op(TensorElem<T1,B1+2,B2  >::get(x), TensorElem<T2,B1+2,B2  >::get(y));
      op(TensorElem<T1,B1  ,B2+1>::get(x), TensorElem<T2,B1  ,B2+1>::get(y));
      op(TensorElem<T1,B1+1,B2+1>::get(x), TensorElem<T2,B1+1,B2+1>::get(y));
      op(TensorElem<T1,B1+2,B2+1>::get(x), TensorElem<T2,B1+2,B2+1>::get(y));
      op(TensorElem<T1,B1  ,B2+2>::get(x), TensorElem<T2,B1  ,B2+2>::get(y));
      op(TensorElem<T1,B1+1,B2+2>::get(x), TensorElem<T2,B1+1,B2+2>::get(y));
      op(TensorElem<T1,B1+2,B2+2>::get(x), TensorElem<T2,B1+2,B2+2>::get(y));
    }
};



#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TensorElements.h,v $   $Author: richi $
// $Revision: 1.18 $   $Date: 2004/11/29 14:19:30 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
