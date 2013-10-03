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

#ifndef POOMA_TINY_VECTOR_ELEMENTS_H
#define POOMA_TINY_VECTOR_ELEMENTS_H

//-----------------------------------------------------------------------------

// Class: 
//    VectorElem
//    VectorAssign
//    VectorBinaryCombine
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// Trait classes for getting elements of Tiny objects at compile time.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class E> class Vector;
template<int D, class T, class E> class VectorEngine;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// The general templates for the class VectorElem.
// VectorElem should be specialized for vector-like classes.
//
// The general definition is for scalars which cannot be subscripted.
//
// The reason for this is so that you can have an object that acts like
// an n-dimensional vector but has some components generated algorithmicly.
// For example, it might have only one nonzero component.  In that case
// VectorElem would be specialized for that class and would return a Zero
// for all but the one nonzero component.
//
//-----------------------------------------------------------------------------

template<class V, int I>
struct VectorElem
{
  typedef V Element_t;
  typedef const V& ConstElementRef_t;
  typedef V& ElementRef_t;
  static const V& get(const V& x) { return x; }
  static       V& get(      V& x) { return x; }
};

template<int D, class T, class E, int I>
struct VectorEngineElem
{
  typedef VectorEngine<D,T,E> V;
  typedef typename V::Element_t Element_t;
  typedef typename V::ConstElementRef_t ConstElementRef_t;
  typedef typename V::ElementRef_t ElementRef_t;
  static ConstElementRef_t get(const V& x) { return x(I); }
  static ElementRef_t      get(      V& x) { return x(I); }
};

template<int D, class T, class E, int I>
struct VectorElem< Vector<D,T,E> , I >
{
  typedef Vector<D,T,E> V;
  typedef VectorEngineElem<D,T,E,I> VE;
  typedef typename VE::Element_t Element_t;
  typedef typename VE::ConstElementRef_t ConstElementRef_t;
  typedef typename VE::ElementRef_t ElementRef_t;
  static ConstElementRef_t get(const V& x) { return VE::get(x.engine()); }
  static ElementRef_t get(V& x) { return VE::get(x.engine()); }
};

//-----------------------------------------------------------------------------
//
// VectorAssign
//
// Template metaprogram for copying out of one vector and into another.
// Input:
//   The vector we're writing into.
//   Something to copy out of.
//
// Evaluate by recursing on the first half and the second half.
// Terminate at either length 1 or length 2.
//
//-----------------------------------------------------------------------------

template<class V1, class V2, class Op, int B, int L>
struct VectorAssign
{
  static void apply(V1& v1, const V2& v2, Op op)
    {
      CTAssert(L>1);
      VectorAssign<V1,V2,Op,B,L/2>::apply(v1,v2,op);
      VectorAssign<V1,V2,Op,B+L/2,L-L/2>::apply(v1,v2,op);
    }
};

template<class V1, class V2, class Op, int B>
struct VectorAssign<V1,V2,Op,B,0>
{
  static void apply(V1&, const V2&, Op) {}
};

template<class V1, class V2, class Op, int B>
struct VectorAssign<V1,V2,Op,B,1>
{
  static void apply(V1& v1, const V2& v2, Op op) 
    {
      op(VectorElem<V1,B>::get(v1),VectorElem<V2,B>::get(v2));
    }
};

template<class V1, class V2, class Op, int B>
struct VectorAssign<V1,V2,Op,B,2>
{
  static void apply(V1& v1, const V2& v2, Op op) 
    {
      op(VectorElem<V1,B  >::get(v1), VectorElem<V2,B  >::get(v2));
      op(VectorElem<V1,B+1>::get(v1), VectorElem<V2,B+1>::get(v2));
    }
};

template<class V1, class V2, class Op, int B>
struct VectorAssign<V1,V2,Op,B,3>
{
  static void apply(V1& v1, const V2& v2, Op op) 
    {
      op(VectorElem<V1,B  >::get(v1), VectorElem<V2,B  >::get(v2));
      op(VectorElem<V1,B+1>::get(v1), VectorElem<V2,B+1>::get(v2));
      op(VectorElem<V1,B+2>::get(v1), VectorElem<V2,B+2>::get(v2));
    }
};



#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: VectorElements.h,v $   $Author: richi $
// $Revision: 1.9 $   $Date: 2004/11/29 14:19:30 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
