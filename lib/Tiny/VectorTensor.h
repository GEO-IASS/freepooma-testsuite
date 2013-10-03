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

#ifndef POOMA_TINY_VECTOR_TENSOR_H
#define POOMA_TINY_VECTOR_TENSOR_H

//-----------------------------------------------------------------------------
// Functions: 
//   vector dot(vector,tensor)
//   vector dot(tensor,vector)
//   tensor outerProduct(vector,vector)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// Dot products between vectors and tensors, both yielding vectors.
// Outer product between vectors, yielding tensor.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------
#include "Pooma/PoomaOperatorTags.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class E> class Tensor;
template<int D, class T, class E> class TensorEngine;
template<int D, class T, class E> class Vector;
template<int D, class T, class E> class VectorEngine;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Subroutines and metaprograms for dotting a vector with a tensor.
//
// As with other things that require template metaprograms, there is
// way more here than you would believe for the seemingly simple job to
// be done.
//
// The chain of operations is:
// 1.  The user calls dot(vector,tensor)
// 2.  That builds a new vector, using a two-arg expression template.
// 3.  The ctor of the new vector evaluates the expr for each element.
// 4.  Each element invokes VectorDotTensor for one vector dot product.
// 5.  VectorDotTensor recurses, split the sum into halves, add results.
// 6.  When the length is one, it multiplies elements and returns that.
// 7.  Elements from the vector and tensor come through VectorElem and
//     TensorElem, so that the type of each one can be different.
//
//-----------------------------------------------------------------------------

//
// General VectorDotTensor
// Takes the dot product of vector of type V1 with column I of 
// a tensor of type T2, and the vector starts with offset B and
// has length L.
//
// Operates by splitting the domain in half, taking the dot product
// of each half, and returning the sum of the results.
//

template<class V1, class T2, int I, int B, int L>
struct VectorDotTensor
{
  typedef typename VectorDotTensor<V1,T2,I,B,L/2>::Type_t E1;
  typedef typename VectorDotTensor<V1,T2,I,B+L/2,L-L/2>::Type_t E2;
  typedef typename BinaryReturn<E1,E2,OpAdd>::Type_t Type_t;
  static Type_t get(const V1& v1, const T2& t2)
    {
      return 
        VectorDotTensor<V1,T2,I,B,L/2>::get(v1,t2) +
        VectorDotTensor<V1,T2,I,B+L/2,L-L/2>::get(v1,t2);
    }
};

//
// Recursion termination for VectorDotTensor
// When the length of the vectors to be gets down to 1, just
// multiply the elements together and return that.
//

template<class V1, class T2, int I, int B>
struct VectorDotTensor<V1,T2,I,B,1>
{
  typedef typename VectorElem<V1,B>::Element_t E1;
  typedef typename TensorElem<T2,B,I>::Element_t E2;
  typedef typename BinaryReturn<E1,E2,OpMultiply>::Type_t Type_t;
  static Type_t get(const V1& v1, const T2& t2)
    {
      return VectorElem<V1,B>::get(v1) * TensorElem<T2,B,I>::get(t2);
    }
};

//
// Extract a single value from a vector engine for vector dot tensor.
// The input vector is of size D and has element type T1 and engine E1
// The input tensor is of size D by D and has element type T2 and engine E2
// The output vector is of size D and has element type T3
// The vector is dotted with column I of the tensor.
//

template<int D, class T1, class T2, class T3, class E1, class E2, int I>
struct VectorEngineElem<D,T3,
  BinaryVectorOp<Vector<D,T1,E1>,Tensor<D,T2,E2>,FnDot>,I>
{
  typedef Vector<D,T1,E1> V1;
  typedef Tensor<D,T2,E2> V2;
  typedef VectorEngine<D,T3,BinaryVectorOp<V1,V2,FnDot> > V;
  typedef typename VectorDotTensor<V1,V2,I,0,D>::Type_t T0;
  typedef T0 Element_t;
  typedef T0 ConstElementRef_t;
  typedef T0 ElementRef_t;
  static T0 get(const V& x) 
    {
      return VectorDotTensor<V1,V2,I,0,D>::get(x.v1_m,x.v2_m); 
    }
};

//
// Define the return type for vector dot tensor.
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D,T1,E1> , Tensor<D,T2,E2> , FnDot >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
  typedef Vector<D,T0,Full> Type_t;
};

//
// Take the dot product of a vector and a tensor returning a vector.
//

template<int D, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn< Vector<D,T1,E1>,Tensor<D,T2,E2> , FnDot >::Type_t
dot( const Vector<D,T1,E1>& v1, const Tensor<D,T2,E2>& v2 )
{
  typedef Vector<D,T1,E1> V1;
  typedef Tensor<D,T2,E2> V2;
  typedef typename BinaryReturn<V1,V2,FnDot>::Type_t Return_t;
  typedef typename Return_t::Element_t T3;
  typedef Vector<D,T3,BinaryVectorOp<V1,V2,FnDot> > Expr_t;
  return Return_t( Expr_t(v1,v2) );
}

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Subroutines and metaprograms for dotting a tensor with a vector.
//
// As with other things that require template metaprograms, there is
// way more here than you would believe for the seemingly simple job to
// be done.
//
// The chain of operations is:
// 1.  The user calls dot(vector,tensor)
// 2.  That builds a new vector, using a two-arg expression template.
// 3.  The ctor of the new vector evaluates the expr for each element.
// 4.  Each element invokes VectorDotTensor for one vector dot product.
// 5.  VectorDotTensor recurses, split the sum into halves, add results.
// 6.  When the length is one, it multiplies elements and returns that.
// 7.  Elements from the vector and tensor come through VectorElem and
//     TensorElem, so that the type of each one can be different.
//
//-----------------------------------------------------------------------------

//
// General TensorDotVector
//
// Much like VectorDotTensor above, this dots a vector with one row
// of a tensor.  It splits that dot product into two halves, adds the
// result and returns that.
//

template<class T1, class V2, int I, int B, int L>
struct TensorDotVector
{
  typedef typename TensorDotVector<T1,V2,I,B,L/2>::Type_t E1;
  typedef typename TensorDotVector<T1,V2,I,B+L/2,L-L/2>::Type_t E2;
  typedef typename BinaryReturn<E1,E2,OpAdd>::Type_t Type_t;
  static Type_t get(const T1& t1, const V2& v2)
    {
      return 
        TensorDotVector<T1,V2,I,B,L/2>::get(t1,v2) +
        TensorDotVector<T1,V2,I,B+L/2,L-L/2>::get(t1,v2);
    }
};

//
// Recursion termination for TensorDotVector
// When you get down to one term each from the vector and tensor, 
// just multiply the terms and return that.
//

template<class T1, class V2, int I, int B>
struct TensorDotVector<T1,V2,I,B,1>
{
  typedef typename TensorElem<T1,I,B>::Element_t E1;
  typedef typename VectorElem<V2,B>::Element_t E2;
  typedef typename BinaryReturn<E1,E2,OpMultiply>::Type_t Type_t;
  static Type_t get(const T1& t1, const V2& v2)
    {
      return TensorElem<T1,I,B>::get(t1) * VectorElem<V2,B>::get(v2);
    }
};

//
// Get one element from a vector of the dot product of a tensor and a vector.
// Gets one value by dotting the tensor and the vector.
//

template<int D, class T1, class T2, class T3, class E1, class E2, int I>
struct VectorEngineElem<D,T3,
  BinaryVectorOp<Tensor<D,T1,E1>,Vector<D,T2,E2>,FnDot>,I>
{
  typedef Tensor<D,T1,E1> V1;
  typedef Vector<D,T2,E2> V2;
  typedef VectorEngine<D,T3,BinaryVectorOp<V1,V2,FnDot> > V;
  typedef typename TensorDotVector<V1,V2,I,0,D>::Type_t T0;
  typedef T0 Element_t;
  typedef T0 ConstElementRef_t;
  typedef T0 ElementRef_t;
  static T0 get(const V& x) 
    {
      return TensorDotVector<V1,V2,I,0,D>::get(x.v1_m,x.v2_m); 
    }
};

//
// Define the return type for dotting a tensor and a vector.
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Tensor<D,T1,E1> , Vector<D,T2,E2> , FnDot >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
  typedef Vector<D,T0,Full> Type_t;
};

//
// Dot a tensor and a vector, returning a vector.
// Build a vector with a simple expression template representing 
// the dot product, and the ctor of the returned vector evaluates
// that expression.
//

template<int D, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn< Tensor<D,T1,E1> , Vector<D,T2,E2>, FnDot >::Type_t
dot( const Tensor<D,T1,E1>& v1 , const Vector<D,T2,E2>& v2 )
{
  typedef Tensor<D,T1,E1> V1;
  typedef Vector<D,T2,E2> V2;
  typedef typename BinaryReturn<V1,V2,FnDot>::Type_t Return_t;
  typedef typename Return_t::Element_t T3;
  typedef Vector<D,T3,BinaryVectorOp<V1,V2,FnDot> > Expr_t;
  return Return_t( Expr_t(v1,v2) );
}


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Subroutine for taking outer product between two vectors, yielding a tensor.
// Takes the outer product of vector of type V1 and vector of type V2
//
// The chain of operations is:
// 1.  The user calls outerProduct(vector,vector)
// 
// TJW: force this to return a Tensor with a Full engine (already so, via the
// the default EngineTag value, really).
//
//-----------------------------------------------------------------------------

//
// Define the return type for outerProduct(vector,vector).
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D,T1,E1> , Vector<D,T2,E2>, FnOuterProduct >
{
  typedef typename BinaryReturn<T1, T2, OpMultiply>::Type_t T0;
  typedef Tensor<D, T0> Type_t;
};


//
// Take the outer product of two vectors, returning a tensor
//

template<int D, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn<Vector<D,T1,E1>, Vector<D,T2,E2>, FnOuterProduct>::Type_t
outerProduct(const Vector<D,T1,E1> &v1, const Vector<D,T2,E2> &v2 )
{
  typedef typename
    BinaryReturn<Vector<D,T1,E1>, Vector<D,T2,E2>, FnOuterProduct>::Type_t
    Return_t;
  Return_t ret;

  for (int i = 0; i < D; ++i)
    for (int j = 0; j < D; ++j)
      ret(i, j) = v1(i) * v2(j);

  return ret;
}

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: VectorTensor.h,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:17:11 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
