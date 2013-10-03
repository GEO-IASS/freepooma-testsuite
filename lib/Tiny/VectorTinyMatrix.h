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

#ifndef POOMA_TINY_VECTOR_TINYMATRIX_H
#define POOMA_TINY_VECTOR_TINYMATRIX_H

//-----------------------------------------------------------------------------

// Functions: 
//   vector dot(vector,TinyMatrix)
//   vector dot(TinyMatrix,vector)
//   TinyMatrix outerProductAsTinyMatrix(vector,vector)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// Dot products between vectors and TinyMatrixs, both yielding vectors.
// Outer product between vectors, yielding TinyMatrix.
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

template<int D1, int D2, class T, class E> class TinyMatrix;
template<int D1, int D2, class T, class E> class TinyMatrixEngine;
template<int D, class T, class E> class Vector;
template<int D, class T, class E> class VectorEngine;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Subroutines and metaprograms for dotting a vector with a TinyMatrix.
//
// As with other things that require template metaprograms, there is
// way more here than you would believe for the seemingly simple job to
// be done.
//
// The chain of operations is:
// 1.  The user calls dot(vector,TinyMatrix)
// 2.  That builds a new vector, using a two-arg expression template.
// 3.  The ctor of the new vector evaluates the expr for each element.
// 4.  Each element invokes VectorDotTinyMatrix for one vector dot product.
// 5.  VectorDotTinyMatrix recurses, split the sum into halves, add results.
// 6.  When the length is one, it multiplies elements and returns that.
// 7.  Elements from the vector and TinyMatrix come through VectorElem and
//     TinyMatrixElem, so that the type of each one can be different.
//
//-----------------------------------------------------------------------------

//
// General VectorDotTinyMatrix
// Takes the dot product of vector of type V1 with column I of 
// a TinyMatrix of type T2, and the vector starts with offset B and
// has length L.
//
// Operates by splitting the domain in half, taking the dot product
// of each half, and returning the sum of the results.
//

template<class V1, class T2, int I, int B, int L>
struct VectorDotTinyMatrix
{
  typedef typename VectorDotTinyMatrix<V1,T2,I,B,L/2>::Type_t E1;
  typedef typename VectorDotTinyMatrix<V1,T2,I,B+L/2,L-L/2>::Type_t E2;
  typedef typename BinaryReturn<E1,E2,OpAdd>::Type_t Type_t;
  static Type_t get(const V1& v1, const T2& t2)
    {
      return 
        VectorDotTinyMatrix<V1,T2,I,B,L/2>::get(v1,t2) +
        VectorDotTinyMatrix<V1,T2,I,B+L/2,L-L/2>::get(v1,t2);
    }
};

//
// Recursion termination for VectorDotTinyMatrix
// When the length of the vectors to be gets down to 1, just
// multiply the elements together and return that.
//

template<class V1, class T2, int I, int B>
struct VectorDotTinyMatrix<V1,T2,I,B,1>
{
  typedef typename VectorElem<V1,B>::Element_t E1;
  typedef typename TinyMatrixElem<T2,B,I>::Element_t E2;
  typedef typename BinaryReturn<E1,E2,OpMultiply>::Type_t Type_t;
  static Type_t get(const V1& v1, const T2& t2)
    {
      return VectorElem<V1,B>::get(v1) * TinyMatrixElem<T2,B,I>::get(t2);
    }
};

//
// Extract a single value from a vector engine for vector dot TinyMatrix.
// The input vector is of size D1 and has element type T1 and engine E1
// The input TinyMatrix is of size D1 by D2 and has element type T2 and engine E2
// The output vector is of size D2 and has element type T3
// The vector is dotted with column I of the TinyMatrix.
//

template<int D1, int D2, class T1, class T2, class T3, class E1, class E2, int I>
struct VectorEngineElem<D2,T3,
  BinaryVectorOp<Vector<D1,T1,E1>,TinyMatrix<D1,D2,T2,E2>,FnDot>,I>
{
  typedef Vector<D1,T1,E1> V1;
  typedef TinyMatrix<D1,D2,T2,E2> V2;
  typedef VectorEngine<D2,T3,BinaryVectorOp<V1,V2,FnDot> > V;
  typedef typename VectorDotTinyMatrix<V1,V2,I,0,D1>::Type_t T0;
  typedef T0 Element_t;
  typedef T0 ConstElementRef_t;
  typedef T0 ElementRef_t;
  static T0 get(const V& x) 
    {
      return VectorDotTinyMatrix<V1,V2,I,0,D1>::get(x.v1_m,x.v2_m); 
    }
};

//
// Define the return type for vector dot TinyMatrix.
//

template<int D1, int D2, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D1,T1,E1> , TinyMatrix<D1,D2,T2,E2> , FnDot >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
  typedef Vector<D2,T0,Full> Type_t;
};

//
// Take the dot product of a vector and a TinyMatrix returning a vector.
//

template<int D1, int D2, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn< Vector<D1,T1,E1>,TinyMatrix<D1,D2,T2,E2> , FnDot >::Type_t
dot( const Vector<D1,T1,E1>& v1, const TinyMatrix<D1,D2,T2,E2>& v2 )
{
  typedef Vector<D1,T1,E1> V1;
  typedef TinyMatrix<D1,D2,T2,E2> V2;
  typedef typename BinaryReturn<V1,V2,FnDot>::Type_t Return_t;
  typedef typename Return_t::Element_t T3;
  typedef Vector<D2,T3,BinaryVectorOp<V1,V2,FnDot> > Expr_t;
  return Return_t( Expr_t(v1,v2) );
}

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Subroutines and metaprograms for dotting a TinyMatrix with a vector.
//
// As with other things that require template metaprograms, there is
// way more here than you would believe for the seemingly simple job to
// be done.
//
// The chain of operations is:
// 1.  The user calls dot(vector,TinyMatrix)
// 2.  That builds a new vector, using a two-arg expression template.
// 3.  The ctor of the new vector evaluates the expr for each element.
// 4.  Each element invokes VectorDotTinyMatrix for one vector dot product.
// 5.  VectorDotTinyMatrix recurses, split the sum into halves, add results.
// 6.  When the length is one, it multiplies elements and returns that.
// 7.  Elements from the vector and TinyMatrix come through VectorElem and
//     TinyMatrixElem, so that the type of each one can be different.
//
//-----------------------------------------------------------------------------

//
// General TinyMatrixDotVector
//
// Much like VectorDotTinyMatrix above, this dots a vector with one row
// of a TinyMatrix.  It splits that dot product into two halves, adds the
// result and returns that.
//

template<class T1, class V2, int I, int B, int L>
struct TinyMatrixDotVector
{
  typedef typename TinyMatrixDotVector<T1,V2,I,B,L/2>::Type_t E1;
  typedef typename TinyMatrixDotVector<T1,V2,I,B+L/2,L-L/2>::Type_t E2;
  typedef typename BinaryReturn<E1,E2,OpAdd>::Type_t Type_t;
  static Type_t get(const T1& t1, const V2& v2)
    {
      return 
        TinyMatrixDotVector<T1,V2,I,B,L/2>::get(t1,v2) +
        TinyMatrixDotVector<T1,V2,I,B+L/2,L-L/2>::get(t1,v2);
    }
};

//
// Recursion termination for TinyMatrixDotVector
// When you get down to one term each from the vector and TinyMatrix, 
// just multiply the terms and return that.
//

template<class T1, class V2, int I, int B>
struct TinyMatrixDotVector<T1,V2,I,B,1>
{
  typedef typename TinyMatrixElem<T1,I,B>::Element_t E1;
  typedef typename VectorElem<V2,B>::Element_t E2;
  typedef typename BinaryReturn<E1,E2,OpMultiply>::Type_t Type_t;
  static Type_t get(const T1& t1, const V2& v2)
    {
      return TinyMatrixElem<T1,I,B>::get(t1) * VectorElem<V2,B>::get(v2);
    }
};

//
// Get one element from a vector of the dot product of a TinyMatrix and a vector.
// Gets one value by dotting the TinyMatrix and the vector.
//

template<int D1, int D2, class T1, class T2, class T3, class E1, class E2, int I>
struct VectorEngineElem<D1,T3,
  BinaryVectorOp<TinyMatrix<D1,D2,T1,E1>,Vector<D2,T2,E2>,FnDot>,I>
{
  typedef TinyMatrix<D1,D2,T1,E1> V1;
  typedef Vector<D2,T2,E2> V2;
  typedef VectorEngine<D1,T3,BinaryVectorOp<V1,V2,FnDot> > V;
  typedef typename TinyMatrixDotVector<V1,V2,I,0,D2>::Type_t T0;
  typedef T0 Element_t;
  typedef T0 ConstElementRef_t;
  typedef T0 ElementRef_t;
  static T0 get(const V& x) 
    {
      return TinyMatrixDotVector<V1,V2,I,0,D2>::get(x.v1_m,x.v2_m); 
    }
};

//
// Define the return type for dotting a TinyMatrix and a vector.
//

template<int D1, int D2, class T1, class T2, class E1, class E2>
struct BinaryReturn< TinyMatrix<D1,D2,T1,E1> , Vector<D2,T2,E2> , FnDot >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
  typedef Vector<D1,T0,Full> Type_t;
};

//
// Dot a TinyMatrix and a vector, returning a vector.
// Build a vector with a simple expression template representing 
// the dot product, and the ctor of the returned vector evaluates
// that expression.
//

template<int D1, int D2, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn< TinyMatrix<D1,D2,T1,E1> , Vector<D2,T2,E2>, FnDot >::Type_t
dot( const TinyMatrix<D1,D2,T1,E1>& v1 , const Vector<D2,T2,E2>& v2 )
{
  typedef TinyMatrix<D1,D2,T1,E1> V1;
  typedef Vector<D2,T2,E2> V2;
  typedef typename BinaryReturn<V1,V2,FnDot>::Type_t Return_t;
  typedef typename Return_t::Element_t T3;
  typedef Vector<D1,T3,BinaryVectorOp<V1,V2,FnDot> > Expr_t;
  return Return_t( Expr_t(v1,v2) );
}


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Subroutine for taking outper product between two vectors, yielding a TinyMatrix.
// Takes the outer product of vector of type V1 and vector of type V2
//
// The chain of operations is:
// 1.  The user calls outerProductAsTinyMatrix(vector,TinyMatrix)
//
//-----------------------------------------------------------------------------

//
// Define the return type for outerProductAsTinyMatrix(vector,vector).
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D,T1,E1> , Vector<D,T2,E2>, FnOuterProductAsTinyMatrix >
{
  typedef typename BinaryReturn<T1, T2, OpMultiply>::Type_t T0;
  typedef TinyMatrix<D, D, T0> Type_t;
};


//
// Take the outer product of two vectors, returning a TinyMatrix
//

template<int D, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn<Vector<D,T1,E1>, Vector<D,T2,E2>, FnOuterProductAsTinyMatrix>::Type_t
outerProductAsTinyMatrix(const Vector<D,T1,E1> &v1, const Vector<D,T2,E2> &v2 )
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

#endif // POOMA_TINY_VECTOR_TINYMATRIX_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: VectorTinyMatrix.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:11 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
