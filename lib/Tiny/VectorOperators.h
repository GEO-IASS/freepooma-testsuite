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

#ifndef POOMA_TINY_VECTOR_OPERATORS_H
#define POOMA_TINY_VECTOR_OPERATORS_H

//-----------------------------------------------------------------------------

// Class: External operators for Vector
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// The various arithmetic operators for vectors are defined here.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "PETE/TypeComputations.h"
#include "Tiny/BinaryVectorOp.h"
#include "Tiny/UnaryVectorOp.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class E> class Vector;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Much like the tensor operators, we define here the external operators
// between vectors.  The general strategy for how these work is to:
//
// 1. Define a function like operator-.
// 2. That returns a type found with UnaryReturn or BinaryReturn
// 3. The return is constructed by passing it a little parse tree for just
//    the unary or binary operator.  The form of that parse tree is
//    a vector with a BinaryVectorOp engine tag.
// 4. That parse tree is evaluated in the ctor for the new object
// 5. The parse tree is evaluated elementwise by evaluating the expression
//    vector at each point using VectorElem to do compile-time lookup.
//
//-----------------------------------------------------------------------------

//
// Unary operators.
// Two things need to be done for each unary operator:
//
// 1. Define the return type by specializing UnaryReturn.
//
// 2. Define the unary function which goes thru the steps above.
//
// These are all pretty much the same, so we define them using a macro.
//

#define POOMA_VECTOR_UNARY_OPERATOR(FUNC,TAG)                                 \
                                                                              \
template <int D, class T, class E>                                            \
struct UnaryReturn< Vector<D,T,E>, TAG >                                      \
{                                                                             \
  typedef Vector< D, typename UnaryReturn<T,TAG>::Type_t, E > Type_t;         \
};                                                                            \
                                                                              \
template <int D, class T, class E>                                            \
inline typename UnaryReturn< Vector<D,T,E>, TAG >::Type_t                     \
FUNC( const Vector<D,T,E>& v1 )                                               \
{                                                                             \
  typedef Vector<D,T,E> V1;                                                   \
  typedef typename UnaryReturn<T,TAG>::Type_t T3;                             \
  typedef Vector< D, T3, UnaryVectorOp<V1,TAG> > Expr_t;                      \
  typedef typename UnaryReturn<V1,TAG>::Type_t Return_t;                      \
  return Return_t( Expr_t(v1) );                                              \
}

POOMA_VECTOR_UNARY_OPERATOR(acos, FnArcCos)
POOMA_VECTOR_UNARY_OPERATOR(asin, FnArcSin)
POOMA_VECTOR_UNARY_OPERATOR(atan, FnArcTan)
POOMA_VECTOR_UNARY_OPERATOR(ceil, FnCeil)
POOMA_VECTOR_UNARY_OPERATOR(cos,  FnCos)
POOMA_VECTOR_UNARY_OPERATOR(cosh, FnHypCos)
POOMA_VECTOR_UNARY_OPERATOR(exp,  FnExp)
POOMA_VECTOR_UNARY_OPERATOR(fabs, FnFabs)
POOMA_VECTOR_UNARY_OPERATOR(floor, FnFloor)
POOMA_VECTOR_UNARY_OPERATOR(log, FnLog)
POOMA_VECTOR_UNARY_OPERATOR(log10, FnLog10)
POOMA_VECTOR_UNARY_OPERATOR(sin, FnSin)
POOMA_VECTOR_UNARY_OPERATOR(sinh, FnHypSin)
POOMA_VECTOR_UNARY_OPERATOR(sqrt, FnSqrt)
POOMA_VECTOR_UNARY_OPERATOR(tan, FnTan)
POOMA_VECTOR_UNARY_OPERATOR(tanh, FnHypTan)
POOMA_VECTOR_UNARY_OPERATOR(operator-, OpUnaryMinus)
POOMA_VECTOR_UNARY_OPERATOR(operator+, OpUnaryPlus)
POOMA_VECTOR_UNARY_OPERATOR(operator~, OpBitwiseNot)

//-----------------------------------------------------------------------------

//
// Binary operators.
// Two things need to be done for each binary operator:
//
// 1. Define the return type by specializing BinaryReturn.
//
// 2. Define the unary function which goes thru the steps above.
//
// We also have to define operations between scalars and vectors.  Any
// type that is not otherwise recognized is considered to be a scalar.
// These are all pretty much the same, so we define them using a macro.
//

#define POOMA_VECTOR_BINARY_OPERATOR(FUNC,TAG)                                \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
struct BinaryReturn< Vector<D,T1,E>, Vector<D,T2,E>, TAG >                    \
{                                                                             \
  typedef Vector< D, typename BinaryReturn<T1,T2,TAG>::Type_t, E > Type_t;    \
};                                                                            \
                                                                              \
template <int D, class T1, class T2, class E1, class E2>                      \
inline                                                                        \
typename BinaryReturn< Vector<D,T1,E1>, Vector<D,T2,E2>, TAG >::Type_t        \
FUNC( const Vector<D,T1,E1>& v1, const Vector<D,T2,E2>& v2 )                  \
{                                                                             \
  typedef Vector<D,T1,E1> V1;                                                 \
  typedef Vector<D,T2,E2> V2;                                                 \
  typedef typename BinaryReturn<V1,V2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef Vector< D, T3, BinaryVectorOp<V1,V2,TAG> > Expr_t;                  \
  return Return_t( Expr_t(v1,v2) );                                           \
}                                                                             \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
struct BinaryReturn< Vector<D,T1,E>, T2, TAG >                                \
{                                                                             \
  typedef Vector< D, typename BinaryReturn<T1,T2,TAG>::Type_t, E > Type_t;    \
};                                                                            \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
struct BinaryReturn< T1, Vector<D,T2,E>, TAG >                                \
{                                                                             \
  typedef Vector< D, typename BinaryReturn<T1,T2,TAG>::Type_t, E > Type_t;    \
};                                                                            \
                                                                              \
template <int D, class T1, class E, class T2>                                 \
inline typename BinaryReturn< Vector<D,T1,E>, T2, TAG >::Type_t               \
FUNC( const Vector<D,T1,E>& v1, const T2& x )                                 \
{                                                                             \
  typedef Vector<D,T1,E> V1;                                                  \
  typedef typename BinaryReturn<V1,T2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef Vector< D, T3, BinaryVectorOp<V1,T2,TAG> > Expr_t;                  \
  return Return_t( Expr_t(v1,x) );                                            \
}                                                                             \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
inline typename BinaryReturn< T1, Vector<D,T2,E>, TAG >::Type_t               \
FUNC( const T1& x, const Vector<D,T2,E>& v2 )                                 \
{                                                                             \
  typedef Vector<D,T2,E> V2;                                                  \
  typedef typename BinaryReturn<T1,V2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef Vector< D, T3, BinaryVectorOp<T1,V2,TAG> > Expr_t;                  \
  return Return_t( Expr_t(x,v2) );                                            \
}

//-----------------------------------------------------------------------------

POOMA_VECTOR_BINARY_OPERATOR(operator+,OpAdd)
POOMA_VECTOR_BINARY_OPERATOR(operator-,OpSubtract)
POOMA_VECTOR_BINARY_OPERATOR(operator*,OpMultiply)
POOMA_VECTOR_BINARY_OPERATOR(operator/,OpDivide)
POOMA_VECTOR_BINARY_OPERATOR(operator%, OpMod)
POOMA_VECTOR_BINARY_OPERATOR(operator&, OpBitwiseAnd)
POOMA_VECTOR_BINARY_OPERATOR(operator|, OpBitwiseOr)
POOMA_VECTOR_BINARY_OPERATOR(operator^, OpBitwiseXor)
POOMA_VECTOR_BINARY_OPERATOR(ldexp, FnLdexp)
POOMA_VECTOR_BINARY_OPERATOR(pow, FnPow)
POOMA_VECTOR_BINARY_OPERATOR(fmod, FnFmod)
POOMA_VECTOR_BINARY_OPERATOR(atan2, FnArcTan2)

//-----------------------------------------------------------------------------
//
// Dot product between two vectors.
//
// Like the functions above, two things need to be done:
//
// 1. Define the return type by specializing BinaryReturn.
// 
// 2. Define the dot product function.
//
// Dot product uses template metaprograms to return a scalar.
//
//-----------------------------------------------------------------------------

//
// General form of VectorDotVector.
// Inputs vectors of type V1 and V2, which element B is the first
// component in the dot product, and L the number of elements in the
// dot product.
//
// Recurses by cutting the interval in half, computing the dot
// product on each and returning the sum of them.
//

template<class V1, class V2, int B, int L>
struct VectorDotVector
{
  typedef typename VectorDotVector<V1,V2,B,L/2>::Type_t Left_t;
  typedef typename VectorDotVector<V1,V2,B+L/2,L-L/2>::Type_t Right_t;
  typedef typename BinaryReturn<Left_t,Right_t,OpAdd>::Type_t Type_t;
  static Type_t get(const V1& x, const V2& y)
    {
      return 
        VectorDotVector<V1,V2,B,L/2>::get(x,y) + 
        VectorDotVector<V1,V2,B+L/2,L-L/2>::get(x,y);
    }
};

//
// Recursion termination for VectorDotVector.
// When the length of the dot product is 1, just multiply those terms
// and return.
//

template<class V1, class V2, int B>
struct VectorDotVector<V1,V2,B,1>
{
  typedef typename VectorElem<V1,B>::Element_t Left_t;
  typedef typename VectorElem<V2,B>::Element_t Right_t;
  typedef typename BinaryReturn<Left_t,Right_t,OpMultiply>::Type_t Type_t;
  static Type_t get(const V1& x, const V2& y)
    {
      return VectorElem<V1,B>::get(x) * VectorElem<V2,B>::get(y);
    }
};

//
// Define the return type for dotting two vectors together.
// The type comes from the result of multiplying T1 and T2 together.
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D,T1,E1> , Vector<D,T2,E2> , FnDot >
{
  typedef Vector<D,T1,E1> V1;
  typedef Vector<D,T2,E2> V2;
  typedef typename VectorDotVector<V1,V2,0,D>::Type_t T0;
  typedef T0 Type_t;
};

//
// The actual dot product function.
// Just call VectorDotVector.
//

template<int D, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< Vector<D,T1,E1>,Vector<D,T2,E2> , FnDot >::Type_t
dot( const Vector<D,T1,E1>& v1, const Vector<D,T2,E2>& v2 )
{
  return VectorDotVector<Vector<D,T1,E1>,Vector<D,T2,E2>,0,D>::get(v1,v2);
}

//-----------------------------------------------------------------------------
// norm, norm2
//
// The norm of any vector is just the square root of the dot product of
// a vector with itself. The norm squared is the dot product with itself.
//
// NOTE: these don't work for Vector<complex<>>!
//-----------------------------------------------------------------------------

template<int D, class T, class E>
struct UnaryReturn< Vector<D,T,E> , FnNorm >
{
  typedef T Type_t;
};

template<int D, class T, class E>
inline T
norm(const Vector<D,T,E>& x)
{
  return sqrt(dot(x,x));
}

template<int D, class T, class E>
inline T
norm2(const Vector<D,T,E>& x)
{
  return dot(x,x);
}

//-----------------------------------------------------------------------------
//
// Equality operator between two vectors.
//
// Like the functions above, two things need to be done:
//
// 1. Define the return type by specializing BinaryReturn.
// 
// 2. Define the equality operation.
//
// Uses template metaprograms to return a bool.
//
//-----------------------------------------------------------------------------

//
// General form of VectorEqualsVector.
// Inputs vectors of type V1 and V2, which element B is the first
// component in the dot product, and L the number of elements in the
// dot product.
//
// Recurses by cutting the interval in half, computing the dot
// product on each and returning the sum of them.
//

template<class V1, class V2, int B, int L>
struct VectorEqualsVector
{
  typedef typename VectorEqualsVector<V1,V2,B,L/2>::Type_t Left_t;
  typedef typename VectorEqualsVector<V1,V2,B+L/2,L-L/2>::Type_t Right_t;
  typedef bool Type_t;
  static Type_t get(const V1& x, const V2& y)
    {
      return 
        VectorEqualsVector<V1,V2,B,L/2>::get(x,y) && 
        VectorEqualsVector<V1,V2,B+L/2,L-L/2>::get(x,y);
    }
};

//
// Recursion termination for VectorEqualsVector.
// When the length of the equality check is 1, just compare those terms
// and return.
//

template<class V1, class V2, int B>
struct VectorEqualsVector<V1,V2,B,1>
{
  typedef typename VectorElem<V1,B>::Element_t Left_t;
  typedef typename VectorElem<V2,B>::Element_t Right_t;
  typedef bool Type_t;
  static Type_t get(const V1& x, const V2& y)
    {
      return VectorElem<V1,B>::get(x) == VectorElem<V2,B>::get(y);
    }
};

//
// Define the return type for checking equality of two vectors.
// Just a bool.
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D,T1,E1> , Vector<D,T2,E2> , OpEQ >
{
  typedef bool Type_t;
};

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Vector<D,T1,E1> , Vector<D,T2,E2> , OpNE >
{
  typedef bool Type_t;
};

//
// The actual equality operators.
// Just call VectorEqualsVector for == and negate result for !=.
//

template<int D, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< Vector<D,T1,E1>,Vector<D,T2,E2> , OpEQ >::Type_t
operator==(const Vector<D,T1,E1>& v1, const Vector<D,T2,E2>& v2)
{
  return VectorEqualsVector<Vector<D,T1,E1>,Vector<D,T2,E2>,0,D>::get(v1,v2);
}

template<int D, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< Vector<D,T1,E1>,Vector<D,T2,E2> , OpNE >::Type_t
operator!=(const Vector<D,T1,E1>& v1, const Vector<D,T2,E2>& v2)
{
  return !(v1 == v2);
}

//
// This version is here to work around compilers that don't hide
// STL internal relops inside a namespace.
// (They define op!=(T,T), so the above version becomes ambiguous.)
//

template<int D, class T, class E>
inline typename 
BinaryReturn<Vector<D, T, E>, Vector<D, T, E> , OpNE>::Type_t
operator!=(const Vector<D, T, E>& v1, const Vector<D, T, E>& v2)
{
  return !(v1 == v2);
}

//-----------------------------------------------------------------------------
//
// General assignment (including things like +=) for vectors.
//
// There are two things that need to be done:
//
// 1.  Define the return type by specializing BinaryReturn.  
//     This is actually done automatically in PETE.
//
// 2.  Define the operator to call VectorAssign.
//
// These are all pretty similar, so we use a macro to define it.
// Types not recognized as vectors are considered to be scalars.
//
//-----------------------------------------------------------------------------

#define POOMA_VECTOR_ACCUM_OPERATOR(FUNC,TAG)                                 \
                                                                              \
template <int D, class T1, class T2, class E1, class E2>                      \
inline Vector<D,T1,E1>&                                                       \
FUNC( Vector<D,T1,E1>& v1, const Vector<D,T2,E2>& v2 )                        \
{                                                                             \
  VectorAssign<Vector<D,T1,E1>,Vector<D,T2,E2>,TAG,0,D>::apply(v1,v2,TAG());  \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2, class E1>                                \
inline Vector<D,T1,E1>&                                                       \
FUNC( Vector<D,T1,E1>& v1, const T2& v2 )                                     \
{                                                                             \
  VectorAssign<Vector<D,T1,E1>,T2,TAG,0,D>::apply(v1,v2,TAG());               \
  return v1;                                                                  \
}

POOMA_VECTOR_ACCUM_OPERATOR(operator+=, OpAddAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator-=, OpSubtractAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator*=, OpMultiplyAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator/=, OpDivideAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator%=, OpModAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator|=, OpBitwiseOrAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator&=, OpBitwiseAndAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator^=, OpBitwiseXorAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator<<=, OpLeftShiftAssign)
POOMA_VECTOR_ACCUM_OPERATOR(operator>>=, OpRightShiftAssign)

#endif // POOMA_TINY_VECTOR_OPERATORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: VectorOperators.h,v $   $Author: richard $
// $Revision: 1.19 $   $Date: 2004/11/01 18:17:11 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
