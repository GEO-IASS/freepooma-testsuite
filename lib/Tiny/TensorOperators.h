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

#ifndef POOMA_TINY_TENSOR_OPERATORS_H
#define POOMA_TINY_TENSOR_OPERATORS_H

//-----------------------------------------------------------------------------
// Class: External operators for Tensor
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// The various arithemetic operators for Tensors.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "PETE/TypeComputations.h"
#include "Tiny/BinaryTensorOp.h"
#include "Tiny/UnaryTensorOp.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class E> class Tensor;
class Full;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// The operators for combinations of Tensors.
//
// These are just overloaded operators.  These don't do
// expression templates, but they do one tricky thing.  In order
// to build the returned temporary object, they build a proxy object
// for the expression and construct the new object using that.  This
// means that the operator is actually evaluated in the ctor of the
// new object.  This makes it easier to avoid unnecessary temporaries.
//
//-----------------------------------------------------------------------------

//
// Unary operators.
// For each operator we need to do two things:
//
// 1. Define the return type by specializing UnaryReturn
// 
// 2. Build the overloaded function which does the operation.
//
// These things are all pretty darn similar, so we can define a macro
// for the definition.
//

#define POOMA_TENSOR_UNARY_OPERATOR(FUNC,TAG)                                 \
                                                                              \
template <int D, class T, class E>                                            \
struct UnaryReturn< Tensor<D,T,E>, TAG >                                      \
{                                                                             \
  typedef Tensor< D, typename UnaryReturn<T,TAG>::Type_t, E > Type_t;         \
};                                                                            \
                                                                              \
template <int D, class T, class E>                                            \
inline typename UnaryReturn< Tensor<D,T,E>, TAG >::Type_t                     \
FUNC( const Tensor<D,T,E>& v1 )                                               \
{                                                                             \
  typedef Tensor<D,T,E> V1;                                                   \
  typedef typename UnaryReturn<T,TAG>::Type_t T3;                             \
  typedef Tensor< D, T3, UnaryTensorOp<V1,TAG> > Expr_t;                      \
  typedef typename UnaryReturn<V1,TAG>::Type_t Return_t;                      \
  return Return_t( Expr_t(v1) );                                              \
}

//-----------------------------------------------------------------------------

POOMA_TENSOR_UNARY_OPERATOR(acos, FnArcCos)
POOMA_TENSOR_UNARY_OPERATOR(asin, FnArcSin)
POOMA_TENSOR_UNARY_OPERATOR(atan, FnArcTan)
POOMA_TENSOR_UNARY_OPERATOR(ceil, FnCeil)
POOMA_TENSOR_UNARY_OPERATOR(cos,  FnCos)
POOMA_TENSOR_UNARY_OPERATOR(cosh, FnHypCos)
POOMA_TENSOR_UNARY_OPERATOR(exp,  FnExp)
POOMA_TENSOR_UNARY_OPERATOR(fabs, FnFabs)
POOMA_TENSOR_UNARY_OPERATOR(floor, FnFloor)
POOMA_TENSOR_UNARY_OPERATOR(log, FnLog)
POOMA_TENSOR_UNARY_OPERATOR(log10, FnLog10)
POOMA_TENSOR_UNARY_OPERATOR(sin, FnSin)
POOMA_TENSOR_UNARY_OPERATOR(sinh, FnHypSin)
POOMA_TENSOR_UNARY_OPERATOR(sqrt, FnSqrt)
POOMA_TENSOR_UNARY_OPERATOR(tan, FnTan)
POOMA_TENSOR_UNARY_OPERATOR(tanh, FnHypTan)
POOMA_TENSOR_UNARY_OPERATOR(operator-, OpUnaryMinus)
POOMA_TENSOR_UNARY_OPERATOR(operator+, OpUnaryPlus)
POOMA_TENSOR_UNARY_OPERATOR(operator~, OpBitwiseNot)

//-----------------------------------------------------------------------------

//
// Elementwise binary operators on Tensors.
//
// Like the unary operators above, two things need to be done:
//
// 1. Define the return type by specializing BinaryReturn.
//
// 2. Define the function.
//
// These are all pretty darn similar, since they just apply an operation
// elementwise.  We define a macro given the name of the function and the
// operator tag, and use that for all the elementwise operators.
//
// Unlike the unary operators above, these need to be able to handle 
// operations with a tensor and a scalar.  Any type that isn't a tensor
// is considered a scalar.
// 
// TJW:
// Things change a bit when you allow the existence of non-Full engines. Now
// we have to account for the tensor engine type of the result of combining
// tensors with different engine types. Basically, you always return something
// as general or more general than the two inputs.
// Tensor<...,Full> + Tensor<...,Symmetric> would return Tensor<...,Full>, 
// for example. Adding Symmetric+Antisymmetric would have to return Full, as
// would adding Symmetric+Diagonal. In fact, adding any heterogeneous types
// should always force return of a Full tensor as far as I (TJW) know, so let's
// make it so. Only adding two non-Full tensors can reliably yield a non-Full.


#define POOMA_TENSOR_BINARY_OPERATOR(FUNC,TAG)                                \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
struct BinaryReturn< Tensor<D,T1,E>, Tensor<D,T2,E>, TAG >                    \
{                                                                             \
  typedef Tensor< D, typename BinaryReturn<T1,T2,TAG>::Type_t, E > Type_t;    \
};                                                                            \
                                                                              \
template <int D, class T1, class T2, class E1, class E2>                      \
struct BinaryReturn< Tensor<D,T1,E1>, Tensor<D,T2,E2>, TAG >                  \
{                                                                             \
  typedef Tensor< D, typename BinaryReturn<T1,T2,TAG>::Type_t, Full > Type_t; \
};                                                                            \
                                                                              \
template <int D, class T1, class T2, class E1, class E2>                      \
inline                                                                        \
typename BinaryReturn< Tensor<D,T1,E1>, Tensor<D,T2,E2>, TAG >::Type_t        \
FUNC( const Tensor<D,T1,E1>& v1, const Tensor<D,T2,E2>& v2 )                  \
{                                                                             \
  typedef Tensor<D,T1,E1> V1;                                                 \
  typedef Tensor<D,T2,E2> V2;                                                 \
  typedef typename BinaryReturn<T1,T2,TAG>::Type_t T3;                        \
  typedef Tensor< D, T3, BinaryTensorOp<V1,V2,TAG> > Expr_t;                  \
  typedef typename BinaryReturn<V1,V2,TAG>::Type_t Return_t;                  \
  return Return_t( Expr_t(v1,v2) );                                           \
}                                                                             \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
struct BinaryReturn< Tensor<D,T1,E>, T2, TAG >                                \
{                                                                             \
  typedef Tensor< D, typename BinaryReturn<T1,T2,TAG>::Type_t, E > Type_t;    \
};                                                                            \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
struct BinaryReturn< T1, Tensor<D,T2,E>, TAG >                                \
{                                                                             \
  typedef Tensor< D, typename BinaryReturn<T1,T2,TAG>::Type_t, E > Type_t;    \
};                                                                            \
                                                                              \
template <int D, class T1, class E, class T2>                                 \
inline typename BinaryReturn< Tensor<D,T1,E>, T2, TAG >::Type_t               \
FUNC( const Tensor<D,T1,E>& v1, const T2& x )                                 \
{                                                                             \
  typedef Tensor<D,T1,E> V1;                                                  \
  typedef typename BinaryReturn<V1,T2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef Tensor< D, T3, BinaryTensorOp<V1,T2,TAG> > Expr_t;                  \
  return Return_t( Expr_t(v1,x) );                                            \
}                                                                             \
                                                                              \
template <int D, class T1, class T2, class E>                                 \
inline typename BinaryReturn< T1, Tensor<D,T2,E>, TAG >::Type_t               \
FUNC( const T1& x, const Tensor<D,T2,E>& v2 )                                 \
{                                                                             \
  typedef Tensor<D,T2,E> V2;                                                  \
  typedef typename BinaryReturn<T1,V2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef Tensor< D, T3, BinaryTensorOp<T1,V2,TAG> > Expr_t;                  \
  return Return_t( Expr_t(x,v2) );                                            \
}

//-----------------------------------------------------------------------------

POOMA_TENSOR_BINARY_OPERATOR(operator+, OpAdd)
POOMA_TENSOR_BINARY_OPERATOR(operator-, OpSubtract)
POOMA_TENSOR_BINARY_OPERATOR(operator*, OpMultiply)
POOMA_TENSOR_BINARY_OPERATOR(operator/, OpDivide)
POOMA_TENSOR_BINARY_OPERATOR(operator%, OpMod)
POOMA_TENSOR_BINARY_OPERATOR(operator&, OpBitwiseAnd)
POOMA_TENSOR_BINARY_OPERATOR(operator|, OpBitwiseOr)
POOMA_TENSOR_BINARY_OPERATOR(operator^, OpBitwiseXor)
POOMA_TENSOR_BINARY_OPERATOR(ldexp,     FnLdexp)
POOMA_TENSOR_BINARY_OPERATOR(pow,       FnPow)
POOMA_TENSOR_BINARY_OPERATOR(fmod,      FnFmod)
POOMA_TENSOR_BINARY_OPERATOR(atan2,     FnArcTan2)

//-----------------------------------------------------------------------------
//
// General form of TensorDotTensor.
//
// Finds one term of the dot product of a tensor dotted with a tensor.
// Inputs a tensor of type T1 and a tensor of type T2 and dots row
// I of T1 with column J of T2, starting with offset K in that row and
// vector length L.
//
// Operates by dividing that row in half, dotting the halves, and 
// returning the sum of those two.
//

template<class T1, class T2, int I, int J, int K, int L>
struct TensorDotTensor
{
  typedef typename TensorDotTensor<T1,T2,I,J,K,L/2>::Type_t Left_t;
  typedef typename TensorDotTensor<T1,T2,I,J,K+L/2,L-L/2>::Type_t Right_t;
  typedef typename BinaryReturn<Left_t,Right_t,OpAdd>::Type_t Type_t;
  static Type_t get(const T1& x, const T2& y)
    {
      return 
        TensorDotTensor<T1,T2,I,J,K,L/2>::get(x,y) + 
        TensorDotTensor<T1,T2,I,J,K+L/2,L-L/2>::get(x,y);
    }
};

//
// Recursion termination for TensorDotTensor.
// When you get down to a vector length of 1, all you need to do
// is multiply the elements and return that.
//

template<class T1, class T2, int I, int J, int K>
struct TensorDotTensor<T1,T2,I,J,K,1>
{
  typedef typename TensorElem<T1,I,K>::Element_t Left_t;
  typedef typename TensorElem<T2,K,J>::Element_t Right_t;
  typedef typename BinaryReturn<Left_t,Right_t,OpMultiply>::Type_t Type_t;
  static Type_t get(const T1& x, const T2& y)
    {
      return TensorElem<T1,I,K>::get(x) * TensorElem<T2,K,J>::get(y);
    }
};

//
// Specialization of TensorElem for getting one value from
// a TensorEngine for dotting two tensors.
// Just calls TensorDotTensor.
//

template<int D, class T, class T1, class T2, int I, int J>
struct TensorEngineElem<D,T,BinaryTensorOp<T1,T2,FnDot>, I, J,
  1> // TJW added true(or, rather, int value 1)
{
  typedef BinaryTensorOp<T1,T2,FnDot> E;
  typedef TensorEngine<D,T,E> T0;
  typedef typename TensorDotTensor<T1,T2,I,J,0,T1::d>::Type_t Element_t;
  typedef Element_t ElementRef_t;
  typedef Element_t ConstElementRef_t;
  static Element_t get(const T0& x) 
    { 
      return TensorDotTensor<T1,T2,I,J,0,T1::d>::get(x.v1_m,x.v2_m);
    }
};

//
// Define the return type for tensor dot tensor.
// Since this has no knowledge of the input tensor types, 
// this returns a full tensor.  Specialized tensor types should
// define specialized versions of this.
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Tensor<D,T1,E1> , Tensor<D,T2,E2> , FnDot >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
  typedef Tensor<D,T0,Full> Type_t;
};

//
// Dot product of two tensors.
// Uses all of the above to build a new tensor.
//

template<int D, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn< Tensor<D,T1,E1>,Tensor<D,T2,E2> , FnDot >::Type_t
dot( const Tensor<D,T1,E1>& v1, const Tensor<D,T2,E2>& v2 )
{
  typedef Tensor<D,T1,E1> V1;
  typedef Tensor<D,T2,E2> V2;
  typedef typename BinaryReturn<V1,V2,FnDot>::Type_t Return_t;
  typedef typename Return_t::Element_t T3;
  typedef Tensor<D,T3,BinaryTensorOp<V1,V2,FnDot> > Expr_t;
  return Return_t( Expr_t(v1,v2) );
}

//-----------------------------------------------------------------------------
// trace()
//
// The sum of the diagonal elements of the Tensor.
//
// Implementation note: these could probably be implemented using template
// metaprograms rather than runtime loops. Metaprogram versions would
// presumably run faster.
//-----------------------------------------------------------------------------

// 1. Define the return type of trace():

template<int D, class T, class E>
struct UnaryReturn< Tensor<D,T,E> , FnTrace >
{
  typedef T Type_t;
};


// 2. Define the trace() functions:

// General template: loop over diagonal elements using (i,j) access:
template<int D, class T, class E>
inline T
trace(const Tensor<D,T,E>& t)
{
  T result(0.0);
  for (int d = 0; d < D; d++) {
    result += t(d,d);
  }
  return result;
}

// Some partial specializations, for efficiency where the trace can be
// computed more quickly:

// Diagonal tensor: loop over stored elements using 1D indexing:
template<int D, class T>
inline T
trace(const Tensor<D,T,Diagonal>& t)
{
  T result(0.0);
  for (int d = 0; d < D; d++) {
    result += t(d);
  }
  return result;
}

// Antisymmetric tensor: trace is zero:
template<int D, class T>
inline T
trace(const Tensor<D,T,Antisymmetric>& t)
{
  return T(0.0);
}


//-----------------------------------------------------------------------------
// det()
//
// The determinant of the Tensor, viewed as a matrix
//
// Implementation note: these could probably be implemented using template
// metaprograms rather than runtime loops. Metaprogram versions would
// presumably run faster.
//-----------------------------------------------------------------------------

// 1. Define the return type of det():

template<int D, class T, class E>
struct UnaryReturn< Tensor<D,T,E> , FnDet >
{
  typedef T Type_t;
};


// 2. Define the det() functions:

// General template: loop over diagonal elements using (i,j) access:
// For now, skip coding up for D>3; should add this later.
template<int D, class T, class E>
inline T
det(const Tensor<D,T,E>& t)
{
  PInsist(D<4, "Tensor det() function not implemented for D>3!");
  return T(-999999.999999);
}
// Partial specializations for D={1,2,3}
template<class T, class E>
inline T
det(const Tensor<3,T,E>& t)
{
  T result;
  result = 
    t(0,0)*(t(1,1)*t(2,2) - t(1,2)*t(2,1)) +
    t(0,1)*(t(1,2)*t(2,0) - t(1,0)*t(2,2)) +
    t(0,2)*(t(1,0)*t(2,1) - t(1,1)*t(2,0));
  return result;
}
template<class T, class E>
inline T
det(const Tensor<2,T,E>& t)
{
  T result;
  result = t(0,0)*t(1,1) - t(0,1)*t(1,0);
  return result;
}
template<class T, class E>
inline T
det(const Tensor<1,T,E>& t)
{
  T result = t(0,0);
  return result;
}

// Some EngineTag partial specializations, for efficiency where the det can be
// computed more quickly:

// Diagonal tensor: loop over stored elements using 1D indexing:
// Partial specializations for D={1,2,3}
template<class T>
inline T
det(const Tensor<3,T,Diagonal>& t)
{
  T result;
  result = t(0)*t(1)*t(2);
  return result;
}
template<class T>
inline T
det(const Tensor<2,T,Diagonal>& t)
{
  T result;
  result = t(0)*t(1);
  return result;
}
template<class T>
inline T
det(const Tensor<1,T,Diagonal>& t)
{
  T result = t(0);
  return result;
}

// Antisymmetric tensor:
// For D=3, det is zero, because diagonal elements are zero:
template<class T>
inline T
det(const Tensor<3,T,Antisymmetric>& t)
{
  return T(0.0);
}
// For D=2, det is nonzero; use linear indexing of only stored element:
template<class T>
inline T
det(const Tensor<2,T,Antisymmetric>& t)
{
  T result;
  result = t(0)*t(0);
  return result;
}
// For D=1, det is zero, because diagonal elements are zero:
template<class T>
inline T
det(const Tensor<1,T,Antisymmetric>& t)
{
  return T(0.0);
}


//-----------------------------------------------------------------------------
// transpose()
//
// The "matrix" transpose of the Tensor.
//
// Implementation note: these could probably be implemented using template
// metaprograms rather than runtime loops. Metaprogram versions would
// presumably run faster.
//-----------------------------------------------------------------------------

// 1. Define the return type of transpose():

template<int D, class T, class E>
struct UnaryReturn< Tensor<D,T,E> , FnTranspose >
{
  typedef Tensor<D,T,E> Type_t;
};


// 2. Define the transpose() functions:

// General template: loop over diagonal elements using (i,j) access:
template<int D, class T, class E>
inline Tensor<D,T,E>
transpose(const Tensor<D,T,E>& t)
{
  Tensor<D,T,E> result;
  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++) {
      result(i,j) = t(j,i);
    }
  }
  return result;
}

// Some EngineTag partial specializations, for efficiency where the transpose
// can be computed more quickly:

// Symmetric tensor: transpose = tensor:
template<int D, class T>
inline Tensor<D,T,Symmetric>
transpose(const Tensor<D,T,Symmetric>& t)
{
  Tensor<D,T,Symmetric> result;
  result = t;
  return result;
}

// Antisymmetric tensor: transpose = -tensor:
template<int D, class T>
inline Tensor<D,T,Antisymmetric>
transpose(const Tensor<D,T,Antisymmetric>& t)
{
  Tensor<D,T,Antisymmetric> result;
  result = -t;
  return result;
}

// Diagonal tensor: transpose = tensor:
template<int D, class T>
inline Tensor<D,T,Diagonal>
transpose(const Tensor<D,T,Diagonal>& t)
{
  Tensor<D,T,Diagonal> result;
  result = t;
  return result;
}


//-----------------------------------------------------------------------------
//
// Equality operator between two tensors.
//
// Like the functions above, two things need to be done:
//
// 1. Define the return type by specializing BinaryReturn.
// 
// 2. Define the equality operation.
//
// Tensors are sufficiently complicated that we just bag the metraprograms
// and write a loop.
//
//-----------------------------------------------------------------------------

//
// Define the return type for checking equality of two tensors.
// Just a bool.
//

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Tensor<D,T1,E1>,Tensor<D,T2,E2> , OpEQ >
{
  typedef bool Type_t;
};

template<int D, class T1, class T2, class E1, class E2>
struct BinaryReturn< Tensor<D,T1,E1>,Tensor<D,T2,E2> , OpNE >
{
  typedef bool Type_t;
};

//
// The actual equality operators.
// Just do the stinkin' loops.
//

template<int D, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< Tensor<D,T1,E1>, Tensor<D,T2,E2>, OpEQ >::Type_t
operator==(const Tensor<D,T1,E1>& t1, const Tensor<D,T2,E2>& t2)
{
  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      if (t1(i,j) != t2(i,j)) return false;
   
  return true;
}

template<int D, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< Tensor<D,T1,E1>, Tensor<D,T2,E2>, OpNE >::Type_t
operator!=(const Tensor<D,T1,E1>& t1, const Tensor<D,T2,E2>& t2)
{
  return !(t1 == t2);
}

//-----------------------------------------------------------------------------
//
// General assignment (including things like +=) for tensors.
//
// There are two things that need to be done:
//
// 1.  Define the return type by specializing BinaryReturn.  
//     This is actually done automatically in PETE.
//
// 2.  Define the operator to call TensorAssign.
//
// These are all pretty similar, so we use a macro to define it.
// Types not recognized as tensors are considered to be scalars.
//
// TJW:
// Things change a bit when you allow the existence of non-Full engines. Now
// we have to account for the tensor engine type of the result of combining
// tensors with different engine types. Basically, you always return something
// as general or more general than the two inputs.  Tensor<...,Full> +
// Tensor<...,Symmetric> would return Tensor<...,Full>, for example. Adding
// Symmetric+Antisymmetric would have to return Full, as would adding
// Symmetric+Diagonal. In fact, combining any heterogeneous types should
// always force return of a Full tensor as far as I (TJW) know, so let's make
// it so. Only adding two non-Full tensors can reliably yield a non-Full. Even
// then, for non-additivie operations like operator*=(), combining two
// non-Full tensors of the same type can result in a tensor with a different
// engine type.
// 
// ... Aaugh! Can't do this, because global reductions on Arrays and Fields of
// these types assume the return value of the reduction is the original
// type. For now, at least, bag everything except operator+= and operator-=,
// for which the return type is the same as the input type. Need to use the
// specialized TensorAssign versions for the non-Full engines. --TJW
//-----------------------------------------------------------------------------

#define POOMA_TENSOR_ACCUM_OPERATOR(FUNC,TAG)                                 \
                                                                              \
template <int D, class T1, class T2, class E1>                                \
inline Tensor<D,T1,E1>&                                                       \
FUNC( Tensor<D,T1,E1>& v1, const Tensor<D,T2,E1>& v2 )                        \
{                                                                             \
  typedef typename Tensor<D,T1,E1>::Engine_t Left_t;                          \
  typedef typename Tensor<D,T2,E1>::Engine_t Right_t;                         \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2>                                          \
inline Tensor<D,T1,Symmetric>&                                                \
FUNC( Tensor<D,T1,Symmetric>& v1, const Tensor<D,T2,Symmetric>& v2 )          \
{                                                                             \
  typedef typename Tensor<D,T1,Symmetric>::Engine_t Left_t;                   \
  typedef typename Tensor<D,T2,Symmetric>::Engine_t Right_t;                  \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2>                                          \
inline Tensor<D,T1,Antisymmetric>&                                            \
FUNC( Tensor<D,T1,Antisymmetric>& v1, const Tensor<D,T2,Antisymmetric>& v2 )  \
{                                                                             \
  typedef typename Tensor<D,T1,Antisymmetric>::Engine_t Left_t;               \
  typedef typename Tensor<D,T2,Antisymmetric>::Engine_t Right_t;              \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2>                                          \
inline Tensor<D,T1,Diagonal>&                                                 \
FUNC( Tensor<D,T1,Diagonal>& v1, const Tensor<D,T2,Diagonal>& v2 )            \
{                                                                             \
  typedef typename Tensor<D,T1,Diagonal>::Engine_t Left_t;                    \
  typedef typename Tensor<D,T2,Diagonal>::Engine_t Right_t;                   \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2, class E1>                                \
inline Tensor<D,T1,E1>&                                                       \
FUNC( Tensor<D,T1,E1>& v1, const T2& v2 )                                     \
{                                                                             \
  typedef typename Tensor<D,T1,E1>::Engine_t Left_t;                          \
  typedef typename T2::Engine_t Right_t;                                      \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2>                                          \
inline Tensor<D,T1,Symmetric>&                                                \
FUNC( Tensor<D,T1,Symmetric>& v1, const T2& v2 )                              \
{                                                                             \
  typedef typename Tensor<D,T1,Symmetric>::Engine_t Left_t;                   \
  typedef typename T2::Engine_t Right_t;                                      \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2>                                          \
inline Tensor<D,T1,Antisymmetric>&                                            \
FUNC( Tensor<D,T1,Antisymmetric>& v1, const T2& v2 )                          \
{                                                                             \
  typedef typename Tensor<D,T1,Antisymmetric>::Engine_t Left_t;               \
  typedef typename T2::Engine_t Right_t;                                      \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D, class T1, class T2>                                          \
inline Tensor<D,T1,Diagonal>&                                                 \
FUNC( Tensor<D,T1,Diagonal>& v1, const T2& v2 )                               \
{                                                                             \
  typedef typename Tensor<D,T1,Diagonal>::Engine_t Left_t;                    \
  typedef typename T2::Engine_t Right_t;                                      \
  TensorAssign<Left_t,Right_t,TAG,0,D,0,D>::                                  \
    apply(v1.engine(),v2.engine(),TAG());                                     \
  return v1;                                                                  \
}

POOMA_TENSOR_ACCUM_OPERATOR(operator+=, OpAddAssign)
POOMA_TENSOR_ACCUM_OPERATOR(operator-=, OpSubtractAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator*=, OpMultiplyAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator/=, OpDivideAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator%=, OpModAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator|=, OpBitwiseOrAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator&=, OpBitwiseAndAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator^=, OpBitwiseXorAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator<<=, OpLeftShiftAssign)
// POOMA_TENSOR_ACCUM_OPERATOR(operator>>=, OpRightShiftAssign)

#endif // POOMA_TINY_TENSOR_OPERATORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TensorOperators.h,v $   $Author: richi $
// $Revision: 1.24 $   $Date: 2004/11/29 14:03:30 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
