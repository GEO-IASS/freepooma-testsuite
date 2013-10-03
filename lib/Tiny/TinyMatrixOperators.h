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

#ifndef POOMA_TINY_TINYMATRIX_OPERATORS_H
#define POOMA_TINY_TINYMATRIX_OPERATORS_H

//-----------------------------------------------------------------------------

// Class: External operators for TinyMatrix
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// The various arithemetic operators for TinyMatrixs.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "PETE/TypeComputations.h"
#include "Tiny/BinaryTinyMatrixOp.h"
#include "Tiny/UnaryTinyMatrixOp.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D1, int D2, class T, class E> class TinyMatrix;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// The operators for combinations of TinyMatrixs.
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

#define POOMA_TINYMATRIX_UNARY_OPERATOR(FUNC,TAG)                             \
                                                                              \
template <int D1, int D2, class T, class E>                                   \
struct UnaryReturn< TinyMatrix<D1,D2,T,E>, TAG >                              \
{                                                                             \
  typedef TinyMatrix< D1, D2, typename UnaryReturn<T,TAG>::Type_t, E >        \
    Type_t;                                                                   \
};                                                                            \
                                                                              \
template <int D1, int D2, class T, class E>                                   \
inline typename UnaryReturn< TinyMatrix<D1,D2,T,E>, TAG >::Type_t             \
FUNC( const TinyMatrix<D1,D2,T,E>& v1 )                                       \
{                                                                             \
  typedef TinyMatrix<D1,D2,T,E> V1;                                           \
  typedef typename UnaryReturn<T,TAG>::Type_t T3;                             \
  typedef TinyMatrix< D1, D2, T3, UnaryTinyMatrixOp<V1,TAG> > Expr_t;         \
  typedef typename UnaryReturn<V1,TAG>::Type_t Return_t;                      \
  return Return_t( Expr_t(v1) );                                              \
}

//-----------------------------------------------------------------------------

POOMA_TINYMATRIX_UNARY_OPERATOR(acos, FnArcCos)
POOMA_TINYMATRIX_UNARY_OPERATOR(asin, FnArcSin)
POOMA_TINYMATRIX_UNARY_OPERATOR(atan, FnArcTan)
POOMA_TINYMATRIX_UNARY_OPERATOR(ceil, FnCeil)
POOMA_TINYMATRIX_UNARY_OPERATOR(cos,  FnCos)
POOMA_TINYMATRIX_UNARY_OPERATOR(cosh, FnHypCos)
POOMA_TINYMATRIX_UNARY_OPERATOR(exp,  FnExp)
POOMA_TINYMATRIX_UNARY_OPERATOR(fabs, FnFabs)
POOMA_TINYMATRIX_UNARY_OPERATOR(floor, FnFloor)
POOMA_TINYMATRIX_UNARY_OPERATOR(log, FnLog)
POOMA_TINYMATRIX_UNARY_OPERATOR(log10, FnLog10)
POOMA_TINYMATRIX_UNARY_OPERATOR(sin, FnSin)
POOMA_TINYMATRIX_UNARY_OPERATOR(sinh, FnHypSin)
POOMA_TINYMATRIX_UNARY_OPERATOR(sqrt, FnSqrt)
POOMA_TINYMATRIX_UNARY_OPERATOR(tan, FnTan)
POOMA_TINYMATRIX_UNARY_OPERATOR(tanh, FnHypTan)
POOMA_TINYMATRIX_UNARY_OPERATOR(operator-, OpUnaryMinus)
POOMA_TINYMATRIX_UNARY_OPERATOR(operator+, OpUnaryPlus)
POOMA_TINYMATRIX_UNARY_OPERATOR(operator~, OpBitwiseNot)

//-----------------------------------------------------------------------------

//
// Elementwise binary operators on TinyMatrixs.
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
// operations with a TinyMatrix and a scalar.  Any type that isn't a TinyMatrix
// is considered a scalar.
//

#define POOMA_TINYMATRIX_BINARY_OPERATOR(FUNC,TAG)                            \
                                                                              \
template <int D1, int D2, class T1, class T2, class E>                        \
struct BinaryReturn< TinyMatrix<D1,D2,T1,E>, TinyMatrix<D1,D2,T2,E>, TAG >    \
{                                                                             \
  typedef TinyMatrix< D1, D2, typename BinaryReturn<T1,T2,TAG>::Type_t, E >   \
    Type_t;                                                                   \
};                                                                            \
                                                                              \
template <int D1, int D2, class T1, class T2, class E1, class E2>             \
inline                                                                        \
typename BinaryReturn< TinyMatrix<D1,D2,T1,E1>, TinyMatrix<D1,D2,T2,E2>,      \
                       TAG >::Type_t                                          \
FUNC( const TinyMatrix<D1,D2,T1,E1>& v1, const TinyMatrix<D1,D2,T2,E2>& v2 )  \
{                                                                             \
  typedef TinyMatrix<D1,D2,T1,E1> V1;                                         \
  typedef TinyMatrix<D1,D2,T2,E2> V2;                                         \
  typedef typename BinaryReturn<T1,T2,TAG>::Type_t T3;                        \
  typedef TinyMatrix< D1, D2, T3, BinaryTinyMatrixOp<V1,V2,TAG> > Expr_t;     \
  typedef typename BinaryReturn<V1,V2,TAG>::Type_t Return_t;                  \
  return Return_t( Expr_t(v1,v2) );                                           \
}                                                                             \
                                                                              \
template <int D1, int D2, class T1, class T2, class E>                        \
struct BinaryReturn< TinyMatrix<D1,D2,T1,E>, T2, TAG >                        \
{                                                                             \
  typedef TinyMatrix< D1, D2, typename BinaryReturn<T1,T2,TAG>::Type_t, E >   \
    Type_t;                                                                   \
};                                                                            \
                                                                              \
template <int D1, int D2, class T1, class T2, class E>                        \
struct BinaryReturn< T1, TinyMatrix<D1,D2,T2,E>, TAG >                        \
{                                                                             \
  typedef TinyMatrix< D1, D2, typename BinaryReturn<T1,T2,TAG>::Type_t, E >   \
    Type_t;                                                                   \
};                                                                            \
                                                                              \
template <int D1, int D2, class T1, class E, class T2>                        \
inline typename BinaryReturn< TinyMatrix<D1,D2,T1,E>, T2, TAG >::Type_t       \
FUNC( const TinyMatrix<D1,D2,T1,E>& v1, const T2& x )                         \
{                                                                             \
  typedef TinyMatrix<D1,D2,T1,E> V1;                                          \
  typedef typename BinaryReturn<V1,T2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef TinyMatrix< D1, D2, T3, BinaryTinyMatrixOp<V1,T2,TAG> > Expr_t;     \
  return Return_t( Expr_t(v1,x) );                                            \
}                                                                             \
                                                                              \
template <int D1, int D2, class T1, class T2, class E>                        \
inline typename BinaryReturn< T1, TinyMatrix<D1,D2,T2,E>, TAG >::Type_t       \
FUNC( const T1& x, const TinyMatrix<D1,D2,T2,E>& v2 )                         \
{                                                                             \
  typedef TinyMatrix<D1,D2,T2,E> V2;                                          \
  typedef typename BinaryReturn<T1,V2,TAG>::Type_t Return_t;                  \
  typedef typename Return_t::Element_t T3;                                    \
  typedef TinyMatrix< D1, D2, T3, BinaryTinyMatrixOp<T1,V2,TAG> > Expr_t;     \
  return Return_t( Expr_t(x,v2) );                                            \
}

//-----------------------------------------------------------------------------

POOMA_TINYMATRIX_BINARY_OPERATOR(operator+, OpAdd)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator-, OpSubtract)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator*, OpMultiply)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator/, OpDivide)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator%, OpMod)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator&, OpBitwiseAnd)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator|, OpBitwiseOr)
POOMA_TINYMATRIX_BINARY_OPERATOR(operator^, OpBitwiseXor)
POOMA_TINYMATRIX_BINARY_OPERATOR(ldexp,     FnLdexp)
POOMA_TINYMATRIX_BINARY_OPERATOR(pow,       FnPow)
POOMA_TINYMATRIX_BINARY_OPERATOR(fmod,      FnFmod)
POOMA_TINYMATRIX_BINARY_OPERATOR(atan2,     FnArcTan2)

//-----------------------------------------------------------------------------
//
// Dot product specialization of TinyMatrixEngine.
// This uses a different specialization because it isn't an elementwise op.
//
//-----------------------------------------------------------------------------

template<int D1, int D2, class T, class V1, class V2>
class TinyMatrixEngine<D1,D2,T,BinaryTinyMatrixOp<V1,V2,FnDot> >
{

public:
 
  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  typedef T Element_t;
  typedef FnDot Op;
  typedef BinaryTinyMatrixOp<V1,V2,Op> EngineTag_t;

  // Return types for accessor functions.
  typedef T ConstElementRef_t;
  typedef T ElementRef_t;

  // Record the type of the current class.
  typedef TinyMatrixEngine<D1,D2,T, BinaryTinyMatrixOp<V1,V2,Op> > This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Construct from two TinyMatrixs and on operator tag.
  TinyMatrixEngine(const V1& v1, const V2& v2, Op op)
    : v1_m(v1), v2_m(v2), op_m(op) {}

  // Construct from two TinyMatrixs and let the op tag contruct itself.
  TinyMatrixEngine(const V1& v1, const V2& v2)
    : v1_m(v1), v2_m(v2) {}

  // Copy ctor just copies the references.
  TinyMatrixEngine(const This_t& x)
    : v1_m(x.v1_m), v2_m(x.v2_m), op_m(x.op_m) {}

  // Let the engine destroy itself.
  ~TinyMatrixEngine() {}

#if !POOMA_NO_TEMPLATE_FRIENDS
  template<int DD1,int DD2, class TT, class EE, int I, int J>
    friend struct TinyMatrixEngineElem;

private:
#endif

  const V1& v1_m;
  const V2& v2_m;
  Op op_m;
};

//-----------------------------------------------------------------------------

//
// General form of TinyMatrixDotTinyMatrix.
//
// Finds one term of the dot product of a TinyMatrix dotted with a TinyMatrix.
// Inputs a TinyMatrix of type T1 and a TinyMatrix of type T2 and dots row
// I of T1 with column J of T2, starting with offset K in that row and
// vector length L.
//
// Operates by dividing that row in half, dotting the halves, and 
// returning the sum of those two.
//

template<class T1, class T2, int I, int J, int K, int L>
struct TinyMatrixDotTinyMatrix
{
  typedef typename TinyMatrixDotTinyMatrix<T1,T2,I,J,K,L/2>::Type_t Left_t;
  typedef typename TinyMatrixDotTinyMatrix<T1,T2,I,J,K+L/2,L-L/2>::Type_t Right_t;
  typedef typename BinaryReturn<Left_t,Right_t,OpAdd>::Type_t Type_t;
  static Type_t get(const T1& x, const T2& y)
    {
      return 
        TinyMatrixDotTinyMatrix<T1,T2,I,J,K,L/2>::get(x,y) + 
        TinyMatrixDotTinyMatrix<T1,T2,I,J,K+L/2,L-L/2>::get(x,y);
    }
};

//
// Recursion termination for TinyMatrixDotTinyMatrix.
// When you get down to a vector length of 1, all you need to do
// is multiply the elements and return that.
//

template<class T1, class T2, int I, int J, int K>
struct TinyMatrixDotTinyMatrix<T1,T2,I,J,K,1>
{
  typedef typename TinyMatrixElem<T1,I,K>::Element_t Left_t;
  typedef typename TinyMatrixElem<T2,K,J>::Element_t Right_t;
  typedef typename BinaryReturn<Left_t,Right_t,OpMultiply>::Type_t Type_t;
  static Type_t get(const T1& x, const T2& y)
    {
      return TinyMatrixElem<T1,I,K>::get(x) * TinyMatrixElem<T2,K,J>::get(y);
    }
};

//
// Specialization of TinyMatrixElem for getting one value from
// a TinyMatrixEngine for dotting two TinyMatrixs.
// Just calls TinyMatrixDotTinyMatrix.
//

template<int D1, int D2, class T, class T1, class T2, int I, int J>
struct TinyMatrixEngineElem<D1,D2,T,BinaryTinyMatrixOp<T1,T2,FnDot>, I, J >
{
  typedef BinaryTinyMatrixOp<T1,T2,FnDot> E;
  typedef TinyMatrixEngine<D1,D2,T,E> T0;
  typedef typename TinyMatrixDotTinyMatrix<T1,T2,I,J,0,T1::d1>::Type_t Element_t;
  typedef Element_t ElementRef_t;
  typedef Element_t ConstElementRef_t;
  static Element_t get(const T0& x) 
    { 
      return TinyMatrixDotTinyMatrix<T1,T2,I,J,0,T1::d2>::get(x.v1_m,x.v2_m);
    }
};

//
// Define the return type for TinyMatrix dot TinyMatrix.
// Since this has no knowledge of the input TinyMatrix types, 
// this returns a full TinyMatrix.  Specialized TinyMatrix types should
// define specialized versions of this.
//

template<int D1, int D2, int D3, class T1, class T2, class E1, class E2>
struct BinaryReturn< TinyMatrix<D1,D2,T1,E1> , TinyMatrix<D2,D3,T2,E2> , FnDot >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
  typedef TinyMatrix<D1,D3,T0,Full> Type_t;
};

//
// Dot product of two TinyMatrixs.
// Uses all of the above to build a new TinyMatrix.
//

template<int D1, int D2, int D3, class T1, class T2, class E1, class E2>
inline
typename BinaryReturn< TinyMatrix<D1,D2,T1,E1>,TinyMatrix<D2,D3,T2,E2> , FnDot >::Type_t
dot( const TinyMatrix<D1,D2,T1,E1>& v1, const TinyMatrix<D2,D3,T2,E2>& v2 )
{
  typedef TinyMatrix<D1,D2,T1,E1> V1;
  typedef TinyMatrix<D2,D3,T2,E2> V2;
  typedef typename BinaryReturn<V1,V2,FnDot>::Type_t Return_t;
  typedef typename Return_t::Element_t T3;
  typedef TinyMatrix<D1,D3,T3,BinaryTinyMatrixOp<V1,V2,FnDot> > Expr_t;
  return Return_t( Expr_t(v1,v2) );
}

//-----------------------------------------------------------------------------
//
// Equality operator between two TinyMatrixs.
//
// Like the functions above, two things need to be done:
//
// 1. Define the return type by specializing BinaryReturn.
// 
// 2. Define the equality operation.
//
// TinyMatrixs are sufficiently complicated that we just bag the metraprograms
// and write a loop.
//
//-----------------------------------------------------------------------------

//
// Define the return type for checking equality of two TinyMatrixs.
// Just a bool.
//

template<int D1, int D2, class T1, class T2, class E1, class E2>
struct BinaryReturn< TinyMatrix<D1,D2,T1,E1>,TinyMatrix<D1,D2,T2,E2> , OpEQ >
{
  typedef bool Type_t;
};

template<int D1, int D2, class T1, class T2, class E1, class E2>
struct BinaryReturn< TinyMatrix<D1,D2,T1,E1>,TinyMatrix<D1,D2,T2,E2> , OpNE >
{
  typedef bool Type_t;
};

//
// The actual equality operators.
// Just do the stinkin' loops.
//

template<int D1, int D2, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< TinyMatrix<D1,D2,T1,E1>,TinyMatrix<D1,D2,T2,E2> , OpEQ >::Type_t
operator==(const TinyMatrix<D1,D2,T1,E1>& t1, const TinyMatrix<D1,D2,T2,E2>& t2)
{
  for (int i = 0; i < D1; i++)
    for (int j = 0; j < D2; j++)
      if (t1(i,j) != t2(i,j)) return false;
   
  return true;
}

template<int D1, int D2, class T1, class T2, class E1, class E2>
inline typename 
BinaryReturn< TinyMatrix<D1,D2,T1,E1>,TinyMatrix<D1,D2,T2,E2> , OpNE >::Type_t
operator!=(const TinyMatrix<D1,D2,T1,E1>& t1, const TinyMatrix<D1,D2,T2,E2>& t2)
{
  return !(t1 == t2);
}

//-----------------------------------------------------------------------------
//
// General assignment (including things like +=) for TinyMatrixs.
//
// There are two things that need to be done:
//
// 1.  Define the return type by specializing BinaryReturn.  
//     This is actually done automatically in PETE.
//
// 2.  Define the operator to call TinyMatrixAssign.
//
// These are all pretty similar, so we use a macro to define it.
// Types not recognized as TinyMatrixs are considered to be scalars.
//
//-----------------------------------------------------------------------------

#define POOMA_TINYMATRIX_ACCUM_OPERATOR(FUNC,TAG)                             \
                                                                              \
template <int D1, int D2, class T1, class T2, class E1, class E2>             \
inline const TinyMatrix<D1,D2,T1,E1>&                                         \
FUNC( TinyMatrix<D1,D2,T1,E1>& v1, const TinyMatrix<D1,D2,T2,E2>& v2 )        \
{                                                                             \
  typedef TinyMatrix<D1,D2,T1,E1> Left_t;                                     \
  typedef TinyMatrix<D1,D2,T2,E2> Right_t;                                    \
  TinyMatrixAssign<Left_t,Right_t,TAG,0,D1,0,D2>::apply(v1,v2,TAG());         \
  return v1;                                                                  \
}                                                                             \
                                                                              \
template <int D1, int D2, class T1, class T2, class E1>                       \
inline const TinyMatrix<D1,D2,T1,E1>&                                         \
FUNC( TinyMatrix<D1,D2,T1,E1>& v1, const T2& v2 )                             \
{                                                                             \
  TinyMatrixAssign<TinyMatrix<D1,D2,T1,E1>,T2,TAG,0,D1,0,D2>::                \
    apply(v1,v2,TAG());                                                       \
  return v1;                                                                  \
}

POOMA_TINYMATRIX_ACCUM_OPERATOR(operator+=, OpAddAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator-=, OpSubtractAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator*=, OpMultiplyAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator/=, OpDivideAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator%=, OpModAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator|=, OpBitwiseOrAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator&=, OpBitwiseAndAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator^=, OpBitwiseXorAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator<<=, OpLeftShiftAssign)
POOMA_TINYMATRIX_ACCUM_OPERATOR(operator>>=, OpRightShiftAssign)

#endif // POOMA_TINY_TINYMATRIX_OPERATORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TinyMatrixOperators.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:11 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
