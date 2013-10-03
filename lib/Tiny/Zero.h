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

#ifndef POOMA_TINY_ZERO_H
#define POOMA_TINY_ZERO_H

/** @file
 * @ingroup Tiny
 * @brief
 * A numeric class for a number that is always zero.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * The point of this class is to be a number of type T which is known at
 * compile time to be zero.  This is reflected in the operations like
 * addition and subtraction that use it.  If you return a Zero object,
 * the compiler should be able to make better optimizations than if you
 * just return zero.
 *
 * Zero is templated on type T, to represent a zero object of type T.
 * Type T is required to be constructable from an integer zero.
 */

template<class T>
class Zero
{
  // If you need to convert to an object of type T,
  // just build one from zero.  This will be used in 
  // the cases where the operators below don't match.

  operator T() { return T(0); }


  // Give it empty ctors and assignment operators
  // to try and keep purify happy.

  Zero() {}
  Zero(const Zero<T>&) {}
  Zero<T>& operator=(const Zero<T>&) { return *this; }
};

//-----------------------------------------------------------------------------
//
// Forward declarations for classes used by TypeComputations for figuring out
// return types in expression templates.
//
//-----------------------------------------------------------------------------

template<class T, class Op> struct UnaryReturn;
template<class T1, class T2, class Op> struct BinaryReturn;

//-----------------------------------------------------------------------------
//
// Operators using Zero.
//
//-----------------------------------------------------------------------------

//
// Binary multiply by a Zero<T> and T returns Zero<T>
//

template<class T> 
  inline Zero<T>  operator*(Zero<T>, const T&)   { return Zero<T>(); }

template<class T> 
  inline Zero<T>  operator*(const T&, Zero<T>)   { return Zero<T>(); }

template<class T> 
  inline Zero<T>  operator*(Zero<T>, Zero<T>)    { return Zero<T>(); }

template<class T> 
  inline Zero<T>  operator/(Zero<T>, const T&)   { return Zero<T>(); }

//
// Trait classes so that expression templates will deal correctly
// with Zeros in multiplicative operations.
//

template<class T>
struct BinaryReturn< Zero<T> , T , OpMultiply >
{
  typedef Zero<T> Type_t;
};

template<class T>
struct BinaryReturn< T, Zero<T> , OpMultiply >
{
  typedef Zero<T> Type_t;
};

template<class T>
struct BinaryReturn< Zero<T>, Zero<T> , OpMultiply >
{
  typedef Zero<T> Type_t;
};

template<class T>
struct BinaryReturn< Zero<T> , T , OpDivide >
{
  typedef Zero<T> Type_t;
};

//
// Adding a Zero<T> to a T returns the T.
//

template<class T> 
  inline const T& operator+(Zero<T>, const T& x) { return  x; }

template<class T> 
  inline const T& operator+(const T& x, Zero<T>) { return  x; }

template<class T> 
  inline Zero<T>  operator+(Zero<T>, Zero<T>)    { return  Zero<T>(); }

template<class T> 
  inline T        operator-(Zero<T>, const T& x) { return -x; }

template<class T> 
  inline const T& operator-(const T& x, Zero<T>) { return  x; }

template<class T> 
  inline Zero<T>  operator-(Zero<T>, Zero<T>)    { return  Zero<T>(); }

//
// Trait class specializations for additive operations.
//

template<class T>
struct BinaryReturn< Zero<T> , T , OpAdd >
{
  typedef T Type_t;
};

template<class T>
struct BinaryReturn< T, Zero<T> , OpAdd >
{
  typedef T Type_t;
};

template<class T>
struct BinaryReturn< Zero<T>, Zero<T> , OpAdd >
{
  typedef Zero<T> Type_t;
};

template<class T>
struct BinaryReturn< Zero<T> , T , OpSubtract >
{
  typedef T Type_t;
};

template<class T>
struct BinaryReturn< T, Zero<T> , OpSubtract >
{
  typedef T Type_t;
};

template<class T>
struct BinaryReturn< Zero<T>, Zero<T> , OpSubtract >
{
  typedef Zero<T> Type_t;
};

//
// Unary Minus of a zero returns a zero.
//

template<class T>
  inline Zero<T>  operator-(Zero<T>) { return Zero<T>(); }

template<class T>
struct UnaryReturn< Zero<T> , OpUnaryMinus >
{
  typedef Zero<T> Type_t;
};


//
// Unary Plus of a zero returns a zero.
//

template<class T>
  inline Zero<T>  operator+(Zero<T>) { return Zero<T>(); }

template<class T>
struct UnaryReturn< Zero<T> , OpUnaryPlus >
{
  typedef Zero<T> Type_t;
};

//
// Trait class definitios for the cases where the Zero<T> gets
// converted to a T.
//

template<class T, class Op>
struct UnaryReturn< Zero<T> , Op >
{
  typedef typename UnaryReturn<T,Op>::Type_t Type_t;
};

template<class T1, class T2, class Op>
struct BinaryReturn< Zero<T1> , T2, Op >
{
  typedef typename BinaryReturn<T1,T2,Op>::Type_t Type_t;
};

template<class T1, class T2, class Op>
struct BinaryReturn< T1 , Zero<T2>, Op >
{
  typedef typename BinaryReturn<T1,T2,Op>::Type_t Type_t;
};

template<class T1, class T2, class Op>
struct BinaryReturn< Zero<T1> , Zero<T2>, Op >
{
  typedef typename BinaryReturn<T1,T2,Op>::Type_t Type_t;
};



#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Zero.h,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:17:11 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
