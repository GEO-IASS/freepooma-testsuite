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

/** @file
 * @ingroup Pooma
 * @brief
 * Pooma Extras for the Portable Expression Template Engine
 */

#ifndef POOMA_POOMA_PETE_EXTRAS_H
#define POOMA_POOMA_PETE_EXTRAS_H

template<int D,class T,class E> class Vector;
template<int DR, int DC, class T, class E> class TinyMatrix;
template<int D, class T, class E> class Tensor;

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "Pooma/PoomaOperatorTags.h"

#include <algorithm>

//-----------------------------------------------------------------------------
// Complex support
//-----------------------------------------------------------------------------

// Attempt to set up macros that tell how to refer to complex numbers. We
// handle templated complex classes in or out of the std::namespace.

#if !POOMA_NO_TEMPLATED_COMPLEX

#if !POOMA_NO_STD_COMPLEX

namespace std {

template<class T> class complex;

}

using std::complex;

#define POOMA_COMPLEX complex<T>
#define POOMA_COMPLEX_CLASS complex

#else

template<class T> class complex;

#define POOMA_COMPLEX complex<T>
#define POOMA_COMPLEX_CLASS complex

#endif

#else

#error "You need to edit Pooma/PETEExtras.h to configure complex numbers."

#endif

//-----------------------------------------------------------------------------
// Complex unary operators
//-----------------------------------------------------------------------------

// The real, imag, abs, arg, and, norm functions need to be defined 
// specially since they don't return complex numbers.

template<class T>
struct UnaryReturn< POOMA_COMPLEX, FnConj >
{
  typedef POOMA_COMPLEX Type_t;
};

template<class T>
struct UnaryReturn<POOMA_COMPLEX, FnReal>
{
  typedef T Type_t;
};

template<class T>
struct UnaryReturn<POOMA_COMPLEX, FnImag>
{
  typedef T Type_t;
};

template<class T>
struct UnaryReturn<POOMA_COMPLEX, FnArg>
{
  typedef T Type_t;
};

template<class T>
struct UnaryReturn<POOMA_COMPLEX, FnNorm>
{
  typedef T Type_t;
};

template<class T>
struct UnaryReturn<T, FnAbs>
{
  typedef T Type_t;
};

template<class T>
struct UnaryReturn<POOMA_COMPLEX, FnAbs>
{
  typedef T Type_t;
};


//-----------------------------------------------------------------------------
// Complex binary operators
//-----------------------------------------------------------------------------

// Define promotions for complex.

// Astoundingly, the binary operators in the standard do not support mixed
// mode arithmetic. 

template<class T>
struct Promote<POOMA_COMPLEX, POOMA_COMPLEX >
{
  typedef POOMA_COMPLEX Type_t;
};

template<class T>
struct Promote<POOMA_COMPLEX, T>
{
  typedef POOMA_COMPLEX Type_t;
};

template<class T>
struct Promote<T, POOMA_COMPLEX>
{
  typedef POOMA_COMPLEX Type_t;
};

// This version is so pow(complex, int) will work correctly.

template<class T>
struct BinaryReturn<POOMA_COMPLEX, int, FnPow>
{
  typedef POOMA_COMPLEX Type_t;
};


//-----------------------------------------------------------------------------
// Special POOMA functions
//-----------------------------------------------------------------------------

template<class T>
struct UnaryReturn<T, FnPow2>
{
  typedef typename BinaryReturn<T, T, OpMultiply>::Type_t Type_t;
};

template<class T>
struct UnaryReturn<T, FnPow3>
{
  typedef typename BinaryReturn<T, T, OpMultiply>::Type_t Type_t;
};

template<class T>
struct UnaryReturn<T, FnPow4>
{
  typedef typename BinaryReturn<T, T, OpMultiply>::Type_t Type_t;
};


#endif // POOMA_POOMA_PETE_EXTRAS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PETEExtras.h,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:17:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
