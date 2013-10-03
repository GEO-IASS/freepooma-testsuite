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

#ifndef POOMA_UTILITIES_UNINITIALIZED_VECTOR_H
#define POOMA_UTILITIES_UNINITIALIZED_VECTOR_H

//-----------------------------------------------------------------------------
// Class:
// UninitializedVector<T,Dim,Elem>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Utilities
 * @brief
 * UninitializedVector stores an array of objects of type 'T' of length 'Dim'
 * in a way that avoids running the default constructors on the objects
 * unless the 'initialize' method is called.
 *
 * It can help avoid unwanted
 * for-loops over array elements which are normally generated in order
 * to run their default constructors, even when the length of the array is
 * known at compile time.  It can be used to unroll the loop over the array
 * elements to run their constructors, in order to initialize a whole array
 * of objects with a non-default constructor in one call.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <new>


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class T, int I> struct InitializeUninitializedVector;


/**
 * UninitializedVector<T,Dim> stores internally an array of object of type
 * 'Elem', of length sizeof(T) * Dim / sizeof(Elem), where T and Elem
 * are some arbitrary types.  It appears
 * to other objects like an array of objects of type T, with a fixed
 * number of elements Dim.  It provides an operator[] for accessing the
 * Nth element.  The type 'Elem' should (but does not have to be) some kind
 * of elementary type, like int or char.
 *
 * UninitializedVector is used primarily by the Domain objects to store
 * N 1-D domain objects, so that T is something like Interval<1>.
 * In that case, we often wish to avoid running the default constructors
 * on the domain array elements, since they will be overwritten with new
 * values when the array is filled in with correct values.  UninitializedVector
 * uses a simple array to provide the storage for the N objects,
 * so that the default constructor for UninitializedVector
 * (the only one publicly available)
 * does not do anything.  If you want to run the default
 * constructors for the domain objects, the 'initialize' method is available,
 * which can take 0, 1, 2, or 3 arguments.  initialize uses simple
 * template meta-programs 'InitializeUninitializedVectorN' to perform the
 * initialization of the elements, by using the placement new operator.
 *
 * When used, UninitializedVector should allow you to have an array of
 * fixed length which is not initialized in any way until you call initialize.
 * And when initialize is called, the loop to call the constructors for
 * the objects in the array should be unrolled and inlined by the template
 * metaprogram used to call the constructors.
 */

template<class T, int Dim, class Elem>
class UninitializedVector
{
public:
  // the constructor does nothing; to initialize the domain, call
  // the initialize method
  UninitializedVector() {
    CTAssert(Dim > 0);
  }

  // destructor: does nothing ... makes assumption domains do not need
  // to have their destructor run.
  ~UninitializedVector() { }

  // initialize the array via placement new and a template metaprogram
  void initialize() {
    InitializeUninitializedVector<T,Dim-1>::initialize(buffer);
  }

  // initialize the array via placement new and a template metaprogram,
  // specifying an argument(s)
  template<class T1>
  void initialize(T1 &a) {
    InitializeUninitializedVector<T,Dim-1>::initialize(buffer, a);
  }

  // initialize the array via placement new and a template metaprogram,
  // specifying an argument(s)
  template<class T1, class T2>
  void initialize(T1 &a, T2 &b) {
    InitializeUninitializedVector<T,Dim-1>::initialize(buffer, a, b);
  }

  // initialize the array via placement new and a template metaprogram,
  // specifying an argument(s)
  template<class T1, class T2, class T3>
  void initialize(T1 &a, T2 &b, T3 &c) {
    InitializeUninitializedVector<T,Dim-1>::initialize(buffer, a, b, c);
  }

  // get the Nth domain
  T &operator[](int n) {
    PAssert(n >= 0 && n < Dim);
    return *(reinterpret_cast<T *>(buffer) + n);
  }

  const T &operator[](int n) const {
    PAssert(n >= 0 && n < Dim);
    return *(reinterpret_cast<const T *>(buffer) + n);
  }

private:
  // a buffer of Elem's used to store bytes for the array of domains
  Elem buffer[(Dim * sizeof(T))/sizeof(Elem)];

  // make the copy constructor and operator= private; they should not be used
  UninitializedVector(const UninitializedVector<T,Dim,Elem> &);
  UninitializedVector<T,Dim,Elem> &
    operator=(const UninitializedVector<T,Dim,Elem> &);
};


/**
 * InitializeUninitializedVector is a simple functor with an N-D and 1-D
 * version used for template meta-programs to initialize the first 'Dim' elems
 * of the data in the given UninitializedVector.  It defines 'initialize
 * which takes from 0 ... 3 arguments (along with the object storage)
 * which are passed on to the objects when placement new is used.  A general
 * N-dimensional version calls the N-1-dimensional version; a specialization
 * for the 0th index terminates the template recursion.
 *
 * The first template parameter is the type of data to initialize; the
 * second is index (0 ... Dim-1) which we are currently initializing.
 */

template<class T, int I>
struct InitializeUninitializedVector
{
  template<class Elem>
  static void initialize(Elem *buffer) {
    CTAssert(I >= 0);
    InitializeUninitializedVector<T,I-1>::initialize(buffer);
    new (static_cast<void *>(reinterpret_cast<char *>(buffer) +
			     (I * sizeof(T)))) T();
  }
  template<class Elem, class T1>
  static void initialize(Elem *buffer, T1 &a) {
    CTAssert(I >= 0);
    InitializeUninitializedVector<T,I-1>::initialize(buffer, a);
    new (static_cast<void *>(reinterpret_cast<char *>(buffer) +
			     (I * sizeof(T)))) T(a);
  }
  template<class Elem, class T1, class T2>
  static void initialize(Elem *buffer, T1 &a, T2 &b) {
    CTAssert(I >= 0);
    InitializeUninitializedVector<T,I-1>::initialize(buffer, a, b);
    new (static_cast<void *>(reinterpret_cast<char *>(buffer) +
			     (I * sizeof(T)))) T(a,b);
  }
  template<class Elem, class T1, class T2, class T3>
  static void initialize(Elem *buffer, T1 &a, T2 &b, T3 &c) {
    CTAssert(I >= 0);
    InitializeUninitializedVector<T,I-1>::initialize(buffer, a, b, c);
    new (static_cast<void *>(reinterpret_cast<char *>(buffer) +
			     (I * sizeof(T)))) T(a,b,c);
  }
};

/// zeroth-index specialization of InitializeUninitializedVector
template<class T>
struct InitializeUninitializedVector<T, 0>
{
  template<class Elem>
  static void initialize(Elem *buffer) {
    new (static_cast<void *>(buffer)) T();
  }
  template<class Elem, class T1>
  static void initialize(Elem *buffer, T1 &a) {
    new (static_cast<void *>(buffer)) T(a);
  }
  template<class Elem, class T1, class T2>
  static void initialize(Elem *buffer, T1 &a, T2 &b) {
    new (static_cast<void *>(buffer)) T(a,b);
  }
  template<class Elem, class T1, class T2, class T3>
  static void initialize(Elem *buffer, T1 &a, T2 &b, T3 &c) {
    new (static_cast<void *>(buffer)) T(a,b,c);
  }
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_UNINITIALIZED_VECTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UninitializedVector.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
