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

//-----------------------------------------------------------------------------
// Functions:
//   reverseBytes<T>(T &)
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 *   Defines a general reverseBytes template and a specialization for
 *   complex. The general template is implemented by the struct
 *   ByteReverser, which is not published. 
 */

#ifndef POOMA_IO_BYTEORDER_H
#define POOMA_IO_BYTEORDER_H

#include "Utilities/ElementProperties.h"

// For reverseBytes specialization

#include <complex>

//-----------------------------------------------------------------------------
// Utility functions to fix byte ordering
//-----------------------------------------------------------------------------

namespace {

  /**
   * The ByteReverser template defines the default behavior for
   * reversing bytes. It is templated on the type and on a bool that
   * indicates whether it is a basic "C" type or not. 
   *   - If it is a basic C type, then the appropriate thing to do is
   *     to reverse the order of the sizeof(T) bytes in which the T
   *     object is stored.
   *   - If it is not a basic C type, then we try to invoke the
   *     "reverseBytes" member funcion. This is appropriate for the
   *     concrete types defined by Pooma (Vector, Tensor,
   *     TinyMatrix). It will not work for complex, which we catch by
   *     directly overloading the global reverseBytes function below.
   *
   * If a user's code fails to compile because he's using a T that
   * does not define T::reverseBytes(), he needs to write a
   * specialization of the global reverseBytes(T&) template for his
   * type. 
   */

  template <class T, bool basicType>
  struct ByteReverser;

  template <class T>
  struct ByteReverser<T,true>
  {
    inline static void reverseBytes(T &t)
    {
      // Might be faster to do this in place, but for now we'll make
      // a temporary. 

      T x;
      char *xb = reinterpret_cast<char*>(&x) + sizeof(T);
      char *b = reinterpret_cast<char*>(&t);
      for (unsigned int i = 0; i < sizeof(T); ++i) 
        {
          *--xb = *b++;
        }

      t = x;
    }
  };

  template <class T>
  struct ByteReverser<T,false>
  {
    inline static void reverseBytes(T &t)
    {
      t.reverseBytes();
    }
  };
}

/// The exported interface to the above template is the global
/// reverseBytes template. The general template delegates directly to
/// the ByteReverser class above. 

template <class T>
inline void reverseBytes(T &t)
{
  typedef ByteReverser<T,ElementProperties<T>::basicType> ByteReverser_t;
  
  ByteReverser_t::reverseBytes(t);
}

/// We overload directly for class complex

template <class FloatType>
inline void reverseBytes(std::complex<FloatType> &t)
{
  // Unfortunately, the C++ std::complex does not give us direct
  // access to the underlying real and imaginary parts - the accessors
  // return by value and not by reference. Obviously we need to access
  // these. In order to do this, we assume that the storage format is:

  struct mycomplex
  { 
    FloatType re;
    FloatType im;
  };

  // and we use a cast to break the encapsulation and get at the
  // underlying types. Clearly this is not guaranteed to work and
  // should be tested for the users specific C++ library.

  FloatType &real = reinterpret_cast<mycomplex*>(&t)->re;
  FloatType &imag = reinterpret_cast<mycomplex*>(&t)->im;

  reverseBytes(real);
  reverseBytes(imag);
}

#endif // POOMA_IO_BYTEORDER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ByteOrder.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
