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

#ifndef POOMA_UTILITIES_METAPROGS_H
#define POOMA_UTILITIES_METAPROGS_H

//-----------------------------------------------------------------------------
// Classes:
//   LoopUtils
// Meta-Programs:
//   LoopUtils<N>::copy
//   LoopUtils<N>::dot
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/** @file
 * @ingroup Utilities
 * @brief
 * LoopUtils: Template meta-programs for carrying out certain operations
 *            on arrays at compile time.
 */


/**
 * Template meta-program loop class with the following functions:
 *  - void LoopUtils<N>::copy(T *dst, const T *src)
 *    copies N contiguous T's from src to dst. Uses assignment
 *    operator to perform the copy.
 *  - T LoopUtils<N>::dot(const T *a, const T *b)
 *    calculates dot product of a[N] and b[N] at compile time.
 */

template <int N>
struct LoopUtils
{
  template <class T>
  static inline void copy(T* dest, const T* src)
  {
    dest[N-1] = src[N-1];
    LoopUtils<N-1>::copy(dest, src);
  }

  template <class T>
  static inline T dot(const T * a, const T * b)
  {
    return *a * *b + LoopUtils<N-1>::dot(a+1,b+1);
  }
};

template <>
struct LoopUtils<1>
{
  template <class T>
  static inline void copy(T* dest, const T* src)
    {
      dest[0] = src[0];
    }

  template <class T>
  static inline T dot(const T * a, const T * b)
    {
      return *a * *b;
    }
};


// } // namespace POOMA
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_METAPROGS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MetaProgs.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
