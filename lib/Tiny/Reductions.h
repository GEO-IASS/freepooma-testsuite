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
 * @ingroup Tiny
 * @brief
 * General global reduction functions.
 *
 * Functions: 
 *   - globalReduction: templated base function called by specific reductions.
 *   - sum: sum all the elements in a tiny object.
 *   - prod: multiply all of the elements in a tiny object.
 *   - max: find the maximum value in a tiny object.
 *   - min: find the minimum value in a tiny object.
 *   - all: returns true if all of the array's elements are != 0.
 *   - any: returns true if any of the array's elements are != 0.
 *   - bitOr: does a bitwise or of all of the elements.
 *   - bitAnd: does a bitwise and of all of the elements.
 *
 * Note: these functions work for reductions that apply pairwise arithmetic 
 * operations on the elements (e.g., sum, prod). This does not work for
 * reductions like all() and any().
 */

#ifndef POOMA_TINY_REDUCTIONS_H
#define POOMA_TINY_REDUCTIONS_H

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "Tiny/Tensor.h"
#include "Tiny/TinyMatrix.h"
#include "Tiny/Vector.h"


template<int Dim, class T, class EngineTag, class Op>
inline T globalReduction(const Vector<Dim, T, EngineTag> &v, 
  const Op &op)
{
  // Prime with the first element.
  
  T val = v(0);
  
  // Loop over the remaining elements.
  
  for (int i = 1; i < Dim; i++)
    op(val, v(i)); 

  return val;
}

template<int Dim, class T, class EngineTag, class Op>
inline T globalReduction(const Tensor<Dim, T, EngineTag> &t, 
  const Op &op)
{
  // Prime with the first element.
  
  T val = t(0, 0);

  // Loop over the remaining elements. 
  
  for (int k = 1; k < Dim; k++)
    op(val,  t(k, 0));
    
  for (int j = 1; j < Dim; j++)
    for (int i = 0; i < Dim; i++)
      op(val, t(i, j));

  return val;
}

template<int Dim, class T, class Op>
inline T globalReduction(const Tensor<Dim, T, Full> &t, 
  const Op &op)
{
  // Prime with the first element.
  
  T val = t(0, 0);

  // Loop over the remaining elements. 
  
  for (int i = 1; i < TensorStorageSize<Dim, Full>::Size; i++)
    op(val,  t(i));

  return val;
}

template<int Dim, class T, class Op>
inline T globalReduction(const Tensor<Dim, T, Antisymmetric> &t, 
  const Op &op)
{
  // The diagonal consists of all zeros. Prime with the first element
  // there.
  
  T val = t(0,0);

  // This loop will be over the off-diagonal elements. The subtraction of
  // 1 / Dim is required because the Dim == 1 case reports a size of 1
  // rather than zero (for dimensioning purposes).
  
  for (int i = 0; i < TensorStorageSize<Dim, Antisymmetric>::Size - 1 / Dim; i++)
    {
      op(val,  t(i));
      op(val, -t(i));
    }

  return val;
}

template<int Dim, class T, class Op>
inline T globalReduction(const Tensor<Dim, T, Diagonal> &t, 
  const Op &op)
{
  // Prime with the first element on the diagonal.
  
  T val = t(0);
  
  // Do the rest of the diagonal.
  
  for (int i = 1; i < TensorStorageSize<Dim, Diagonal>::Size; i++)
    op(val, t(i));

  // Handle the off-diagonal elements in one whack. This assumes
  // that the operator is such that repeated applications with
  // zero does not change the answer.
  
  if (Dim > 1)
    op(val, T(0));
  
  return val;
}

template<int Dim1, int Dim2, class T, class EngineTag, class Op>
inline T globalReduction(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m, 
  const Op &op)
{
  // Prime with the first element.
  
  T val = m(0, 0);

  // Loop over the remaining elements. 
  
  for (int k = 1; k < Dim1; k++)
    op(val, m(k, 0));
    
  for (int j = 1; j < Dim2; j++)
    for (int i = 0; i < Dim1; i++)
      op(val, m(i, j));

  return val;
}

template<int Dim1, int Dim2, class T, class Op>
inline T globalReduction(const TinyMatrix<Dim1, Dim2, T, Full> &m, 
  const Op &op)
{
  // Prime with the first element on the diagonal.
  
  T val = m(0, 0);

  // Loop over all of the other elements.
  
  for (int i = 1; i < Dim1 * Dim2; i++)
    op(val, m(i));

  return val;
}


//-----------------------------------------------------------------------------
// Specific global reduction functions.
//-----------------------------------------------------------------------------

// Vectors.

/// Sum up the elements of a Vector.

template<int Dim, class T, class EngineTag>
T sum(const Vector<Dim, T, EngineTag> &v)
{
  return globalReduction(v, OpAddAssign());
}

/// Compute the product of the elements of a Vector.

template<int Dim, class T, class EngineTag>
T prod(const Vector<Dim, T, EngineTag> &v)
{
  return globalReduction(v, OpMultiplyAssign());
}

/// Find the smallest element of a Vector.

template<int Dim, class T, class EngineTag>
T min(const Vector<Dim, T, EngineTag> &v)
{
  return globalReduction(v, FnMinAssign());
}

/// Find the largest element of a Vector.

template<int Dim, class T, class EngineTag>
T max(const Vector<Dim, T, EngineTag> &v)
{
  return globalReduction(v, FnMaxAssign());
}

/// Report if all of the elements of a Vector are true.

template<int Dim, class T, class EngineTag>
bool all(const Vector<Dim, T, EngineTag> &v)
{
  for (int i = 0; i < Dim; i++)
    if (!v(i))
      return false;
  return true;
}

/// Report if some of the elments of a Vector are true.

template<int Dim, class T, class EngineTag>
bool any(const Vector<Dim, T, EngineTag> &v)
{
  for (int i = 0; i < Dim; i++)
    if (v(i))
      return true;
  return false;
}

/// Bitwise-or all of the elements together.

template<int Dim, class T, class EngineTag>
T bitOr(const Vector<Dim, T, EngineTag> &v)
{
  return globalReduction(v, OpBitwiseOrAssign());
}

/// Bitwise-and all of the elements together.

template<int Dim, class T, class EngineTag>
T bitAnd(const Vector<Dim, T, EngineTag> &v)
{
  return globalReduction(v, OpBitwiseAndAssign());
}

// Tensors.

/// Sum up the elements of a Tensor.

template<int Dim, class T, class EngineTag>
T sum(const Tensor<Dim, T, EngineTag> &t)
{
  return globalReduction(t, OpAddAssign());
}

/// Trivial case (elements must sum to zero).

template<int Dim, class T>
T sum(const Tensor<Dim, T, Antisymmetric> &t)
{
  return T(0);
}

/// Compute the product of the elements of a Tensor.

template<int Dim, class T, class EngineTag>
T prod(const Tensor<Dim, T, EngineTag> &t)
{
  return globalReduction(t, OpMultiplyAssign());
}

/// Trivial case (diagonal is zero).

template<int Dim, class T>
T prod(const Tensor<Dim, T, Antisymmetric> &t)
{
  return T(0);
}

/// Find the smallest element of a Tensor.

template<int Dim, class T, class EngineTag>
T min(const Tensor<Dim, T, EngineTag> &t)
{
  return globalReduction(t, FnMinAssign());
}

/// Find the largest element of a Tensor.

template<int Dim, class T, class EngineTag>
T max(const Tensor<Dim, T, EngineTag> &t)
{
  return globalReduction(t, FnMaxAssign());
}

/// Report if all of the elements of a Tensor are true.

template<int Dim, class T, class EngineTag>
bool all(const Tensor<Dim, T, EngineTag> &t)
{
  for (int j = 0; j < Dim; j++)
    for (int i = 0; i < Dim; i++)
      if (!t(i, j))
        return false;
  return true;
}

/// Trivial case (diagonal is zero).

template<int Dim, class T>
bool all(const Tensor<Dim, T, Antisymmetric> &t)
{
  return false;
}

/// Report if some of the elments of a Tensor are true.

template<int Dim, class T, class EngineTag>
bool any(const Tensor<Dim, T, EngineTag> &t)
{
  for (int j = 0; j < Dim; j++)
    for (int i = 0; i < Dim; i++)
      if (t(i, j))
        return true;
  return false;
}

/// Bitwise-or all of the elements together.

template<int Dim, class T, class EngineTag>
T bitOr(const Tensor<Dim, T, EngineTag> &t)
{
  return globalReduction(t, OpBitwiseOrAssign());
}

/// Bitwise-and all of the elements together.

template<int Dim, class T, class EngineTag>
T bitAnd(const Tensor<Dim, T, EngineTag> &t)
{
  return globalReduction(t, OpBitwiseAndAssign());
}

// TinyMatrices.

/// Sum up the elements of a TinyMatrix.

template<int Dim1, int Dim2, class T, class EngineTag>
T sum(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  return globalReduction(m, OpAddAssign());
}

/// Compute the product of the elements of a TinyMatrix.

template<int Dim1, int Dim2, class T, class EngineTag>
T prod(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  return globalReduction(m, OpMultiplyAssign());
}

/// Find the smallest element of a TinyMatrix.

template<int Dim1, int Dim2, class T, class EngineTag>
T min(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  return globalReduction(m, FnMinAssign());
}

/// Find the largest element of a TinyMatrix.

template<int Dim1, int Dim2, class T, class EngineTag>
T max(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  return globalReduction(m, FnMaxAssign());
}

/// Report if all of the elements of a TinyMatrix are true.

template<int Dim1, int Dim2, class T, class EngineTag>
bool all(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  for (int j = 0; j < Dim2; j++)
    for (int i = 0; i < Dim1; i++)
      if (!m(i, j))
        return false;
  return true;
}

/// Report if some of the elments of a TinyMatrix are true.

template<int Dim1, int Dim2, class T, class EngineTag>
bool any(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  for (int j = 0; j < Dim2; j++)
    for (int i = 0; i < Dim1; i++)
      if (m(i, j))
        return true;
  return false;
}

/// Bitwise-or all of the elements together.

template<int Dim1, int Dim2, class T, class EngineTag>
T bitOr(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  return globalReduction(m, OpBitwiseOrAssign());
}

/// Bitwise-and all of the elements together.

template<int Dim1, int Dim2, class T, class EngineTag>
T bitAnd(const TinyMatrix<Dim1, Dim2, T, EngineTag> &m)
{
  return globalReduction(m, OpBitwiseAndAssign());
}

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Reductions.h,v $   $Author: richard $
// $Retision: 1.10 $   $Date: 2004/11/01 18:17:11 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
