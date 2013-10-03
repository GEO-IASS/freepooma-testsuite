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
//   sum             - sum all the elements.
//   prod            - multiply all of the elements.
//   max             - find the maximum value.
//   min             - find the minimum value.
//   all             - returns true if all of the elements are != 0.
//   any             - returns true if any of the elements are != 0.
//   bitOr           - does a bitwise or of all of the elements.
//   bitAnd          - does a bitwise and of all of the elements.
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Functions
 * @brief
 * Reduction functions for Fields and Arrays.
 */

#ifndef POOMA_FUNCTIONS_REDUCTIONS_H
#define POOMA_FUNCTIONS_REDUCTIONS_H

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "Evaluator/Reduction.h"
#include "Utilities/WrappedInt.h"


//-----------------------------------------------------------------------------
// Specific global reduction functions for Fields.
//-----------------------------------------------------------------------------

/// Sum up the elements.

template<class Subject>
typename Subject::Element_t sum(const Subject &s)
{
  typename Subject::Element_t ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), s);
  return ret;
}

/// Compute the product of the elements.

template<class Subject>
typename Subject::Element_t prod(const Subject &s)
{
  typename Subject::Element_t ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, OpMultiplyAssign(), s);
  return ret;
}

/// Find the smallest element.

template<class Subject>
typename Subject::Element_t min(const Subject &s)
{
  typename Subject::Element_t ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, FnMinAssign(), s);
  return ret;
}

/// Find the largest element.

template<class Subject>
typename Subject::Element_t max(const Subject &s)
{
  typename Subject::Element_t ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, FnMaxAssign(), s);
  return ret;
}

/// Report if all of the elements are true.

template<class Subject>
bool all(const Subject &s)
{
  bool ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, FnAndAssign(), s);
  return ret;
}

/// Report if some of the elements are true.

template<class Subject>
bool any(const Subject &s)
{
  bool ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, FnOrAssign(), s);
  return ret;
}

/// Bitwise-or all of the elements together.

template<class Subject>
typename Subject::Element_t bitOr(const Subject &s)
{
  typename Subject::Element_t ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, OpBitwiseOrAssign(), s);
  return ret;
}

/// Bitwise-and all of the elements together.

template<class Subject>
typename Subject::Element_t bitAnd(const Subject &s)
{
  typename Subject::Element_t ret;
  Reduction<MainEvaluatorTag>().evaluate(ret, OpBitwiseAndAssign(), s);
  return ret;
}

#endif // POOMA_FUNCTIONS_REDUCTIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Reductions.h,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:49 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
