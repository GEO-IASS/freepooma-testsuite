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
// Class:
// LoopApplyEvaluator
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_LOOPAPPLY_H
#define POOMA_EVALUATOR_LOOPAPPLY_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * LoopApplyEvaluator is a simple class that wraps a set of 7 functions that
 * provide efficient loops over Interval<Dim> type domains and call operator()
 * with the integers on a user provided functor.
 *
 * For example, calling:
 *
 * LoopApplyEvaluator::evaluate(op, Interval<2>(2, 2));
 *
 * would cause the following statements to be performed:
 *
 * <PRE>
 * op(0, 0);
 * op(0, 1);
 * op(1, 0);
 * op(1, 1);
 * </PRE>
 *
 * This class is sufficiently general that we could rewrite the
 * InlineEvaluator class using the LoopApplyEvaluator by writing a functor
 * op(i, j) that evaluates equationOp(lhs(i, j), rhs(i, j));
 * For now, we are using this evaluator to serialize engines and in the
 * ExtendedPatchEvaluator which evaluates user's scalar code inside the
 * op(i, j);
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/WrappedInt.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// LoopApplyEvaluator:
//
//-----------------------------------------------------------------------------

struct LoopApplyEvaluator
{
  // The main interface takes an operator and a domain object and dispatches
  // to a specific function with the right number of loops for the dimension
  // of the domain.

  template<class Op, class Dom>
  inline static
  void evaluate(const Op &op, const Dom &domain)
  {
    CTAssert(Dom::unitStride);
    evaluate(op, domain, WrappedInt<Dom::dimensions>());
  }

  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<1>)
  {
    int f0 = domain[0].first();
    int e0 = domain[0].last();
#pragma omp parallel for if (e0-f0 > 512)
    for (int i0 = f0; i0 <= e0; ++i0)
      op(i0);
  }
  
  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<2>)
  {
    int f0 = domain[0].first();
    int f1 = domain[1].first();
    int e0 = domain[0].last();
    int e1 = domain[1].last();
#pragma omp parallel for
    for (int i1 = f1; i1 <= e1; ++i1)
      for (int i0 = f0; i0 <= e0; ++i0)
	op(i0, i1);
  }
  
  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<3>)
  {
    int f0 = domain[0].first();
    int f1 = domain[1].first();
    int f2 = domain[2].first();
    int e0 = domain[0].last();
    int e1 = domain[1].last();
    int e2 = domain[2].last();
#pragma omp parallel for
    for (int i2 = f2; i2 <= e2; ++i2)
      for (int i1 = f1; i1 <= e1; ++i1)
	for (int i0 = f0; i0 <= e0; ++i0)
	  op(i0,i1,i2);
  }
  
  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<4>)
  {
    int f0 = domain[0].first();
    int f1 = domain[1].first();
    int f2 = domain[2].first();
    int f3 = domain[3].first();
    int e0 = domain[0].last();
    int e1 = domain[1].last();
    int e2 = domain[2].last();
    int e3 = domain[3].last();
#pragma omp parallel for
    for (int i3 = f3; i3 <= e3; ++i3)
      for (int i2 = f2; i2 <= e2; ++i2)
	for (int i1 = f1; i1 <= e1; ++i1)
	  for (int i0 = f0; i0 <= e0; ++i0)
	    op(i0,i1,i2,i3);
  }

  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<5>)
  {
    int f0 = domain[0].first();
    int f1 = domain[1].first();
    int f2 = domain[2].first();
    int f3 = domain[3].first();
    int f4 = domain[4].first();
    int e0 = domain[0].last();
    int e1 = domain[1].last();
    int e2 = domain[2].last();
    int e3 = domain[3].last();
    int e4 = domain[4].last();
#pragma omp parallel for
    for (int i4 = f4; i4 <= e4; ++i4)
      for (int i3 = f3; i3 <= e3; ++i3)
	for (int i2 = f2; i2 <= e2; ++i2)
	  for (int i1 = f1; i1 <= e1; ++i1)
	    for (int i0 = f0; i0 <= e0; ++i0)
	      op(i0,i1,i2,i3,i4);
  }

  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<6>)
  {
    int f0 = domain[0].first();
    int f1 = domain[1].first();
    int f2 = domain[2].first();
    int f3 = domain[3].first();
    int f4 = domain[4].first();
    int f5 = domain[5].first();
    int e0 = domain[0].last();
    int e1 = domain[1].last();
    int e2 = domain[2].last();
    int e3 = domain[3].last();
    int e4 = domain[4].last();
    int e5 = domain[5].last();
#pragma omp parallel for
    for (int i5 = f5; i5 <= e5; ++i5)
      for (int i4 = f4; i4 <= e4; ++i4)
	for (int i3 = f3; i3 <= e3; ++i3)
	  for (int i2 = f2; i2 <= e2; ++i2)
	    for (int i1 = f1; i1 <= e1; ++i1)
	      for (int i0 = f0; i0 <= e0; ++i0)
		op(i0,i1,i2,i3,i4,i5);
  }

  template<class Op, class Domain>
  inline static void evaluate(const Op &op, const Domain &domain, WrappedInt<7>)
  {
    int f0 = domain[0].first();
    int f1 = domain[1].first();
    int f2 = domain[2].first();
    int f3 = domain[3].first();
    int f4 = domain[4].first();
    int f5 = domain[5].first();
    int f6 = domain[6].first();
    int e0 = domain[0].last();
    int e1 = domain[1].last();
    int e2 = domain[2].last();
    int e3 = domain[3].last();
    int e4 = domain[4].last();
    int e5 = domain[5].last();
    int e6 = domain[6].last();
#pragma omp parallel for
    for (int i6 = f6; i6 <= e6; ++i6)
      for (int i5 = f5; i5 <= e5; ++i5)
	for (int i4 = f4; i4 <= e4; ++i4)
	  for (int i3 = f3; i3 <= e3; ++i3)
	    for (int i2 = f2; i2 <= e2; ++i2)
	      for (int i1 = f1; i1 <= e1; ++i1)
		for (int i0 = f0; i0 <= e0; ++i0)
		  op(i0,i1,i2,i3,i4,i5,i6);
  }

};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_LOOPAPPLY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LoopApply.h,v $   $Author: richi $
// $Revision: 1.10 $   $Date: 2004/11/04 20:07:26 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
