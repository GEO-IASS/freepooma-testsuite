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
// Class: Evaluator
//----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_EVALUATOR_H
#define POOMA_EVALUATOR_EVALUATOR_H

/** @file
 * @ingroup Evaluator
 * @brief
 * Evaluator evaluates expressions by examining the engines that are
 * participating in the expression and dispatching to custom code.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Evaluator/CompressibleEval.h"
#include "Evaluator/ExpressionKernel.h"
#include "Evaluator/EvaluatorTags.h"
#include "Engine/Intersector.h"
#include "Engine/IntersectEngine.h"
#include "Engine/NotifyEngineWrite.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * The point of this class is to input an expression with the
 * 'evaluate' member function and evaluate it by breaking it up into
 * appropriate sub-blocks, looping over the whole domain, and
 * evaluating the expression at each point.
 */

template <class EvalTag>
struct Evaluator
{
};


/**
 * This evaluator is the one that gets called for a data-parallel expression.
 * It just determines the appropriate evaluator from the types of the LHS and
 * RHS.  Also, we block if appropriate.
 */

template <>
struct Evaluator<MainEvaluatorTag>
{
  /// Default ctor.

  Evaluator() { }

  /// Destructor

  ~Evaluator() { }

  /// evaluate(expression)
  /// Input an expression and cause it to be evaluated.
  /// We just pass the buck to a special evaluator.

  template <class LHS, class RHS, class Op>
  void evaluate(const LHS& lhs, const Op& op, const RHS& rhs) const
  {
    typedef typename EvaluatorTag<LHS, RHS>::Evaluator_t Eval_t;
    Evaluator<Eval_t> evaluator;

    Pooma::beginExpression();
    evaluator.evaluate(lhs(), op, rhs());
    notifyEngineWrite(lhs.engine());
    Pooma::endExpression();
    
    POOMA_INCREMENT_STATISTIC(NumExpressions)
  }

  /// evaluateZeroBased(expression)
  /// Input an expression and cause it to be evaluated.
  /// We just pass the buck to a special evaluator.
  /// This version does not bother to take views of the expression
  /// since the caller is assuring us they are already zero-based.

  template <class LHS, class RHS, class Op>
  void evaluateZeroBased(const LHS& lhs, const Op& op, const RHS& rhs) const
  {
    typedef typename EvaluatorTag<LHS, RHS>::Evaluator_t Eval_t;
    Evaluator<Eval_t> evaluator;

    Pooma::beginExpression();
    evaluator.evaluate(lhs, op, rhs);
    notifyEngineWrite(lhs.engine());
    Pooma::endExpression();
    
    POOMA_INCREMENT_STATISTIC(NumZBExpressions)
  }
};


/**
 * The single patch version just passes the tag on to generate
 * an expression kernel.
 */

template <>
struct Evaluator<SinglePatchEvaluatorTag>
{
  /// Default ctor.

  Evaluator() { }

  /// Destructor

  ~Evaluator() { }

  /// evaluate(expression)
  /// Input an expression and cause it to be evaluated.
  /// We just pass the buck to a special evaluator.

  template <class LHS, class RHS, class Op>
  void evaluate(const LHS& lhs, const Op& op, const RHS& rhs) const
  {
    typedef typename KernelTag<LHS,RHS>::Kernel_t Kernel_t;
    Pooma::Iterate_t *iterate = ::generateKernel(lhs, op, rhs, Kernel_t());
    Pooma::scheduler().handOff(iterate);
  }
};


/**
 * The multiple patch version makes patches and sends them out to
 * the single patch evaluator.
 */

template <>
struct Evaluator<MultiPatchEvaluatorTag>
{
  /// Default ctor.

  Evaluator() { }

  /// Destructor

  ~Evaluator() { }

  /// evaluate(expression)
  /// Input an expression and cause it to be evaluated.
  /// We just pass the buck to a special evaluator.

  template <class LHS, class RHS, class Op>
  void evaluate(const LHS& lhs, const Op& op, const RHS& rhs) const
  {
    typedef Intersector<LHS::dimensions> Inter_t;
    Inter_t inter;

    expressionApply(lhs, IntersectorTag<Inter_t>(inter));
    expressionApply(rhs, IntersectorTag<Inter_t>(inter));
  
    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
    {
      Evaluator<SinglePatchEvaluatorTag>().evaluate(lhs(*i), op, rhs(*i));
      ++i;
    }
    
    POOMA_INCREMENT_STATISTIC(NumMultiPatchExpressions)
    POOMA_INCREMENT_STATISTIC_BY(NumLocalPatchesEvaluated, inter.size())
  }
};


#endif // POOMA_EVALUATOR_EVALUATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Evaluator.h,v $   $Author: richard $
// $Revision: 1.61 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
