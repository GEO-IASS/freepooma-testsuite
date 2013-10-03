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
// Classes: 
//   Reduction base template
//   Reduction<SinglePatchEvaluatorTag>
//   Reduction<MuliPatchEvaluatorTag>
//----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_REDUCTION_H
#define POOMA_EVALUATOR_REDUCTION_H

/** @file
 * @ingroup Evaluator
 * @brief
 * Reduction performs global reductions on expressions by examining the 
 * engines that are participating in the expression and dispatching to custom 
 * code.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/NullDomain.h"
#include "Engine/Intersector.h"
#include "Engine/IntersectEngine.h"
#include "Evaluator/ReductionKernel.h"
#include "Evaluator/EvaluatorTags.h"
#include "Evaluator/WhereProxy.h"
#include "Threads/PoomaCSem.h"
#include "Utilities/PerformUpdate.h"

#include <vector>
#include <iterator>


/**
 * The point of this class is to input an expression with the
 * 'evaluate' member function and reduce it by breaking it up into
 * appropriate sub-blocks, looping over the whole domain, and
 * evaluating the expression at each point.
 */

template <class EvalTag>
struct Reduction
{ };


/**
 * This reduction is the one that gets called for a data-parallel expression.
 * It just determines the appropriate reduction from the types of the LHS and
 * RHS.  We don't need to do a blockAndEvaluate() because  all reductions
 * naturally involve some sort of blocking using counting semaphores. 
 * This approach is superior to blockAndEvaluate() because iterates not 
 * related to the reduction can continue to execute out-of-order.
 */

template <>
struct Reduction<MainEvaluatorTag>
{
  //---------------------------------------------------------------------------
  // Default ctor.

  Reduction() { }

  //---------------------------------------------------------------------------
  // Destructor

  ~Reduction() { }

  /// Helper to check validity of the expression, general version.

  template <class Expr>
  static inline bool checkValidity(const Expr &e, WrappedInt<false>)
  {
    return true;
  }

  /// Helper to check validity of the expression, version for fields.

  template <class Expr>
  static inline bool checkValidity(const Expr &e, WrappedInt<true>)
  {
    return e.centeringSize() == 1 && e.numMaterials() == 1;
  }

  /// Un-wrap where() expression operation and pass on to generic evaluator.

  template<class T, class Op, class Cond, class Expr>
  void evaluate(T &ret, const Op &op, const WhereProxy<Cond, Expr> &w) const
  {
    evaluate(ret, w.opMask(op), w.whereMask());
  }

  /// Input an expression and cause it to be reduced.
  /// We just pass the buck to a special reduction after updating
  /// the expression leafs and checking its validity (we can handle
  /// one subfield only).

  template<class T, class Op, class Expr>
  void evaluate(T &ret, const Op &op, const Expr &e) const
  {
    typedef typename EvaluatorTag1<Expr>::Evaluator_t Evaluator_t;
    PAssert(checkValidity(e, WrappedInt<Expr::hasRelations>()));
    forEach(e, PerformUpdateTag(), NullCombine());
    Reduction<Evaluator_t>().evaluate(ret, op, e());
 
    POOMA_INCREMENT_STATISTIC(NumReductions)
  }
};


/** Single-patch Reduction:
 *
 * The single patch version just passes the tag on to generate
 * a reduction kernel.
 */

template <>
struct Reduction<SinglePatchEvaluatorTag>
{
  //---------------------------------------------------------------------------
  // Default ctor.

  Reduction() { }

  //---------------------------------------------------------------------------
  // Destructor

  ~Reduction() { }

  //---------------------------------------------------------------------------
  // Input an expression and cause it to be reduced.
  // We just pass the buck to a special reduction.

  // Include versions expecting and not expecting counting semaphores.
  
  template<class T, class Op, class Expr>
  void evaluate(T &ret, const Op &op, const Expr &e,
		Pooma::CountingSemaphore &csem) const
  {
    typedef typename KernelTag1<Expr>::Kernel_t Kernel_t;

    Pooma::Iterate_t *iterate = 
      new ReductionKernel<T, Op, Expr, Kernel_t>(ret, op, e, csem);
    Pooma::scheduler().handOff(iterate);
  }
  
  template<class T, class Op, class Expr>
  void evaluate(T &ret, const Op &op, const Expr &e) const
  {
    Pooma::CountingSemaphore csem;
    csem.height(1);

    Pooma::scheduler().beginGeneration();

    evaluate(ret, op, e, csem);

    Pooma::scheduler().endGeneration();

    csem.wait();
  }
};


/** Multiple-patch Reduction:
 *
 * The multiple patch version makes patches and sends them out to
 * the single patch reduction.
 */

template <>
struct Reduction<MultiPatchEvaluatorTag>
{
  //---------------------------------------------------------------------------
  // Default ctor.

  Reduction() { }

  //---------------------------------------------------------------------------
  // Destructor

  ~Reduction() { }

  //---------------------------------------------------------------------------
  // Input an expression and cause it to be reduced according to the 
  // computational scheme:
  //   1. Perform the intersection calculation to deduce the patches that 
  //      computation will proceed on.
  //   2. Construct a counting sempahore with a height equal to the number 
  //      of patches.
  //   3. Construct a vector vals to hold the results from each patch 
  //      reduction.
  //   4. For each patch, take a view over the patch and do a reduction of 
  //      the resulting array. Increment the semaphore and store the result in 
  //      the appropriate slot of the vals vector.
  //   5. Wait for all reductions to finish.
  //   6. Finish by doing an immediate reduction of the vals array.

  template<class T, class Op, class Expr>
  void evaluate(T &ret, const Op &op, const Expr &e) const
  {
    typedef Intersector<Expr::dimensions> Inter_t;
    Inter_t inter;

    expressionApply(e, IntersectorTag<Inter_t>(inter));
  
    const int n = std::distance(inter.begin(), inter.end());
    Pooma::CountingSemaphore csem;
    csem.height(n);
    T *vals = new T[n];

    Pooma::scheduler().beginGeneration();
    
    typename Inter_t::const_iterator i = inter.begin();
    int j = 0;
    while (j < n)
      {
        Reduction<SinglePatchEvaluatorTag>().
          evaluate(vals[j], op, e(*i), csem);
        ++i; ++j;
      }

    Pooma::scheduler().endGeneration();

    csem.wait();

    ret = vals[0];
    for (j = 1; j < n; j++)
      op(ret, vals[j]);
    delete [] vals;
  }
};


#endif // POOMA_EVALUATOR_REDUCTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Reduction.h,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
