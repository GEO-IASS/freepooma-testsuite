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

#ifndef POOMA_EVALUATOR_REDUCTIONKERNEL_H
#define POOMA_EVALUATOR_REDUCTIONKERNEL_H

//-----------------------------------------------------------------------------
// Class:
//   ReductionKernel
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Evaluator
 * @brief
 * A ReductionKernel encapsulates reducing an expression on a domain.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Evaluator/EvaluatorTags.h"
#include "Evaluator/ReductionEvaluator.h"
#include "Evaluator/RequestLocks.h"
#include "Threads/PoomaCSem.h"
#include "Threads/PoomaSmarts.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * A ReductionKernel is a specific kind of iterate which
 * reduces a particular expression with over a given domain.
 *
 * A ReductionKernel IS-AN Interate. That means that it
 * has three primary functions:
 *   -# Construct given an expression and a domain. This will acquire
 *      locks on the data referenced by the expression.
 *   -# Destruct. This releases the locks on the data.
 *   -# Run the kernel by calling the member function run.
 */

template<class T, class Op, class Expr, class KernelTag>
class ReductionKernel : public Pooma::Iterate_t
{
public:

  typedef ReductionKernel<T, Op, Expr, KernelTag> This_t;

  //---------------------------------------------------------------------------
  // Construct from an Expr.
  // Build the kernel that will reduce the expression on the given domain.
  // Acquire locks on the data referred to by the expression.

  ReductionKernel(T &ret, const Op &op, const Expr &e,
		  Pooma::CountingSemaphore &csem);

  //---------------------------------------------------------------------------
  // Virtual Destructor.
  // Release locks on the data referred to by the expression.

  virtual ~ReductionKernel();

  //---------------------------------------------------------------------------
  // Do the reduction.

  virtual void run();

private:

  // The expression we will reduce.

  T &ret_m;
  Op op_m;
  Expr expr_m;
  Pooma::CountingSemaphore &csem_m;
};

//-----------------------------------------------------------------------------
//
// Implementation of the functions for ReductionKernel.
// Since ReductionKernel is templated on an expression type,
// there is little hope of instantiating by hand, so we put
// the definition in the header file.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Constructor
// Input an expression and record it for later use.

template<class T, class Op, class Expr, class KernelTag>
ReductionKernel<T, Op, Expr, KernelTag>::
ReductionKernel(T &ret, const Op &op, const Expr &e,
		Pooma::CountingSemaphore &csem)
  : Pooma::Iterate_t(Pooma::scheduler()),
    ret_m(ret), op_m(op), expr_m(e), csem_m(csem)
{
  // Request read lock.

  DataObjectRequest<ReadRequest> readReq(*this);
  engineFunctor(expr_m, readReq);
}

//---------------------------------------------------------------------------
// Destroy the kernel.
// Just let the expression destroy itself.
// Release locks and increment the semaphore.

template<class T, class Op, class Expr, class KernelTag>
ReductionKernel<T, Op, Expr, KernelTag>::~ReductionKernel()
{
  // Release read lock.

  DataObjectRequest<ReadRelease> readRelease;
  engineFunctor(expr_m, readRelease);
  
  // Increment semaphore.
  
  csem_m.incr();
}

//---------------------------------------------------------------------------
// Evaluate the kernel

template<class T, class Op, class Expr, class KernelTag>
void ReductionKernel<T, Op, Expr, KernelTag>::run()
{
  // Just evaluate the expression.
  
  ReductionEvaluator<KernelTag>::evaluate(ret_m, op_m, expr_m);
}

#endif     // POOMA_EVALUATOR_REDUCTIONKERNEL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ReductionKernel.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
