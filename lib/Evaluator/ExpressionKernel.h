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

#ifndef POOMA_EVALUATOR_EXPRESSIONKERNEL_H
#define POOMA_EVALUATOR_EXPRESSIONKERNEL_H

//-----------------------------------------------------------------------------
// Class:
// ExpressionKernel
//-----------------------------------------------------------------------------


//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * An ExpressionKernel encapsulates evaluating an expression on a domain.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/PoomaSmarts.h"
#include "Evaluator/InlineEvaluator.h"
#include "Evaluator/EvaluatorTags.h"
#include "Evaluator/RequestLocks.h"
#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Pooma/Configuration.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * An ExpressionKernel is a specific kind of iterate which
 * evaluates a particular expression with a loop over a given domain.
 *
 * An ExpressionKernel IS-AN ExpressionTask. That means that it
 * has three primary functions:
 * -# Construct given an expression and a domain. This will acquire
 *    locks on the data referenced by the expression.
 * -# Destruct. This releases the locks on the data.
 * -# Run the kernel by calling the member function run.
 */

template<class LHS,class Op,class RHS,class EvalTag>
class ExpressionKernel : public Pooma::Iterate_t
{
public:

  typedef ExpressionKernel<LHS,Op,RHS,EvalTag> This_t;
  //
  // Construct from an Expr.
  // Build the kernel that will evaluate the expression on the given domain.
  // Acquire locks on the data referred to by the expression.
  //
  ExpressionKernel(const LHS&,const Op&,const RHS&);

  //
  // Virtual Destructor.
  // Release locks on the data referred to by the expression.
  //
  virtual ~ExpressionKernel();

  //
  // Do the loop.
  //
  virtual void run();

private:

  // The expression we will evaluate.
  LHS lhs_m;
  Op  op_m;
  RHS rhs_m;
};

//////////////////////////////////////////////////////////////////////
//
// Implementation of the functions for ExpressionKernel.
// Since ExpressionKernel is templated on an expression type,
// there is little hope of instantiating by hand, so we put
// the definition in the header file.
//
//////////////////////////////////////////////////////////////////////

//
// Constructor
// Input an expression and record it for later use.
//
template<class LHS,class Op,class RHS,class EvalTag>
ExpressionKernel<LHS,Op,RHS,EvalTag>::
ExpressionKernel(const LHS& lhs,const Op& op,const RHS& rhs)
  : Pooma::Iterate_t(Pooma::scheduler()),
    lhs_m(lhs), op_m(op), rhs_m(rhs)
{
  hintAffinity(engineFunctor(lhs, DataObjectRequest<BlockAffinity>()));

  // Request locks

  // First make the write request.
  // The write request tag remembers the data block
  // for the left hand side so we can check if there is a stencil
  // on the right.

  DataObjectRequest<WriteRequest> writeReq(*this);
  engineFunctor(lhs_m, writeReq);

  // Now make the read request.
  // Use the remembered write request block to check and see
  // if that block is used on the right.  If so, do a notify to the
  // iterated instead of a request of the data object. 

  DataObjectRequest<ReadRequest> readReq(writeReq);
  engineFunctor(rhs_m, readReq);
}

//
// Destroy the kernel.
// Just let the expression destroy itself.
// Release locks.
//

template<class LHS,class Op,class RHS,class EvalTag>
ExpressionKernel<LHS,Op,RHS,EvalTag>::~ExpressionKernel()
{
  // The write request remembers the data block it sees on the left
  // so that it can check and see if it finds it on the right.

  DataObjectRequest<WriteRelease> writeReq;
  engineFunctor(lhs_m, writeReq);

  // The read request checks to see if the data object for the left
  // appears on the right.  If it does, it doesn't do a release for it
  // since a request wasn't generated above.

  DataObjectRequest<ReadRelease> readReq(writeReq);
  engineFunctor(rhs_m, readReq);
}

//
// Evaluate the kernel
// Just tell an InlineEvaluator to do it.
//

template<class LHS,class Op,class RHS,class EvalTag>
void
ExpressionKernel<LHS,Op,RHS,EvalTag>::run()
{
  // Just evaluate the expression.
  KernelEvaluator<EvalTag>::evaluate(lhs_m,op_m,rhs_m);

  // we could release the locks here instead of in the dtor.
}

template<class LHS,class Op,class RHS,class EvalTag>
inline static
ExpressionKernel<LHS,Op,RHS,EvalTag>*
generateKernel(const LHS& lhs, const Op& op, const RHS& rhs, const EvalTag&)
{
  return new ExpressionKernel<LHS,Op,RHS,EvalTag>(lhs, op, rhs);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_EXPRESSIONKERNEL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ExpressionKernel.h,v $   $Author: richard $
// $Revision: 1.47 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
