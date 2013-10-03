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
// KernelEvaluator<CompressibleViewKernelTag>
// KernelEvaluator<CompressibleKernelTag>
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_COMPRESSIBLEEVAL_H
#define POOMA_EVALUATOR_COMPRESSIBLEEVAL_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * These evaluators are used in the iterates that perform evaluation on
 * expressions with CompressibleBricks. 
 *
 * There are two versions here.
 * If there are any Bricks or BrickViews on the RHS, then it doesn't make
 * sense to do compressed assignment, so we have CompressibleViewEvalTag
 * which just takes a BrickView of the LHS and calls the inline evaluator.
 * If the expression is completely compressible, then we invoke the evaluator
 * with CompressibleEvalTag, which checks the compression status and perhaps
 * performs a compressed assign.
 *
 * This file should really be called CompressibleKernel.h.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "Engine/CompressibleBrick.h"
#include "Evaluator/KernelTags.h"
#include "Evaluator/CompressibleEngines.h"
#include "Evaluator/InlineEvaluator.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<>
struct KernelEvaluator<CompressibleViewKernelTag>
{
  template<class LHS, class Op, class RHS>
  static void evaluate(const LHS &lhs, const Op &op, const RHS &rhs)
  {
    // Try this for now.  We just pass the uncompressed engine
    // to the inline evaluator.  If you want to pack this in an
    // array go ahead, but then you'll have to include Array.h.
    KernelEvaluator<InlineKernelTag>().
      evaluate(engineFunctor(lhs, UnCompressedViewEngine()),
	       op, rhs);
    
    POOMA_INCREMENT_STATISTIC(NumAssignsRequiringUnCompression)
  }
};

template<>
struct KernelEvaluator<CompressibleKernelTag>
{
  template<class LHS, class Op, class RHS>
  static void evaluate(const LHS &lhs, const Op &op, const RHS &rhs)
  {
    // If everybody is compressed, then we do a compressed assign,
    // provided the left-hand-side is viewing the entire compressed
    // block or the value being assigned is the same as the compressed
    // value on the left-hand-side.
    // If either side of the expression is uncompressed, then we perform
    // an assign to a BrickView.  When the BrickView goes away, the block
    // will try to compress itself.

    // Note: the CompressibleBlockController is NOT locked at this
    // point, so asking "are you compressed" only gives you the answer
    // at a particular point in time. The current evaluation mechanism
    // does not allow multiple iterates to be writing to sub-blocks
    // of a CompressibleBrick simultaneously, and the parse thread
    // should never be changing the LHS while iterates are outstanding,
    // so this should be safe. If we later change to allowing all writes
    // within a single generation to occur in parallel, we'll need to
    // make sure that this remains thread-safe. 

    typedef typename LHS::Element_t LHST_t;
    typedef typename RHS::Element_t RHST_t;
    
    if (engineFunctor(lhs, Compressed()) &&
	engineFunctor(rhs, Compressed()))
      {
        // Get the compressed values on the LHS and RHS. Make a copy of the
        // LHS value and apply the operation to this copy.
        
        LHST_t &l   = engineFunctor(lhs, CompressedReadWrite());
        LHST_t test = l;
        RHST_t r    = engineFunctor(rhs, CompressedRead());
        op(test, r);
        
        // If the test value has not changed, we're done. If it has and the LHS
        // represents the entire view, we just need to assign the test value to
        // the LHS. Otherwise, we need to uncompress and do the operation.
        
        if (test != l)
          {
	    if (engineFunctor(lhs, CompressedBrickIsWholeView()))
	      {
	        l = test;
	        POOMA_INCREMENT_STATISTIC(NumCompressedAssigns)
	      }
	    else
	      {
	        KernelEvaluator<CompressibleViewKernelTag>::
                  evaluate(lhs, op, rhs);
	      }
          }
        else
          {
            POOMA_INCREMENT_STATISTIC(NumCompressedAssigns)
          }
      }
    else
      {
        KernelEvaluator<CompressibleViewKernelTag>::evaluate(lhs, op, rhs);
      }
  }
};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CompressibleEval.h,v $   $Author: richard $
// $Revision: 1.27 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo


