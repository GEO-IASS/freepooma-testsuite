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
// Structs:
//   EvaluatorTag1<Expr>
//   EvaluatorTag<LHS, RHS>
//   EvaluatorCombine<LHSTag, RHSTag>
// LeafFunctors:
//   LeafFunctor<Scalar<T>, EvaluatorTypeTag>
//   LeafFunctor<A, EvaluatorTypeTag>
// Combiners:
//   Combine2<Eval1, Eval2, Op, EvaluatorCombineTag>
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_EVALUATORTAGS_H
#define POOMA_EVALUATOR_EVALUATORTAGS_H

/** @file
 * @ingroup Evaluator
 * @brief
 * Evaluator Tags are used for picking the appropriate evaluator given the
 * engines in an expression.
 *
 * The external interface for EvaluatorTags are the traits
 * - EvaluatorTag1<Expr>::Evaluator_t
 * - EvaluatorTag<LHS, RHS>::Evaluator_t
 *
 * which yields an evaluator tag, given the expression type or the types for 
 * the left hand side and right hand side.  To add new evaluators or new
 * engines, specialize the EvaluatorEngineTraits struct to give the evaluator
 * tag for each engine, and specialize the EvaluatorCombine
 * struct to determine how to chose a new evaluator given two evaluators.
 *
 * Evaluator Tags are used for picking the appropriate evaluator given the
 * engines in an expression.  Each evaluator tag represents a set of engines
 * that it is capable of dealing with.
 *
 * This file provides the following
 * interface:
 * - EvaluatorTag1<Expr>::Evaluator_t
 *        this is the main function that produces an evaluator tag for the
 *        expression
 * - EvaluatorTag<LHS,RHS>::Evaluator_t
 *        this is the main function that produces an evaluator tag for the
 *        expression given the types on the LHS and RHS
 * - EvaluatorCombine<LHSEval, RHSEval>::Evaluator_t
 *        combines the evaluator for the left hand side and the right hand
 *        side to yield an evalutor type for the expression.
 *
 * To add new evaluator types, the user must specialize EvaluatorCombine for
 * the new tags and combinations of the new tags with old ones.
 *
 * To add new engines tags, you must specialize the struct
 * EvaluatorEngineTraits<EngineTag> to return the appropriate evaluator
 * for that engine.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Evaluator/EngineTraits.h"
#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * EvaluatorCombine<LHSEval, RHSEval>
 *
 * This struct computes an evaluator that can evaluate an expression
 * given evaluator tags for the left and right hand sides.
 *
 * The current rules are:
 *   -# RemoteMultiPatchEvaluatorTag and anything -> 
 *      RemoteMultiPatchEvaluatorTag
 *   -# MultiPatchEvaluatorTag and SinglePatchEvaluatorTag -> 
 *      MultiPatchEvaluatorTag
 *   -# MultiPatchEvaluatorTag and RemoteSinglePatchEvaluatorTag ->
 *      RemoteMultiPatchEvaluatorTag
 *   -# RemoteSinglePatchEvaluatorTag and SinglePatchEvaluatorTag ->
 *      RemoteSinglePatchEvaluatorTag
 *   -# Combining anything with itself is a no-op.
 */

// Implementation of rule #1 and rule #5 for RemoteMultiPatchEvaluatorTag.
// Note: if anyone adds a new evaluator type, this probably needs to 
// be changed.

template<class Eval1, class Eval2>
struct EvaluatorCombine
{
 EvaluatorCombine() {}
 ~EvaluatorCombine() {}
 typedef RemoteMultiPatchEvaluatorTag Evaluator_t;
};

// Implementation of rule #2.

template<>
struct EvaluatorCombine<SinglePatchEvaluatorTag, 
                        MultiPatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef MultiPatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorCombine<MultiPatchEvaluatorTag, 
                        SinglePatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef MultiPatchEvaluatorTag Evaluator_t;
};

// Implementation of rule #3.

template<>
struct EvaluatorCombine<RemoteSinglePatchEvaluatorTag,
                        MultiPatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef RemoteMultiPatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorCombine<MultiPatchEvaluatorTag, 
                        RemoteSinglePatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef RemoteMultiPatchEvaluatorTag Evaluator_t;
};

// Implementation of rule #4.

template<>
struct EvaluatorCombine<RemoteSinglePatchEvaluatorTag, 
                        SinglePatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef RemoteSinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorCombine<SinglePatchEvaluatorTag, 
                        RemoteSinglePatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef RemoteSinglePatchEvaluatorTag Evaluator_t;
};

// Implementation of rule #5 for everything except 
// RemoteMultiPatchEvaluatorTag.

template<>
struct EvaluatorCombine<MultiPatchEvaluatorTag, 
                        MultiPatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef MultiPatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorCombine<RemoteSinglePatchEvaluatorTag, 
                        RemoteSinglePatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef RemoteSinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorCombine<SinglePatchEvaluatorTag, 
                        SinglePatchEvaluatorTag>
{
  EvaluatorCombine() {}
  ~EvaluatorCombine() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};


//-----------------------------------------------------------------------------
// LeafFunctor constructs the determine the type of evaluator associated with
// a particular leaf. We handle scalars specially and then assume that
// everything else has an Engine_t yypedef that return that Engine's tag, which
// is then passed to EvaluatorEngineTraits.
//-----------------------------------------------------------------------------

template<class T>
struct LeafFunctor<Scalar<T>, EvaluatorTypeTag>
{
 LeafFunctor() {}
 ~LeafFunctor() {}
 typedef typename EvaluatorEngineTraits<ScalarEngineTag>::Evaluator_t Type_t;
};

template<class A>
struct LeafFunctor<A, EvaluatorTypeTag>
{
  LeafFunctor() {}
  ~LeafFunctor() {}
  typedef typename
    EvaluatorEngineTraits<typename A::Engine_t::Tag_t>::Evaluator_t Type_t;
};

template<class Eval1,class Eval2,class Op>
struct Combine2<Eval1, Eval2, Op, EvaluatorCombineTag>
{
  Combine2() {}
  ~Combine2() {}
  typedef typename EvaluatorCombine<Eval1, Eval2>::Evaluator_t Type_t;
};


//-----------------------------------------------------------------------------
// EvaluatorTag<Expr>
//
// This struct computes the evaluator tag for a single expression..
//-----------------------------------------------------------------------------

template<class Expr>
struct EvaluatorTag1
{
  EvaluatorTag1() {}
  ~EvaluatorTag1() {}

  typedef typename LeafFunctor<Expr, EvaluatorTypeTag>::Type_t Evaluator_t;
};


//-----------------------------------------------------------------------------
// EvaluatorTag<LHS, RHS>
//
// Finally, this struct computes the evaluator tag for the whole expression
// given the types of the arrays on the left and right hand sides.
//-----------------------------------------------------------------------------

template<class LHS, class RHS>
struct EvaluatorTag
{
  EvaluatorTag() {}
  ~EvaluatorTag() {}
  typedef typename LeafFunctor<LHS, EvaluatorTypeTag>::Type_t LHSEval_t;
  typedef typename LeafFunctor<RHS, EvaluatorTypeTag>::Type_t RHSEval_t;

  typedef typename EvaluatorCombine<LHSEval_t, RHSEval_t>::Evaluator_t 
    Evaluator_t;
};


#endif     // POOMA_EVALUATOR_EVALUATORTAGS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EvaluatorTags.h,v $   $Author: richard $
// $Revision: 1.19 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
