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

#ifndef POOMA_EVALUATOR_REDUCTIONEVALUATOR_H
#define POOMA_EVALUATOR_REDUCTIONEVALUATOR_H

//-----------------------------------------------------------------------------
// Class: 
//   ReductionEvaluator<InlineKernelTag>
//   ReductionEvaluator<CompressibleKernelTag>
//   CompressibleReduce<T, Op>
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Evaluator
 * @brief
 * ReductionEvaluator<InlineKernelTag> reduces expressions by inlining a 
 * simple loop. ReductionEvaluator<CompressibleKernelTag> can optionally take
 * advantage of compression.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/EngineFunctor.h"
#include "Evaluator/CompressibleEngines.h"
#include "Evaluator/KernelTags.h"
#include "PETE/OperatorTags.h"
#include "Pooma/PoomaOperatorTags.h"
#include "Utilities/WrappedInt.h"
#include "Utilities/PAssert.h"
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif


/**
 * Traits class defining identity element for type T under
 * operation Op.  Needs to be specialized for Op and possibly T.
 */

template<class Op, class T>
struct ReductionTraits {
};

template<class T>
struct ReductionTraits<OpAddAssign, T> {
  static inline T identity() { return T(0); }
};

template<class T>
struct ReductionTraits<OpMultiplyAssign, T> {
  static inline T identity() { return T(1); }
};

template<class T>
struct ReductionTraits<FnMinAssign, T> {
  static inline T identity() { return std::numeric_limits<T>::max(); }
};

template<class T>
struct ReductionTraits<FnMaxAssign, T> {
  static inline T identity() { return std::numeric_limits<T>::min(); }
};

template<class T>
struct ReductionTraits<FnOrAssign, T> {
  static inline T identity() { return T(false); }
};

template<class T>
struct ReductionTraits<FnAndAssign, T> {
  static inline T identity() { return T(true); }
};

template<class T>
struct ReductionTraits<OpBitwiseOrAssign, T> {
  static inline T identity() { return T(); }
};

template<class T>
struct ReductionTraits<OpBitwiseAndAssign, T> {
  static inline T identity() { return ~T(); }
};


/**
 * Class to hold static array for partial reduction results
 * and routine for final reduction. Two versions, one dummy
 * for non-OpenMP, one for OpenMP operation.
 */

#ifndef _OPENMP
template<class T>
struct PartialReduction {
  static inline void init() {}
  inline void storePartialResult(const T& result)
  {
    answer = result;
  }
  template <class Op>
  inline void reduce(T& ret, const Op&)
  {
    ret = answer;
  }
  T answer;
};
#else
template<class T>
struct PartialReduction {
  static inline void init()
  {
    if (!answer)
      answer = new T[omp_get_max_threads()];
  }
  inline void storePartialResult(const T& result)
  {
    int n = omp_get_thread_num();
    answer[n] = result;
    if (n == 0)
      num_threads = omp_get_num_threads();
  }
  template <class Op>
  inline void reduce(T& ret, const Op& op)
  {
    T res = answer[0];
    for (int i = 1; i<num_threads; ++i)
      op(res, answer[i]);
    ret = res;
  }
  int num_threads;
  static T *answer;
};
template <class T>
T *PartialReduction<T>::answer = NULL;
#endif


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class KernelTag>
struct ReductionEvaluator;

/**
 * The point of this class is to input an expression with the
 * 'evaluate' member function and reduce it by looping over the
 * whole domain.
 *
 * This is the simplest possible reduction. It makes the simplifying
 * assumption that expression passed in can handle random access to 
 * all of its elements efficiently.
 */

template<>
struct ReductionEvaluator<InlineKernelTag>
{

  //---------------------------------------------------------------------------
  // Input an expression and cause it to be evaluated.
  // All this template function does is extract the domain
  // from the expression and call evaluate on that.

  template<class T, class Op, class Expr>
  static void evaluate(T &ret, const Op &op, const Expr &e)
  {
    typedef typename Expr::Domain_t Domain_t;

    // The reduction evaluators assume unit-stride zero-based domains.
    CTAssert(Domain_t::unitStride);
    for (int i=0; i<Domain_t::dimensions; ++i)
      PAssert(e.domain()[i].first() == 0);

    PartialReduction<T>::init();
    evaluate(ret, op, e, e.domain(),
      WrappedInt<Domain_t::dimensions>());
  }

  //---------------------------------------------------------------------------
  // This is the function both of the above functions call.
  // It adds a third argument which is a tag class templated on
  // the dimension of the domain.
  //
  // This parameter lets us specialize the function based on
  // that dimension.
  //
  // Some day, we will figure out how to specialize template 
  // member functions outside the class declaration...
  //
  // These functions are all inline for efficiency. That means that if
  // they are being used at the user level we will get the optimization
  // of recognizing multiple uses of a single Array on the right hand
  // side.
  //
  // There are seven specializations here, for dimension 1 through 7.
  // Rather than use template metaprograms for these seven cases we
  // simply enumerate them explicitly.  This is done to reduce the
  // burden on the compiler, which would otherwise have to jump through
  // a bunch of hoops to get the code that is here.
  //
  // For each of the specializations it builds a nested loop for each
  // dimension. Each loop is constructed last() from the
  // appropriate dimension of the zero-based domain.
  //
  // NOTE: These loops assume that the domain passed in is a unit-stride
  // domain starting at 0.  Assertions are made to make sure this is true
  // in the dispatching evaluate.

  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<1>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();

    PartialReduction<T> reduction;
#pragma omp parallel if (e0 > 512)
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i0 = 0; i0 < e0; ++i0)
        op(answer, localExpr.read(i0));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }

  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<2>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();
    int e1 = domain[1].length();

    PartialReduction<T> reduction;
#pragma omp parallel
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i1 = 0; i1 < e1; ++i1)
	for (int i0 = 0; i0 < e0; ++i0)
	  op(answer, localExpr.read(i0, i1));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }
  
  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<3>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    
    PartialReduction<T> reduction;
#pragma omp parallel
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i2 = 0; i2 < e2; ++i2)
	for (int i1 = 0; i1 < e1; ++i1)
	  for (int i0 = 0; i0 < e0; ++i0)
	    op(answer, localExpr.read(i0, i1, i2));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }

  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<4>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    
    PartialReduction<T> reduction;
#pragma omp parallel
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i3 = 0; i3 < e3; ++i3)
	for (int i2 = 0; i2 < e2; ++i2)
	  for (int i1 = 0; i1 < e1; ++i1)
	    for (int i0 = 0; i0 < e0; ++i0)
	      op(answer, localExpr.read(i0, i1, i2, i3));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }

  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<5>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    int e4 = domain[4].length();
    
    PartialReduction<T> reduction;
#pragma omp parallel
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i4 = 0; i4 < e4; ++i4)
	for (int i3 = 0; i3 < e3; ++i3)
	  for (int i2 = 0; i2 < e2; ++i2)
	    for (int i1 = 0; i1 < e1; ++i1)
	      for (int i0 = 0; i0 < e0; ++i0)
		op(answer, localExpr.read(i0, i1, i2, i3, i4));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }

  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<6>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    int e4 = domain[4].length();
    int e5 = domain[5].length();
    
    PartialReduction<T> reduction;
#pragma omp parallel
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i5 = 0; i5 < e5; ++i5)
	for (int i4 = 0; i4 < e4; ++i4)
	  for (int i3 = 0; i3 < e3; ++i3)
	    for (int i2 = 0; i2 < e2; ++i2)
	      for (int i1 = 0; i1 < e1; ++i1)
		for (int i0 = 0; i0 < e0; ++i0)
		  op(answer, localExpr.read(i0, i1, i2, i3, i4, i5));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }

  template<class T, class Op, class Expr, class Domain>
  inline static void evaluate(T &ret, const Op &op, const Expr &e,
    const Domain &domain, WrappedInt<7>)
  {
    Expr localExpr(e);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    int e4 = domain[4].length();
    int e5 = domain[5].length();
    int e6 = domain[6].length();
    
    PartialReduction<T> reduction;
#pragma omp parallel
    {
      T answer = ReductionTraits<Op, T>::identity();
#pragma omp for nowait
      for (int i6 = 0; i6 < e6; ++i6)
	for (int i5 = 0; i5 < e5; ++i5)
	  for (int i4 = 0; i4 < e4; ++i4)
	    for (int i3 = 0; i3 < e3; ++i3)
	      for (int i2 = 0; i2 < e2; ++i2)
		for (int i1 = 0; i1 < e1; ++i1)
		  for (int i0 = 0; i0 < e0; ++i0)
		    op(answer, localExpr.read(i0, i1, i2, i3, i4, i5, i6));
      reduction.storePartialResult(answer);
    }
    reduction.reduce(ret, op);
  }
};


//-----------------------------------------------------------------------------
// This class handles the evaluation of a reduction from a single compressed
// value. The current possibilies are:
//   o sum:    N * val
//   o prod:   val^N
//   o min:    val
//   o max:    val
//   o any:    val
//   o all:    val
//   o bitOr:  val
//   o bitAnd: val
//-----------------------------------------------------------------------------

template<class T, class Op>
struct CompressibleReduce
{
  template<class T1>
  inline static void evaluate(T &ret, const Op &, const T1 &val, int)
  {
    ret = static_cast<T>(val);
  }
};

template<class T>
struct CompressibleReduce<T, OpAddAssign>
{
  template<class T1>
  inline static void evaluate(T &ret, const OpAddAssign &, const T1 &val, 
    int n)
  {
    ret = static_cast<T>(n * val);
  }
};

template<class T>
struct CompressibleReduce<T, OpMultiplyAssign>
{
  template<class T1>
  inline static void evaluate(T &ret, const OpMultiplyAssign &, const T1 &val,
    int n)
  {
    ret = static_cast<T>(val);
    while (--n > 0)
      ret *= static_cast<T>(val);
  }
};


//-----------------------------------------------------------------------------
// The point of this class is to input an expression with the
// 'evaluate' member function and reduce it, optionally taking advantage of
// compression.
//-----------------------------------------------------------------------------

template<>
struct ReductionEvaluator<CompressibleKernelTag>
{
  //---------------------------------------------------------------------------
  // Input an expression and cause it to be reduced.
  // This class relies on another class, CompressibleReduce<T, Op> to
  // perform the correct reduction based on the operator if the expression
  // is compressed. If it is not, we simply use 
  // ReductionEvaluator<InlineKernelTag>.

  template<class T, class Op, class Expr>
  inline static void evaluate(T &ret, const Op &op, const Expr &e)
  {
    if (engineFunctor(e, Compressed()))
      {
        CompressibleReduce<T, Op>::
          evaluate(ret, op, engineFunctor(e, CompressedRead()), 
            e.domain().size());
      }
    else
      {
        ReductionEvaluator<InlineKernelTag>::evaluate(ret, op, e);
      }
  }
};

#endif // POOMA_EVALUATOR_REDUCTIONEVALUATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ReductionEvaluator.h,v $   $Author: richi $
// $Revision: 1.12 $   $Date: 2004/11/04 20:07:26 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
