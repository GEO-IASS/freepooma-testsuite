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
// ScalarCode
// ApplyMultiArgLoc
// EvaluateLocLoop
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Evaluator
 * @brief
 * Undocumented.
 */

#ifndef POOMA_EVALUATOR_SCALARCODE_H
#define POOMA_EVALUATOR_SCALARCODE_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Functions/MultiArg.h"
#include "Evaluator/MultiArgEvaluator.h"
#include "Evaluator/ScalarCodeInfo.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * ApplyMultiArgLoc<MultiArg, Function>
 *
 * ApplyMultiArgLoc is a helper class that permits generic application of
 * a function to a set of arguments with a Loc<>.  The user can provide a
 * functor which takes a set of arguments:
 *
 * <PRE>
 * struct UserFunctor {
 *   operator()(F1 &f1, F2 &f2, F3 &f3, Loc<4> loc);
 * }
 * </PRE>
 *
 * By wrapping the arguments and the functor in this object, we get a
 * common interface that takes integers and can be used by the evaluator
 * in the inner loop.
 *
 * <PRE>
 * ApplyMultiArgLoc<MultiArg3<F1, F2, F3>, UserFunctor>
 *    op(MultiArg3<F1, F2, F3>(f1, f2, f3), function);
 *
 * op(i, j, k, l); 
 * </PRE>
 *
 * translates to UserFunctor::operator()(f1, f2, f3, Loc<4>(i, j, k, l));
 *
 * This helper class is not for external use.  It gets used inside
 * EvaluateLocLoop which generates the loops for the MultiArgKernel.
 *
 * WARNING: This object is intended to have a very short lifetime.  It
 * stores non-const references to the function and multi-arg, so you should
 * not store one of these objects for later use.  EvaluateLocLoop creates
 * an ApplyMultiArgLoc that just lives for the length of the function call.
 */

template<class MA, class Function>
struct ApplyMultiArgLoc;

template<class A1, class Function>
struct ApplyMultiArgLoc<MultiArg1<A1>, Function>
{
  ApplyMultiArgLoc(const MultiArg1<A1> &multiArg, const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg1<A1> &multiArg_m;
  const Function &function_m;
};

template<class A1, class A2, class Function>
struct ApplyMultiArgLoc<MultiArg2<A1, A2>, Function>
{
  ApplyMultiArgLoc(const MultiArg2<A1, A2> &multiArg,const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg2<A1, A2> &multiArg_m;
  const Function &function_m;
};

template<class A1, class A2, class A3, class Function>
struct ApplyMultiArgLoc<MultiArg3<A1, A2, A3>, Function>
{
  ApplyMultiArgLoc(const MultiArg3<A1, A2, A3> &multiArg,
		   const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
	       Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
	       Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
	       Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
	       Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg3<A1, A2, A3> &multiArg_m;
  const Function &function_m;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class Function>
struct ApplyMultiArgLoc<MultiArg7<A1, A2, A3, A4, A5, A6, A7>, Function>
{
  ApplyMultiArgLoc(const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &multiArg,
		   const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               multiArg_m.a7_m, Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               multiArg_m.a7_m, Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               multiArg_m.a7_m, Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               multiArg_m.a7_m, Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &multiArg_m;
  const Function &function_m;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class Function>
struct ApplyMultiArgLoc<MultiArg6<A1, A2, A3, A4, A5, A6>, Function>
{
  ApplyMultiArgLoc(const MultiArg6<A1, A2, A3, A4, A5, A6> &multiArg,
		   const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
                Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m, multiArg_m.a6_m,
               Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg6<A1, A2, A3, A4, A5, A6> &multiArg_m;
  const Function &function_m;
};

template<class A1, class A2, class A3, class A4, class Function>
struct ApplyMultiArgLoc<MultiArg4<A1, A2, A3, A4>, Function>
{
  ApplyMultiArgLoc(const MultiArg4<A1, A2, A3, A4> &multiArg,
		   const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m,
               Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m,
                Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m,
               Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m,
               Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg4<A1, A2, A3, A4> &multiArg_m;
  const Function &function_m;
};

template<class A1, class A2, class A3, class A4, class A5, class Function>
struct ApplyMultiArgLoc<MultiArg5<A1, A2, A3, A4, A5>, Function>
{
  ApplyMultiArgLoc(const MultiArg5<A1, A2, A3, A4, A5> &multiArg,
		   const Function &function)
    : multiArg_m(multiArg), function_m(function)
  {
  }

  void operator()(int i0) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m,
               Loc<1>(i0));
  }

  void operator()(int i0, int i1) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m,
                Loc<2>(i0, i1));
  }

  void operator()(int i0, int i1, int i2) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m,
               Loc<3>(i0, i1, i2));
  }

  void operator()(int i0, int i1, int i2, int i3) const
  {
    function_m(multiArg_m.a1_m, multiArg_m.a2_m, multiArg_m.a3_m,
               multiArg_m.a4_m, multiArg_m.a5_m,
               Loc<4>(i0, i1, i2, i3));
  }

  const MultiArg5<A1, A2, A3, A4, A5> &multiArg_m;
  const Function &function_m;
};

template<class Function, int Dim>
struct EvaluateLocLoop
{
  EvaluateLocLoop()
  {
  }

  EvaluateLocLoop(const Function &function, const Interval<Dim> &domain)
    : function_m(function), domain_m(domain)
  {
  }

  template<class MultiArg>
  void operator()(MultiArg &multiArg) const
  {
    ApplyMultiArgLoc<MultiArg, Function> op(multiArg, function_m);
    LoopApplyEvaluator::evaluate(op, domain_m);
  }

  Function function_m;
  Interval<Dim> domain_m;
};


/**
 * ScalarCode is a Stencil like operation that allows for more than one
 * field to be operated on. Generally the functor is a local (set of)
 * function(s) which could be described as
 *
 *   (f1..fM) = op(fM+1..fN)
 *
 * where fM+1 to fN are input fields read from and f1 to fM are output
 * fields written to (this distinction nor its ordering is strictly
 * required, but both will result in the least possible surprises).
 */

template<class Function>
struct ScalarCode
{
  ScalarCode()
  {
  }

  ScalarCode(const Function &function)
    : function_m(function)
  {
  }

  /// Constructor to allow ScalarCode being used as RelationFunctor

  template <class LHS>
  ScalarCode(const ScalarCode<Function>& sc, const LHS&)
    : function_m(sc.function_m)
  {
  }

  template<class F>
  static inline bool checkValidity(const F& f, WrappedInt<false>)
  {
    // Arrays are ok always.
    return true;
  }

  template<class F>
  static inline bool checkValidity(const F& f, WrappedInt<true>)
  {
    // We cannot handle fields with arbitrary centering size
    // due to physicalDomain() giving potentially wrong answers.
    // We don't handle materials either.
    return f.centeringSize() == 1 && f.numMaterials() == 1;
  }

  /// @name Evaluators
  /// Evaluate the ScalarCode functor on the fields f1 to fN using the
  /// specified evaluation domain. Note that views of the evaluation domain
  /// are taken of every field, so domains of the fields should be strictly
  /// conforming (in fact, passing views to these operators is a bug unless
  /// you really know what you are doing).
  ///
  /// The evaluation domain defaults to the physical domain of
  /// the first field which should usually be (on of) the left hand side(s).
  /// If you want the functor to operate on a different domain use the
  /// operators with the explicit specified evaluation domain.
  //@{

  template<class F1>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg1<F1> multiArg(f1);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1>
  inline void operator()(const F1 &f1) const
  {
    (*this)(f1, f1.physicalDomain());
  }


  template<class F1, class F2>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom,
                  const F2 &f2) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg2<F1, F2> multiArg(f1, f2);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1, class F2>
  inline void operator()(const F1 &f1, const F2 &f2) const
  {
    (*this)(f1, f1.physicalDomain(), f2);
  }


  template<class F1, class F2, class F3>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom,
                  const F2 &f2, const F3 &f3) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg3<F1, F2, F3> multiArg(f1, f2, f3);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1, class F2, class F3>
  inline void operator()(const F1 &f1, const F2 &f2, const F3 &f3) const
  {
    (*this)(f1, f1.physicalDomain(), f2, f3);
  }


  template<class F1, class F2, class F3, class F4>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom,
		  const F2 &f2, const F3 &f3, const F4 &f4) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg4<F1, F2, F3, F4> multiArg(f1, f2, f3, f4);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1, class F2, class F3, class F4>
  inline void operator()(const F1 &f1, const F2 &f2, const F3 &f3, const F4 &f4) const
  {
    (*this)(f1, f1.physicalDomain(), f2, f3, f4);
  }


  template<class F1, class F2, class F3, class F4, class F5>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom,
		  const F2 &f2, const F3 &f3, const F4 &f4, const F5 &f5) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg5<F1, F2, F3, F4, F5> multiArg(f1, f2, f3, f4, f5);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1, class F2, class F3, class F4, class F5>
  inline void operator()(const F1 &f1, const F2 &f2, const F3 &f3, const F4 &f4,
			 const F5 &f5) const
  {
    (*this)(f1, f1.physicalDomain(), f2, f3, f4, f5);
  }


  template<class F1, class F2, class F3, class F4, class F5, class F6>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom,
		  const F2 &f2, const F3 &f3, const F4 &f4, const F5 &f5,
		  const F6 &f6) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg6<F1, F2, F3, F4, F5, F6> multiArg(f1, f2, f3, f4, f5, f6);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1, class F2, class F3, class F4, class F5, class F6>
  inline void operator()(const F1 &f1, const F2 &f2, const F3 &f3, const F4 &f4,
			 const F5 &f5, const F6 &f6) const
  {
    (*this)(f1, f1.physicalDomain(), f2, f3, f4, f5, f6);
  }


  template<class F1, class F2, class F3, class F4, class F5, class F6, class F7>
  void operator()(const F1 &f1, const Interval<F1::dimensions> &evalDom,
		  const F2 &f2, const F3 &f3, const F4 &f4,
		  const F5 &f5, const F6 &f6, const F7 &f7) const
  {
    PAssert(checkValidity(f1, WrappedInt<F1::hasRelations>()));
    MultiArg7<F1, F2, F3, F4, F5, F6, F7> multiArg(f1, f2, f3, f4, f5, f6, f7);
    EvaluateLocLoop<Function, F1::dimensions> kernel(function_m, evalDom);
    MultiArgEvaluator<MainEvaluatorTag>::
      evaluate(multiArg, function_m, evalDom, kernel);
  }

  template<class F1, class F2, class F3, class F4, class F5, class F6, class F7>
  inline void operator()(const F1 &f1, const F2 &f2, const F3 &f3, const F4 &f4,
			 const F5 &f5, const F6 &f6, const F7 &f7) const
  {
    (*this)(f1, f1.physicalDomain(), f2, f3, f4, f5, f6, f7);
  }

  //@}

  Function function_m;
};



//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_SCALARCODE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ScalarCode.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
