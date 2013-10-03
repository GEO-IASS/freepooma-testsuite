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
// MultiArgEvaluator
// MultiArgEvaluatorTag
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_MULTIARGEVALUATOR_H
#define POOMA_EVALUATOR_MULTIARGEVALUATOR_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * MultiArgEvaluator is an evaluator that takes a MultiArg object.
 *
 * Unlike
 * previous evaluators, this evaluator is not limited to a set number of
 * arguments.  All you need to do is write a MultiArg object that packs the
 * number of arguments you need.  This evaluator takes a templated function
 * object, so it is not restricted to the standard expression evaluation over
 * a loop.  Currently this evaluator is used to implement the ScalarCode
 * evaluation, which allows users to write functions that operate on a point
 * in a field and its neighbors and have the function applied to all points
 * in the field.
 *
 * This evaluator currently only works with conforming layouts.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/Intersector.h"
#include "Evaluator/MultiArgKernel.h"
#include "Evaluator/SimpleIntersector.h"
#include "Evaluator/ScalarCodeInfo.h"
#include "Utilities/PerformUpdate.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class MultiArg> struct MultiArgEvaluatorTag;
template<class MeshTag, class T, class EngineTag> class Field;
template<int Dim, class T, class EngineTag> class Array;

/**
 * Implements: MultiArgEvaluator<MainEvaluatorTag>::evaluate
 * (multiArg, function, domain, kernel)
 * - multiArg: a MultiArgX object containing several fields
 * - function: a function that is used to construct the function applied to a
 *            patch (see kernel)
 * - domain: the domain over which the evaluation is performed
 * - kernel: this object is just used to provide a type, which is the actual
 *            function applied to the patch.  That kernel function is
 *            constructed by Kernel_t(function, domain') where function is the
 *            function from the input, and domain' is the domain of the patch.
 */

template<class EvalTag>
struct MultiArgEvaluator
{
};

struct EngineWriteNotifier
{
  EngineWriteNotifier()
  {
  }

  template<class A>
  inline void dirtyRelations(const A &a, const WrappedInt<true>&) const
  {
    a.setDirty();
  }

  template<class A>
  inline void dirtyRelations(const A &, const WrappedInt<false>&) const
  {
  }

  template<class A>
  void operator()(const A &a) const
  {
    // This isn't quite what we want here, because we may want to
    // write to a field containing multiple centering engines.
    // Need to rewrite notifyEngineWrite as an ExpressionApply,
    // and create a version of ExpressionApply that goes through
    // all the engines in a field.

    notifyEngineWrite(a.engine());
    dirtyRelations(a, WrappedInt<A::hasRelations>());
  }

  // overload for ExpressionTag engines to not fall on our faces compile time
  template<class MeshTag, class T, class Expr>
  void operator()(const Field<MeshTag, T, ExpressionTag<Expr> >&) const
  {
    // we must be able to compile this, but never execute
    PInsist(false, "writing to expression engine?");
  }
  template<int Dim, class T, class Expr>
  void operator()(const Array<Dim, T, ExpressionTag<Expr> >&) const
  {
    // we must be able to compile this, but never execute
    PInsist(false, "writing to expression engine?");
  }
};

struct UpdateNotifier
{
  UpdateNotifier()
  {
  }

  template<class A>
  void operator()(const A &a) const
  {
    forEach(a, PerformUpdateTag(), NullCombine());
  }
};

template<>
struct MultiArgEvaluator<MainEvaluatorTag>
{
public:

  // Default ctor.

  MultiArgEvaluator() {}

  // Destructor

  ~MultiArgEvaluator() {}

  template<class MultiArg, class Function, int Dim, class Kernel>
  static void
  evaluate(const MultiArg &multiArg,
	   const Function &function,
	   const Interval<Dim> &domain,
	   const Kernel &kernel)
  {
    typedef typename MultiArgEvaluatorTag<MultiArg>::Evaluator_t Evaluator_t;

    ScalarCodeInfo info;
    function.scalarCodeInfo(info);

    Pooma::beginExpression();

    applyMultiArg(multiArg, UpdateNotifier());

    MultiArgEvaluator<Evaluator_t>::evaluate(multiArg, function,
					     domain, info, kernel);

    applyMultiArgIf(multiArg, EngineWriteNotifier(), info.writers());

    Pooma::endExpression();
  }

  template<class A1, class Function, int Dim, class Kernel>
  static void
  createIterate(const A1& a1,
		const Function& function,
		const Interval<Dim> &domain,
		ScalarCodeInfo &info,
		const Kernel &)
  {
    Kernel kernelf(function, domain);

    Pooma::Iterate_t *iterate =
      new MultiArgKernel<A1, Kernel>(a1, kernelf,
				     info.writers(), info.readers());
    Pooma::scheduler().handOff(iterate);
  }

};

// The single patch version just passes the tag on to generate
// an expression kernel.

template<>
struct MultiArgEvaluator<SinglePatchEvaluatorTag>
{
public:

  //
  // Default ctor.
  // The only member data can construct itself, so we
  // don't need to specify anything.
  //
  MultiArgEvaluator() {}

  //
  // Destructor
  //
  ~MultiArgEvaluator() {}

  template<class MultiArg, class Function, int Dim, class Kernel>
  static void
  evaluate(const MultiArg& multiArg,
	   const Function& function,
	   Interval<Dim> domain,
	   ScalarCodeInfo &info,
	   const Kernel &kernel)
  {
    Interval<Dim> newDom = info.extendDomain(domain);
    Interval<Dim> evalDom = info.evaluationDomain(domain);
    MultiArgEvaluator<MainEvaluatorTag>::createIterate(multiArg(newDom),
						       function,
						       evalDom, info,
						       kernel);

  }

};



// The multiple patch version makes patches and sends them out to
// the single patch evaluator.

template<>
struct MultiArgEvaluator<MultiPatchEvaluatorTag>
{
public:

  //
  // Default ctor.
  // The only member data can construct itself, so we
  // don't need to specify anything.
  //
  MultiArgEvaluator() {}

  //
  // Destructor
  //
  ~MultiArgEvaluator() {}

  template<class MultiArg, class Function, int Dim, class Kernel>
  static void
  evaluate(const MultiArg &multiArg,
	   const Function &function,
	   const Interval<Dim> &domain,
	   ScalarCodeInfo &info,
	   const Kernel &kernel)
  {
    typedef SimpleIntersector<Dim> Inter_t;
    GuardLayers<Dim> extent;
    for (int i=0; i<Dim; ++i) {
      extent.lower(i) = info.lowerExtent(i);
      extent.upper(i) = info.upperExtent(i);
    }
    Inter_t inter(domain, extent);

    applyMultiArg(multiArg, inter, info.useGuards());
 
    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
    {
      INode<Dim> inode = info.extendDomain(*i);
      Interval<Dim> evalDom = info.evaluationDomain((*i).domain());

      MultiArgEvaluator<MainEvaluatorTag>::createIterate(multiArg(inode),
							 function,
							 evalDom, info,
							 kernel);
      ++i;
    }
  }

};

//-----------------------------------------------------------------------------
// Single-patch Evaluator involving remote engines:
//
// This evaluator handles a single patch involving engines that may be remote.
//-----------------------------------------------------------------------------

template <>
struct MultiArgEvaluator<RemoteSinglePatchEvaluatorTag>
{
  // Default ctor.

  MultiArgEvaluator() { }

  // Destructor.

  ~MultiArgEvaluator() { }

  // evaluate(expression)
  // Input an expression and cause it to be evaluated.
  // We just pass the buck to a special evaluator.

  template<class MultiArg, class Function, int Dim, class Kernel>
  static void
  evaluate(const MultiArg &multiArg,
	   const Function &function,
	   const Interval<Dim> &domain,
	   ScalarCodeInfo &info,
	   const Kernel &kernel)
  {
    // This code is still untested.  Field doesn't
    // support remote engines yet.

    GatherContexts gtag;
    engineFunctor(multiArg.a1_m.engine(), gtag);
    int lhsContext = gtag.mostCommonContext();

    expressionApply(multiArg, RemoteSend(lhsContext));

    EngineView<RemoteView> view;

    if (lhsContext == -1 || Pooma::context() == lhsContext)
    {
      MultiArgEvaluator<SinglePatchEvaluatorTag> speval;
      speval.evaluate(
		      forEach(multiArg, view, TreeCombine()),
		      function, domain, info, kernel
		      );
    }
  }
};


//-----------------------------------------------------------------------------
// Multiple-patch MultiArgEvaluator involving remote engines:
//
// The remote multiple patch version makes patches and sends them out to
// the remote single patch evaluator.
//-----------------------------------------------------------------------------

template <>
struct MultiArgEvaluator<RemoteMultiPatchEvaluatorTag>
{
  // Default ctor.

  MultiArgEvaluator() { }

  // Destructor.

  ~MultiArgEvaluator() { }

  // evaluate(expression)
  // Input an expression and cause it to be evaluated.
  // We just pass the buck to a special evaluator.

  template<class MultiArg, class Function, int Dim, class Kernel>
  static void
  evaluate(const MultiArg &multiArg,
	   const Function &function,
	   const Interval<Dim> &domain,
	   ScalarCodeInfo &info,
	   const Kernel &kernel)
  {
    typedef SimpleIntersector<Dim> Inter_t;
    GuardLayers<Dim> extent;
    for (int i=0; i<Dim; ++i) {
      extent.lower(i) = info.lowerExtent(i);
      extent.upper(i) = info.upperExtent(i);
    }
    Inter_t inter(domain, extent);

    applyMultiArg(multiArg, inter, info.useGuards());
 
    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
    {
      INode<Dim> inode = info.extendDomain(*i);
      Interval<Dim> evalDom = info.evaluationDomain((*i).domain());

      MultiArgEvaluator<RemoteSinglePatchEvaluatorTag>().
	evaluate(multiArg(inode),
		 function, evalDom, info, kernel);
      ++i;
    }
  }
};

template<class A1>
struct MultiArgEvaluatorTag<MultiArg1<A1> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Evaluator_t;
};

template<class A1, class A2>
struct MultiArgEvaluatorTag<MultiArg2<A1, A2> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
  typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
  typedef typename EvaluatorCombine<Eval1_t, Eval2_t>::Evaluator_t Evaluator_t;
};

template<class A1, class A2, class A3>
struct MultiArgEvaluatorTag<MultiArg3<A1, A2, A3> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
  typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
  typedef typename EvaluatorTag1<A3>::Evaluator_t Eval3_t;
  typedef typename EvaluatorCombine<Eval1_t, Eval2_t>::Evaluator_t Eval12_t;
  typedef typename EvaluatorCombine<Eval3_t, Eval12_t>::Evaluator_t
  Evaluator_t;
};

template<class A1, class A2, class A3, class A4>
struct MultiArgEvaluatorTag<MultiArg4<A1, A2, A3, A4> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
  typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
  typedef typename EvaluatorTag1<A3>::Evaluator_t Eval3_t;
  typedef typename EvaluatorTag1<A4>::Evaluator_t Eval4_t;
  typedef typename EvaluatorCombine<Eval1_t, Eval2_t>::Evaluator_t Eval12_t;
  typedef typename EvaluatorCombine<Eval3_t, Eval4_t>::Evaluator_t  Eval34_t;
  typedef typename EvaluatorCombine<Eval12_t, Eval34_t>::Evaluator_t  Evaluator_t;
};

template<class A1, class A2, class A3, class A4, class A5>
struct MultiArgEvaluatorTag<MultiArg5<A1, A2, A3, A4, A5> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
  typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
  typedef typename EvaluatorTag1<A3>::Evaluator_t Eval3_t;
  typedef typename EvaluatorTag1<A4>::Evaluator_t Eval4_t;
  typedef typename EvaluatorTag1<A5>::Evaluator_t Eval5_t;
  typedef typename EvaluatorCombine<Eval1_t, Eval2_t>::Evaluator_t Eval12_t;
  typedef typename EvaluatorCombine<Eval3_t, Eval4_t>::Evaluator_t  Eval34_t;
  typedef typename EvaluatorCombine<Eval12_t, Eval34_t>::Evaluator_t  Eval1234_t;
  typedef typename EvaluatorCombine<Eval1234_t, Eval5_t>::Evaluator_t  Evaluator_t;
};

template<class A1, class A2, class A3, class A4, class A5, class A6>
struct MultiArgEvaluatorTag<MultiArg6<A1, A2, A3, A4, A5, A6> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
  typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
  typedef typename EvaluatorTag1<A3>::Evaluator_t Eval3_t;
  typedef typename EvaluatorTag1<A4>::Evaluator_t Eval4_t;
  typedef typename EvaluatorTag1<A5>::Evaluator_t Eval5_t;
  typedef typename EvaluatorTag1<A6>::Evaluator_t Eval6_t;
  typedef typename EvaluatorCombine<Eval1_t, Eval2_t>::Evaluator_t Eval12_t;
  typedef typename EvaluatorCombine<Eval3_t, Eval4_t>::Evaluator_t  Eval34_t;
  typedef typename EvaluatorCombine<Eval5_t, Eval6_t>::Evaluator_t  Eval56_t;
  typedef typename EvaluatorCombine<Eval12_t, Eval34_t>::Evaluator_t  Eval1234_t;
  typedef typename EvaluatorCombine<Eval1234_t, Eval56_t>::Evaluator_t  Evaluator_t;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7>
struct MultiArgEvaluatorTag<MultiArg7<A1, A2, A3, A4, A5, A6, A7> >
{
  typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
  typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
  typedef typename EvaluatorTag1<A3>::Evaluator_t Eval3_t;
  typedef typename EvaluatorTag1<A4>::Evaluator_t Eval4_t;
  typedef typename EvaluatorTag1<A5>::Evaluator_t Eval5_t;
  typedef typename EvaluatorTag1<A6>::Evaluator_t Eval6_t;
  typedef typename EvaluatorTag1<A7>::Evaluator_t Eval7_t;
  typedef typename EvaluatorCombine<Eval1_t, Eval2_t>::Evaluator_t Eval12_t;
  typedef typename EvaluatorCombine<Eval3_t, Eval12_t>::Evaluator_t  Eval123_t;
  typedef typename EvaluatorCombine<Eval4_t, Eval123_t>::Evaluator_t  Eval1234_t;
  typedef typename EvaluatorCombine<Eval5_t, Eval1234_t>::Evaluator_t  Eval12345_t;
  typedef typename EvaluatorCombine<Eval6_t, Eval12345_t>::Evaluator_t  Eval123456_t;
  typedef typename EvaluatorCombine<Eval7_t, Eval123456_t>::Evaluator_t
  Evaluator_t;
};



//----------------------------------------------------------------------------
// LeafFunctor Specializations for ExpressionApply and EngineView
//


template<class A1, class Tag>
struct LeafFunctor<MultiArg1<A1>, ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg1<A1> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    return 0;
  }
};

template<class A1, class Tag>
struct LeafFunctor<MultiArg1<A1>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef MultiArg1<Type1_t> Type_t;

  inline static
  Type_t apply(const MultiArg1<A1> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag)
		  );
  }
};

template<class A1, class A2, class Tag>
struct LeafFunctor<MultiArg2<A1, A2>, ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg2<A1, A2> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    leafFunctor(multiarg.a2_m, tag);
    return 0;
  }
};

template<class A1, class A2, class Tag>
struct LeafFunctor<MultiArg2<A1, A2>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef typename LeafFunctor<A2, EngineView<Tag> >::Type_t Type2_t;
  typedef MultiArg2<Type1_t, Type2_t> Type_t;

  inline static
  Type_t apply(const MultiArg2<A1, A2> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag),
		  leafFunctor(multiarg.a2_m, tag)
		  );
  }
};

template<class A1, class A2, class A3, class Tag>
struct LeafFunctor<MultiArg3<A1, A2, A3>,
  ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg3<A1, A2, A3> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    leafFunctor(multiarg.a2_m, tag);
    leafFunctor(multiarg.a3_m, tag);
    return 0;
  }
};

template<class A1, class A2, class A3, class Tag>
struct LeafFunctor<MultiArg3<A1, A2, A3>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef typename LeafFunctor<A2, EngineView<Tag> >::Type_t Type2_t;
  typedef typename LeafFunctor<A3, EngineView<Tag> >::Type_t Type3_t;
  typedef MultiArg3<Type1_t, Type2_t, Type3_t> Type_t;

  inline static
  Type_t apply(const MultiArg3<A1, A2, A3> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag),
		  leafFunctor(multiarg.a2_m, tag),
		  leafFunctor(multiarg.a3_m, tag)
		  );
  }
};

template<class A1, class A2, class A3, class A4, class Tag>
struct LeafFunctor<MultiArg4<A1, A2, A3, A4>,
  ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg4<A1, A2, A3, A4> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    leafFunctor(multiarg.a2_m, tag);
    leafFunctor(multiarg.a3_m, tag);
    leafFunctor(multiarg.a4_m, tag);
    return 0;
  }
};

template<class A1, class A2, class A3, class A4, class Tag>
struct LeafFunctor<MultiArg4<A1, A2, A3, A4>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef typename LeafFunctor<A2, EngineView<Tag> >::Type_t Type2_t;
  typedef typename LeafFunctor<A3, EngineView<Tag> >::Type_t Type3_t;
  typedef typename LeafFunctor<A4, EngineView<Tag> >::Type_t Type4_t;
  typedef MultiArg4<Type1_t, Type2_t, Type3_t, Type4_t> Type_t;

  inline static
  Type_t apply(const MultiArg4<A1, A2, A3, A4> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag),
		  leafFunctor(multiarg.a2_m, tag),
		  leafFunctor(multiarg.a3_m, tag),
		  leafFunctor(multiarg.a4_m, tag)
		  );
  }
};

template<class A1, class A2, class A3, class A4, class A5, class Tag>
struct LeafFunctor<MultiArg5<A1, A2, A3, A4, A5>,
  ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg5<A1, A2, A3, A4, A5> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    leafFunctor(multiarg.a2_m, tag);
    leafFunctor(multiarg.a3_m, tag);
    leafFunctor(multiarg.a4_m, tag);
    leafFunctor(multiarg.a5_m, tag);
    return 0;
  }
};

template<class A1, class A2, class A3, class A4, class A5, class Tag>
struct LeafFunctor<MultiArg5<A1, A2, A3, A4, A5>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef typename LeafFunctor<A2, EngineView<Tag> >::Type_t Type2_t;
  typedef typename LeafFunctor<A3, EngineView<Tag> >::Type_t Type3_t;
  typedef typename LeafFunctor<A4, EngineView<Tag> >::Type_t Type4_t;
  typedef typename LeafFunctor<A5, EngineView<Tag> >::Type_t Type5_t;
  typedef MultiArg5<Type1_t, Type2_t, Type3_t, Type4_t, Type5_t> Type_t;

  inline static
  Type_t apply(const MultiArg5<A1, A2, A3, A4, A5> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag),
		  leafFunctor(multiarg.a2_m, tag),
		  leafFunctor(multiarg.a3_m, tag),
		  leafFunctor(multiarg.a4_m, tag),
		  leafFunctor(multiarg.a5_m, tag)
		  );
  }
};

template<class A1, class A2, class A3, class A4, class A5,
  class A6, class Tag>
struct LeafFunctor<MultiArg6<A1, A2, A3, A4, A5, A6>,
  ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg6<A1, A2, A3, A4, A5, A6> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    leafFunctor(multiarg.a2_m, tag);
    leafFunctor(multiarg.a3_m, tag);
    leafFunctor(multiarg.a4_m, tag);
    leafFunctor(multiarg.a5_m, tag);
    leafFunctor(multiarg.a6_m, tag);
    return 0;
  }
};

template<class A1, class A2, class A3, class A4, class A5,
  class A6, class Tag>
struct LeafFunctor<MultiArg6<A1, A2, A3, A4, A5, A6>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef typename LeafFunctor<A2, EngineView<Tag> >::Type_t Type2_t;
  typedef typename LeafFunctor<A3, EngineView<Tag> >::Type_t Type3_t;
  typedef typename LeafFunctor<A4, EngineView<Tag> >::Type_t Type4_t;
  typedef typename LeafFunctor<A5, EngineView<Tag> >::Type_t Type5_t;
  typedef typename LeafFunctor<A6, EngineView<Tag> >::Type_t Type6_t;
  typedef MultiArg6<Type1_t, Type2_t, Type3_t, Type4_t, Type5_t,
    Type6_t> Type_t;

  inline static
  Type_t apply(const MultiArg6<A1, A2, A3, A4, A5, A6> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag),
		  leafFunctor(multiarg.a2_m, tag),
		  leafFunctor(multiarg.a3_m, tag),
		  leafFunctor(multiarg.a4_m, tag),
		  leafFunctor(multiarg.a5_m, tag),
		  leafFunctor(multiarg.a6_m, tag)
		  );
  }
};

template<class A1, class A2, class A3, class A4, class A5,
  class A6, class A7, class Tag>
struct LeafFunctor<MultiArg7<A1, A2, A3, A4, A5, A6, A7>,
  ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &multiarg,
	       const ExpressionApply<Tag> &tag)
  {
    leafFunctor(multiarg.a1_m, tag);
    leafFunctor(multiarg.a2_m, tag);
    leafFunctor(multiarg.a3_m, tag);
    leafFunctor(multiarg.a4_m, tag);
    leafFunctor(multiarg.a5_m, tag);
    leafFunctor(multiarg.a6_m, tag);
    leafFunctor(multiarg.a7_m, tag);
    return 0;
  }
};

template<class A1, class A2, class A3, class A4, class A5,
  class A6, class A7, class Tag>
struct LeafFunctor<MultiArg7<A1, A2, A3, A4, A5, A6, A7>, EngineView<Tag> >
{
  typedef typename LeafFunctor<A1, EngineView<Tag> >::Type_t Type1_t;
  typedef typename LeafFunctor<A2, EngineView<Tag> >::Type_t Type2_t;
  typedef typename LeafFunctor<A3, EngineView<Tag> >::Type_t Type3_t;
  typedef typename LeafFunctor<A4, EngineView<Tag> >::Type_t Type4_t;
  typedef typename LeafFunctor<A5, EngineView<Tag> >::Type_t Type5_t;
  typedef typename LeafFunctor<A6, EngineView<Tag> >::Type_t Type6_t;
  typedef typename LeafFunctor<A7, EngineView<Tag> >::Type_t Type7_t;
  typedef MultiArg7<Type1_t, Type2_t, Type3_t, Type4_t, Type5_t,
    Type6_t, Type7_t> Type_t;

  inline static
  Type_t apply(const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &multiarg,
	       const EngineView<Tag> &tag)
  {
    return Type_t(
		  leafFunctor(multiarg.a1_m, tag),
		  leafFunctor(multiarg.a2_m, tag),
		  leafFunctor(multiarg.a3_m, tag),
		  leafFunctor(multiarg.a4_m, tag),
		  leafFunctor(multiarg.a5_m, tag),
		  leafFunctor(multiarg.a6_m, tag),
		  leafFunctor(multiarg.a7_m, tag)
		  );
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_MULTIARGEVALUATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MultiArgEvaluator.h,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
