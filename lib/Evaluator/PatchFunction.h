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
// PatchFunction
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_PATCHFUNCTION_H
#define POOMA_EVALUATOR_PATCHFUNCTION_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * PatchFunction is mix-in class that encapsulates evaluation of patch-based
 * functors in parallel.
 *
 * PatchFunctions are tools that allow you to apply a functor to the patches
 * in an array in parallel.  For example, you could write a functor:
 *
 * <PRE>
 * struct MyFunction
 * {
 *   template<class ArrayPatch>
 *   void apply(const ArrayPatch& a) const
 *   {
 *     for (i=0;i<a.domain().size();++i)
 *     {
 *       a(i) += 2;
 *     }
 *   }
 * };
 * </PRE>
 *
 * and apply it to an array with the PatchFunction:
 *
 * <PRE>
 * PatchFunction<MyFunction,PatchTag1> myFunc;
 * myFunc(array);
 * </PRE>
 *
 * Iterates will be spawned for each patch in array, and MyFuntion::apply()
 * will be called for each patch.  The general form of PatchFunction is
 * PatchFunction<Functor,Tag> func(c1,c2...);
 * Constructor arguments are passed to the constructor of a functor object
 * of type Functor.  Tag is a policy tag specifying the type of action
 * performed.  Currently the following tags are supported:
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "Pooma/PETE/AssertEquals.h"
#include "Evaluator/EvaluatorTags.h"
#include "Evaluator/Evaluator.h"
#include "Evaluator/PatchKernel.h"
#include "Engine/EnginePatch.h"
#include "Threads/PoomaCSem.h"


//-----------------------------------------------------------------------------
// Policy tags for patch functions.
//
// PatchParticleN - recommended for operations on particles, bypasses the
//                  intersection process and just loops through the patches.
//-----------------------------------------------------------------------------

struct PatchTag1 { };
struct PatchReadTag1 { };
struct PatchTag2 { };
struct PatchTag3 { };

template<bool Write1>
struct PatchParticle1 { };

template<bool Write1, bool Write2>
struct PatchParticle2 { };

template<bool Write1, bool Write2, bool Write3>
struct PatchParticle3 { };

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<class EvalTag>
class PatchEvaluator
{
};

template<>
class PatchEvaluator<MainEvaluatorTag>
{
public:

  // Default ctor.

  PatchEvaluator() {}

  // Destructor

  ~PatchEvaluator() {}

  template<class A1, class Function>
  void evaluate(const A1& a1, const Function& function) const
  {
    typedef typename EvaluatorTag1<A1>::Evaluator_t Evaluator_t;
    PatchEvaluator<Evaluator_t> evaluator;
    Pooma::beginExpression();
    evaluator.evaluate(a1(), function);
    notifyEngineWrite(a1.engine());
    Pooma::endExpression();
  }

  template<class A1, class Function>
  void evaluateRead(const A1& a1, const Function& function) const
  {
    typedef typename EvaluatorTag1<A1>::Evaluator_t Evaluator_t;
    PatchEvaluator<Evaluator_t> evaluator;
    Pooma::beginExpression();
    evaluator.evaluateRead(a1(), function);
    Pooma::endExpression();
  }

  template<class A1,class A2,class Function>
  void evaluate2(const A1& a1, const A2& a2,
		const Function& function) const
  {
    typedef typename EvaluatorTag<A1,A2>::Evaluator_t Eval_t;
    PatchEvaluator<Eval_t> evaluator;
    Pooma::beginExpression();
    evaluator.evaluate2(a1(), a2(), function);
    notifyEngineWrite(a1.engine());
    Pooma::endExpression();
  }

  template<class A1, class A2, class A3, class Function>
  void evaluate3(const A1& a1, const A2& a2, const A3 &a3,
		const Function& function) const
  {
    typedef typename EvaluatorTag1<A2>::Evaluator_t Eval2_t;
    typedef typename EvaluatorTag1<A3>::Evaluator_t Eval3_t;
    typedef typename EvaluatorCombine<Eval2_t,Eval3_t>::Evaluator_t Eval23_t;
    typedef typename EvaluatorTag1<A1>::Evaluator_t Eval1_t;
    typedef typename EvaluatorCombine<Eval1_t,Eval23_t>::Evaluator_t Eval_t;

    PatchEvaluator<Eval_t> evaluator;
    Pooma::beginExpression();
    evaluator.evaluate3(a1(), a2(), a3(), function);
    notifyEngineWrite(a1.engine());
    Pooma::endExpression();
  }

private:
};

// The single patch version just passes the tag on to generate
// an expression kernel.

template<>
class PatchEvaluator<SinglePatchEvaluatorTag>
{
public:

  //
  // Default ctor.
  // The only member data can construct itself, so we
  // don't need to specify anything.
  //
  PatchEvaluator() {}

  //
  // Destructor
  //
  ~PatchEvaluator() {}

  template<class A1, class Function>
  void evaluate(const A1& a1, const Function& function) const
  {
    Pooma::Iterate_t *iterate = new PatchKernel<A1,Function>(a1,function);
    Pooma::scheduler().handOff(iterate);
  }

  template<class A1, class Function>
  void evaluateRead(const A1& a1, const Function& function) const
  {
    Pooma::Iterate_t *iterate = new PatchKernel<A1,Function>(a1,function);
    Pooma::scheduler().handOff(iterate);
  }

  template<class A1,class A2,class Function>
  void evaluate2(const A1 &a1, const A2 &a2,
		const Function &function) const
  {
    Pooma::Iterate_t *iterate =
      new PatchKernel2<A1,A2,Function>(a1,a2,function);

    Pooma::scheduler().handOff(iterate);
  }

  template<class A1, class A2, class A3, class Function>
  void evaluate3(const A1 &a1, const A2 &a2, const A3 &a3,
		const Function &function) const
  {
    Pooma::Iterate_t *iterate =
      new PatchKernel3<A1,A2,A3,Function>(a1,a2,a3,function);

    Pooma::scheduler().handOff(iterate);
  }

private:
};



// The multiple patch version makes patches and sends them out to
// the single patch evaluator.

template<>
class PatchEvaluator<MultiPatchEvaluatorTag>
{
public:

  //
  // Default ctor.
  // The only member data can construct itself, so we
  // don't need to specify anything.
  //
  PatchEvaluator() {}

  //
  // Destructor
  //
  ~PatchEvaluator() {}

  template<class A1,class Function>
  void evaluate(const A1& a1,const Function& function) const
  {
    typedef Intersector<A1::dimensions> Inter_t;
    Inter_t inter;

    expressionApply(a1, IntersectorTag<Inter_t>(inter));

    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
    {
      PatchEvaluator<SinglePatchEvaluatorTag>().evaluate(a1(*i),function);
      ++i;
    }
  }

  template<class A1,class Function>
  inline
  void evaluateRead(const A1& a1,const Function& function) const
  {
    evaluate(a1,function);
  }

  template<class A1,class A2,class Function>
  void evaluate2(const A1& a1, const A2& a2,
		const Function& function) const
  {
    typedef Intersector<A1::dimensions> Inter_t;
    Inter_t inter;
    
    expressionApply(a1, IntersectorTag<Inter_t>(inter));
    expressionApply(a2, IntersectorTag<Inter_t>(inter));
  
    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
    {
      PatchEvaluator<SinglePatchEvaluatorTag>().evaluate2(a1(*i), a2(*i),
							  function
							  );
      ++i;
    }
  }

  template<class A1, class A2, class A3, class Function>
  void evaluate3(const A1 &a1, const A2 &a2, const A3 &a3,
		 const Function& function) const
  {
    typedef Intersector<A1::dimensions> Inter_t;
    Inter_t inter;
    
    expressionApply(a1, IntersectorTag<Inter_t>(inter));
    expressionApply(a2, IntersectorTag<Inter_t>(inter));
    expressionApply(a3, IntersectorTag<Inter_t>(inter));
  
    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
    {
      PatchEvaluator<SinglePatchEvaluatorTag>().
        evaluate3(a1(*i), a2(*i), a3(*i), function);
      ++i;
    }
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class ParticleEvaluator
{
public:
  // Default ctor.

  ParticleEvaluator() {}

  // Destructor

  ~ParticleEvaluator() {}

  template<class A1, class Function, bool Write1>
  void evaluate(const A1 &a1, const Function &function,
		const PatchParticle1<Write1> &) const
  {
    Pooma::Scheduler_t &scheduler = Pooma::scheduler();
    Pooma::beginExpression();

    int n = a1.numPatchesLocal();
    int i;

    typedef typename PatchView<A1>::Type_t PatchA1_t;

    for (i = 0; i < n; ++i)
    {
      Pooma::Iterate_t *iterate =
	new ParticleKernel<PatchA1_t,
	Function>(a1.patchLocal(i)(), function, i, Write1);
      scheduler.handOff(iterate);
    }

    notifyEngineWrite(a1.engine(), WrappedInt<Write1>());

    Pooma::endExpression();
  }

  template<class A1, class Function, bool Write1>
  void evaluateBlock(const A1 &a1, const Function &function,
		     const PatchParticle1<Write1> &) const
  {
    Pooma::Scheduler_t &scheduler = Pooma::scheduler();

    int n = a1.numPatchesLocal();
    int i;

    Pooma::CountingSemaphore *csem = new Pooma::CountingSemaphore;
    csem->height(n);

    scheduler.beginGeneration();

    typedef typename PatchView<A1>::Type_t PatchA1_t;

    for (i = 0; i < n; ++i)
    {
      Pooma::Iterate_t *iterate =
        new ParticleKernelBlock<PatchA1_t,Function>(a1.patchLocal(i)(),
        function, i, Write1, csem);
      scheduler.handOff(iterate);
    }

    notifyEngineWrite(a1.engine(), WrappedInt<Write1>());

    scheduler.endGeneration();

    csem->wait();
    delete csem;

    if (Pooma::blockingExpressions()) 
    {
      Pooma::blockAndEvaluate();
    }
  }

  struct NoOp { };

  template<class A1, class A2, class Function, bool Write1, bool Write2>
  void evaluate2(const A1& a1, const A2 &a2, const Function& function,
		 const PatchParticle2<Write1, Write2> &) const
  {
    Pooma::Scheduler_t &scheduler = Pooma::scheduler();
    Pooma::beginExpression();

    int n1 = a1.numPatchesLocal();
    int n2 = a2.numPatchesLocal();
    int i;

    int n = peteCombine(n1, n2, NoOp(), AssertEquals());

    typedef typename PatchView<A1>::Type_t PatchA1_t;
    typedef typename PatchView<A2>::Type_t PatchA2_t;

    for (i = 0; i < n; ++i)
    {
      Pooma::Iterate_t *iterate =
        new ParticleKernel2<PatchA1_t,PatchA2_t,Function>(a1.patchLocal(i)(),
        a2.patchLocal(i)(), function, i, Write1, Write2);
      scheduler.handOff(iterate);
    }

    notifyEngineWrite(a1.engine(), WrappedInt<Write1>());
    notifyEngineWrite(a2.engine(), WrappedInt<Write2>());

    Pooma::endExpression();
  }

  template<class A1, class A2, class Function, bool Write1, bool Write2>
  void evaluate2Block(const A1& a1, const A2 &a2, const Function& function,
		 const PatchParticle2<Write1, Write2> &) const
  {
    Pooma::Scheduler_t &scheduler = Pooma::scheduler();
    scheduler.beginGeneration();

    int n1 = a1.numPatchesLocal();
    int n2 = a2.numPatchesLocal();
    int i;

    int n = peteCombine(n1, n2, NoOp(), AssertEquals());

    Pooma::CountingSemaphore *csem = new Pooma::CountingSemaphore;
    csem->height(n);

    typedef typename PatchView<A1>::Type_t PatchA1_t;
    typedef typename PatchView<A2>::Type_t PatchA2_t;

    for (i = 0; i < n; ++i)
    {
      Pooma::Iterate_t *iterate =
        new ParticleKernel2Block<PatchA1_t,PatchA2_t,Function>(
        a1.patchLocal(i)(), a2.patchLocal(i)(), function, i, Write1, Write2, 
        csem);
      scheduler.handOff(iterate);
    }

    notifyEngineWrite(a1.engine(), WrappedInt<Write1>());
    notifyEngineWrite(a2.engine(), WrappedInt<Write2>());

    scheduler.endGeneration();

    csem->wait();
    delete csem;

    if (Pooma::blockingExpressions()) 
    {
      Pooma::blockAndEvaluate();
    }
  }

  template<class A1, class A2, class A3, class Function,
    bool Write1, bool Write2, bool Write3>
  void evaluate3(const A1& a1, const A2& a2, const A3 &a3,
		 const Function& function,
		 const PatchParticle3<Write1, Write2, Write3> &) const
  {
    Pooma::Scheduler_t &scheduler = Pooma::scheduler();
    Pooma::beginExpression();

    int n1 = a1.numPatchesLocal();
    int n2 = a2.numPatchesLocal();
    int n3 = a3.numPatchesLocal();
    int i;

    int n = peteCombine(n1, n2, n3, NoOp(), AssertEquals());

    typedef typename PatchView<A1>::Type_t PatchA1_t;
    typedef typename PatchView<A2>::Type_t PatchA2_t;
    typedef typename PatchView<A3>::Type_t PatchA3_t;

    for (i = 0; i < n; ++i)
    {
      Pooma::Iterate_t *iterate =
        new ParticleKernel3<PatchA1_t,PatchA2_t,PatchA3_t,Function>(
        a1.patchLocal(i)(), a2.patchLocal(i)(), a3.patchLocal(i)(), function, 
        i, Write1, Write2, Write3);
      scheduler.handOff(iterate);
    }

    notifyEngineWrite(a1.engine(), WrappedInt<Write1>());
    notifyEngineWrite(a2.engine(), WrappedInt<Write2>());
    notifyEngineWrite(a3.engine(), WrappedInt<Write3>());

    Pooma::endExpression();
  }

  template<class A1, class A2, class A3, class Function,
    bool Write1, bool Write2, bool Write3>
  void evaluate3Block(const A1& a1, const A2& a2, const A3 &a3,
		 const Function& function,
		 const PatchParticle3<Write1, Write2, Write3> &) const
  {
    Pooma::Scheduler_t &scheduler = Pooma::scheduler();
    scheduler.beginGeneration();

    int n1 = a1.numPatchesLocal();
    int n2 = a2.numPatchesLocal();
    int n3 = a3.numPatchesLocal();
    int i;

    int n = peteCombine(n1, n2, n3, NoOp(), AssertEquals());

    Pooma::CountingSemaphore *csem = new Pooma::CountingSemaphore;
    csem->height(n);

    typedef typename PatchView<A1>::Type_t PatchA1_t;
    typedef typename PatchView<A2>::Type_t PatchA2_t;
    typedef typename PatchView<A3>::Type_t PatchA3_t;

    for (i = 0; i < n; ++i)
    {
      Pooma::Iterate_t *iterate =
        new ParticleKernel3Block<PatchA1_t,PatchA2_t,PatchA3_t,Function>(
        a1.patchLocal(i)(), a2.patchLocal(i)(), a3.patchLocal(i)(), function, 
        i, Write1, Write2, Write3, csem);
      scheduler.handOff(iterate);
    }

    notifyEngineWrite(a1.engine(), WrappedInt<Write1>());
    notifyEngineWrite(a2.engine(), WrappedInt<Write2>());
    notifyEngineWrite(a3.engine(), WrappedInt<Write3>());

    scheduler.endGeneration();

    csem->wait();
    delete csem;

    if (Pooma::blockingExpressions()) 
    {
      Pooma::blockAndEvaluate();
    }
  }

private:
};


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

#define POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(CLASS,MEMBER)            \
template <class I1>                                                        \
CLASS(const I1 &i1)                                                        \
  : MEMBER(i1)                                                             \
{ }                                                                        \
                                                                           \
template <class I1, class I2>                                              \
CLASS(const I1 &i1, const I2 &i2)                                          \
  : MEMBER(i1,i2)                                                          \
{ }                                                                        \
                                                                           \
template <class I1, class I2, class I3>                                    \
CLASS(const I1 &i1, const I2 &i2, const I3 &i3)                            \
  : MEMBER(i1,i2,i3)                                                       \
{ }                                                                        \
                                                                           \
template <class I1, class I2, class I3, class I4>                          \
CLASS(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4)              \
  : MEMBER(i1,i2,i3,i4)                                                    \
{ }                                                                        \
                                                                           \
template <class I1, class I2, class I3, class I4, class I5>                \
CLASS(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4,              \
      const I5 &i5)                                                        \
  : MEMBER(i1,i2,i3,i4,i5)                                                 \
{ }                                                                        \
                                                                           \
template <class I1, class I2, class I3, class I4, class I5, class I6>      \
CLASS(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4,              \
      const I5 &i5, const I6 &i6)                                          \
  : MEMBER(i1,i2,i3,i4,i5,i6)                                              \
{ }                                                                        \
                                                                           \
template <class I1, class I2, class I3, class I4, class I5, class I6,      \
          class I7>                                                        \
CLASS(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4,              \
      const I5 &i5, const I6 &i6, const I7 &i7)                            \
  : MEMBER(i1,i2,i3,i4,i5,i6,i7)                                           \
{ }


template<class Function, class Patch>
class PatchFunction
{
};

template<class Function>
class PatchFunction<Function,PatchTag1>
{
public:

  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array>
  inline void
  operator()(const Array& a) const
  {
    PatchEvaluator<MainEvaluatorTag>().evaluate(a,function());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

template<class Function>
class PatchFunction<Function,PatchReadTag1>
{
public:

  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array>
  inline void
  operator()(const Array& a) const
  {
    PatchEvaluator<MainEvaluatorTag>().evaluateRead(a,function());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

template<class Function>
class PatchFunction<Function,PatchTag2>
{
public:
  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array1, class Array2>
  inline void
  operator()(const Array1 &a1, const Array2 &a2) const
  {
    PatchEvaluator<MainEvaluatorTag>().evaluate2(a1,a2,function());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

template<class Function>
class PatchFunction<Function,PatchTag3>
{
public:
  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array1, class Array2, class Array3>
  inline void
  operator()(const Array1 &a1, const Array2 &a2, const Array3 &a3) const
  {
    PatchEvaluator<MainEvaluatorTag>().evaluate3(a1,a2,a3,function());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

template<class Function, bool Write1>
class PatchFunction<Function, PatchParticle1<Write1> >
{
public:

  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array>
  inline void
  operator()(const Array& a) const
  {
    ParticleEvaluator().evaluate(a, function(), PatchParticle1<Write1>());
  }

  template<class Array>
  inline void
  block(const Array& a) const
  {
    ParticleEvaluator().
      evaluateBlock(a, function(), PatchParticle1<Write1>());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

template<class Function, bool Write1, bool Write2>
class PatchFunction<Function, PatchParticle2<Write1, Write2> >
{
public:

  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array1, class Array2>
  inline void
  operator()(const Array1 &a1, const Array2 &a2) const
  {
    ParticleEvaluator().evaluate2(a1, a2, function(),
				  PatchParticle2<Write1, Write2>());
  }

  template<class Array1, class Array2>
  inline void
  block(const Array1 &a1, const Array2 &a2) const
  {
    ParticleEvaluator().evaluate2Block(a1, a2, function(),
				       PatchParticle2<Write1, Write2>());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

template<class Function, bool Write1, bool Write2, bool Write3>
class PatchFunction<Function, PatchParticle3<Write1, Write2, Write3> >
{
public:

  PatchFunction() { }
  PatchFunction(const Function &function) : function_m(function) { }

  POOMA_PATCHFUNCTION_ARGUMENT_CONSTRUCTORS(PatchFunction,function_m)

  template<class Array1, class Array2, class Array3>
  inline void
  operator()(const Array1 &a1, const Array2 &a2, const Array3 &a3) const
  {
    ParticleEvaluator().evaluate3(
      a1, a2, a3, function(), PatchParticle3<Write1, Write2, Write3>());
  }

  template<class Array1, class Array2, class Array3>
  inline void
  block(const Array1 &a1, const Array2 &a2, const Array3 &a3) const
  {
    ParticleEvaluator().evaluate3Block(
      a1, a2, a3, function(), PatchParticle3<Write1, Write2, Write3>());
  }

  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_PATCHFUNCTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PatchFunction.h,v $   $Author: richard $
// $Revision: 1.33 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
