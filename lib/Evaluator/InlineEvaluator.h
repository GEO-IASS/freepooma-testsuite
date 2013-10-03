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

#ifndef POOMA_EVALUATOR_INLINEEVALUATOR_H
#define POOMA_EVALUATOR_INLINEEVALUATOR_H

//-----------------------------------------------------------------------------
// Class: InlineEvaluator
//-----------------------------------------------------------------------------


/** @file
 * @ingroup Evaluator
 * @brief
 * InlineEvaluator evaluates expressions by inlining a simple loop.
 * It does no dependency checking, locking, where blocks, etc.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Evaluator/KernelTags.h"
#include "Utilities/WrappedInt.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class KernelTag>
struct KernelEvaluator;

/**
 * The point of this class is to input an expression with the
 * 'evaluate' member function and evaluate it by looping over the
 * whole domain.
 *
 * This is the simplest possible evaluator. It makes a number of
 * simplifying assumptions about the expressions it tries to evaluate
 * and the context in which they are evaluated.  These assumptions let
 * it do some things very efficiently, but limit the contexts in which
 * it can be used.
 *
 * These assumptions are:
 * -# There are no where blocks. That means that the InlineEvaluator
 *    does not need to have any state.
 * -# The expression passed in can handle random access to all of its
 *    elements efficiently.  That basically means that it can only be
 *    used with BrickEngine or its equivalent.
 */

template<>
struct KernelEvaluator<InlineKernelTag>
{
  /// Evaluate an expression on a given domain.  This function must be
  /// specialized for particular domain types.  The expectation is that
  /// it will just loop over the domain and use random access in the
  /// expression to evaluate it.

  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain)
  {
    // All the evaluators assume unit-stride, zero-based domains.
    CTAssert(Domain::unitStride);
    for (int i=0; i<Domain::dimensions; ++i)
      PAssert(domain[i].first() == 0);

    evaluate(lhs,op,rhs,domain,
	     WrappedInt<Domain::dimensions>());

    POOMA_INCREMENT_STATISTIC(NumInlineEvaluations)
  }

  /// Input an expression and cause it to be evaluated.
  /// All this template function does is extract the domain
  /// from the expression and call evaluate on that.

  template<class LHS,class Op,class RHS>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs)
  {
    evaluate(lhs,op,rhs,lhs.domain());
  }


  ///@name evaluate(expression,domain,domain_dimension)
  //@{
  /// This is the function both of the above functions call.
  /// It adds a third argument which is a tag class templated on
  /// the dimension of the domain.
  ///
  /// This parameter lets us specialize the function based on
  /// that dimension.
  ///
  /// Some day, we will figure out how to specialize template 
  /// member functions outside the class declaration...
  ///
  /// These functions are all inline for efficiency. That means that if
  /// they are being used at the user level we will get the optimization
  /// of recognizing multiple uses of a single Array on the right hand
  /// side.
  ///
  /// There are seven specializations here, for dimension 1 through 7.
  /// Rather than use template metaprograms for these seven cases we
  /// simply enumerate them explicitly.  This is done to reduce the
  /// burden on the compiler, which would otherwise have to jump through
  /// a bunch of hoops to get the code that is here.
  ///
  /// For each of the specializations it builds a nested loop for each
  /// dimension. Each loop is constructed with last() from the
  /// appropriate dimension of the zero-based domain.
  ///
  /// NOTE: These loops assume that the domain passed in is a unit-stride
  /// domain starting at 0.  Assertions are made to make sure this is true
  /// in the dispatching evaluate above.

  
  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<1>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
#pragma omp parallel for if (e0 > 512)
    for (int i0=0; i0<e0; ++i0)
      op(localLHS(i0),localRHS.read(i0));
  }

  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<2>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
#pragma omp parallel for
    for (int i1=0; i1<e1; ++i1)
      for (int i0=0; i0<e0; ++i0)
	op(localLHS(i0,i1),localRHS.read(i0,i1));
  }
  
  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<3>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
#pragma omp parallel for
    for (int i2=0; i2<e2; ++i2)
      for (int i1=0; i1<e1; ++i1)
	for (int i0=0; i0<e0; ++i0)
	  op(localLHS(i0,i1,i2),localRHS.read(i0,i1,i2));
  }

  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<4>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
#pragma omp parallel for
    for (int i3=0; i3<e3; ++i3)
      for (int i2=0; i2<e2; ++i2)
	for (int i1=0; i1<e1; ++i1)
	  for (int i0=0; i0<e0; ++i0)
	    op(localLHS(i0,i1,i2,i3),localRHS.read(i0,i1,i2,i3));
  }

  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<5>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    int e4 = domain[4].length();
#pragma omp parallel for
    for (int i4=0; i4<e4; ++i4)
      for (int i3=0; i3<e3; ++i3)
	for (int i2=0; i2<e2; ++i2)
	  for (int i1=0; i1<e1; ++i1)
	    for (int i0=0; i0<e0; ++i0)
	      op(localLHS(i0,i1,i2,i3,i4),localRHS.read(i0,i1,i2,i3,i4));
  }

  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<6>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    int e4 = domain[4].length();
    int e5 = domain[5].length();
#pragma omp parallel for
    for (int i5=0; i5<e5; ++i5)
      for (int i4=0; i4<e4; ++i4)
	for (int i3=0; i3<e3; ++i3)
	  for (int i2=0; i2<e2; ++i2)
	    for (int i1=0; i1<e1; ++i1)
	      for (int i0=0; i0<e0; ++i0)
		op(localLHS(i0,i1,i2,i3,i4,i5),
		   localRHS.read(i0,i1,i2,i3,i4,i5));
  }

  template<class LHS,class Op,class RHS,class Domain>
  inline static void evaluate(const LHS& lhs,const Op& op,const RHS& rhs,
			      const Domain& domain,WrappedInt<7>)
  {
    LHS localLHS(lhs);
    RHS localRHS(rhs);
    int e0 = domain[0].length();
    int e1 = domain[1].length();
    int e2 = domain[2].length();
    int e3 = domain[3].length();
    int e4 = domain[4].length();
    int e5 = domain[5].length();
    int e6 = domain[6].length();
#pragma omp parallel for
    for (int i6=0; i6<e6; ++i6)
      for (int i5=0; i5<e5; ++i5)
	for (int i4=0; i4<e4; ++i4)
	  for (int i3=0; i3<e3; ++i3)
	    for (int i2=0; i2<e2; ++i2)
	      for (int i1=0; i1<e1; ++i1)
		for (int i0=0; i0<e0; ++i0)
		  op(localLHS(i0,i1,i2,i3,i4,i5,i6),
		     localRHS.read(i0,i1,i2,i3,i4,i5,i6));
  }

  //@}

private:

};

//-----------------------------------------------------------------------------

#endif // POOMA_EVALUATOR_INLINEEVALUATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: InlineEvaluator.h,v $   $Author: richi $
// $Revision: 1.31 $   $Date: 2004/11/04 20:07:26 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
