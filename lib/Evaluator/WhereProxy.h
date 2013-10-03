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
// WhereProxy<F,B>
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_WHEREPROXY_H
#define POOMA_EVALUATOR_WHEREPROXY_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * WhereProxy is used to implement 2 argument where().
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Evaluator/OpMask.h"
#include "Pooma/PETE/ExpressionTraits.h"
#include "Engine/ExpressionEngine.h"

//-----------------------------------------------------------------------------
/// We need the tools to convert a WhereProxy into an Array or Field or
/// whatever, so users need to specialize
//-----------------------------------------------------------------------------

template<class ETrait, class Tree>
struct ConvertWhereProxy
{ };

/**
 * The only legal use of where(a,b) is in an expression like:
 *
 * A = where(f,B);
 *
 * Rather than have where(f,B) return an array that could be combined in
 * an expression, we return a WhereProxy that is recognized by assignment
 * operators.
 *
 * The WhereProxy is also necessary because the elements returned by
 * where are MaskAssign<T> objects, so a special assignment operator,
 * OpMask<Op> must be used.
 */

template<class F, class B>
struct WhereProxy
{
  template <class Cond, class Val, class F1, class B1>
  struct WhereProxyTraits {
    enum { dimensions = F1::dimensions };
    typedef typename ForEach<Val, EvalLeaf<dimensions>, OpCombine>::Type_t Element_t;
  };
  template <class Cond, class T, class F1, class B1>
  struct WhereProxyTraits<Cond, Scalar<T>, F1, B1> {
    enum { dimensions = F1::dimensions };
    typedef T Element_t;
  };
  template <class Val, class T, class F1, class B1>
  struct WhereProxyTraits<Scalar<T>, Val, F1, B1> {
    enum { dimensions = B1::dimensions };
    typedef typename ForEach<Val, EvalLeaf<dimensions>, OpCombine>::Type_t Element_t;
  };
  template <class T1, class T2, class F1, class B1>
  struct WhereProxyTraits<Scalar<T1>, Scalar<T2>, F1, B1> {
    // We open a can of worms, if we try to support this strange case.
    // Just use the simpler
    // if (cond)
    //   lhs = val;
  };

  WhereProxy(const F& f, const B& b) : f_m(f), b_m(b) { }

  typedef BinaryNode<WhereMask,
    typename CreateLeaf<F>::Leaf_t,
    typename CreateLeaf<B>::Leaf_t> Tree_t;

  typedef typename ExpressionTraits<Tree_t>::Type_t           ETrait_t;
  typedef typename ConvertWhereProxy<ETrait_t,Tree_t>::Make_t MakeFromTree_t;
  typedef typename MakeFromTree_t::Expression_t               WhereMask_t;
  typedef typename WhereProxyTraits<typename CreateLeaf<F>::Leaf_t,
	typename CreateLeaf<B>::Leaf_t, F, B>::Element_t      Element_t;

  inline WhereMask_t
  whereMask() const
  {
    return MakeFromTree_t::make(Tree_t(CreateLeaf<F>::make(f_m),
				       CreateLeaf<B>::make(b_m)));
  }

  template<class Op>
  inline OpMask<Op>
  opMask(const Op &op) const
  {
    return OpMask<Op>(op);
  }

  inline const F &flag() { return f_m; }
  inline const B &value() { return b_m; }
  const F &f_m;
  const B &b_m;
};

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

template<class F, class B>
inline WhereProxy<F,B>
where(const F &f, const B &b)
{
  return WhereProxy<F,B>(f,b);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_WHEREPROXY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: WhereProxy.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
