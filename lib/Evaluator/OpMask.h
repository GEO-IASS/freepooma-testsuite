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
// OpMask<Op>
// WhereMask
// MaskAssign<T>
// ForEach<BinaryNode<OpMask>,FTag,OpCombine>
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_OPMASK_H
#define POOMA_EVALUATOR_OPMASK_H

/** @file
 * @ingroup Evaluator
 * @brief
 * These classes implement the two argument where.
 *
 * The expression,
 *
 * a += where(f, b);
 *
 * implements the following code:
 *
 * <PRE>
 * for (loc in domain)
 * {
 *   if (f(loc))
 *     a(loc) += b(loc);
 * }
 * </PRE>
 *
 * To implement this behaviour, the expression is translated to the tree:
 *
 * <PRE>
 *        OpMask<OpAddAssign>
 *           /        |
 *          A        WhereMask
 *                     /    |
 *                    F     B
 * </PRE>
 *
 * ForEach is specialized for WhereMask to evaluate B only if F is true.
 * The result is returned in a MaskAssign<B::Element_t>, which contains
 * the bool from F and the value from B if F is true.  OpMask applies the
 * operator Op to the result from A and B if F is true.
 *
 * This design has the advantage that A is on the LHS and F and B are on the
 * RHS, so we apply write locks to A and read locks to F and B.  Another
 * approach is the tree OpMask(F, OpAddAssign(A,B)), which would lead to
 * complicated interpretation of expressions.
 * Pooma r1 used the implementation OpAddAssign(A, OpWhere(F, B)), which meant
 * that operators needed specializations to deal with the MaskAssign<T>
 * object. 
 *
 * WhereProxy wraps f and b so that the assignment operator can convert it 
 * to the WhereMask object that goes on the RHS and convert the operator
 * into OpMask.
 */ 

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
//
// Conditional Assign wraps the result of WhereMask(F,B);
// It contains the value of B and a bool which is true if the value is defined.
//
//-----------------------------------------------------------------------------

template<class T>
struct MaskAssign
{
  MaskAssign() { }
  MaskAssign(bool q) : cond_m(q) { }
  MaskAssign(bool q, const T& v) : cond_m(q), value_m(v) { }
  ~MaskAssign() { }

  inline bool defined() const { return cond_m; }
  inline const T &value() const { return value_m; }

  inline bool operator!=(const MaskAssign<T> &other) const
  {
    if (defined())
    {
      return ((other.defined() != defined()) || (other.value() != value()));
    }
    else
    {
      return other.defined();
    }
  }

  // To make purify happy:

  MaskAssign(const MaskAssign<T> &) { }
  MaskAssign<T> &operator=(const MaskAssign<T> &) { return *this; }

  bool cond_m;
  T value_m;
};

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

template<class Op>
struct OpMask
{
  OpMask() { }
  OpMask(const Op &op) : op_m(op) { }
  ~OpMask() { }

  /// WhereProxy Op, embed a conditional operation.
  template<class T1, class T2>
  inline void
  operator()(T1 &a, const MaskAssign<T2> &b) const
  {
    if (b.defined())
    {
      op_m(a, b.value());
    }
  }

  /// Fall back to native operation.
  template<class T1, class T2>
  inline void
  operator()(T1 &a, const T2 &b) const
  {
    op_m(a, b);
  }

  Op op_m;
};

template<class T1, class T2, class Op>
struct BinaryReturn<T1, T2, OpMask<Op> >
{
  typedef T1 &Type_t;
};

template <class Op, class T>
struct ReductionTraits;

template <class Op, class T>
struct ReductionTraits<OpMask<Op>, T>
{
  static T identity() { return ReductionTraits<Op, T>::identity(); }
};


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

struct WhereMask
{
  WhereMask() { }
  WhereMask(const WhereMask &) { }
  WhereMask &operator=(const WhereMask &) { return *this; }

  ~WhereMask() { }
  // Should never explitly evaluate WhereMask
};

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

template<class T1, class T2>
struct BinaryReturn<T1, T2, WhereMask>
{
  typedef MaskAssign<T2> Type_t;
};

//-----------------------------------------------------------------------------
//
// ForEach
//
// To evaluate a WhereMask, we have to specialize ForEach.  The right hand
// expression should not be evaluated at all if the mask isn't true, so
// the stucture of the member function apply() has to be different from
// the generic ForEach.
//
//-----------------------------------------------------------------------------

template<class A, class B, class FTag>
struct ForEach< BinaryNode<WhereMask, A, B>, FTag, OpCombine >
{
  // The return type for the left tree.
  // really should be a bool
  typedef typename ForEach<A,FTag,OpCombine>::Type_t TypeA_t;

  // The return type for the right tree.
  typedef typename ForEach<B,FTag,OpCombine>::Type_t TypeB_t;

  // The return type for the expression.
  typedef MaskAssign<TypeB_t> Type_t;

  // How to evaluate the expression.
  inline
  static Type_t
  apply(const BinaryNode<WhereMask,A,B>& expr, 
	const FTag &f, const OpCombine &c) 
  {
    // Evaluate the left.
    bool mask = forEach(expr.left(), f, c);

    if ( mask )
    {
      return Type_t(mask, forEach(expr.right(), f, c));
    }
    else
    {
      return Type_t(mask);
    }
  }
};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: OpMask.h,v $   $Author: richard $
// $Revision: 1.23 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
