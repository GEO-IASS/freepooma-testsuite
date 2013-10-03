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

#ifndef POOMA_ARRAY_CREATELEAF_H
#define POOMA_ARRAY_CREATELEAF_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Classes:
// CreateLeaf<Expr>
// MakeReturn<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Array
 * @brief
 * These are the external traits classes that are used to build trees.
 *
 * CreateLeaf is used to convert arbitrary classes into expression objects.
 * MakeReturn are used to combine expressions together with operators.
 */

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Traits classes
//-----------------------------------------------------------------------------

/** @name CreateLeaf
 * CreateLeaf - external functor class used to convert objects into the
 * leaves of the expression tree.
 *
 * For the current version of pooma, our objects are arrays and we actually
 * stick arrays at the leaves of the expression.  Also, we store expressions
 * in arrays and want to extract the expression if we use the array in a
 * bigger expression.
 * (Other users of expression templates might wish to store iterates in the
 * expression objects.)
 *
 * CreateLeaf<T> converts objects of type T to leaf objects and requires
 * the following interface:
 *  - typedef ... Leaf_t;       The leaf object
 *  - typedef ... Return_t;     Type returned by make()
 *  - Return_t make(const T&);  make the leaf object from the T object
 *
 * Return_t should equivalent to Leaf_t. (Leaf_t needs to be able
 * be constructed with a Return_t.)  We avoid making extra copies by building
 * expression trees from references, so define Return_t to be a const ref to
 * an Leaf_t.  (Returning by value would be bad, since we would create
 * a temporary copy that won't survive until the whole expression is put
 * together.)
 *
 * CreateLeaf is used to construct expression trees.  It should also be
 * used when performing operations on the expression tree, such as forEach,
 * in order to extract the expression.  For example:
 * template<int D,class T,class E>
 * void func(const Array<D,T,E>& array)
 * {
 *   forEach(CreateLeaf<Array<D,T,E> >::make(array),...,...);
 * }
 */

//@{

/// Arrays are leaf objects, we just pass them through unless they have
/// Expression engines. 

template<int Dim, class T, class EngineTag>
struct CreateLeaf<Array<Dim, T, EngineTag> >
{
  typedef Array<Dim, T, EngineTag> Input_t;
  typedef Reference<Input_t> Leaf_t;
  typedef Leaf_t Return_t;
  inline static
  Return_t make(const Input_t &a)
    {
      return Leaf_t(a);
    }
};

template<int Dim, class T, class Expr>
struct CreateLeaf<Array<Dim, T, ExpressionTag<Expr> > >
{
  typedef Array<Dim, T, ExpressionTag<Expr> > Input_t;
  typedef Expr Leaf_t;
  typedef const Leaf_t &Return_t;
  inline static
  Return_t make(const Input_t &a)
    {
      return a.engine().expression();
    }
};

/// Special case for Scalar<Array> returns ErrorType to avoid
/// hairy type computations.

template<int Dim, class T, class EngineTag>
struct CreateLeaf<Scalar<Array<Dim, T, EngineTag> > >
{
  typedef Scalar<Array<Dim, T, EngineTag> > Input_t;
  typedef Scalar<ErrorType> Leaf_t;
  typedef Leaf_t Return_t;
  inline static
  Return_t make(const Input_t &)
    {
      return ErrorType();
    }
};

//@}


/** @name MakeReturn
 * MakeReturn is tool used by operator functions to construct the
 * expression tree representing that function.
 */

//@{

/// Unary node version

template<class Op,class Leaf>
struct MakeReturn<UnaryNode<Op,Leaf> >
{
  typedef UnaryNode<Op,Leaf> Tree_t;

  typedef typename ForEach<Tree_t,
    DomainFunctorTag, DomainFunctorTag>::Type_t Domain_t;
  enum { dim = Domain_t::dimensions };
  typedef typename UnaryReturn<typename ForEach<Leaf,
    EvalLeaf<dim>,OpCombine>::Type_t,
    Op>::Type_t T_t;
  typedef Engine<dim,T_t,ExpressionTag<Tree_t> > Engine_t;
  typedef Array<dim,T_t,ExpressionTag<Tree_t > > Expression_t;

  inline static
  Expression_t make(const Tree_t &tree)
    {
      return Expression_t(Engine_t(tree));
    }
};

/// Binary node version

template<class Op,class Left,class Right>
struct MakeReturn<BinaryNode<Op,Left,Right> >
{
  typedef BinaryNode<Op,Left,Right> Tree_t;

  typedef typename ForEach<Tree_t,
    DomainFunctorTag, DomainFunctorTag>::Type_t Domain_t;
  enum { dim = Domain_t::dimensions };
  typedef typename BinaryReturn<typename ForEach<Left,
    EvalLeaf<dim>,OpCombine>::Type_t,
    typename ForEach<Right,EvalLeaf<dim>,OpCombine>::Type_t,
    Op>::Type_t T_t;
  typedef Engine<dim,T_t,ExpressionTag<Tree_t> > Engine_t;
  typedef Array<dim,T_t,ExpressionTag<Tree_t > > Expression_t;

  inline static
  Expression_t make(const Tree_t &tree)
    {
      return Expression_t(Engine_t(tree));
    }
};

/// Trinary node version

template<class Op,class Cl,class Tr,class Fl>
struct MakeReturn<TrinaryNode<Op,Cl,Tr,Fl> >
{
  typedef TrinaryNode<Op,Cl,Tr,Fl> Tree_t;

  typedef typename ForEach<Tree_t,
    DomainFunctorTag, DomainFunctorTag>::Type_t Domain_t;
  enum { dim = Domain_t::dimensions };
  typedef typename TrinaryReturn<typename ForEach<Cl,
    EvalLeaf<dim>,OpCombine>::Type_t,
    typename ForEach<Tr,EvalLeaf<dim>,OpCombine>::Type_t,
    typename ForEach<Fl,EvalLeaf<dim>,OpCombine>::Type_t,
    Op>::Type_t T_t;
  typedef Engine<dim,T_t,ExpressionTag<Tree_t> > Engine_t;
  typedef Array<dim,T_t,ExpressionTag<Tree_t > > Expression_t;

  inline static
  Expression_t make(const Tree_t &tree)
    {
      return Expression_t(Engine_t(tree));
    }
};

//@}

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CreateLeaf.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:13 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
