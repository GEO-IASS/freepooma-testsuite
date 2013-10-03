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
// MakeFieldReturn<T>
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_FIELDMAKERETURN_H
#define POOMA_FIELD_FIELDMAKERETURN_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

struct FarLeftTag;

template<class G, class T, class E> class Field;


/** @file
 * @ingroup Field
 * @brief
 * MakeFieldReturn is used to combine expressions together with operators.
 */

/**
 * MakeFieldReturn is a tool used by operator functions to construct the
 * expression tree representing that function.  Each function needs to define
 * a corresponding operator functor Op which is used to compute the return
 * type.  The required interface for MakeFieldReturn is:
 *  - typedef ... Expression_t;    // type of the expression (UnaryNode<...>)
 *  - Expression_t make(const T&); // construct the tree
 *
 * These versions are a little more complicated than those for Array because we
 * want to preserve Geometry information to the largest extent possible.
 */

template<class Expr> 
struct MakeFieldReturn;

//-----------------------------------------------------------------------------
// op(Expression)
//
// This includes all combinations of scalars, arrays, and fields.
//-----------------------------------------------------------------------------

template<class Op, class Leaf>
struct MakeFieldReturn<UnaryNode<Op, Leaf> >
{
  // The node type.
  
  typedef UnaryNode<Op, Leaf> Tree_t;

  // Deduce the template parameters of the expression engine we're building.
  
  typedef typename 
    ForEach<Tree_t, DomainFunctorTag, DomainFunctorTag>::Type_t Domain_t;
  enum { dim = Domain_t::dimensions };
  typedef typename ForEach<Tree_t, EvalLeaf<dim>, OpCombine>::Type_t T_t;
  typedef Engine<dim, T_t, ExpressionTag<Tree_t> > Engine_t;
  typedef typename ForEach<Tree_t, FarLeftTag, FarLeftTag>::
    Type_t::MeshTag_t MeshTag_t;

  // Construct the type of the Field we're making.

  typedef Field<MeshTag_t, T_t, ExpressionTag<Tree_t> > Expression_t;

  // This function turns the tree node into a Field.
  
  inline static
  Expression_t make(const Tree_t &tree)
    {
      return Expression_t(Engine_t(tree));
    }
};


//-----------------------------------------------------------------------------
// Expression op Expression
//
// This includes all combinations of scalars, arrays, and fields.
//-----------------------------------------------------------------------------

template<class Op, class Left, class Right>
struct MakeFieldReturn<BinaryNode<Op, Left, Right> >
{
  // The node type.
  
  typedef BinaryNode<Op, Left, Right> Tree_t;

  // Deduce the template parameters of the expression engine we're building.
  
  typedef typename 
    ForEach<Tree_t, DomainFunctorTag, DomainFunctorTag>::Type_t Domain_t;
  enum { dim = Domain_t::dimensions };
  typedef typename ForEach<Tree_t, EvalLeaf<dim>, OpCombine>::Type_t T_t;
  typedef Engine<dim, T_t, ExpressionTag<Tree_t> > Engine_t;
  typedef typename ForEach<Tree_t, FarLeftTag, FarLeftTag>::
    Type_t::MeshTag_t MeshTag_t;
  
  // Construct the type of the Field we're making.

  typedef Field<MeshTag_t, T_t, ExpressionTag<Tree_t> > 
    Expression_t;

  // This function turns the tree node into a Field.

  inline static
  Expression_t make(const Tree_t &tree)
    {
      return Expression_t(Engine_t(tree));
    }
};


//-----------------------------------------------------------------------------
// Expression "?" Expression ":" Expression
//
// This includes all combinations of scalars, arrays, and fields.
//-----------------------------------------------------------------------------

template<class Op, class Left, class Middle, class Right>
struct MakeFieldReturn<TrinaryNode<Op, Left, Middle, Right> >
{
  // The node type.
  
  typedef TrinaryNode<Op, Left, Middle, Right> Tree_t;

  // Deduce the template parameters of the expression engine we're building.
  
  typedef typename 
    ForEach<Tree_t, DomainFunctorTag, DomainFunctorTag>::Type_t Domain_t;
  enum { dim = Domain_t::dimensions };
  typedef typename ForEach<Tree_t, EvalLeaf<dim>, OpCombine>::Type_t T_t;
  typedef Engine<dim, T_t, ExpressionTag<Tree_t> > Engine_t;
  typedef typename ForEach<Tree_t, FarLeftTag, FarLeftTag>::
    Type_t::MeshTag_t MeshTag_t;
  
  // Construct the type of the Field we're making.

  typedef Field<MeshTag_t, T_t, ExpressionTag<Tree_t> > Expression_t;

  // This function turns the tree node into a Field.

  inline static
  Expression_t make(const Tree_t &tree)
    {
      return Expression_t(Engine_t(tree));
    }
};


#endif // POOMA_FIELD_FIELDMAKERETURN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldMakeReturn.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:43 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
