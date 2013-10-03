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
// ExpressionTraits
// ExpressionIsScalar
//-----------------------------------------------------------------------------

/** @file
 * @ingroup PETE
 * @brief
 * Undocumented.
 */

#ifndef POOMA_PETE_EXPRESSIONTRAITS_H
#define POOMA_PETE_EXPRESSIONTRAITS_H

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

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * ExpressionTraits<Expr>::Type_t is a traits class for expression objects.
 * We use it to determine if an expression is an array expression or a field
 * expression to decide what kind of object to wrap the expression inside.
 */

template<class T>
struct ExpressionTraits
{
  // Some default here seems to be necessary.
  typedef void Type_t;
};

/**
 * ExpressionIsScalar is the return type for scalar expression objects.
 */

struct ExpressionIsScalar { };

template<class T>
struct ExpressionTraits<Scalar<T> >
{
  typedef ExpressionIsScalar Type_t;
};

/**
 * CombineExpressionTraits<Trait1,Trait2>::Type_t is used to describe
 * the result of various expression types.  If we ever decide to something
 * wacky like sin(field) is a field but field+field is not, then we will
 * have to replace this concept with a general ForEach computation.
 */

template<class A, class B>
struct CombineExpressionTraits
{ };

/**
 * Determine the ExpressionTraits for all the expression objects:
 */

template<class T>
struct ExpressionTraits<Reference<T> >
{
  typedef typename ExpressionTraits<T>::Type_t Type_t;
};

template<class Op, class Child>
struct ExpressionTraits<UnaryNode<Op, Child> >
{
  typedef typename ExpressionTraits<Child>::Type_t Type_t;
};

template<class Op, class Left, class Right>
struct ExpressionTraits<BinaryNode<Op, Left, Right> >
{
  typedef typename ExpressionTraits<Left>::Type_t  Left_t;
  typedef typename ExpressionTraits<Right>::Type_t Right_t;
  typedef typename CombineExpressionTraits<Left_t, Right_t>::Type_t Type_t;
};

template<class Op, class Left, class Middle, class Right>
struct ExpressionTraits<TrinaryNode<Op, Left, Middle, Right> >
{
  typedef typename ExpressionTraits<Left>::Type_t    Left_t;
  typedef typename ExpressionTraits<Middle>::Type_t  Middle_t;
  typedef typename ExpressionTraits<Right>::Type_t   Right_t;
  typedef typename CombineExpressionTraits<Left_t, Right_t>::Type_t Temp_t;
  typedef typename CombineExpressionTraits<Temp_t, Middle_t>::Type_t Type_t;
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PETE_EXPRESSIONTRAITS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ExpressionTraits.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:05 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
