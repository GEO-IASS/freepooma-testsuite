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
// CreateLeaf<Expr>
// MakeFieldReturn<T>
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_FIELDCREATELEAF_H
#define POOMA_FIELD_FIELDCREATELEAF_H

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
 * CreateLeaf is used to convert arbitrary classes into expression objects.
 */

//-----------------------------------------------------------------------------
// Traits classes
//-----------------------------------------------------------------------------

/**
 * CreateLeaf is an external functor class used to convert objects into the
 * leaves of the expression tree.
 *
 * CreateLeaf<T> converts objects of type T to leaf objects and requires
 * the following interface:
 *  - typedef ... Leaf_t;        // The leaf object
 *  - typedef ... Return_t;      // Type returned by make()
 *  - Return_t make(const T&);   // make the leaf object from the T object
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
 * <PRE>
 * template<class G, class T, class E>
 * void func(const Field<G,T,E>& f)
 * {
 *   forEach(CreateLeaf<Field<G,T,E >::make(f),...,...);
 * }
 * </PRE>
 */

//-----------------------------------------------------------------------------
// Fields are leaf objects, we just pass them through unless they have
// Expression engines. Then, we remove the expression from the field and form
// leaves with tree-nodes. 

template<class GeometryTag, class T, class Expr>
struct CreateLeaf<Field<GeometryTag, T, ExpressionTag<Expr> > >
{
  typedef Field<GeometryTag, T, ExpressionTag<Expr> > Input_t;
  typedef Expr Leaf_t;
  typedef const Leaf_t &Return_t;
  inline static
  Return_t make(const Input_t &f)
    {
      return f.engine().expression();
    }
};
  
template<class GeometryTag, class T, class EngineTag>
struct CreateLeaf<Field<GeometryTag, T, EngineTag> >
{
  typedef Field<GeometryTag, T, EngineTag> Input_t;
  typedef Reference<Input_t> Leaf_t;
  typedef Leaf_t Return_t;
  inline static
  Return_t make(const Input_t &f)
    {
      return Leaf_t(f);
    }
};

//-----------------------------------------------------------------------------
// Special case for Scalar<Field> returns ErrorType to avoid
// hairy type computations.

template<class GeometryTag, class T, class EngineTag>
struct CreateLeaf<Scalar<Field<GeometryTag, T, EngineTag> > >
{
  typedef Scalar<Field<GeometryTag, T, EngineTag> > Input_t;
  typedef Scalar<ErrorType> Leaf_t;
  typedef Leaf_t Return_t;
  inline static
  Return_t make(const Input_t &)
    {
      return ErrorType();
    }
};


#endif // POOMA_FIELD_FIELDCREATELEAF_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldCreateLeaf.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:42 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
