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
// ErrorType
// ForEachRef
//-----------------------------------------------------------------------------

#ifndef POOMA_PETE_POOMAPETE_H
#define POOMA_PETE_POOMAPETE_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//
// This file contains POOMA specific extensions to PETE.  Please refrain from
// changing other PETE files in this directory without making corresponding
// changes to the PETE repository.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/ErrorType.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Reference versions of ForEach.
//-----------------------------------------------------------------------------

template<class Expr, class FTag, class CTag>
struct ForEachRef
{
  typedef typename LeafFunctor<Expr, FTag>::Type_t Type_t;
  inline static
  const Type_t &apply(const Expr &expr, const FTag &f, const CTag &)
    {
      return LeafFunctor<Expr, FTag>::apply(expr, f);
    }
};

template<class Expr, class FTag, class CTag>
inline const typename ForEachRef<Expr,FTag,CTag>::Type_t &
forEachRef(const Expr &e, const FTag &f, const CTag &c)
{
  return ForEachRef<Expr, FTag, CTag>::apply(e, f, c);
}

template<class Op, class A, class FTag, class CTag>
struct ForEachRef<UnaryNode<Op, A>, FTag, CTag>
{
  typedef typename ForEachRef<A, FTag, CTag>::Type_t TypeA_t;
  typedef typename Combine1<TypeA_t, Op, CTag>::Type_t Type_t;
  inline static
  const Type_t &apply(const UnaryNode<Op, A> &expr, const FTag &f, 
    const CTag &c) 
    {
      return Combine1<TypeA_t, Op, CTag>::
        combine(ForEachRef<A, FTag, CTag>::apply(expr.child(), f, c),
                c);
    }
};

template<class Op, class A, class B, class FTag, class CTag>
struct ForEachRef<BinaryNode<Op, A, B>, FTag, CTag >
{
  typedef typename ForEachRef<A, FTag, CTag>::Type_t TypeA_t;
  typedef typename ForEachRef<B, FTag, CTag>::Type_t TypeB_t;
  typedef typename Combine2<TypeA_t, TypeB_t, Op, CTag>::Type_t Type_t;
  inline static
  const Type_t &apply(const BinaryNode<Op, A, B> &expr, const FTag &f,
	                  const CTag &c) 
    {
      return Combine2<TypeA_t, TypeB_t, Op, CTag>::
        combine(ForEachRef<A, FTag, CTag>::apply(expr.left(), f, c),
                ForEachRef<B, FTag, CTag>::apply(expr.right(), f, c),
	        c);
    }
};

template<class Op, class A, class B, class C, class FTag, class CTag>
struct ForEachRef<TrinaryNode<Op, A, B, C>, FTag, CTag >
{
  typedef typename ForEachRef<A, FTag, CTag>::Type_t TypeA_t;
  typedef typename ForEachRef<B, FTag, CTag>::Type_t TypeB_t;
  typedef typename ForEachRef<C, FTag, CTag>::Type_t TypeC_t;
  typedef typename Combine3<TypeA_t, TypeB_t, TypeC_t, Op, CTag>::Type_t 
    Type_t;
  inline static
  const Type_t &apply(const TrinaryNode<Op, A, B, C> &expr, const FTag &f,
	                  const CTag &c) 
    {
      return Combine3<TypeA_t, TypeB_t, TypeC_t, Op, CTag>::
        combine(ForEachRef<A, FTag, CTag>::apply(expr.left(), f, c),
	        ForEachRef<B, FTag, CTag>::apply(expr.middle(), f, c),
	        ForEachRef<C, FTag, CTag>::apply(expr.right(), f, c),
	        c);
    }
};

template<class T, class FTag, class CTag>
struct ForEachRef<Expression<T>, FTag, CTag>
{
  typedef typename ForEachRef<T, FTag, CTag>::Type_t Type_t;
  inline static
  const Type_t &apply(const Expression<T> &expr, const FTag &f, 
	                  const CTag &c) 
    {
      return ForEachRef<T, FTag, CTag>::apply(expr.expression(), f, c);
    }
};

template<class T, class FTag, class CTag>
struct ForEachRef<Reference<T>, FTag, CTag>
{
  typedef typename ForEachRef<T, FTag, CTag>::Type_t Type_t;
  inline static
  const Type_t &apply(const Reference<T> &ref, const FTag &f, 
	                  const CTag &c) 
    {
      return ForEachRef<T, FTag, CTag>::apply(ref.reference(), f, c);
    }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PETE_POOMAPETE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PoomaPETE.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:56 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
