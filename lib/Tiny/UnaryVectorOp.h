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

#ifndef POOMA_TINY_UNARY_VECTOR_OP_H
#define POOMA_TINY_UNARY_VECTOR_OP_H

//-----------------------------------------------------------------------------
// Class: UnaryVectorEngine
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// An engine class for representing the sum, product etc between two
// Vectors.  This is used as part of evaluating expressions.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<class V1, class Op>
class UnaryVectorOp;


template<int D, class T, class V1, class Op>
class VectorEngine<D,T,UnaryVectorOp<V1,Op> >
{

public:
 
  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=1 };
  enum { d1 = D };
  typedef T Element_t;
  typedef UnaryVectorOp<V1,Op> EngineTag_t;

  // Return types for accessor functions.
  typedef T ConstElementRef_t;
  typedef T ElementRef_t;

  // Record the type of the current class.
  typedef VectorEngine<D,T, UnaryVectorOp<V1,Op> > This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Construct from one Vector.
  explicit VectorEngine(const V1& v1)
    : v1_m(v1) {}

  Element_t operator()(int i) const
  {
    return Op()(v1_m(i));
  }

#if !POOMA_NO_TEMPLATE_FRIENDS
  template<int DD,class TT, class EE, int I>
    friend struct VectorEngineElem;

private:
#endif

  const V1& v1_m;
};

//-----------------------------------------------------------------------------
//
// Specializaton of TensorElem for UnaryTensorOp
// Compile time lookup for unary operators.
//
//-----------------------------------------------------------------------------

template<int D, class T, class V1, class Op, int I>
struct VectorEngineElem<D,T,UnaryVectorOp<V1,Op>, I>
{
  typedef VectorEngine<D,T,UnaryVectorOp<V1,Op> > V;
  typedef typename VectorElem<V1,I>::Element_t T1;
  typedef typename UnaryReturn<T1,Op>::Type_t Element_t;
  typedef Element_t ElementRef_t;
  typedef Element_t ConstElementRef_t;
  static Element_t get(const V& x) 
    { 
      return Op()(VectorElem<V1,I>::get(x.v1_m));
    }
};


#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UnaryVectorOp.h,v $   $Author: richi $
// $Revision: 1.11 $   $Date: 2004/11/26 12:42:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
