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

#ifndef POOMA_TINY_BINARY_VECTOR_OP_H
#define POOMA_TINY_BINARY_VECTOR_OP_H

//-----------------------------------------------------------------------------
// Class: BinaryVectorEngine
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Tiny
 * @brief
 * An engine class for representing the sum, product etc between two
 * vectors.  This is used as part of evaluating expressions.
 *
 * This specialization of VectorEngine implements an optimization
 * for expression evaluation.  Binary operations are done by constructing
 * a Vector with a BinaryVectorOp engine tag, then the result of the 
 * operation is constructed by "copying out of" that one.  Copying each 
 * element out of one of these causes the binary expression to be evaluated 
 * for that point.
 */

template<class V1, class V2, class Op>
class BinaryVectorOp;


template<int D, class T, class V1, class V2, class Op>
class VectorEngine<D,T, BinaryVectorOp<V1,V2,Op> >
{

public:
 
  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=1 };
  enum { d1=D };
  typedef T Element_t;
  typedef BinaryVectorOp<V1,V2,Op> EngineTag_t;

  // Return types for accessor functions.
  typedef T ConstElementRef_t;
  typedef T ElementRef_t;

  // Record the type of the current class.
  typedef VectorEngine<D,T, BinaryVectorOp<V1,V2,Op> > This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Construct from two vectors and let the op tag contruct itself.
  VectorEngine(const V1& v1, const V2& v2)
    : v1_m(v1), v2_m(v2) {}

  Element_t operator()(int i) const
  {
    return Op()(v1_m(i), v2_m(i));
  }

  //----------------------------------------------------------------------
  // Element access

#if !POOMA_NO_TEMPLATE_FRIENDS
  template<int DD,class TT, class EE, int I>
    friend struct VectorEngineElem;

private:
#endif

  const V1& v1_m;
  const V2& v2_m;
};


/**
 * Specialization of VectorElem for BinaryVectorOp.
 * Compile time element lookup.
 */

template<int D, class T, class V1, class V2, class Op, int I>
struct VectorEngineElem<D,T,BinaryVectorOp<V1,V2,Op>, I >
{
  typedef VectorEngine<D,T,BinaryVectorOp<V1,V2,Op> > V;
  typedef typename VectorElem<V1,I>::Element_t T1;
  typedef typename VectorElem<V2,I>::Element_t T2;
  typedef typename BinaryReturn<T1,T2,Op>::Type_t Element_t;
  typedef Element_t ElementRef_t;
  typedef Element_t ConstElementRef_t;
  static Element_t get(const V& x) 
    { 
      return Op()(
        VectorElem<V1,I>::get(x.v1_m),
        VectorElem<V2,I>::get(x.v2_m));
    }
};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BinaryVectorOp.h,v $   $Author: richi $
// $Revision: 1.14 $   $Date: 2004/11/26 12:42:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
