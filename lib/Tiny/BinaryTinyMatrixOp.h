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

#ifndef POOMA_TINY_BINARY_TINYMATRIX_OP_H
#define POOMA_TINY_BINARY_TINYMATRIX_OP_H

//-----------------------------------------------------------------------------
// Class: BinaryTinyMatrixEngine
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Tiny
 * @brief
 * An engine class for representing the sum, product, etc. between two
 * TinyMatrixs.  This is used as part of evaluating expressions.
 *
 * This specialization of TinyMatrixEngine implements an optimization
 * for expression evaluation.  Binary operations are done by constructing
 * TinyMatrix with a BinaryTinyMatrixOp engine tag, then the result of the 
 * operation is constructed by "copying out of" that one.  Copying each 
 * element out of one of these causes the binary expression to be evaluated 
 * for that point.
 *
 * This strategy reduces the number of copies that get made duing expression
 * evaluations.
 */

template<class V1, class V2, class Op>
class BinaryTinyMatrixOp;


/**
 * Specialization of TinyMatrixEngine for BinaryTinyMatrixOp.
 * This has the same interface as any other TinyMatrixEngine of course,
 * but the implementation has a const ref to a V1 and a V2.  Getting
 * an element applies the operator to V1 and V2, returning the result.
 *
 * Since it stores references to V1 and V2, don't try to keep one
 * of these things around!  This should only be used as part of
 * expression evaluation.
 */

template<int D1, int D2, class T, class V1, class V2, class Op>
class TinyMatrixEngine<D1,D2,T,BinaryTinyMatrixOp<V1,V2,Op> >
{

public:
 
  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  typedef T Element_t;
  typedef BinaryTinyMatrixOp<V1,V2,Op> EngineTag_t;

  // Return types for accessor functions.
  typedef T ConstElementRef_t;
  typedef T ElementRef_t;

  // Record the type of the current class.
  typedef TinyMatrixEngine<D1,D2,T, BinaryTinyMatrixOp<V1,V2,Op> > This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Construct from two TinyMatrixs and let the op tag contruct itself.
  TinyMatrixEngine(const V1& v1, const V2& v2)
    : v1_m(v1), v2_m(v2) {}

  Element_t operator()(int i, int j) const
  {
    return Op()(v1_m(i,j), v2_m(i,j));
  }

#if !POOMA_NO_TEMPLATE_FRIENDS
  template<int DD1,int DD2, class TT, class EE, int I, int J>
    friend struct TinyMatrixEngineElem;

private:
#endif

  // Just store const refs to objects of type V1 and V2.
  const V1& v1_m;
  const V2& v2_m;
};


/**
 * Specialization of TinyMatrixElem for BinaryTinyMatrixOp.
 *
 * The type of the return can be different for each element of the 
 * tinymatrix.
 */

template<int D1, int D2, class T, class V1, class V2, class Op, int I, int J>
struct TinyMatrixEngineElem<D1,D2,T,BinaryTinyMatrixOp<V1,V2,Op>, I, J >
{
  typedef TinyMatrixEngine<D1,D2,T,BinaryTinyMatrixOp<V1,V2,Op> > V;
  typedef typename TinyMatrixElem<V1,I,J>::Element_t T1;
  typedef typename TinyMatrixElem<V2,I,J>::Element_t T2;
  typedef typename BinaryReturn<T1,T2,Op>::Type_t Element_t;
  typedef Element_t ElementRef_t;
  typedef Element_t ConstElementRef_t;
  static Element_t get(const V& x) 
    { 
      return Op()(
        TinyMatrixElem<V1,I,J>::get(x.v1_m),
        TinyMatrixElem<V2,I,J>::get(x.v2_m));
    }
};


#endif // POOMA_TINY_BINARY_TINYMATRIX_OP_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BinaryTinyMatrixOp.h,v $   $Author: richi $
// $Revision: 1.8 $   $Date: 2004/11/26 12:42:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
