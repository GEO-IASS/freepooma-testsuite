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

#ifndef POOMA_TINY_TINYMATRIX_ELEMENTS_H
#define POOMA_TINY_TINYMATRIX_ELEMENTS_H

//-----------------------------------------------------------------------------

// Class: 
//    TinyMatrixElem
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
// Trait classes for getting elements of Tiny objects at compile time.
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

template<int D1, int D2, class T, class E> class TinyMatrix;
template<int D1, int D2, class T, class E> class TinyMatrixEngine;

//-----------------------------------------------------------------------------
//
// Full Description:
//
// The general templates for the class TinyMatrixElem.
// VectorElem should be specialized for TinyMatrix-like classes.
//
// The general definition is for scalars which cannot be subscripted.
// We also have specializations for TinyMatrixs with arbitrary engines
// which just use operator() with both integers.  This is the 
// fallback if a given engine type doesn't specify anything else.
//
//-----------------------------------------------------------------------------

template<class V, int I, int J>
struct TinyMatrixElem
{
  typedef       V  Element_t;
  typedef const V& ConstElementRef_t;
  typedef       V& ElementRef_t;
  static ConstElementRef_t get(const V& x) { return x; }
  static      ElementRef_t get(      V& x) { return x; }
};

template<int D1, int D2, class T, class E, int I, int J>
struct TinyMatrixEngineElem
{
  typedef TinyMatrixEngine<D1,D2,T,E> V;
  typedef typename V::Element_t         Element_t;
  typedef typename V::ConstElementRef_t ConstElementRef_t;
  typedef typename V::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return x(I,J); }
  static ElementRef_t      get(      V& x) { return x(I,J); }
};

template<int D1, int D2, class T, class E, int I, int J>
struct TinyMatrixElem< TinyMatrix<D1,D2,T,E> , I , J>
{
  typedef TinyMatrix<D1,D2,T,E> V;
  typedef TinyMatrixEngineElem<D1,D2,T,E,I,J> TE;
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x.engine()); }
  static      ElementRef_t get(V& x)       { return TE::get(x.engine()); }
};

//-----------------------------------------------------------------------------
//
// TinyMatrixAssign
//
// Template metaprogram for copying out of one TinyMatrix and into another.
// Input:
//   The TinyMatrix we're writing into.
//   Something to copy out of.
//
// Evaluate by recursing on the quadrants of the TinyMatrix.
//
//-----------------------------------------------------------------------------

//
// The general case of copying divides the TinyMatrix into quadrants
// and calls copy on each quadrant.
// This will be applied if each axis of the TinyMatrix is larger than 1.
//

template<class T1, class T2, class Op, int B1, int L1, int B2, int L2>
struct TinyMatrixAssign
{
  enum { B11=B1 , L11=L1/2 , B12=B1+L1/2 , L12 = L1-L1/2 };
  enum { B21=B2 , L21=L2/2 , B22=B2+L2/2 , L22 = L2-L2/2 };
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      TinyMatrixAssign<T1,T2,Op,B11,L11,B21,L21>::apply(x,y,op);
      TinyMatrixAssign<T1,T2,Op,B12,L12,B21,L21>::apply(x,y,op);
      TinyMatrixAssign<T1,T2,Op,B11,L11,B22,L22>::apply(x,y,op);
      TinyMatrixAssign<T1,T2,Op,B12,L12,B22,L22>::apply(x,y,op);
    }
};

//
// The case for a column.
// Divide the column in two, and recurse.
//

template<class T1, class T2, class Op, int B1, int L1, int B2>
struct TinyMatrixAssign<T1,T2,Op,B1,L1,B2,1>
{
  enum { B11=B1 , L11=L1/2 , B12=B1+L1/2 , L12 = L1-L1/2 };
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      TinyMatrixAssign<T1,T2,Op,B11,L11,B2,1>::apply(x,y,op);
      TinyMatrixAssign<T1,T2,Op,B12,L12,B2,1>::apply(x,y,op);
    }
};

//
// The case for a row.
// Divide the row in two, and recurse.
//

template<class T1, class T2, class Op, int B1, int B2, int L2>
struct TinyMatrixAssign<T1,T2,Op,B1,1,B2,L2>
{
  enum { B21=B2 , L21=L2/2 , B22=B2+L2/2 , L22 = L2-L2/2 };
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      TinyMatrixAssign<T1,T2,Op,B1,1,B21,L21>::apply(x,y,op);
      TinyMatrixAssign<T1,T2,Op,B1,1,B22,L22>::apply(x,y,op);
    }
};

//
// The case for a single element.
// Just do it.
//

template<class T1, class T2, class Op, int B1, int B2>
struct TinyMatrixAssign<T1,T2,Op,B1,1,B2,1>
{
  static void apply(T1& x, const T2& y,Op op=Op())
    {
      op(TinyMatrixElem<T1,B1,B2>::get(x), TinyMatrixElem<T2,B1,B2>::get(y));
    }
};

//
// The case for a two by two block.
// Just do it.
//

template<class T1, class T2, class Op, int B1, int B2>
struct TinyMatrixAssign<T1,T2,Op,B1,2,B2,2>
{
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      op(TinyMatrixElem<T1,B1  ,B2  >::get(x), TinyMatrixElem<T2,B1  ,B2  >::get(y));
      op(TinyMatrixElem<T1,B1+1,B2  >::get(x), TinyMatrixElem<T2,B1+1,B2  >::get(y));
      op(TinyMatrixElem<T1,B1  ,B2+1>::get(x), TinyMatrixElem<T2,B1  ,B2+1>::get(y));
      op(TinyMatrixElem<T1,B1+1,B2+1>::get(x), TinyMatrixElem<T2,B1+1,B2+1>::get(y));
    }
};

//
// The case for a three by three block.
// Just do it.
//

template<class T1, class T2, class Op, int B1, int B2>
struct TinyMatrixAssign<T1,T2,Op,B1,3,B2,3>
{
  static void apply(T1& x, const T2& y, Op op=Op())
    {
      op(TinyMatrixElem<T1,B1  ,B2  >::get(x), TinyMatrixElem<T2,B1  ,B2  >::get(y));
      op(TinyMatrixElem<T1,B1+1,B2  >::get(x), TinyMatrixElem<T2,B1+1,B2  >::get(y));
      op(TinyMatrixElem<T1,B1+2,B2  >::get(x), TinyMatrixElem<T2,B1+2,B2  >::get(y));
      op(TinyMatrixElem<T1,B1  ,B2+1>::get(x), TinyMatrixElem<T2,B1  ,B2+1>::get(y));
      op(TinyMatrixElem<T1,B1+1,B2+1>::get(x), TinyMatrixElem<T2,B1+1,B2+1>::get(y));
      op(TinyMatrixElem<T1,B1+2,B2+1>::get(x), TinyMatrixElem<T2,B1+2,B2+1>::get(y));
      op(TinyMatrixElem<T1,B1  ,B2+2>::get(x), TinyMatrixElem<T2,B1  ,B2+2>::get(y));
      op(TinyMatrixElem<T1,B1+1,B2+2>::get(x), TinyMatrixElem<T2,B1+1,B2+2>::get(y));
      op(TinyMatrixElem<T1,B1+2,B2+2>::get(x), TinyMatrixElem<T2,B1+2,B2+2>::get(y));
    }
};



#endif // POOMA_TINY_TINYMATRIX_ELEMENTS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TinyMatrixElements.h,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/29 14:19:30 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
