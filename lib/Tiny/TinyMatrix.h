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
// Class: TinyMatrix, Full, TinyMatrixEngine<D,T,Full>
//-----------------------------------------------------------------------------

#ifndef POOMA_TINY_TINYMATRIX_H
#define POOMA_TINY_TINYMATRIX_H

/** @file
 * @ingroup Tiny
 * @brief
 * An interface class for an N-dimensional TinyMatrix of numeric objects,
 * and an engine class for defining a general TinyMatrix.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

class Full;

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Utilities/ElementProperties.h"
#include "PETE/PETE.h"
#include "Pooma/PoomaOperatorTags.h"
#include "Domain/Loc.h"
#include "Tiny/TinyMatrixElements.h"
#include "Tiny/TinyMatrixOperators.h"
#if POOMA_NO_IOS_HEADER
#include <iostream>
#else
#include <ios>
#endif

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D1, int D2, class T, class E> class TinyMatrix;
template<int D1, int D2, class T, class E> class TinyMatrixEngine;

template <class T>
void reverseBytes(T&);

/**
 * TinyMatrix is an interface class that takes three template parameters:
 * - int D1, int D2: The number of components in each rank of the TinyMatrix.
 * - class T: The type of the components.  
 * - class EngineTag: A policy parameter for the storage type.
 */

#if !POOMA_NO_DEPENDENT_TEMPLATE_ARGS
template<int D1, int D2=D1, class T=double, class EngineTag=Full>
#else
template<int D1, int D2, class T=double, class EngineTag=Full>
#endif
class TinyMatrix
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d1=D1, d2=D2 };
  typedef T Element_t;
  typedef EngineTag EngineTag_t;

  // Deduce the engine type from the tamplate parameters.
  typedef TinyMatrixEngine<D1,D2,T,EngineTag> Engine_t;

  // Return types for accessor functions.
  typedef typename Engine_t::ElementRef_t       ElementRef_t;
  typedef typename Engine_t::ConstElementRef_t  ConstElementRef_t;

  // Record the type of the current class.
  typedef TinyMatrix<D1,D2,T,EngineTag> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor uses the engine's null ctor.
  TinyMatrix() {}

  // Copy ctor is deep.
  TinyMatrix(const This_t& x) : engine_m(x.engine_m) {}

  // Construct from TinyMatrix with different type and engine.
  // Used as automatic expression expander.
  template<int D3, int D4, class T2, class EngineTag2>
  TinyMatrix(const TinyMatrix<D3, D4, T2, EngineTag2>& x) : engine_m(x) {}

  // Construct from an arbitrary single object.
  // The object must be indexable using TinyMatrixElem.
  template<class X>
  explicit TinyMatrix(const X& x) : engine_m(x) {}

  // Construct from two, three, ..., nine objects.
  template<class X1, class X2>
    TinyMatrix(const X1& x1, const X2& x2) : engine_m(x1,x2) {}
  template<class X1, class X2, class X3>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3) 
    : engine_m(x1,x2,x3) {}
  template<class X1, class X2, class X3, class X4>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3, const X4& x4) 
    : engine_m(x1,x2,x3,x4) {}
  template<class X1, class X2, class X3, class X4, class X5>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5) 
    : engine_m(x1,x2,x3,x4,x5) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6) 
    : engine_m(x1,x2,x3,x4,x5,x6) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6,
           class X7>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6, const X7& x7) 
    : engine_m(x1,x2,x3,x4,x5,x6,x7) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6,
           class X7, class X8>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6, const X7& x7, const X8& x8) 
    : engine_m(x1,x2,x3,x4,x5,x6,x7,x8) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6,
           class X7, class X8, class X9>
    TinyMatrix(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6, const X7& x7, const X8& x8, const X9& x9) 
    : engine_m(x1,x2,x3,x4,x5,x6,x7,x8,x9) {}

  // Let the engine destroy itself.
  ~TinyMatrix() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of TinyMatrix.
  This_t& operator=(const This_t& x) 
    { 
      if ( this != &x )
        engine() = x.engine();
      return *this;
    }

  // Assign from an arbitrary type.
  template<class V>
    This_t& operator=(const V& x)
        {
          engine() = x;
          return *this;
        }

  //----------------------------------------------------------------------
  // Element access

  ConstElementRef_t operator()(int i,int j) const 
    {
      return engine()(i,j); 
    }
  ElementRef_t operator()(int i,int j)
    {
      return engine()(i,j); 
    }

  ConstElementRef_t  operator()(int i) const { return engine()(i); }
  ElementRef_t       operator()(int i)       { return engine()(i); }


  // An accessor to get the engine.
  const Engine_t& engine() const { return engine_m; }
  Engine_t& engine() { return engine_m; }


  template<class Out>
  void print(Out &out) const;

  // This is only used when reading and writing data to disk
  inline void reverseBytes() { engine_m.reverseBytes(); }

private:

  // The only data is the engine itself.
  Engine_t engine_m;

};


/// Output to a stream.
/// The format is: ((t(0,0) t(0,1),... ) (t(1,0) t(1,1) ... ) ... ))

template<int D1, int D2, class T, class EngineTag>
template<class Out>
void TinyMatrix<D1, D2, T, EngineTag>::print(Out &out) const
{
  // Maintain the input formatting state through the multiple output
  // statements following:
  std::ios::fmtflags incomingFormatFlags = out.flags();
  long width = out.width();
  long precision = out.precision();
  out.width(0);
  out << "(";
  for (int i = 0; i < D1; i++) {
    out << "(";
    out.flags(incomingFormatFlags);
    out.width(width);
    out.precision(precision);
    out << (*this)(i, 0);
    for (int j = 1; j < D2; j++) {
      out << " ";
      out.flags(incomingFormatFlags);
      out.width(width);
      out.precision(precision);
      out << (*this)(i, j);
    }
    out << ")";
  }
  out << ")";
}


/// Output to a stream.
/// The format is: ( ( t(0,0) t(0,1),... ) ( t(1,0) t(1,1) ... ) ... )

template<int D1, int D2, class T, class E>
std::ostream &operator<<(std::ostream &out, const TinyMatrix<D1,D2,T,E> &t)
{
  t.print(out);
  return out;
}


//-----------------------------------------------------------------------------
// Specialization of ElementProperties struct for TinyMatrix.
//-----------------------------------------------------------------------------

template <int D1, int D2, class T, class E>
struct ElementProperties< TinyMatrix<D1,D2,T,E> > 
  : public TrivialElementProperties< TinyMatrix<D1,D2,T,E> >
{ };


/**
 * Definitions for a Full TinyMatrix.
 */

template<int D1, int D2, class T>
class TinyMatrixEngine<D1,D2,T,Full>
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d1=D1, d2=D2 };
  typedef T Element_t;
  typedef Full EngineTag_t;

  // Return types for accessor functions.
  typedef       T&  ElementRef_t;
  typedef const T&  ConstElementRef_t;

  // Record the type of the current class.
  typedef TinyMatrixEngine<D1,D2,T,Full> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor takes no action.
  TinyMatrixEngine()
  {
    CTAssert(ElementProperties<T>::hasTrivialDefaultConstructor
	     && ElementProperties<T>::hasTrivialDestructor
	     && ElementProperties<T>::concrete);
  }

  // Copy ctor is deep.
  TinyMatrixEngine(const TinyMatrixEngine<D1,D2,T,Full>& x)
    {
      TinyMatrixAssign<This_t,This_t,OpAssign,0,D1,0,D2>
	::apply(*this,x,OpAssign());
    }

  // Construct from an argument of arbitrary type.
  // The arg must be indexable using TinyMatrixElem.
  template<class X>
  explicit TinyMatrixEngine(const X& x)
      {
        TinyMatrixAssign<This_t,X,OpAssign,0,D1,0,D2>::apply(*this,x,OpAssign());
      }

  // Construct from two, three, ... nine objects.
  // Just fill in the 2D array using them.
  template<class X1, class X2>
    TinyMatrixEngine(const X1& x1, const X2& x2)
      {
        CTAssert( D1*D2 == 2 );
        (*this)(0) = x1;
        (*this)(1) = x2;
      }
  template<class X1, class X2, class X3>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3)
      {
        CTAssert( D1*D2 == 3 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
      }

  template<class X1, class X2, class X3, class X4>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4) 
      {
        CTAssert( D1*D2 == 4 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
        (*this)(3) = x4;
      }
    
  template<class X1, class X2, class X3, class X4, class X5>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5) 
      {
        CTAssert( D1*D2 == 5 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
        (*this)(3) = x4;
        (*this)(4) = x5;
      }
    
  template<class X1, class X2, class X3, class X4, class X5, class X6>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6) 
      {
        CTAssert( D1*D2 == 6 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
        (*this)(3) = x4;
        (*this)(4) = x5;
        (*this)(5) = x6;
      }
    
  template<class X1, class X2, class X3, class X4, class X5, class X6,
           class X7>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6, const X7& x7) 
      {
        CTAssert( D1*D2 == 7 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
        (*this)(3) = x4;
        (*this)(4) = x5;
        (*this)(5) = x6;
        (*this)(6) = x7;
      }
    
  template<class X1, class X2, class X3, class X4, class X5, class X6,
           class X7, class X8>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6, const X7& x7, const X8& x8) 
      {
        CTAssert( D1*D2 == 8 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
        (*this)(3) = x4;
        (*this)(4) = x5;
        (*this)(5) = x6;
        (*this)(6) = x7;
        (*this)(7) = x8;
      }
    
  template<class X1, class X2, class X3, class X4, class X5, class X6,
           class X7, class X8, class X9>
    TinyMatrixEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
           const X5& x5, const X6& x6, const X7& x7, const X8& x8, const X9& x9) 
      {
        CTAssert( D1*D2 == 9 );
        (*this)(0) = x1;
        (*this)(1) = x2;
        (*this)(2) = x3;
        (*this)(3) = x4;
        (*this)(4) = x5;
        (*this)(5) = x6;
        (*this)(6) = x7;
        (*this)(7) = x8;
        (*this)(8) = x9;
      }
    
  // Let the engine destroy itself.
  ~TinyMatrixEngine() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of TinyMatrix.
  const This_t&
    operator=(const This_t& x)
      {
        if ( this != &x )
          TinyMatrixAssign<This_t,This_t,OpAssign,0,D1,0,D2>
	    ::apply(*this,x,OpAssign());
        return *this;
      }

  // Assign from an arbitrary type.
  template<class V>
    const This_t& 
      operator=(const V& x)
        {
          TinyMatrixAssign<This_t,V,OpAssign,0,D1,0,D2>::apply(*this,x,OpAssign());
          return *this;
        }

  //----------------------------------------------------------------------
  // Element access

  ConstElementRef_t  operator()(int i,int j) const
    {
      PBoundAssert((i>=0)&&(i<D1)
		   && (j>=0)&&(j<D2));
      return x_m[i+D1*j]; 
    }
  ElementRef_t operator()(int i,int j)
    {
      PBoundAssert((i>=0)&&(i<D1)
		   && (j>=0)&&(j<D2));
      return x_m[i+D1*j]; 
    }

  ConstElementRef_t  operator()(int i) const 
    {
      PBoundAssert((i>=0)&&(i<D1*D2));
      return x_m[i]; 
    }
  ElementRef_t operator()(int i)
    {
      PBoundAssert((i>=0)&&(i<D1*D2));
      return x_m[i]; 
    }

  // This is only supposed to be used in the IO stuff, which should
  // make the inline function reverseBytes visable. 

  inline void reverseBytes() 
  { 
    for (int d = 0; d < D1*D2; ++d) ::reverseBytes(x_m[d]);
  }

private:

  // The actual data is just an array of T's.
  // Store in fortran order.
  T x_m[D1*D2];
};


//-----------------------------------------------------------------------------
// ComponentAccess is an interface class that is used to provide an API for
// accessing components of a composite type. This version works with TinyMatrixs.
//-----------------------------------------------------------------------------

template<class T, class Components> struct ComponentAccess;

template<int D1, int D2, class T, class E, int N>
struct ComponentAccess< TinyMatrix<D1, D2, T, E>, Loc<N> >
{
  typedef TinyMatrix<D1, D2, T, E> V;
  typedef typename V::Element_t Element_t;
  typedef typename V::ElementRef_t ElementRef_t;
  
  static inline ElementRef_t indexRef(V &t, const Loc<N> &l)
    {
      CTAssert(N==2);
      return t(l[0].first(), l[1].first());
    }
  
  static inline Element_t index(const V &t, const Loc<N> &l)
    {
      CTAssert(N==2);
      return t(l[0].first(), l[1].first());
    }
};


//-----------------------------------------------------------------------------
// TinyMatrixElem specialization for Full TinyMatrix engines:
//-----------------------------------------------------------------------------

template<int D1, int D2, class T, int I, int J>
struct TinyMatrixElem< TinyMatrixEngine<D1,D2,T,Full> , I , J>
{
  typedef TinyMatrixEngine<D1,D2,T,Full> V;
  typedef TinyMatrixEngineElem<D1,D2,T,Full,I,J> TE;
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x); }
  static      ElementRef_t get(V& x)       { return TE::get(x); }
};

#endif // POOMA_TINY_TINYMATRIX_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TinyMatrix.h,v $   $Author: richi $
// $Revision: 1.23 $   $Date: 2004/11/26 13:00:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
