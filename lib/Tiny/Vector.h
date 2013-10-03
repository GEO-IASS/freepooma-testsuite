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
// Class: Vector, Full, VectorEngine<D,T,Full>
//-----------------------------------------------------------------------------

#ifndef POOMA_TINY_VECTOR_H
#define POOMA_TINY_VECTOR_H

/** @file
 * @ingroup Tiny
 * @brief
 * An interface class for an N-dimensional vector of numeric objects,
 * and an engine class for defining a general vector.
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
#include "Tiny/VectorElements.h"
#include "Tiny/VectorOperators.h"
#if POOMA_NO_IOS_HEADER
#include <iostream>
#else
#include <ios>
#endif

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class E> class Vector;
template<int D, class T, class E> class VectorEngine;

template <class T>
void reverseBytes(T&);

/**
 * Vector is an interface class that takes three template parameters:
 * - int Dim: The number of components in the vector.
 * - class T: The type of the components.  
 * - class EngineTag: A policy parameter for the storage type.
 */

template<int Dim, class T=double, class EngineTag=Full>
class Vector
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=1 };
  enum { d1 = Dim };
  typedef T Element_t;
  typedef EngineTag EngineTag_t;

  // Deduce the engine type from the tamplate parameters.
  typedef VectorEngine<Dim,T,EngineTag> Engine_t;

  // Return types for accessor functions.
  typedef typename Engine_t::ElementRef_t       ElementRef_t;
  typedef typename Engine_t::ConstElementRef_t  ConstElementRef_t;

  // Record the type of the current class.
  typedef Vector<Dim,T,EngineTag> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor uses the engine's null ctor.
  Vector() {}

  // Copy ctor is deep.
  Vector(const This_t& x) : engine_m(x.engine_m) {}

  // Construct from Vector with different type and engine.
  // Used as automatic expression expander.
  template<int D2, class T2, class EngineTag2>
  Vector(const Vector<D2, T2, EngineTag2>& x) : engine_m(x) {}

  // Construct from an arbitrary single object.
  // The object must be indexable using VectorElem.
  template<class X>
  explicit Vector(const X& x) : engine_m(x) {}

  // Construct from 2-7 scalars.
  template<class X1, class X2>
  Vector(const X1& x, const X2& y)
    : engine_m(x,y) {}
  template<class X1, class X2, class X3>
  Vector(const X1& x, const X2& y, const X3& z)
    : engine_m(x,y,z) {}
  template<class X1, class X2, class X3, class X4>
  Vector(const X1& x, const X2& y, const X3& z, const X4& a)
    : engine_m(x,y,z,a) {}
  template<class X1, class X2, class X3, class X4, class X5>
  Vector(const X1& x, const X2& y, const X3& z, const X4& a,
	 const X5& b)
    : engine_m(x,y,z,a,b) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  Vector(const X1& x, const X2& y, const X3& z, const X4& a,
	 const X5& b, const X6& c)
    : engine_m(x,y,z,a,b,c) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7>
  Vector(const X1& x, const X2& y, const X3& z, const X4& a,
	 const X5& b, const X6& c, const X7& d)
    : engine_m(x,y,z,a,b,c,d) {}

  // Let the engine destroy itself.
  ~Vector() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of vector.
  This_t& operator=(const This_t& x) 
    { 
      if ( this != &const_cast<This_t &>(x) )
        engine() = x.engine();
      return *this;
    }

  // Assign from an arbitrary type.
  // The engine should recognize if it is given a vector versus a scalar.
  template<class V>
    This_t&
      operator=(const V& x)
        {
          engine() = x;
          return *this;
        }

  //----------------------------------------------------------------------
  // Element access, bounds checking is done in engine.

  ConstElementRef_t  operator()(int i) const 
    {
      return engine()(i); 
    }
  ElementRef_t operator()(int i) 
    { 
      return engine()(i); 
    }

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
/// The format is: (v(0),v(1),...,v(D-1))

template<int Dim, class T, class EngineTag>
template<class Out>
void Vector<Dim, T, EngineTag>::print(Out &out) const
{
  // Maintain the input formatting state through the multiple output
  // statements following:
  std::ios::fmtflags incomingFormatFlags = out.flags();
  long width = out.width();
  long precision = out.precision();
  out.width(0);
  out << "(";
  out.flags(incomingFormatFlags);
  out.width(width);
  out.precision(precision);
  out << (*this)(0) ;
  for (int i = 1; i < Dim; i++) {
    out << ",";
    out.flags(incomingFormatFlags);
    out.width(width);
    out.precision(precision);
    out << (*this)(i);
  }
  out << ")";
}


/// Output to a stream.
/// The format is: (v(0),v(1),...,v(D-1))

template<int D, class T, class E>
std::ostream &operator<<(std::ostream &out, const Vector<D,T,E> &v)
{
  v.print(out);
  return out;
}

//-----------------------------------------------------------------------------
// Specialization of ElementProperties struct for Vector.
//-----------------------------------------------------------------------------

template <int D, class T, class E>
struct ElementProperties< Vector<D,T,E> > 
  : public TrivialElementProperties< Vector<D,T,E> >
{ };

/**
 * Definitions for a Full vector.
 */

template<int D, class T>
class VectorEngine<D,T,Full>
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=1 };
  enum { d1 = D };
  typedef T Element_t;
  typedef Full EngineTag_t;

  // Return types for accessor functions.
  typedef       T&  ElementRef_t;
  typedef const T&  ConstElementRef_t;

  // Record the type of the current class.
  typedef VectorEngine<D,T,Full> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor takes no action.
  VectorEngine() 
  {
    CTAssert(ElementProperties<T>::hasTrivialDefaultConstructor
	     && ElementProperties<T>::hasTrivialDestructor
	     && ElementProperties<T>::concrete);

    for (int i = 0; i < D; ++i)
      {
	ElementProperties<T>::construct(&x_m[i]);
      }
  }

  // Copy ctor is deep.
  inline VectorEngine(const VectorEngine<D,T,Full>&);

  // Copy from an argument of arbitrary type.
  // The arg must be indexable using VectorElem.
  template<class X>
    inline explicit VectorEngine(const X& x)
      {
        VectorAssign< VectorEngine<D,T,Full> , X , OpAssign, 0, D >
          ::apply(*this,x,OpAssign());
      }
    
  // Construct from two or three scalars.
  template<class X1, class X2>
  inline VectorEngine(const X1& x, const X2& y)
  {
    CTAssert(D == 2);
    x_m[0] = x;
    x_m[1] = y;
  }
  template<class X1, class X2, class X3>
  inline VectorEngine(const X1& x, const X2& y, const X3& z)
  {
    CTAssert(D == 3);
    x_m[0] = x;
    x_m[1] = y;
    x_m[2] = z;
  }
  template<class X1, class X2, class X3, class X4>
  inline VectorEngine(const X1& x, const X2& y, const X3& z, 
		      const X4& a)
  {
    CTAssert(D == 4);
    x_m[0] = x;
    x_m[1] = y;
    x_m[2] = z;
    x_m[3] = a;
  }
  template<class X1, class X2, class X3, class X4, class X5>
  inline VectorEngine(const X1& x, const X2& y, const X3& z, 
		      const X4& a, const X5& b)
  {
    CTAssert(D == 5);
    x_m[0] = x;
    x_m[1] = y;
    x_m[2] = z;
    x_m[3] = a;
    x_m[4] = b;
  }
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  inline VectorEngine(const X1& x, const X2& y, const X3& z, 
		      const X4& a, const X5& b, const X6& c)
  {
    CTAssert(D == 6);
    x_m[0] = x;
    x_m[1] = y;
    x_m[2] = z;
    x_m[3] = a;
    x_m[4] = b;
    x_m[5] = c;
  }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7>
  inline VectorEngine(const X1& x, const X2& y, const X3& z, 
		      const X4& a, const X5& b, const X6& c, const X7& d)
  {
    CTAssert(D == 7);
    x_m[0] = x;
    x_m[1] = y;
    x_m[2] = z;
    x_m[3] = a; 
    x_m[4] = b;
    x_m[5] = c;
    x_m[6] = d;
  }

  // Let the engine destroy itself.
  ~VectorEngine() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of vector.
  This_t&
    operator=(const This_t& x)
      {
        if ( this != &x )
          VectorAssign<This_t,This_t,OpAssign,0,D>::apply(*this,x,OpAssign());
        return *this;
      }

  // Assign from an arbitrary type.
  // The engine should recognize if it is given a vector versus a scalar.
  template<class V>
    This_t& 
      operator=(const V& x)
        {
          VectorAssign<This_t,V,OpAssign,0,D>::apply(*this,x,OpAssign());
          return *this;
        }

  //----------------------------------------------------------------------
  // Element access

  ConstElementRef_t operator()(int i) const 
    {
      PBoundAssert((i>=0)&&(i<D));
      return x_m[i]; 
    }
  ElementRef_t operator()(int i)
    {
      PBoundAssert((i>=0)&&(i<D));
      return x_m[i]; 
    }

  // This is only supposed to be used in the IO stuff, which should
  // make the inline function reverseBytes visible. 

  inline void reverseBytes() 
  { 
    for (int d = 0; d < D; ++d) ::reverseBytes(x_m[d]);
  }

private:

  // The actual data is just an array of T's of length D.
  T x_m[D];
};


//----------------------------------------------------------------------
// VectorElem class for the Full engine.
//----------------------------------------------------------------------

template<int D, class T, int I>
struct VectorElem< VectorEngine<D,T,Full> , I >
{
  typedef VectorEngine<D,T,Full> V;
  typedef VectorEngineElem<D,T,Full,I> VE;
  typedef typename VE::Element_t Element_t;
  typedef typename VE::ConstElementRef_t ConstElementRef_t;
  typedef typename VE::ElementRef_t ElementRef_t;
  static ConstElementRef_t get(const V& x) { return VE::get(x); }
  static ElementRef_t get(V& x) { return VE::get(x); }
};


//-----------------------------------------------------------------------------
/// ComponentAccess is an interface class that is used to provide an API for
/// accessing components of a composite type. This version works with Vectors.
//-----------------------------------------------------------------------------

template<class T, class Components> struct ComponentAccess;

template<int D, class T, class E, int N>
struct ComponentAccess< Vector<D, T, E>, Loc<N> >
{
  typedef Vector<D, T, E> V;
  typedef typename V::Element_t Element_t;
  typedef typename V::ElementRef_t ElementRef_t;
  
  static inline ElementRef_t indexRef(V &v, const Loc<N> &l)
    {
      CTAssert(N==1);
      return v(l[0].first());
    }
  
  static inline Element_t index(const V &v, const Loc<N> &l)
    {
      CTAssert(N==1);
      return v(l[0].first());
    }
};


//----------------------------------------------------------------------
// Copy ctor uses VectorAssign to copy the elements of x into this.
//----------------------------------------------------------------------

template<int D, class T>
inline 
VectorEngine<D,T,Full>::VectorEngine(const VectorEngine<D,T,Full>& x)
{
  VectorAssign<This_t,This_t,OpAssign,0,D>::apply(*this,x,OpAssign());
}

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Vector.h,v $   $Author: richi $
// $Revision: 1.40 $   $Date: 2004/11/26 13:00:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
