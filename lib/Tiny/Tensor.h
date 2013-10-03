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
// 
//  Tensor, 
//  Full, Antisymmetric, Symmetric, Diagonal, 
//  TensorEngine<D,T,Full>, TensorEngine<D,T,Antisymmetric>
//  TensorEngine<D,T,Symmetric>, TensorEngine<D,T,Diagonal> 
// 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 
// Global Function Templates:
// 
// symmetrize() : Takes a Tensor as input, and converts it to a Tensor with
//                the requested output symmetry.
//-----------------------------------------------------------------------------

#ifndef POOMA_TINY_TENSOR_H
#define POOMA_TINY_TENSOR_H

/** @file
 * @ingroup Tiny
 * @brief
 * An interface class for an N-dimensional tensor of numeric objects,
 * and engines class for defining a general tensor, using Full and 
 * Antisymmetric engine tag classes
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Utilities/ElementProperties.h"
#include "PETE/PETE.h"
#include "Pooma/PoomaOperatorTags.h"
#include "Domain/Loc.h"
#include "Tiny/TensorElements.h"
#include "Tiny/TensorOperators.h"
#include <iosfwd>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int D, class T, class EngineTag> class Tensor;
template<int D, class T, class EngineTag> class TensorEngine;

template <class T>
void reverseBytes(T&);

//-----------------------------------------------------------------------------
// Classes to be used for the "EngineTag" template parameter of Tensor:
//-----------------------------------------------------------------------------

/** @class Full 
 * Full storage (general tensor):
 * 
 * Like this, for the 3D case; just chop off pieces for 1D and 2D:
 * <PRE>
 *      Tensor Structure       (i,j) Indices      array storage of elements
 *      -----------------    -----------------    --------------------------
 *      | x00  x01  x02 |    | 0,0  0,1  0,2 |    | x_m[0]  x_m[3]  x_m[6] |
 *      | x10  x11  x12 |    | 1,0  1,1  1,2 |    | x_m[1]  x_m[4]  x_m[7] |
 *      | x20  x21  x22 |    | 2,0  2,1  2,2 |    | x_m[2]  x_m[5]  x_m[8] |
 * </PRE>
 */
class Full;

/** @class Antisymmetric
 * Like this, for the 3D case; just chop off pieces for 1D and 2D:
 * <PRE>
 *      Tensor Structure       (i,j) Indices      array storage of elements 
 *      -----------------    -----------------    --------------------------
 *      |  O  -x10 -x20 |    | 0,0  0,1  0,2 |    |                        |
 *      | x10   0  -x21 |    | 1,0  1,1  1,2 |    | x_m[0]                 |
 *      | x20  x21   0  |    | 2,0  2,1  2,2 |    | x_m[1]  x_m[2]         |
 * </PRE>
 */
class Antisymmetric;

/** @class Symmetric
 * Like this, for the 3D case:
 * <PRE>
 *      Tensor Structure       (i,j) Indices      array storage of elements
 *      -----------------    -----------------    --------------------------
 *      | xOO  x10  x20 |    | 0,0  0,1  0,2 |    | x_m[0]                 |
 *      | x10  x11  x21 |    | 1,0  1,1  1,2 |    | x_m[1]  x_m[2]         |
 *      | x20  x21  x22 |    | 2,0  2,1  2,2 |    | x_m[3]  x_m[4]  x_m[5] |
 * </PRE>
 */
class Symmetric;

/** @class Diagonal
 * Like this, for the 3D case; just chop off pieces for 1D and 2D:
 * <PRE>
 *      Tensor Structure       (i,j) Indices      array storage of elements 
 *      -----------------    -----------------    --------------------------
 *      | xOO   0    0  |    | 0,0  0,1  0,2 |    | x_m[0]                 |
 *      |  0   x11   0  |    | 1,0  1,1  1,2 |    |         x_m[1]         |
 *      |  0    0   x22 |    | 2,0  2,1  2,2 |    |                 x_m[2] |
 * </PRE>
 */
class Diagonal;


// ----------------------------------------------------------------------------
// Sizes of linear arrays used to store the elements, for the various
// EngineTag types:
// ----------------------------------------------------------------------------

// General template: leave empty to trigger compile errors for unsupported
// EngineTag types and/or dimensionalities:

template<int D, class EngineTag>
class TensorStorageSize {};

// These partial specializations spell out everything for the supported
// EngineTag types and dimensionalities (arbitrary dimensionality, note):

// Full:
template<int D>
class TensorStorageSize<D, Full>
{
public:
  enum { Size = D*D };
};

// Antisymmetric:
template<int D>
class TensorStorageSize<D, Antisymmetric>
{
public:
  enum { Size = (D*D - D)/2 + 1/D };
};

// Symmetric:
template<int D>
class TensorStorageSize<D, Symmetric>
{
public:
  enum { Size = (D*D - D)/2 + D };
};

// Diagonal:
template<int D>
class TensorStorageSize<D, Diagonal>
{
public:
  enum { Size = D };
};


/**
 * Tensor is an interface class that takes three template parameters:
 *   - int D: The number of components in each rank (row or col) of the Tensor.
 *          That is, the tensor is a DxD object. For a non-square tiny object
 *          of this nature, which wouldn't be a tensor mathematically, use
 *          the TinyMatrix class.
 *   - class T: The type of the components.  
 *   - class EngineTag: A policy parameter for the storage type.
 */

template<int D, class T=double, class EngineTag=Full>
class Tensor
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d=D }; // needed, still??? (TJW)  
  typedef T Element_t;
  typedef EngineTag EngineTag_t;

  // Deduce the engine type from the tamplate parameters.
  typedef TensorEngine<D,T,EngineTag> Engine_t;

  // Return types for accessor functions.
  typedef typename Engine_t::ElementRef_t       ElementRef_t;
  typedef typename Engine_t::ConstElementRef_t  ConstElementRef_t;

  // Record the type of the current class.
  typedef Tensor<D,T,EngineTag> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor uses the engine's null ctor.
  Tensor() {}

  // Copy ctor is deep.
  Tensor(const This_t& x) : engine_m(x.engine_m) {}

  // Construct from Tensor with different type and engine.
  // Used as automatic expression expander.
  template<int D2, class T2, class EngineTag2>
  Tensor(const Tensor<D2, T2, EngineTag2>& x) : engine_m(x) {}

  // Construct from an arbitrary single object.
  // The object must be indexable using TensorElem.
  template<class X>
  explicit Tensor(const X& x) : engine_m(x) {}

  // Construct from two, three, ..., nine objects.
  template<class X1, class X2>
  Tensor(const X1& x, const X2& y) : engine_m(x,y) {}
  template<class X1, class X2, class X3>
  Tensor(const X1& x, const X2& y, const X3& z) : engine_m(x,y,z) {}
  template<class X1, class X2, class X3, class X4>
  Tensor(const X1& x1, const X2& x2, const X3& x3, const X4& x4) 
    : engine_m(x1,x2,x3,x4) {}
  template<class X1, class X2, class X3, class X4, class X5>
  Tensor(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
         const X5& x5) 
    : engine_m(x1,x2,x3,x4,x5) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  Tensor(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
         const X5& x5, const X6& x6) 
    : engine_m(x1,x2,x3,x4,x5,x6) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6,
    class X7>
  Tensor(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
         const X5& x5, const X6& x6, const X7& x7) 
    : engine_m(x1,x2,x3,x4,x5,x6,x7) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6,
    class X7, class X8>
  Tensor(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
         const X5& x5, const X6& x6, const X7& x7, const X8& x8) 
    : engine_m(x1,x2,x3,x4,x5,x6,x7,x8) {}
  template<class X1, class X2, class X3, class X4, class X5, class X6,
    class X7, class X8, class X9>
  Tensor(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
         const X5& x5, const X6& x6, const X7& x7, const X8& x8, const X9& x9) 
    : engine_m(x1,x2,x3,x4,x5,x6,x7,x8,x9) {}

  // Let the engine destroy itself.
  ~Tensor() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of Tensor.
  This_t& operator=(const This_t& x) 
  { 
    if ( this != &x )
      engine() = x.engine();
    return *this;
  }

  // Assign from the another kind of Tensor.
  template<class T1, class EngineTag1>
  This_t& operator=(const Tensor<d,T1,EngineTag1> &x) 
  { 
    engine() = x.engine();
    return *this;
  }

  // Assign from an arbitrary type.
  template<class V>
  This_t&
  operator=(const V& x)
  {
    engine() = x;
    return *this;
  }

  //----------------------------------------------------------------------
  // Element access

  // Compile-time indices: only have this in partial specializations, as some
  // engines won't define this functionality.

  // Runtime indices:

  ConstElementRef_t operator()(int i, int j) const 
  {
    return engine()(i,j); 
  }
  ElementRef_t operator()(int i, int j)
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
/// The format is: ((t(0,0) t(0,1),... ) ( t(1,0) t(1,1) ... ) ... ))

template<int D, class T, class EngineTag>
template<class Out>
void Tensor<D, T, EngineTag>::print(Out &out) const
{
  // Maintain the input formatting state through the multiple output
  // statements following:
  std::ios::fmtflags incomingFormatFlags = out.flags();
  long width = out.width();
  long precision = out.precision();
  out.width(0);
  out << "(";
  for (int i = 0; i < D; i++) {
    out << "(";
    out.flags(incomingFormatFlags);
    out.width(width);
    out.precision(precision);
    out << (*this)(i,0);
    for (int j = 1; j < D; j++) {
      out << " ";
      out.flags(incomingFormatFlags);
      out.width(width);
      out.precision(precision);
      out << (*this)(i,j);
    }
    out << ")";
  }
  out << ")";
}


/// Output to a stream.
/// The format is: ( ( t(0,0) t(0,1),... ) ( t(1,0) t(1,1) ... ) ... )

template<int D, class T, class E>
std::ostream &operator<<(std::ostream &out, const Tensor<D,T,E> &t)
{
  t.print(out);
  return out;
}


//-----------------------------------------------------------------------------
// Specialization of ElementProperties struct for Tensor.
//-----------------------------------------------------------------------------

template <int D, class T, class E>
struct ElementProperties< Tensor<D,T,E> > 
  : public TrivialElementProperties< Tensor<D,T,E> >
{ };


/**
 * TensorEngine definitions for a Full Tensor.
 */

template<int D, class T>
class TensorEngine<D, T, Full>
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d=D };
  typedef T Element_t;
  typedef Full EngineTag_t;

  // Return types for accessor functions.

  // Runtime indices:
  typedef       T&  ElementRef_t;
  typedef const T&  ConstElementRef_t;

  // Compile-time indices:
  typedef        T&  CTElementRef_t;
  typedef  const T&  CTConstElementRef_t;

  // Record the type of the current class.
  typedef TensorEngine<D,T,Full> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor takes no action.
  TensorEngine()
  {
    CTAssert(ElementProperties<T>::hasTrivialDefaultConstructor
	     && ElementProperties<T>::hasTrivialDestructor
	     && ElementProperties<T>::concrete);
  }

  // Copy ctor is deep.
  TensorEngine(const TensorEngine<D,T,Full>& x)
  {
    TensorAssign<This_t,This_t,OpAssign,0,D,0,D>
      ::apply(*this,x,OpAssign());
  }

  // Construct from an argument of arbitrary type.
  // The arg must be indexable using TensorElem.
  template<class X>
  explicit TensorEngine(const X& x)
  {
    TensorAssign<This_t,X,OpAssign,0,D,0,D>::apply(*this,x,OpAssign());
  }

  // Construct a Tensor from 1 object of type T:
  explicit TensorEngine(const T& x)
  {
    for (int i = 0 ; i < D*D ; i++) {
      (*this)(i) = x;
    }
  }

  // Construct a 2D Tensor from 4 general objects:
  template<class X1, class X2, class X3, class X4>
  TensorEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4) 
  {
    CTAssert( D == 2 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
  }
    
  // Construct a 3D Tensor from 9 general objects:
  template<class X1, class X2, class X3, class X4, class X5, class X6,
    class X7, class X8, class X9>
  TensorEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
               const X5& x5, const X6& x6, const X7& x7, const X8& x8, 
               const X9& x9) 
  {
    CTAssert( D == 3 );
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
  ~TensorEngine() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of Tensor.
  This_t&
  operator=(const This_t& x)
  {
    if (this != &x) {
      TensorAssign<This_t,This_t,OpAssign,0,D,0,D>::apply(*this,x,OpAssign());
    }
    return *this;
  }

  // Assign from an arbitrary type.
  template<class V>
  This_t& 
  operator=(const V& x)
  {
    TensorAssign<This_t,V,OpAssign,0,D,0,D>::apply(*this,x,OpAssign());
    return *this;
  }

  //----------------------------------------------------------------------
  // Element access

  // Compile-time indices, used by things like TensorAssign:

  // Const:
  template<int I, int J>
  CTConstElementRef_t getIJ() const
  {
    CTAssert((I >= 0) && (I < D)
	     && (J >= 0) && (J < D));
    return x_m[I + D*J]; 
  }
    
  // Non-const:
  template<int I, int J>
  CTElementRef_t getIJ()
  {
    CTAssert((I >= 0) && (I < D)
	     && (J >= 0) && (J < D));
    return x_m[I + D*J]; 
  }

  // Runtime indices:

  // (i,j) indexing:

  ConstElementRef_t  operator()(int i,int j) const
  {
    PBoundAssert((i >= 0) && (i < D)
		 && (j >= 0) && (j < D));
    return x_m[i + D*j]; 
  }
  ElementRef_t operator()(int i,int j)
  {
    PBoundAssert((i >= 0) && (i < D)
		 && (j >= 0) && (j < D));
    return x_m[i + D*j]; 
  }

  // (i) indexing (directly into linear array x_m):
  ConstElementRef_t  operator()(int i) const 
  {
    PBoundAssert((i >= 0) && (i < D*D));
    return x_m[i]; 
  }
  ElementRef_t operator()(int i)
  {
    PBoundAssert((i >= 0) && (i < D*D));
    return x_m[i]; 
  }

//-----------------------------------------------------------------------------
// For direct access to the data array.

  const T* data() const{
    return (const T*) &x_m[0];
  }

  T* data(){
    return (T*) &x_m[0];
  }

  // This is only used when reading and writing data to disk
  inline void reverseBytes() 
  {
    const int sz = TensorStorageSize<d, EngineTag_t>::Size;
    for (int i = 0; i < sz; ++i)
      ::reverseBytes(x_m[i]);
  }

private:

  // The actual data is just an array of T's. Store in fortran order.
  T x_m[TensorStorageSize<d, EngineTag_t>::Size];
};


/**
 * TensorEngine definitions for an Antisymmetric Tensor.
 */

template<int D, class T>
class TensorEngine<D, T, Antisymmetric>
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d=D };
  typedef T Element_t;
  typedef Antisymmetric EngineTag_t;

  // ----------------------------------------------------------------------
  // The AssignProxy class is needed for non-const Tensor<Antisymmetric>
  // operator(i,j) implementation:

  class AssignProxy {
  public:
    AssignProxy(Element_t &elem, int where)
      : elem_m(elem), where_m(where) { }
    AssignProxy(const AssignProxy &model)
      : elem_m(model.elem_m), where_m(where_m) { }
    AssignProxy &operator=(const AssignProxy &a) const
    {
      PAssert(where_m != 0 || a.elem_m == -a.elem_m);
      elem_m = where_m < 0 ? -a.elem_m : a.elem_m;
      return const_cast<AssignProxy &>(*this);
    }
    AssignProxy &operator=(const Element_t &e) const
    {
      PAssert(where_m != 0 || e == -e);
      elem_m = where_m < 0 ? -e : e;
      return const_cast<AssignProxy &>(*this);
    }
    operator Element_t() const
    {
      return (where_m < 0 ? -elem_m : elem_m);
    }
  private:
    mutable Element_t &elem_m;
    mutable int where_m;
  };

  // Return types for accessor functions.

  // Runtime indices:
  // For reading/writing via runtime operator(i,j), need to use proxy object:
  typedef AssignProxy  ElementRef_t;
  // Really return-by-value for Antisym; const qualifier no longer appropriate:
  typedef T ConstElementRef_t;

  // Compile-time indices:
  typedef T& CTElementRef_t;
  // Really return-by-value for Antisym; const qualifier no longer appropriate:
  typedef T  CTConstElementRef_t;

  // Record the type of the current class.
  typedef TensorEngine<d,T,Antisymmetric> This_t;

  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor takes no action.
  TensorEngine()
  {
    CTAssert(ElementProperties<T>::hasTrivialDefaultConstructor
	     && ElementProperties<T>::hasTrivialDestructor
	     && ElementProperties<T>::concrete);
  }

  // Copy ctor is deep.
  TensorEngine(const TensorEngine<d,T,Antisymmetric> &x) 
  {
    TensorAssign<This_t,This_t,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
  }

  // Construct from an argument of arbitrary type.
  // The arg must be indexable using TensorElem.
  template<class X>
  explicit TensorEngine(const X& x)
  {
    TensorAssign<This_t,X,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
  }

  // Construct from one argument of type T. Assign this to all the actually
  // stored (not computed) values in the antisymmetric tensor:
  explicit TensorEngine(const T &x) { 
    for (int i = 0; i < TensorStorageSize<d, EngineTag_t>::Size; i++) {
      x_m[i] = x; // N.B.: may have to move this out to allow for D=1 spec.
    }
  }

  // Construct a 3D Antisymmetric Tensor from 3 general objects:
  template<class X1, class X2, class X3>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3) 
  {
    CTAssert( D == 3 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
  }
    
  // Construct a 4D Antisymmetric Tensor from 6 general objects:
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  TensorEngine(const X1& x1, const X2& x2, const X3& x3, const X4& x4,
               const X5& x5, const X6& x6)
  {
    CTAssert( D == 4 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
    (*this)(4) = x5;
    (*this)(5) = x6;
  }

  // Let the engine destroy itself.
  ~TensorEngine() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of Tensor.
  This_t&
  operator=(const This_t& x)
  {
    if (this != &x) {
      TensorAssign<This_t,This_t,OpAssign,0,d,0,d>::
        apply(*this,x,OpAssign());
    }
    return *this;
  }

  // Assign from an arbitrary type.
  template<class V>
  This_t& 
  operator=(const V& x)
  {
    TensorAssign<This_t,V,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
    return *this;
  }

  //----------------------------------------------------------------------
  // Element access; use the antisymmetry to compute any value except the 
  // stored value(s):

  // Compile-time indices, used by things like TensorAssign:

  // Const; really returning by value:
  template<int I, int J>
  CTConstElementRef_t getIJ() const
  {
    CTAssert((I >= 0) && (I < d)
	     && (J >= 0) && (J < d));
    if (I == J) {
      return This_t::Zero;
    } else {
      int lo = I < J ? I : J;
      int hi = I > J ? I : J;
      int symmetrySign  = I < J ? -1 : 1;
      return symmetrySign * x_m[((hi - 1)*hi/2) + lo];
    }
    // N.B.: the following code would be faster, but generates spurious
    // compiler warnings from KCC:
    //     } else if (I < J) {
    //       return -x_m[((J - 1)*J/2) + I];
    //     } else {
    //       return  x_m[((I - 1)*I/2) + J];
    //     }
  }
    
  // Non-const. Assume only stored (I,J) values, not computed values, will
  // ever be requested:
  template<int I, int J>
  CTElementRef_t getIJ()
  {
    CTAssert((I >= 0) && (I < d)
	     && (J >= 0) && (J < d)
	     && (I > J));
    return x_m[((I - 1)*I/2) + J];
  }

  // Runtime indices:

  // (i,j) indexing:
  ConstElementRef_t  operator()(int i,int j) const
  {
    PBoundAssert((i >= 0) && (i < d)
		 && (j >= 0) && (j < d));
    if (i == j) {
      return This_t::Zero;
    } else if (i < j) {
      return -x_m[((j - 1)*j/2) + i];
    } else {
      return  x_m[((i - 1)*i/2) + j];
    }
  }
  ElementRef_t operator()(int i,int j)
  {
    PBoundAssert((i >= 0) && (i < d)
		 && (j >= 0) && (j < d));
    if (i == j) {
      return AssignProxy(This_t::Zero, 0);
    } else {
      int lo = i < j ? i : j;
      int hi = i > j ? i : j;
      return AssignProxy(x_m[((hi-1)*hi/2) + lo], i - j);
    }
  }

  // (i) indexing (directly into linear array x_m):
  ConstElementRef_t  operator()(int i) const 
  {
    PBoundAssert((i >= 0) && (i < TensorStorageSize<d, EngineTag_t>::Size));
    return x_m[i]; 
  }
  ElementRef_t operator()(int i)
  {
    PBoundAssert((i >= 0) && (i < TensorStorageSize<d, EngineTag_t>::Size));
    return AssignProxy(x_m[i], 1); 
  }

//-----------------------------------------------------------------------------
// For direct access to the data array.

  const T* data() const{
    return (const T*) &x_m[0];
  }

  T* data(){
    return (T*) &x_m[0];
  }

  // This is only used when reading and writing data to disk
  inline void reverseBytes() 
  {
    const int sz = TensorStorageSize<d, EngineTag_t>::Size;
    for (int i = 0; i < sz; ++i)
      ::reverseBytes(x_m[i]);
  }

private:

  // The actual data is just an array of T's. Store in fortran order.
  T x_m[TensorStorageSize<d, EngineTag_t>::Size];

  // A place to store a zero element.  We need to return a reference to this
  // for the diagonal element:
  static T Zero;
};


// Assign the static zero element value:
template<int D, class T>
T TensorEngine<D,T,Antisymmetric>::Zero = 0.0;


/**
 * Special antisymmetric assignment class: Has specializations for different
 * dimensionalities (for 1, 2, and 3, so far).  This may ultimately be
 * replaceable with a dimensionality-independent equivalent that uses template
 * metaprogramming.
 */

// 1D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<1,T,Antisymmetric>,T2,Op,0,1,0,1>
{
  static void apply(TensorEngine<1,T,Antisymmetric> &x, const T2 &y, 
                    Op op=Op())
  { }
};
// 2D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<2,T,Antisymmetric>,T2,Op,0,2,0,2>
{
  static void apply(TensorEngine<2,T,Antisymmetric> &x, const T2 &y, 
                    Op op=Op())
  {
    TensorAssign<TensorEngine<2,T,Antisymmetric>,T2,Op,1,1,0,1>::apply(x,y,op);
  }
};
// 3D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<3,T,Antisymmetric>,T2,Op,0,3,0,3>
{
  static void apply(TensorEngine<3,T,Antisymmetric> &x, const T2 &y,
		    Op op=Op())
  {
    TensorAssign<TensorEngine<3,T,Antisymmetric>,T2,Op,1,1,0,1>::apply(x,y,op);
    TensorAssign<TensorEngine<3,T,Antisymmetric>,T2,Op,2,1,0,2>::apply(x,y,op);
  }
};



/**
 * TensorEngine definitions for a Symmetric Tensor.
 */

template<int D, class T>
class TensorEngine<D, T, Symmetric>
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d=D };
  typedef T Element_t;
  typedef Symmetric EngineTag_t;

  // Return types for accessor functions.

  // Runtime indices:
  typedef       T&  ElementRef_t;
  typedef const T&  ConstElementRef_t;

  // Compile-time indices:
  typedef        T&  CTElementRef_t;
  typedef  const T&  CTConstElementRef_t;

  // Record the type of the current class.
  typedef TensorEngine<D,T,Symmetric> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor takes no action.
  TensorEngine()
  {
    CTAssert(ElementProperties<T>::hasTrivialDefaultConstructor
	     && ElementProperties<T>::hasTrivialDestructor
	     && ElementProperties<T>::concrete);
  }

  // Copy ctor is deep.
  TensorEngine(const TensorEngine<D,T,Symmetric> &x)
  {
    TensorAssign<This_t,This_t,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
  }

  // Construct from an argument of arbitrary type.
  // The arg must be indexable using TensorElem.
  template<class X>
  explicit TensorEngine(const X &x)
  {
    *this = x; // Use operator=(); see definition below
  }

  // Construct from one argument of type T. Assign this to all the actually
  // stored (not computed) values in the symmetric tensor:
  explicit TensorEngine(const T &x) { 
    for (int i = 0; i < TensorStorageSize<d, EngineTag_t>::Size; i++) {
      x_m[i] = x; // N.B.: may have to move this out to allow for D=1 spec.
    }
  }

  // Construct a 2D Symmetric Tensor from 3 general objects:
  template<class X1, class X2, class X3>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3) 
  {
    CTAssert( D == 2 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
  }
    
  // Construct a 3D Symmetric Tensor from 6 general objects:
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3, const X4 &x4,
               const X5 &x5, const X6 &x6)
  {
    CTAssert( D == 3 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
    (*this)(4) = x5;
    (*this)(5) = x6;
  }
    
  // Let the engine destroy itself.
  ~TensorEngine() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of Tensor.
  This_t&
  operator=(const This_t &x)
  {
    if (this != &x) {
      TensorAssign<This_t,This_t,OpAssign,0,d,0,d>::
        apply(*this,x,OpAssign());
    }
    return *this;
  }

  // Assign from an arbitrary type.
  template<class V>
  This_t& 
  operator=(const V &x)
  {
    TensorAssign<This_t,V,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
    return *this;
  }

  //----------------------------------------------------------------------
  // Element access

  // Compile-time indices, used by things like TensorAssign:

  // Const:
  template<int I, int J>
  CTConstElementRef_t getIJ() const
  {
    CTAssert((I >= 0) && (I < D)
	     && (J >= 0) && (J < D));
    int lo = I < J ? I : J;
    int hi = I > J ? I : J;
    return x_m[((hi + 1)*hi/2) + lo];
  }
    
  // Non-const:
  template<int I, int J>
  CTElementRef_t getIJ()
  {
    CTAssert((I >= 0) && (I < D)
	     && (J >= 0) && (J < D));
    int lo = I < J ? I : J;
    int hi = I > J ? I : J;
    return x_m[((hi + 1)*hi/2) + lo];
  }

  // Runtime indices:

  // (i,j) indexing:

  ConstElementRef_t  operator()(int i,int j) const
  {
    PBoundAssert((i >= 0) && (i < D)
		 && (j >= 0) && (j < D));
    int lo = i < j ? i : j;
    int hi = i > j ? i : j;
    return x_m[((hi + 1)*hi/2) + lo];
  }
  ElementRef_t operator()(int i,int j)
  {
    PBoundAssert((i >= 0) && (i < D)
		 && (j >= 0) && (j < D));
    int lo = i < j ? i : j;
    int hi = i > j ? i : j;
    return x_m[((hi + 1)*hi/2) + lo];
  }

  // (i) indexing (directly into linear array x_m):
  ConstElementRef_t  operator()(int i) const 
  {
    PBoundAssert((i >= 0) && (i < TensorStorageSize<d, EngineTag_t>::Size));
    return x_m[i]; 
  }
  ElementRef_t operator()(int i)
  {
    PBoundAssert((i >= 0) && (i < TensorStorageSize<d, EngineTag_t>::Size));
    return x_m[i]; 
  }

//-----------------------------------------------------------------------------
// For direct access to the data array.

  const T* data() const{
    return (const T*) &x_m[0];
  }

  T* data(){
    return (T*) &x_m[0];
  }

  // This is only used when reading and writing data to disk
  inline void reverseBytes() 
  {
    const int sz = TensorStorageSize<d, EngineTag_t>::Size;
    for (int i = 0; i < sz; ++i)
      ::reverseBytes(x_m[i]);
  }

private:

  // The actual data is just an array of T's. Store in fortran order.
  T x_m[TensorStorageSize<d, EngineTag_t>::Size];
};


/**
 * Special symmetric assignment class: Has specializations for different
 * dimensionalities (for 2, and 3, so far).  This may ultimately be
 * replaceable with a dimensionality-independent equivalent that uses template
 * metaprogramming.
 */

// 2D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<2,T,Symmetric>,T2,Op,0,2,0,2>
{
  static void apply(TensorEngine<2,T,Symmetric> &x, const T2 &y, 
                    Op op=Op())
  {
    TensorAssign<TensorEngine<2,T,Symmetric>,T2,Op,0,1,0,1>::apply(x,y,op);
    TensorAssign<TensorEngine<2,T,Symmetric>,T2,Op,1,1,0,2>::apply(x,y,op);
  }
};
// 3D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<3,T,Symmetric>,T2,Op,0,3,0,3>
{
  static void apply(TensorEngine<3,T,Symmetric> &x, const T2 &y, 
                    Op op=Op())
  {
    TensorAssign<TensorEngine<3,T,Symmetric>,T2,Op,0,1,0,1>::apply(x,y,op);
    TensorAssign<TensorEngine<3,T,Symmetric>,T2,Op,1,1,0,2>::apply(x,y,op);
    TensorAssign<TensorEngine<3,T,Symmetric>,T2,Op,2,1,0,3>::apply(x,y,op);
  }
};



/**
 * TensorEngine definitions for a Diagonal Tensor.
 */

template<int D, class T>
class TensorEngine<D, T, Diagonal>
{
public:

  //----------------------------------------------------------------------
  // Typedefs
  //----------------------------------------------------------------------

  // Export the input types.
  enum { dimensions=2 };
  enum { d=D };
  typedef T Element_t;
  typedef Diagonal EngineTag_t;

  // ----------------------------------------------------------------------
  // The AssignProxy class is needed for non-const Tensor<Diagonal>
  // operator(i,j) implementation:

  class AssignProxy {
  public:
    AssignProxy(Element_t &elem, int where)
      : elem_m(elem), where_m(where) { }
    AssignProxy(const AssignProxy &model)
      : elem_m(model.elem_m), where_m(where_m) { }
    AssignProxy &operator=(const AssignProxy &a) const
    {
      PAssert(where_m != 0);
      elem_m = a.elem_m;
      return const_cast<AssignProxy &>(*this);
    }
    AssignProxy &operator=(const Element_t &e) const
    {
      PAssert(where_m != 0);
      elem_m = e;
      return const_cast<AssignProxy &>(*this);
    }
    operator Element_t() const
    {
      return (elem_m);
    }
  private:
    mutable Element_t &elem_m;
    mutable int where_m;
  };

  // Return types for accessor functions.

  // Runtime indices:
  typedef AssignProxy ElementRef_t;
  // Really return-by-value for Diagonal;const qualifier no longer appropriate:
  typedef T ConstElementRef_t;

  // Compile-time indices:
  typedef T& CTElementRef_t;
  // Really return-by-value for Diagonal;const qualifier no longer appropriate:
  typedef T  CTConstElementRef_t;

  // Record the type of the current class.
  typedef TensorEngine<D,T,Diagonal> This_t;


  //----------------------------------------------------------------------
  // Constructors and Destructor

  // Null ctor takes no action.
  TensorEngine()
  {
    CTAssert(ElementProperties<T>::hasTrivialDefaultConstructor
	     && ElementProperties<T>::hasTrivialDestructor
	     && ElementProperties<T>::concrete);
  }

  // Copy ctor is deep.
  TensorEngine(const TensorEngine<D,T,Diagonal> &x)
  {
    TensorAssign<This_t,This_t,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
  }

  // Construct from an argument of arbitrary type.
  // The arg must be indexable using TensorElem.
  template<class X>
  explicit TensorEngine(const X &x)
  {
    *this = x; // Use operator=(); see definition below
  }

  // Construct from one argument of type T. Assign this to all the actually
  // stored (not computed) values in the symmetric tensor:
  explicit TensorEngine(const T &x) { 
    for (int i = 0; i < TensorStorageSize<d, EngineTag_t>::Size; i++) {
      x_m[i] = x; // N.B.: may have to move this out to allow for D=1 spec.
    }
  }

  // Construct a 2D Diagonal Tensor from 2 general objects:
  template<class X1, class X2>
  TensorEngine(const X1 &x1, const X2 &x2) 
  {
    CTAssert( D == 2 );
    (*this)(0) = x1;
    (*this)(1) = x2;
  }
    
  // Construct a 3D Diagonal Tensor from 3 general objects:
  template<class X1, class X2, class X3>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3) 
  {
    CTAssert( D == 3 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
  }
    
  // Construct a 4D Diagonal Tensor from 4 general objects:
  template<class X1, class X2, class X3, class X4>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3, const X4 &x4)
  {
    CTAssert( D == 4 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
  }
    
  // Construct a 5D Diagonal Tensor from 5 general objects:
  template<class X1, class X2, class X3, class X4, class X5>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3, const X4 &x4,
               const X5 &x5)
  {
    CTAssert( D == 5 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
    (*this)(4) = x5;
  }
    
  // Construct a 6D Diagonal Tensor from 6 general objects:
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3, const X4 &x4,
               const X5 &x5, const X6 &x6)
  {
    CTAssert( D == 6 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
    (*this)(4) = x5;
    (*this)(5) = x6;
  }
    
  // Construct a 7D Diagonal Tensor from 7 general objects:
  template<class X1, class X2, class X3, class X4, class X5, class X6,
    class X7>
  TensorEngine(const X1 &x1, const X2 &x2, const X3 &x3, const X4 &x4,
               const X5 &x5, const X6 &x6, const X7 &x7)
  {
    CTAssert( D == 7 );
    (*this)(0) = x1;
    (*this)(1) = x2;
    (*this)(2) = x3;
    (*this)(3) = x4;
    (*this)(4) = x5;
    (*this)(5) = x6;
    (*this)(6) = x7;
  }
    
  // Let the engine destroy itself.
  ~TensorEngine() {}


  //----------------------------------------------------------------------
  // Assignment

  // Assign from the same kind of Tensor.
  This_t&
  operator=(const This_t &x)
  {
    if (this != &x) {
      TensorAssign<This_t,This_t,OpAssign,0,d,0,d>::
        apply(*this,x,OpAssign());
    }
    return *this;
  }

  // Assign from an arbitrary type.
  template<class V>
  This_t& 
  operator=(const V &x)
  {
    TensorAssign<This_t,V,OpAssign,0,d,0,d>::
      apply(*this,x,OpAssign());
    return *this;
  }

  // ----------------------------------------------------------------------
  // Element access; use the diagonality to compute any value except the
  // stored value(s); all computed values are zero:

  // Compile-time indices, used by things like TensorAssign:

  // Const:
  template<int I, int J>
  CTConstElementRef_t getIJ() const
  {
    CTAssert((I >= 0) && (I < d)
	     && (J >= 0) && (J < d));
    if (I != J) {
      return This_t::Zero;
    } else {
      return x_m[I];
    }
  }
    
  // Non-const. Assume only stored (I,J) values, not computed values, will
  // ever be requested:
  template<int I, int J>
  CTElementRef_t getIJ()
  {
    CTAssert((I >= 0) && (I < d)
	     && (J >= 0) && (J < d));
    CTAssert(I == J);
    return x_m[I];
  }

  // Runtime indices:

  // (i,j) indexing:

  ConstElementRef_t  operator()(int i, int j) const
  {
    PBoundAssert((i >= 0) && (i < D)
		 && (j >= 0) && (j < D));
    if (i != j) {
      return This_t::Zero;
    } else {
      return x_m[i];
    }
  }
  ElementRef_t operator()(int i, int j)
  {
    PBoundAssert((i >= 0) && (i < D)
		 && (j >= 0) && (j < D));
    if (i != j) {
      return AssignProxy(This_t::Zero, 0);
    } else {
      return AssignProxy(x_m[i], 1);
    }
  }

  // (i) indexing (directly into linear array x_m):
  ConstElementRef_t  operator()(int i) const 
  {
    PBoundAssert((i >= 0) && (i < TensorStorageSize<d, EngineTag_t>::Size));
    return x_m[i]; 
  }
  ElementRef_t operator()(int i)
  {
    PBoundAssert((i >= 0) && (i < TensorStorageSize<d, EngineTag_t>::Size));
    return AssignProxy(x_m[i], 1); 
  }

//-----------------------------------------------------------------------------
// For direct access to the data array.

  const T* data() const{
    return (const T*) &x_m[0];
  }

  T* data(){
    return (T*) &x_m[0];
  }

  // This is only used when reading and writing data to disk
  inline void reverseBytes() 
  {
    const int sz = TensorStorageSize<d, EngineTag_t>::Size;
    for (int i = 0; i < sz; ++i)
      ::reverseBytes(x_m[i]);
  }

private:

  // The actual data is just an array of T's. Store in fortran order.
  T x_m[TensorStorageSize<d, EngineTag_t>::Size];

  // A place to store a zero element.  We need to return a reference to this
  // for the diagonal element:
  static T Zero;
};


// Assign the static zero element value:
template<int D, class T>
T TensorEngine<D,T,Diagonal>::Zero = 0.0;


/**
 * Special diagonal assignment class: Has specializations for different
 * dimensionalities (for 2, and 3, so far).  This may ultimately be
 * replaceable with a dimensionality-independent equivalent that uses template
 * metaprogramming.
 */

// 2D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<2,T,Diagonal>,T2,Op,0,2,0,2>
{
  static void apply(TensorEngine<2,T,Diagonal> &x, const T2 &y, Op op=Op())
  {
    TensorAssign<TensorEngine<2,T,Diagonal>,T2,Op,0,1,0,1>::apply(x,y,op);
    TensorAssign<TensorEngine<2,T,Diagonal>,T2,Op,1,1,1,1>::apply(x,y,op);
  }
};
// 3D partial specialization:
template<class T, class T2, class Op>
struct TensorAssign<TensorEngine<3,T,Diagonal>,T2,Op,0,3,0,3>
{
  static void apply(TensorEngine<3,T,Diagonal> &x, const T2 &y, Op op=Op())
  {
    TensorAssign<TensorEngine<3,T,Diagonal>,T2,Op,0,1,0,1>::apply(x,y,op);
    TensorAssign<TensorEngine<3,T,Diagonal>,T2,Op,1,1,1,1>::apply(x,y,op);
    TensorAssign<TensorEngine<3,T,Diagonal>,T2,Op,2,1,2,1>::apply(x,y,op);
  }
};



/**
 * ComponentAccess is an interface class that is used to provide an API for
 * accessing components of a composite type. This version works with Tensors.
 */

template<class T, class Components> struct ComponentAccess;

template<int D, class T, class E, int N>
struct ComponentAccess< Tensor<D, T, E>, Loc<N> >
{
  typedef Tensor<D, T, E> V;
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
// TensorElem specialization for Full tensor engines:
//-----------------------------------------------------------------------------

template<int D, class T, int I, int J>
struct TensorElem< TensorEngine<D,T,Full> , I , J>
{
  typedef TensorEngine<D,T,Full> V;
  typedef TensorEngineElem<D,T,Full,I,J,
    Writable<D,Full,I,J>::value> TE; // TJW added writable...
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x); }
  static      ElementRef_t get(V& x)       { return TE::get(x); }
};


//-----------------------------------------------------------------------------
// TensorElem specialization for Antisymmetric tensor engines:
//-----------------------------------------------------------------------------

template<int D, class T, int I, int J>
struct TensorElem< TensorEngine<D,T,Antisymmetric> , I , J>
{
  typedef TensorEngine<D,T,Antisymmetric> V;
  typedef TensorEngineElem<D,T,Antisymmetric,I,J,
    Writable<D,Antisymmetric,I,J>::value> TE; // TJW added writable...
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x); }
  static      ElementRef_t get(V& x)       { return TE::get(x); }
};



//-----------------------------------------------------------------------------
// TensorElem specialization for Symmetric tensor engines:
//-----------------------------------------------------------------------------

template<int D, class T, int I, int J>
struct TensorElem< TensorEngine<D,T,Symmetric> , I , J>
{
  typedef TensorEngine<D,T,Symmetric> V;
  typedef TensorEngineElem<D,T,Symmetric,I,J,
    Writable<D,Symmetric,I,J>::value> TE; // TJW added writable...
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x); }
  static      ElementRef_t get(V& x)       { return TE::get(x); }
};

//-----------------------------------------------------------------------------
// TensorElem specialization for Diagonal tensor engines:
//-----------------------------------------------------------------------------

template<int D, class T, int I, int J>
struct TensorElem< TensorEngine<D,T,Diagonal> , I , J>
{
  typedef TensorEngine<D,T,Diagonal> V;
  typedef TensorEngineElem<D,T,Diagonal,I,J,
    Writable<D,Diagonal,I,J>::value> TE; // TJW added writable...
  typedef typename TE::Element_t         Element_t;
  typedef typename TE::ConstElementRef_t ConstElementRef_t;
  typedef typename TE::ElementRef_t      ElementRef_t;
  static ConstElementRef_t get(const V& x) { return TE::get(x); }
  static      ElementRef_t get(V& x)       { return TE::get(x); }
};



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



// ----------------------------------------------------------------------------
// 
// Global Function Templates:
//
// ----------------------------------------------------------------------------

template<class OutputEngineTag, int D, class T, class EngineTag>
class Symmetrize;

/// The actual symmetrize() global function template. To get around problems
/// with partially specializating on return type, introduce the Symmetrize
/// functor class, which the symmetrize() function uses. The partial
/// specializations are of the Symmetrize functor. Forward declare the functor
/// before the function template:

template<class OutputEngineTag, int D, class T, class EngineTag>
Tensor<D, T, OutputEngineTag>
symmetrize(const Tensor<D, T, EngineTag> &x)
{
  Symmetrize<OutputEngineTag, D, T, EngineTag> s;
  return s.apply(x);
}

// Needed for application of symmetrize<> to Arrays or Fields of Tensors:
template<class OutputEngineTag, int D, class T, class E>
struct UnaryReturn< Tensor<D,T,E> , FnSymmetrize<OutputEngineTag> >
{
  typedef Tensor<D,T,OutputEngineTag> Type_t;
};



// Now define the Symmetrize functors.

// General template; return a red-flag value (should never be calling this one)
template<class OutputEngineTag, int D, class T, class EngineTag>
class Symmetrize
{
public:
  Symmetrize() {}; // Ctor does nothing.
  // The "guts;" the apply() function of the functor:
  Tensor<D, T, OutputEngineTag>
  apply(const Tensor<D, T, EngineTag> &x) {
    Tensor<D, T, OutputEngineTag> y(-99.99);
    return y;
  }
};

// Partial specializations for specific output symmetries (OutputEngineTag),
// and (where appropriate), specific input symmetries (EngineTag):

// --------------------------------------------------
// Symmetric output.

// Symmetric output; general input:
template<int D, class T, class EngineTag>
class Symmetrize<Symmetric, D, T, EngineTag>
{
public:
  Tensor<D, T, Symmetric>
  apply(const Tensor<D, T, EngineTag> &x)
  {
    Tensor<D, T, Symmetric> y;
    for (int i = 0; i < D; i++) {
      y(i,i) = x(i,i);
      for (int j = i + 1; j < D; j++) {
        y(i,j) = (x(i,j) + x(j,i))*0.5;
      }
    }
    return y;
  }
};
// Symmetric output; antisymmetric input:
template<int D, class T>
class Symmetrize<Symmetric, D, T, Antisymmetric>
{
public:
  Tensor<D, T, Symmetric>
  apply(const Tensor<D, T, Antisymmetric> &x)
  {
    // Result is the zero tensor:
    Tensor<D, T, Symmetric> y(0.0);
    return y;
  }
};
// Symmetric output; Diagonal input:
template<int D, class T>
class Symmetrize<Symmetric, D, T, Diagonal>
{
public:
  Tensor<D, T, Symmetric>
  apply(const Tensor<D, T, Diagonal> &x)
  {
    Tensor<D, T, Symmetric> y(0.0);
    for (int i = 0; i < D; i++) {
      y(i,i) = x(i,i);
    }
    return y;
  }
};


// --------------------------------------------------
// Antisymmetric output.

// Antisymmetric output; general input:
template<int D, class T, class EngineTag>
class Symmetrize<Antisymmetric, D, T, EngineTag>
{
public:
  Tensor<D, T, Antisymmetric>
  apply(const Tensor<D, T, EngineTag> &x)
  {
    Tensor<D, T, Antisymmetric> y;
    // Take the symmetric part of the input:
    for (int i = 1; i < D; i++) {
      for (int j = 0; j < i; j++) {
        y(((i - 1)*i/2)+j) = (x(i,j) - x(j,i))*0.5;
      }
    }
    return y;
  }
};
// Antisymmetric output; Symmetric input:
template<int D, class T>
class Symmetrize<Antisymmetric, D, T, Symmetric>
{
public:
  Tensor<D, T, Antisymmetric>
  apply(const Tensor<D, T, Symmetric> &x)
  {
    // Result is the zero tensor:
    Tensor<D, T, Antisymmetric> y(0.0);
    return y;
  }
};
// Antisymmetric output; Diagonal input:
template<int D, class T>
class Symmetrize<Antisymmetric, D, T, Diagonal>
{
public:
  Tensor<D, T, Antisymmetric>
  apply(const Tensor<D, T, Diagonal> &x)
  {
    // Result is the zero tensor:
    Tensor<D, T, Antisymmetric> y(0.0);
    return y;
  }
};

// --------------------------------------------------
// Diagonal output.

// Diagonal output; general input:
template<int D, class T, class EngineTag>
class Symmetrize<Diagonal, D, T, EngineTag>
{
public:
  Tensor<D, T, Diagonal>
  apply(const Tensor<D, T, EngineTag> &x)
  {
    Tensor<D, T, Diagonal> y;
    for (int i = 0; i < D; i++) {
      y(i) = x(i,i);
    }
    return y;
  }
};
// Diagonal output; Antisymmetric input:
template<int D, class T>
class Symmetrize<Diagonal, D, T, Antisymmetric>
{
public:
  Tensor<D, T, Diagonal>
  apply(const Tensor<D, T, Antisymmetric> &x)
  {
    // Result is the zero tensor:
    Tensor<D, T, Diagonal> y(0.0);
    return y;
  }
};
// Diagonal output; Symmetric input:
template<int D, class T>
class Symmetrize<Diagonal, D, T, Symmetric>
{
public:
  Tensor<D, T, Diagonal>
  apply(const Tensor<D, T, Symmetric> &x)
  {
    Tensor<D, T, Diagonal> y(0.0);
    for (int i = 0; i < D; i++) {
      y(i) = x(i,i);
    }
    return y;
  }
};

// --------------------------------------------------
// Full output.

// Full output; general input:
template<int D, class T, class EngineTag>
class Symmetrize<Full, D, T, EngineTag>
{
public:
  Tensor<D, T, Full>
  apply(const Tensor<D, T, EngineTag> &x)
  {
    Tensor<D, T, Full> y;
    for (int i = 0; i < D; i++) {
      for (int j = 0; j < D; j++) {
        y(i,j) = x(i,j);
      }
    }
    return y;
  }
};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Tensor.h,v $   $Author: richi $
// $Revision: 1.54 $   $Date: 2004/11/29 14:03:30 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
