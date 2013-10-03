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
// Div
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_DIFFOPS_DIV_UR_H
#define POOMA_FIELD_DIFFOPS_DIV_UR_H

/** @file
 * @ingroup DiffOps
 * @brief
 * Divergence operator on Fields, using 2nd-order centered differences
 * These are used by the div() template function.
 *
 * See Div.h for
 * details, and the div() function template definition.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tiny/Vector.h"
#include "Field/DiffOps/FieldStencil.h"
#include "Field/Mesh/UniformRectilinearMesh.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Partial specializations of the generic Div (divergence) functor. See Div.h
 * for general comments. These are for UniformRectilinear-based 
 * DiscreteGeometry.
 * 
 * Div is a functor class serving as the "Functor" template parameter for
 * FieldStencil<Functor>. Its operator() functions implement 2nd
 * centered differences on an input Field and return output elements of an
 * output Field.  The types of the input and output Field differ in two ways:
 *	-# The input centering is (possibly) different than the output 
 *	   centering.
 *	-# The input element type is Vector<Dim,T2> (or Tensor<Dim,T2>) and
 *         the output type is a scalar type T2 (or Vector<Dim,T2>).
 * Partial specializations implement various combinations of input and output
 * centerings, for specific coordinate systems.
 * 
 * Exported typedefs:
 *  - OutputElement_t: Type of the elements in the output ConstField; 
 *                     restricted to a scalar type (vector input) or vector
 *                     (tensor input)
 * 
 * Accessors:
 *  - inputCentering(): Returns the centering of the input field.  This
 *                      function is just provided as a sanity check for when
 *                      the stencil is created.
 *  - outputCentering(): The centering of the output field. This centering is
 *                      used to construct the return value of the stencil.
 *  - lowerExtent(int d): Returns the stencil width in direction d, at the "low"
 *                      end of the (logically) rectilinear mesh. This is the
 *                      maximum positive integer offset from the element 
 *                      indexed by integer i in the input Field's index space
 *                      along dimension d used in outputting the element
 *                      indexed by integer i in the output Field's index space
 *                      along dimension d
 *  - upperExtent(int d): Same as lowerExtent(), but for the "high" end of the 
 *                      mesh. That is, the maximum (magnitude) *negative*
 *                      offset from i in direction d.
 * 
 * Other methods:
 *  - operator()(...): The actual implementation for the stencil. This acts on a 
 *                   set of scalar-indexed values in the input Field's index
 *                   space making up the stencil, as offset from the fixed
 *                   index point specified by the function's input arguments
 *                   (list of scalar index values).  The stencil must be
 *                   written so that the same fixed index point specified by
 *                   the input arguments where the values are to be assigned in
 *                   the index space of the output Field. This means, for
 *                   example, that if the operator is going from one centering
 *                   to a different output centering, the index bookkeeping
 *                   must be done correctly by this operator()() function
 *                   implementation.
 */

// ----------------------------------------------------------------------------
// Partial specializations of DivVertToCell
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Divergence Vector/Vert -> Scalar/Cell
// ----------------------------------------------------------------------------

template<class T2, class Mesh>
class DivVertToCell;

template<class T2, int Dim, class TM>
class DivVertToCell<Vector<Dim, T2>, UniformRectilinearMesh<Dim, TM> >
{
public:

  typedef T2   OutputElement_t;

  Centering<Dim> outputCentering() const
  {
    return canonicalCentering<Dim>(CellType, Continuous, AllDim);
  }

  Centering<Dim> inputCentering() const
  {
    return canonicalCentering<Dim>(VertexType, Continuous, AllDim);
  }

  // 
  // Constructors.
  // 

  // default version is required by default stencil engine constructor.

  DivVertToCell()
  {
    for (int d = 0; d < Dim; ++d)
    {
      fact_m(d) = 1.0;
    }
  }

  template<class FE>
  DivVertToCell(const FE &fieldEngine)
  {
    for (int d = 0; d < Dim; ++d)
    {
      fact_m(d) = 1 / fieldEngine.mesh().spacings()(d);
    }
  }

  //
  // Methods.
  //

  int lowerExtent(int d) const { return 0; }
  int upperExtent(int d) const { return 1; }
      
  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1) const
  {
    return OutputElement_t
      (fact_m(0)*(f.read(i1+1)(0) - f.read(i1)(0)));
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2) const
  {
    return OutputElement_t
      (0.5*(fact_m(0)*(f.read(i1+1,i2)(0)   - f.read(i1,i2)(0)
		     + f.read(i1+1,i2+1)(0) - f.read(i1,i2+1)(0))
	  + fact_m(1)*(f.read(i1,  i2+1)(1) - f.read(i1,  i2)(1)
		     + f.read(i1+1,i2+1)(1) - f.read(i1+1,i2)(1))));
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2, int i3) const
  {
    return OutputElement_t
      (0.25*(fact_m(0)*(f.read(i1+1,i2,  i3)(0)   - f.read(i1,i2,  i3)(0)
		      + f.read(i1+1,i2+1,i3)(0)   - f.read(i1,i2+1,i3)(0)
		      + f.read(i1+1,i2,  i3+1)(0) - f.read(i1,i2,  i3+1)(0)
		      + f.read(i1+1,i2+1,i3+1)(0) - f.read(i1,i2+1,i3+1)(0))
	   + fact_m(1)*(f.read(i1,  i2+1,i3)(1)   - f.read(i1,  i2,i3)(1)
		      + f.read(i1+1,i2+1,i3)(1)   - f.read(i1+1,i2,i3)(1)
		      + f.read(i1,  i2+1,i3+1)(1) - f.read(i1,  i2,i3+1)(1)
		      + f.read(i1+1,i2+1,i3+1)(1) - f.read(i1+1,i2,i3+1)(1))
	   + fact_m(2)*(f.read(i1,  i2,  i3+1)(2) - f.read(i1,  i2,  i3)(2)
		      + f.read(i1+1,i2,  i3+1)(2) - f.read(i1+1,i2,  i3)(2)
		      + f.read(i1,  i2+1,i3+1)(2) - f.read(i1,  i2+1,i3)(2)
		      + f.read(i1+1,i2+1,i3+1)(2) - f.read(i1+1,i2+1,i3)(2))));
  }

private:

  Vector<Dim, TM> fact_m;
};


template<class T2, class Mesh>
class DivCellToVert;

template<class T2, int Dim, class TM>
class DivCellToVert<Vector<Dim, T2>, UniformRectilinearMesh<Dim, TM> >
{
public:

  typedef T2   OutputElement_t;

  Centering<Dim> outputCentering() const
  {
    return canonicalCentering<Dim>(CellType, Continuous, AllDim);
  }

  Centering<Dim> inputCentering() const
  {
    return canonicalCentering<Dim>(VertexType, Continuous, AllDim);
  }

  // 
  // Constructors.
  // 

  // default version is required by default stencil engine constructor.

  DivCellToVert()
  {
    for (int d = 0; d < Dim; ++d)
    {
      fact_m(d) = 1.0;
    }
  }

  template<class FE>
  DivCellToVert(const FE &fieldEngine)
  {
    for (int d = 0; d < Dim; ++d)
    {
      fact_m(d) = 1 / fieldEngine.mesh().spacings()(d);
    }
  }

  //
  // Methods.
  //

  int lowerExtent(int d) const { return 0; }
  int upperExtent(int d) const { return 1; }
      
  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1) const
  {
    return OutputElement_t
      (fact_m(0)*(f.read(i1)(0) - f.read(i1-1)(0)));
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2) const
  {
    return OutputElement_t
      (0.5*(fact_m(0)*(f.read(i1,i2-1)(0) - f.read(i1-1,i2-1)(0)
		     + f.read(i1,i2)(0)   - f.read(i1-1,i2)(0))
	  + fact_m(1)*(f.read(i1-1,i2)(1) - f.read(i1-1,i2-1)(1)
		     + f.read(i1,  i2)(1) - f.read(i1,  i2-1)(1))));
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2, int i3) const
  {
    return OutputElement_t
      (0.25*(fact_m(0)*(f.read(i1,i2-1,i3-1)(0) - f.read(i1-1,i2-1,i3-1)(0)
		      + f.read(i1,i2,  i3-1)(0) - f.read(i1-1,i2,  i3-1)(0)
		      + f.read(i1,i2-1,i3)(0)   - f.read(i1-1,i2-1,i3)(0)
		      + f.read(i1,i2,  i3)(0)   - f.read(i1-1,i2,  i3)(0))
	   + fact_m(1)*(f.read(i1-1,i2,i3-1)(1) - f.read(i1-1,i2-1,i3-1)(1)
		      + f.read(i1,  i2,i3-1)(1) - f.read(i1,  i2-1,i3-1)(1)
		      + f.read(i1-1,i2,i3)(1)   - f.read(i1-1,i2-1,i3)(1)
		      + f.read(i1,  i2,i3)(1)   - f.read(i1,  i2-1,i3)(1))
	   + fact_m(2)*(f.read(i1-1,i2-1,i3)(2) - f.read(i1-1,i2-1,i3-1)(2)
		      + f.read(i1,  i2-1,i3)(2) - f.read(i1,  i2-1,i3-1)(2)
		      + f.read(i1-1,i2,  i3)(2) - f.read(i1-1,i2,  i3-1)(2)
		      + f.read(i1,  i2,  i3)(2) - f.read(i1,  i2,  i3-1)(2))));
  }

private:

  Vector<Dim, TM> fact_m;
};


template<class T2, class Mesh, CenteringType OC>
class DivSameToSame;

template<class T2, int Dim, class TM, CenteringType OC>
class DivSameToSame<Vector<Dim, T2>, UniformRectilinearMesh<Dim, TM>, OC>
{
public:

  typedef T2   OutputElement_t;

  Centering<Dim> outputCentering() const
  {
    return canonicalCentering<Dim>(OC, Continuous);
  }

  Centering<Dim> inputCentering() const
  {
    return canonicalCentering<Dim>(OC, Continuous);
  }

  // 
  // Constructors.
  // 

  // default version is required by default stencil engine constructor.

  DivSameToSame()
  {
    for (int d = 0; d < Dim; ++d)
    {
      fact_m(d) = 0.5;
    }
  }

  template<class FE>
  DivSameToSame(const FE &fieldEngine)
  {
    for (int d = 0; d < Dim; ++d)
    {
      fact_m(d) = 0.5 / fieldEngine.mesh().spacings()(d);
    }
  }

  //
  // Methods.
  //

  int lowerExtent(int d) const { return 1; }
  int upperExtent(int d) const { return 1; }
      
  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1) const
  {
    return OutputElement_t
      (fact_m(0)*(f.read(i1+1)(0) - f.read(i1-1)(0)));
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2) const
  {
    return OutputElement_t
      (fact_m(0)*(f.read(i1+1,i2)(0)   - f.read(i1-1,i2)(0))
     + fact_m(1)*(f.read(i1,  i2+1)(1) - f.read(i1,  i2-1)(1)));
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2, int i3) const
  {
    return OutputElement_t
      (fact_m(0)*(f.read(i1+1,i2,  i3)(0)   - f.read(i1-1,i2,  i3)(0))
     + fact_m(1)*(f.read(i1,  i2+1,i3)(1)   - f.read(i1,  i2-1,i3)(1))
     + fact_m(2)*(f.read(i1,  i2,  i3+1)(2) - f.read(i1,  i2,  i3-1)(2)));
  }

private:

  Vector<Dim, TM> fact_m;
};

#endif     // POOMA_FIELD_DIFFOPS_DIV_UR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Div.UR.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:44 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
