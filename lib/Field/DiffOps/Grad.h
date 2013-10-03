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
//   Div
// Global Function Templates:
//   div
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_DIFFOPS_GRAD_H
#define POOMA_FIELD_DIFFOPS_GRAD_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup DiffOps
 * @brief
 * Gradient operator (functor) on discrete Fields.
 *
 * Wrapper function around FieldStencil<Grad>::operator() . The Div
 * functors actually used are partial specializations of the generic
 * Grad that come from Grad.UR for example.
 *
 * Grad is a functor class serving as the "Functor" template parameter for
 * FieldStencil<Functor,Expression>, which implements a discrete 
 * gradient operator.
 * Partial specializations implement various combinations of input and output
 * centerings, for specific coordinate systems, and different finite-difference
 * orders, are defined in other headers like Grad.[URM,RM].h .
 * 
 * grad(): Gradient. Takes a scalar Field a 
 * discrete geometry with one centering and returns a Field of
 * vectors on a geometry that's the same except
 * (possibly) for the centering. All the work happens in the embedded
 * Grad functor partial specialization, in its operator() methods.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Field/Field.h"
#include "Field/FieldCentering.h"
#include "Field/DiffOps/FieldStencil.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// General Grad template
// ----------------------------------------------------------------------------

template<class T2, class Mesh>
class GradCellToVert;

template<class T2, class Mesh>
class GradVertToCell;

template<class T2, class Mesh, CenteringType OC>
class GradSameToSame;


// ----------------------------------------------------------------------------
// 
// Global Function Templates:
//
// ----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<GradVertToCell<T, Mesh>,
  Field<Mesh, T, EngineTag> >::Type_t
gradVertToCell(const Field<Mesh, T, EngineTag> &f)
{
  typedef GradVertToCell<T, Mesh> Grad_t;
  typedef FieldStencilSimple<Grad_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Grad_t(f.fieldEngine()), f);
}

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<GradCellToVert<T, Mesh>,
  Field<Mesh, T, EngineTag> >::Type_t
gradCellToVert(const Field<Mesh, T, EngineTag> &f)
{
  typedef GradCellToVert<T, Mesh> Grad_t;
  typedef FieldStencilSimple<Grad_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Grad_t(f.fieldEngine()), f);
}

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<GradSameToSame<T, Mesh, CellType>,
  Field<Mesh, T, EngineTag> >::Type_t
gradCellToCell(const Field<Mesh, T, EngineTag> &f)
{
  typedef GradSameToSame<T, Mesh, CellType> Grad_t;
  typedef FieldStencilSimple<Grad_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Grad_t(f.fieldEngine()), f);
}

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<GradSameToSame<T, Mesh, VertexType>,
  Field<Mesh, T, EngineTag> >::Type_t
gradVertToVert(const Field<Mesh, T, EngineTag> &f)
{
  typedef GradSameToSame<T, Mesh, VertexType> Grad_t;
  typedef FieldStencilSimple<Grad_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Grad_t(f.fieldEngine()), f);
}


#endif     // POOMA_FIELD_DIFFOPS_GRAD_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Grad.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:44 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
