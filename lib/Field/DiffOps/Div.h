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

#ifndef POOMA_FIELD_DIFFOPS_DIV_H
#define POOMA_FIELD_DIFFOPS_DIV_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup DiffOps
 * @brief
 * Divergence operator (functor) on discrete Fields.
 *
 * Wrapper function around FieldStencil<Div>::operator() . The Div
 * functors actually used are partial specializations of the generic
 * Div that come from Div.UR for example.
 *
 * Div is a functor class serving as the "Functor" template parameter for
 * FieldStencil<Functor,Expression>, which implements a discrete 
 * divergence operator.
 * Partial specializations implement various combinations of input and output
 * centerings, for specific coordinate systems, and different finite-difference
 * orders, are defined in other headers like Div.[URM,RM].h .
 * 
 * div(): Divergence. Takes a Field of Vectors (or Tensors) on a 
 * discrete geometry with one centering and returns a Field of
 * scalars (or Vectors) on a geometry that's the same except
 * (possibly) for the centering. All the work happens in the embedded
 * Div functor partial specialization, in its operator() methods.
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
// General Div template
// ----------------------------------------------------------------------------

template<class T2, class Mesh>
class DivCellToVert;

template<class T2, class Mesh>
class DivVertToCell;

template<class T2, class Mesh, CenteringType OC>
class DivSameToSame;


// ----------------------------------------------------------------------------
// 
// Global Function Templates:
//
// ----------------------------------------------------------------------------

// Divergence.

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<DivSameToSame<T, Mesh, CellType>,
  Field<Mesh, T, EngineTag> >::Type_t
divCellToCell(const Field<Mesh, T, EngineTag> &f)
{
  typedef DivSameToSame<T, Mesh, CellType> Div_t;
  typedef FieldStencilSimple<Div_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Div_t(f.fieldEngine()), f);
}

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<DivVertToCell<T, Mesh>,
  Field<Mesh, T, EngineTag> >::Type_t
divVertToCell(const Field<Mesh, T, EngineTag> &f)
{
  typedef DivVertToCell<T, Mesh> Div_t;
  typedef FieldStencilSimple<Div_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Div_t(f.fieldEngine()), f);
}

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<DivCellToVert<T, Mesh>,
  Field<Mesh, T, EngineTag> >::Type_t
divCellToVert(const Field<Mesh, T, EngineTag> &f)
{
  typedef DivCellToVert<T, Mesh> Div_t;
  typedef FieldStencilSimple<Div_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Div_t(f.fieldEngine()), f);
}

template<class Mesh, class T, class EngineTag>
typename
FieldStencilSimple<DivSameToSame<T, Mesh, VertexType>,
  Field<Mesh, T, EngineTag> >::Type_t
divVertToVert(const Field<Mesh, T, EngineTag> &f)
{
  typedef DivSameToSame<T, Mesh, VertexType> Div_t;
  typedef FieldStencilSimple<Div_t, Field<Mesh, T, EngineTag> > Ret_t;
  return Ret_t::make(Div_t(f.fieldEngine()), f);
}

#endif     // POOMA_FIELD_DIFFOPS_DIV_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Div.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:44 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
