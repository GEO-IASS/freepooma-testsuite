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
// Class:
// XField
// Functions:
// xField
//-----------------------------------------------------------------------------

/** @file
 * @ingroup XUnused
 * @brief
 * Computes position locations for Uniform Rectilinear meshes.
 */

#ifndef POOMA_FIELD_MESH_POSITIONFUNCTIONS_H
#define POOMA_FIELD_MESH_POSITIONFUNCTIONS_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tiny/Vector.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Computes position locations for Uniform Rectilinear meshes.
 * The constructor takes the origin and spacings, the physical domain over
 * which the mesh is defined, and the Loc describing the centering position.
 */


template<int Dim, class TM>
struct PositionFunctorUR
{
  typedef Vector<Dim, TM> PointType_t;

  PositionFunctorUR()
  {
  }

  PositionFunctorUR(const PointType_t &origin,
                    const PointType_t &spacings,
                    const Interval<Dim> &physicalCellDomain,
                    const Vector<Dim, double> &position)
    : origin_m(origin), spacings_m(spacings)
  {
    for (int i = 0; i < Dim; i++)
    {
      origin_m(i) += position(i) * spacings_m(i);
      origin_m(i) -= physicalCellDomain[i].first() * spacings_m(i);
    }
  }

  PointType_t operator()(int i0) const
  {
    return origin_m + PointType_t(i0) * spacings_m;
  }
      
  PointType_t operator()(int i0, int i1) const
  {
    return origin_m + PointType_t(i0, i1) * spacings_m;
  }

  PointType_t operator()(int i0, int i1, int i2) const
  {
    return origin_m + PointType_t(i0, i1, i2) * spacings_m;
  }

  PointType_t origin_m, spacings_m;
};

template<int Dim, class TM>
struct FixPositionFunctorUR
{
  typedef Vector<Dim, TM> PointType_t;

  FixPositionFunctorUR(const PointType_t &origin,
		       const PointType_t &spacings)
    : origin_m(origin), spacings_m(spacings)
  {
  }

  template<class FEB>
  void operator()(FEB &fieldEngineBase) const
  {
    typedef PositionFunctorUR<Dim, TM> Functor_t;

    for (int m = 0; m < fieldEngineBase.numMaterials(); ++m)
    {
      for (int c = 0; c < fieldEngineBase.centering().size(); ++ c)
      {
        fieldEngineBase.data(m, c).engine().
          setFunctor(Functor_t(origin_m, spacings_m,
                               fieldEngineBase.physicalCellDomain(),
                               fieldEngineBase.centering().position(c))
                     );
      }
    }
  }

  PointType_t origin_m, spacings_m;
};

template<int Dim, class TM, class T, class EngineTag>
FixPositionFunctorUR<Dim, TM>
fixPositionFunctor(const Field<UniformRectilinearMesh<Dim, TM>,
		   T, EngineTag> &field)
{
  return FixPositionFunctorUR<Dim, TM>(field.mesh().origin(),
				       field.mesh().spacings());
}

template<class Geom>
struct XField
{
};

template<int Dim, class TM>
struct XField<UniformRectilinearMesh<Dim, TM> >
{
  typedef UniformRectilinearMesh<Dim, TM> Mesh_t;
  typedef Vector<Dim, TM> PointType_t;
  typedef IndexFunction<PositionFunctorUR<Dim, TM> > PositionEngine_t;
  typedef Field<Mesh_t, PointType_t, PositionEngine_t> Type_t;
};

template<class F>
void setXField(F &f)
{
  fixPositionFunctor(f)(f.fieldEngine());
}

template<class F, class Init>
typename XField<typename F::Mesh_t>::Type_t
xField(const F &f, const Init &centering)
{
  typedef typename XField<typename F::Mesh_t>::Type_t Field_t;
  Field_t ret(centering, f.layout(), f.mesh());
  setXField(ret);
  return ret;
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_FIELD_MESH_POSITIONFUNCTIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PositionFunctions.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:46 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
