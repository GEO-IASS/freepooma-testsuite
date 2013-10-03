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

/** @file
 * @ingroup Field
 * @brief
 * TypeInfo<> specializations for Field things.
 */

#ifndef POOMA_FIELD_FIELDTYPEINFO_H
#define POOMA_FIELD_FIELDTYPEINFO_H

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

#include "Utilities/TypeInfo.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class GeometryTag, class T, class EngineTag> class Field;
template<class GeometryTag, class T, class EngineTag> class FieldEngine;
template<int Dim, class T, class EngineTag>           class FieldEngineBase;

template<int Dim> struct NoGeometry;

template<int Dim, class T, class CoordinateSystem>
struct UniformRectilinear;

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Compile-time TypeInfo for Field.
//-----------------------------------------------------------------------------

template<class GeometryTag, class T, class EngineTag>
struct TypeInfo<Field<GeometryTag, T, EngineTag> >
{
  static inline std::string name()
  {
    return "Field<"
      + TypeInfo<GeometryTag>::name() + ","
      + TypeInfo<T>::name() + ","
      + TypeInfo<EngineTag>::name() + " >";
  }
};

template<int Dim>
struct TypeInfo<NoGeometry<Dim> >
{
  static inline std::string name()
  {
    return "NoGeometry<" + TypeInfoInt<Dim>::name() + ">";
  }
};

template<int Dim, class T, class CoordinateSystem>
struct TypeInfo<UniformRectilinear<Dim, T, CoordinateSystem> >
{
  static inline std::string name()
  {
    return "UniformRectilinear<"
      + TypeInfoInt<Dim>::name() + ","
      + TypeInfo<T>::name() + ","
      + TypeInfo<CoordinateSystem>::name() + " >";
  }
};

template<class GeometryTag, class T, class EngineTag>
struct TypeInfo<FieldEngine<GeometryTag, T, EngineTag> >
{
  static inline std::string name()
  {
    return "FieldEngine<"
      + TypeInfo<GeometryTag>::name() + ","
      + TypeInfo<T>::name() + ","
      + TypeInfo<EngineTag>::name() + " >";
  }
};

template<int Dim, class T, class EngineTag>
struct TypeInfo<FieldEngineBase<Dim, T, EngineTag> >
{
  static inline std::string name()
  {
    return "FieldEngineBase<"
      + TypeInfoInt<Dim>::name() + ","
      + TypeInfo<T>::name() + ","
      + TypeInfo<EngineTag>::name() + " >";
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_FIELD_FIELDTYPEINFO_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldTypeInfo.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:43 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
