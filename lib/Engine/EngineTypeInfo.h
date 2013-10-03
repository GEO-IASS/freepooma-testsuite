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
// TypeInfo<> specializations for Engine things.
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_ENGINETYPEINFO_H
#define POOMA_ENGINE_ENGINETYPEINFO_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Engine
 * @brief
 * Undocumented.
 */

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

struct Brick;
struct BrickView;
struct CompressibleBrick;
struct ConstantFunction;
template <class LayoutTag, class PatchTag> struct MultiPatch;
template <class LayoutTag, class PatchTag, int Dim2> struct MultiPatchView;
struct GridTag;

template<class Expr> struct ExpressionTag;

template<int D, class T, class E> class Array;

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Compile-time TypeInfo for Engines.
//-----------------------------------------------------------------------------

template<>
struct TypeInfo<Brick>
{
  static inline std::string name()
  {
    return "Brick";
  }
};

template<>
struct TypeInfo<BrickView>
{
  static inline std::string name()
  {
    return "BrickView";
  }
};

template<>
struct TypeInfo<CompressibleBrick>
{
  static inline std::string name()
  {
    return "CompressibleBrick";
  }
};

template<>
struct TypeInfo<ConstantFunction>
{
  static inline std::string name()
  {
    return "ConstantFunction";
  }
};

template<class Expr>
struct TypeInfo<ExpressionTag<Expr> >
{
  static inline std::string name()
  {
    return "ExpressionTag<" + TypeInfo<Expr>::name() + " >";
  }
};

template<>
struct TypeInfo<GridTag>
{
  static inline std::string name()
  {
    return "GridTag";
  }
};

template<class LayoutTag, class PatchTag>
struct TypeInfo<MultiPatch<LayoutTag, PatchTag> >
{
  static inline std::string name()
  {
    return "MultiPatch<"
      + TypeInfo<LayoutTag>::name() + ","
      + TypeInfo<PatchTag>::name()
      + " >";
  }
};

template<class LayoutTag, class PatchTag, int D2>
struct TypeInfo<MultiPatchView<LayoutTag, PatchTag, D2> >
{
  static inline std::string name()
  {
    return "MultiPatchView<"
      + TypeInfo<LayoutTag>::name() + ","
      + TypeInfo<PatchTag>::name() + ","
      + TypeInfoInt<D2>::name()
      + " >";
  }
};

template<int D, class T, class E>
struct TypeInfo<Array<D, T, E> >
{
  static inline std::string name()
  {
    return "Array<"
      + TypeInfoInt<D>::name() + ","
      + TypeInfo<T>::name() + ","
      + TypeInfo<E>::name()
      + " >";
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_ENGINETYPEINFO_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EngineTypeInfo.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
