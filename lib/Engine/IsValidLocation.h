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

#ifndef POOMA_ENGINE_ISVALIDLOCATION_H
#define POOMA_ENGINE_ISVALIDLOCATION_H

/** @file
 * @ingroup Engine
 * @brief
 * helper functions to determine if a particular location
 * or region of a object is defined.
 *
 * These helper functions are used to determine if a particular location
 * or region of a object is defined. For all objects _not_ based
 * on SparseTileLayout, simply return true. 
 * For STL based objects, do a touches, and if anything is found, return true
 * otherwise return false. 
 *
 * For Expression Engines, use the EngineFunctor to search the expression
 * tree and do logic on the result. 
 *
 * Used in PrintArray....
 */

template <class LayoutTag, class PatchTag>
struct MultiPatch;

template <class LayoutTag, class PatchTag, int Dim2>
struct MultiPatchView;

template<class Expr>
struct ExpressionTag;

#include "Engine/Engine.h"
#include "PETE/PETE.h"
#include "Utilities/WrappedInt.h"
#include "Engine/EngineFunctor.h"
#include "Layout/SparseTileLayout.h"

template<int Dim >
struct IsValid
{
  IsValid(Loc<Dim> loc) : loc_m(loc) { }
  typedef AndCombine Combine_t;
  Loc<Dim> loc_m;
};

//-----------------------------------------------------------------------------
// Scalars are always valid.
//-----------------------------------------------------------------------------

template<class T, int Dim>
struct EngineFunctorScalar<T, IsValid<Dim> >
{
  typedef bool Type_t;
  static inline
  Type_t apply(const T &, const IsValid<Dim> &)
  {
    return true;
  }
};

template<class Engine, int Dim>
struct EngineFunctorDefault<Engine, IsValid<Dim> >
{
  typedef bool  Type_t;
  static inline
  Type_t apply(const Engine &, const IsValid<Dim> &)
  {
    return true;
  }
};

template<int Dim,class T,class ptag>
struct EngineFunctor<Engine<Dim, T,MultiPatch<SparseTileTag,ptag> >, IsValid<Dim> >
{
  typedef Engine<Dim,T,MultiPatch<SparseTileTag,ptag> > Engine_t;

  typedef bool Type_t;

  static inline
  Type_t apply(const Engine_t &e, const IsValid<Dim> &f)
  {
    typedef typename Engine_t::Domain_t domain_t;
    typedef Node<domain_t,domain_t> node_t;
    std::vector<node_t> v;
    int count = e.layout().touches(f.loc_m,std::back_inserter(v));
    return (count!=0);
  }
};

template<class Object,class Dom,class tag>
inline bool isValidLocation(const Object &,
			    const Dom &,
			    const tag &)
{
 return true;
}

#endif
