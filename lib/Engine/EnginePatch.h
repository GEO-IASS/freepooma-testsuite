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
 * @ingroup Engine
 * @brief
 * Patch handling with engines:
 *  - EnginePatch, functor for getting the nth patch from an engine.
 *  - EngineNumPatches, gives the number of patches in an engine.
 *  - EnginePatchEmpty, tells you if the nth patch is empty.
 *  - PatchView<Container>, traits class for the type of object returned
 *    by the patch() member function.
 */

#ifndef POOMA_ENGINE_ENGINEPATCH_H
#define POOMA_ENGINE_ENGINEPATCH_H

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

#include "PETE/PETE.h"
#include "Pooma/PETE/AssertEquals.h"
#include "Pooma/View.h"
#include "Engine/EngineFunctor.h"


/**
 * EnginePatch is a tag used with engineFunctor to return the ith patch
 * from a multi-patch engine in a generic way.  Non-multipatch engines
 * are defined to have 1 patch, so you can use EnginePatch on them as well.
 * The syntax looks like:
 *
 * patchEngine = engineFunctor(engine, EnginePatch(i));
 */

struct EnginePatch
{
  // Expression engines need to be combine the patch engines into
  // a new tree:
  typedef TreeCombine Combine_t;
  typedef int PatchID_t;

  explicit EnginePatch(PatchID_t patch) : patch_m(patch) { }
  PatchID_t patch_m;
};


template<class Eng>
struct EngineFunctorDefault<Eng, EnginePatch>
{
  typedef Eng Type_t;

  inline static
  const Type_t &apply(const Eng &e, const EnginePatch &)
  {
    // Engines that are multipatch must specialize this functor
    // to access the patch
    CTAssert(!(Eng::multiPatch));
    return e;
  }
};

template<class T>
struct LeafFunctor<Scalar<T>, EnginePatch>
{
  typedef Scalar<T> Type_t;

  inline static
  Type_t apply(const Scalar<T> &scalar, const EnginePatch &)
  {
    return scalar;
  }
};


/**
 * PatchView is the traits class used to describe the result of
 * .patch(i) for a container class.  PatchView is specialized for
 * Arrays etc. to give the Array with the engine given by
 * EngineFunctor<Engine,EnginePatch>::Type_t;
 */

template<class Container>
struct Patch
{
};


template<class Container>
struct PatchView
{
  typedef typename Patch<Container>::Type_t Patch_t;
  typedef typename Patch_t::Domain_t Dom_t;
  typedef typename View1<Patch_t, Dom_t>::Type_t Type_t;

  inline static
  Type_t make(const Container &subject, int i)
  {
    return subject.patchLocal(i)();
  }
};

template<class Node>
struct LeafFunctor<Node, EnginePatch>
{
  typedef typename PatchView<Node>::Type_t Type_t;

  inline static
  Type_t apply(const Node &node, const EnginePatch &tag)
  {
    return node.patchLocal(tag.patch_m)();
  }
};


/**
 * EngineNumPatches is used to find out how many patches an engine has.
 * (Or throw an assertion if the answer is ambiguous.)  Typical use would
 * look like:
 *
 * <PRE>
 * int n = engineFunctor(a.engine(), EngineNumPatches());
 * for (i = 0; i < n; ++i)
 *     calculate(a.patch(i));
 * </PRE>
 */

struct EngineNumPatches
{
  // Throw an assertion if not all the engines have the same number of
  // patches:
  typedef AssertEquals Combine_t;
};

// Generic engines have 1 patch.

template<class Eng>
struct EngineFunctorDefault<Eng, EngineNumPatches>
{
  typedef int Type_t;

  inline static
  Type_t apply(const Eng &, const EngineNumPatches &)
  {
    // Engines that are multipatch must specialize this functor
    // to access the patch
    CTAssert(!(Eng::multiPatch));
    return 1;
  }
};

// Scalars have 0 patches.

template<class T>
struct EngineFunctorScalar<T, EngineNumPatches>
{
  typedef int Type_t;

  inline static
  Type_t apply(const T &, const EngineNumPatches &)
  {
    return 0;
  }
};

/////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_ENGINEPATCH_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EnginePatch.h,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
