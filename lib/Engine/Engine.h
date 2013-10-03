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

#ifndef POOMA_ENGINE_ENGINE_H
#define POOMA_ENGINE_ENGINE_H

//-----------------------------------------------------------------------------
// Classes:
//   Engine<Dim,T,EngineTag>
//   NewEngine<Engine,SubDomain>
//   NewEngineEngine<Engine,SubDomain>
//   NewEngineDomain<Engine,SubDomain>
//   EngineConstructTag
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

/** @file
 * @ingroup Engine
 * @brief
 * Engine
 *  - Engine<Dim,T,EngineTag>
 *    General Engine template.
 *  - NewEngine<Engine,SubDomain>
 *    Traits class used for taking views.
 *  - NewEngineEngine<Engine,SubDomain>
 *  - NewEngineDomain<Engine,SubDomain>
 *    Optional functors that can be used to simplify view construction.
 *  - EngineConstructTag
 *    tag class to disambiguate certain constructors.
 */


/**
 * General Engine template. All concrete Engine classes specialize this
 * template for a particular Engine tag.
 */

template <int Dim, class T, class EngineTag>
class Engine;


/**
 * Traits class for determining the type obtained by subsetting a 
 * particular Engine with a particular SubDomain.
 *
 * Concrete Engine classes will make specializations of this class
 * for the pairs that can result in that particular Engine being created.
 */

template <class Engine, class SubDomain>
struct NewEngine
{
};


/**
 * NewEngineEngine<Engine,SubDomain>
 * NewEngineDomain<Engine,SubDomain>
 *
 * These two traits classes allow you to modify the engine and domain
 * that get passed to the view engine that results from a subset
 * operation.  This indirection allows you to define view operations
 * without requiring that the view engine knows about the engine you
 * start with (for example, BrickView shouldn't have to know about
 * patch engines that contain it).
 *
 * The natural location for these functors is inside the NewEngine
 * traits class, but defining them separately allows us to provide
 * the default behaviour: just forward the engine and domain through.
 */

template <class Engine, class SubDomain>
struct NewEngineEngine
{
  typedef Engine Type_t;
  static inline
  const Engine &apply(const Engine &e, const SubDomain &)
  {
    return e;
  }
} ;

template <class Engine, class SubDomain>
struct NewEngineDomain
{
  typedef SubDomain Type_t;
  typedef const SubDomain &Return_t;
  static inline
  Return_t apply(const Engine &, const SubDomain &i)
  {
    return i;
  }
} ;


template<class Eng, class Dom>
inline typename NewEngineEngine<Eng, Dom>::Type_t
newEngineEngine(const Eng &e, const Dom &dom)
{
  return NewEngineEngine<Eng, Dom>::apply(e, dom);
}


template<class Eng, class Dom>
inline typename NewEngineDomain<Eng, Dom>::Type_t
newEngineDomain(const Eng &e, const Dom &dom)
{
  return NewEngineDomain<Eng, Dom>::apply(e, dom);
}


/**
 * EngineConstructTag is used by Array to disambiguate
 * engine-based constructor calls.  (also used by some engines)
 */

struct EngineConstructTag 
{
  EngineConstructTag() { };
  ~EngineConstructTag() { };
  EngineConstructTag(const EngineConstructTag &) { };
  EngineConstructTag &operator=(const EngineConstructTag &) { return *this; }
};

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Engine.h,v $   $Author: richard $
// $Revision: 1.27 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
