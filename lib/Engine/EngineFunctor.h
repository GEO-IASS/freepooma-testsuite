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
// EngineFunctor<Engine,Tag>
// EngineFunctorDefault<Engine,Tag>
// EngineFunctorScalar<T,Tag>
// EngineView<Tag>
// ExpressionApply<Tag>
//
// function:
// engineFunctor(eng,tag)
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_ENGINEFUNCTOR_H
#define POOMA_ENGINE_ENGINEFUNCTOR_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Engine
 * @brief
 * EngineFunctor provides a common interface to a variety of engine queries,
 * like "are you compressed" or "are you shifted with respect to another
 * engine". 
 *
 * By providing a common interface, we minimize the number of changes
 * you need to make to add a new capability that makes non-standard queries of
 * engines.
 *
 * This approach replaces the previous message() function which made queries
 * of engines.  The current version doesn't require new member functions to be
 * added to engines to support new capabilities.  Also this approach allows for
 * simple default cases, since you use partial specialization on EngineFunctor.
 *
 * WARNING: If you use a default action, you should probably have some
 * verification mechanism to ensure that the engine doesn't need to have a
 * special action defined.  For example, the DataObject default action checks
 * the dataObject enum in the engine to make sure it doesn't have a data
 * object.
 *
 * Default actions for a given functor are specified using
 * EngineFunctorDefault.  We define EngineFunctor for Expression Engines and
 * general tags, so you cannot specify the action of a specific functor
 * for general engines.  (The compiler would not be able to choose between
 * EngineFunctor<ExpressionEngine,T> and EngineFunctor<T,YourFunctor>)
 * The generic version of EngineFunctor<> calls DefaultEngineFunctor<>, so
 * you may define EngineFunctorDefault for a general engine and tag.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag> class Engine;


/**
 * The most generic version of EngineFunctor falls through to
 * EngineFunctorDefault.
 *
 * (We define
 * template<class Tag> EngineFunctor<ExpressionEngine, Tag>
 * so you cannot define
 * template<class Engine> EngineFunctor<Engine, SpecificTag>
 * without ambiguities.
 *
 * By using EngineFunctorDefault, you can define general operations
 * for a given tag that will be caught by all non-expression engines.
 */

template<class Eng, class Tag>
struct EngineFunctorDefault
{ };


/**
 * EngineFunctor<Eng,Tag> defines the action corresponding to Tag on the
 * engine of type Eng.  Requires the action return type and apply function:
 *
 * <PRE>
 *   typedef return_type   Type_t;     // return type of the action
 *   static Type_t apply(const Eng &engine, const Tag &tag);
 * </PRE>
 *   
 * Tag is a policy tag for the action.  Requires one typedef:
 *
 *   typedef combiner_type Combine_t;  // PETE combiner for the return
 *
 * EngineFunctorScalar<T,Tag> computes a value for a scalar leaf. Requires:
 *
 *   typedef return_type   Type_t;     // return type of the action
 *   static Type_t apply(const Scalar<T> &scalar, const Tag &tag);
 *
 * EngineFunctorDefault<Eng,Tag> allows you to define generic actions
 * for a specific tag.
 */

template<class Eng, class Tag>
struct EngineFunctor
{
  typedef typename EngineFunctorDefault<Eng,Tag>::Type_t Type_t;
  inline static
  Type_t apply(const Eng &e, const Tag &t)
  {
    return EngineFunctorDefault<Eng,Tag>::apply(e,t);
  }
};


/// This helper function is shorthand for:
///
/// EngineFunctor<Engine, Tag>::apply(engine, tag)

template<class Eng, class Tag>
inline typename EngineFunctor<Eng,Tag>::Type_t
engineFunctor(const Eng &e, const Tag &tag)
{
  return EngineFunctor<Eng,Tag>::apply(e,tag);
}

/**
 * Users must specialize this struct for all tags.
 * The specialization needs to contain Type_t and an apply method:
 *
 * <PRE>
 * template<class T>
 * struct EngineFunctorScalar<T, MyTag >
 * {
 *   typedef ... Type_t;
 *
 *   inline static
 *   Type_t apply(const T &s, const MyTag &t)
 *   {
 *     ...
 *   }
 * };
 * </PRE>
 */

template<class T, class Tag>
struct EngineFunctorScalar
{ };


/** 
 * EngineView<Tag> and ExpressionApply<Tag> are replacements for EngineFunctor.
 * EngineFunctor applied to an expression uses forEach, which means there are
 * two levels of indirection at the leaves.  EngineView and EngineTag are
 * forEach functors, reducing the number of levels of indirection.
 */

template<class Tag>
struct EngineView;

/**
 * LeafFunctor specializations for EngineView.
 *
 * Applying EngineView to a general node is an error.
 * Applying EngineView to a scalar just returns the scalar.
 */

template<class Node, class Tag>
struct LeafFunctor<Node, EngineView<Tag> >
{
};

template<class T, class Tag>
struct LeafFunctor<Scalar<T>, EngineView<Tag> >
{
  typedef Scalar<T> Type_t;
  static inline
  Type_t apply(const Scalar<T> &s, const EngineView<Tag> &)
  {
    return s;
  }
};


/**
 * For a given type of engine view, you must either specialize LeafFunctor
 * for all engines or provide a specialization of DefaultEngineView.
 *
 * This level of indirection is necessary to avoid the ambiguity that would
 * result from attempting to provide the specializations:
 *
 * LeafFunctor<ExpressionEngine, EngineView<GeneralTag>>
 * LeafFunctor<GeneralEngine, EngineView<SpecificTag>>
 */

template<class Engine, class Tag>
struct DefaultEngineView;

template<int Dim, class T, class E, class Tag>
struct LeafFunctor<Engine<Dim, T, E>, EngineView<Tag> >
{
  typedef Engine<Dim, T, E> Subject_t;
  typedef DefaultEngineView<Subject_t, Tag> EngineView_t;
  typedef typename EngineView_t::Type_t Type_t;

  static inline
  Type_t apply(const Subject_t &engine,
	       const EngineView<Tag> &tag)
  {
    return EngineView_t::apply(engine, tag);
  }
};


/**
 * The default version of ExpressionApply contains a reference to the tag
 * class, which is typically empty.  Users can store information in the tag
 * needed in the LeafFunctors.
 */

template<class Tag>
struct ExpressionApply
{
  inline
  ExpressionApply(const Tag &tag)
    : tag_m(tag)
  {
  }

  template<class A>
  void operator()(const A &a) const
  {
    forEach(a, *this, NullCombine());
  }

  const Tag &tag() const { return tag_m; }
  const Tag &tag_m;
};

template<class A, class Tag>
inline void
expressionApply(const A &a, const Tag &tag)
{
  forEach(a, ExpressionApply<Tag>(tag), NullCombine());
}


/**
 * LeafFunctor specializations for ExpressionApply.
 *
 * Applying EngineView to a general node is an error.
 */

template<class Node, class Tag>
struct LeafFunctor<Node, ExpressionApply<Tag> >
{
};

template<class T, class Tag>
struct LeafFunctor<Scalar<T>, ExpressionApply<Tag> >
{
  typedef int Type_t;
  static inline
  Type_t apply(const Scalar<T> &, const ExpressionApply<Tag> &)
  {
    return 0;
  }
};


/**
 * For a given type of engine view, you must either specialize LeafFunctor
 * for all engines or provide a specialization of DefaultEngineView.
 *
 * This level of indirection is necessary to avoid the ambiguity that would
 * result from attempting to provide the specializations:
 *
 * LeafFunctor<ExpressionEngine, EngineView<GeneralTag>>
 * LeafFunctor<GeneralEngine, EngineView<SpecificTag>>
 */

template<class Engine, class Tag>
struct DefaultExpressionApply;

template<int Dim, class T, class E, class Tag>
struct LeafFunctor<Engine<Dim, T, E>, ExpressionApply<Tag> >
{
  typedef Engine<Dim, T, E> Subject_t;
  typedef DefaultExpressionApply<Subject_t, Tag> ExpressionApply_t;
  typedef int Type_t;

  static inline
  Type_t apply(const Subject_t &engine,
	       const ExpressionApply<Tag> &tag)
  {
    return ExpressionApply_t::apply(engine, tag);
  }
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_ENGINEFUNCTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EngineFunctor.h,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
