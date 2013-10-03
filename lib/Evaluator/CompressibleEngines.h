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
// Tags:
// Compressible
// Compressed
// CompressedRead
// CompressedReadWrite
// CompressedBrickIsWholeView
// UnCompressedViewEngine
//
// Classes:
// EngineFunctor<Engine,Tag>  for above tags.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Work to do:
//   -test compression with constant function
//   -figure out interaction of stencil engine with compressed eval.
//    (right now, they are just viewed as uncompressed which is probably
//    inefficient)
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_COMPRESSIBLEENGINES_H
#define POOMA_EVALUATOR_COMPRESSIBLEENGINES_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * We define EngineFunctor<Engine,Tag> for tags that express the
 * functionality that compressible bricks have.
 *
 * Since compressible bricks can
 * appear inside other engines, we need to provide some mechanism for those
 * other engines to provide the same interface.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/Engine.h"
#include "PETE/PETE.h"
#include "Utilities/WrappedInt.h"
#include "Engine/EngineFunctor.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

struct Brick;
struct BrickView;
struct CompressibleBrick;
struct CompressibleBrickView;
template<class Eng, class Components> struct CompFwd;
template<class A1,class A2> struct IndirectionTag;
struct ConstantFunction;
template<class Expr> struct ExpressionTag;
template<class Stencil, class Expression> struct StencilEngine;


/**
 * Compressible<Engine>
 * CompressedRead<Engine>
 * CompressedReadWrite<Engine>
 *
 * These functors wrap the functions necessary to perform compressed
 * evaluation with compressible bricks.  This layer is necessary, because there
 * are engines that can contain compressible brick engines, which should
 * therefore work with compressed evaluation.  We don't want to force those
 * engines to support the compressible brick member functions, however, so we
 * perform all the queries through this functor class.
 * 
 * These functors are separate entities, because not all three are implemented
 * for all engines.  You should query Compressible<Engine>::compressible before
 * even attempting to instantiate the other two classes.
 *
 * The required interface of a specialization of Compressible<Engine> is:
 *  - enum { compressible = ? }        - True if the engine is compressible.
 *  - bool compressed(const Engine &)  - True if the engine is compressed.
 *
 * If the engine supports compressed reads then CompressedRead<Engine> has:
 *
 * T compressedRead(const Engine &) - returns the compressed value
 *
 * If the engine can be written to, then it must support the following
 * interface in CompressedReadWrite<Engine>:
 *
 * T& compressedReadWrite(const Engine &) - a reference to the compressed value
 * bool compressedBrickIsWholeView(const Engine &)
 * typedef ... ViewEngine_t; - Tag for efficient uncompressed view engine.
 * ViewEngine_t viewEngine(const Engine &) - an uncompressed view.
 */

struct Compressible
{
  typedef AndCombine Combine_t;
};

struct Compressed
{
  typedef AndCombine Combine_t;
};

struct CompressedRead
{
  typedef OpCombine Combine_t;
};

struct CompressedReadWrite
{
  //  typedef NullCombine Combine_t;
};

struct CompressedBrickIsWholeView
{
  //  typedef NullCombine Combine_t;
};

struct UnCompressedViewEngine
{
  //  typedef NullCombine Combine_t;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<int A, int B, class Op>
struct Combine2<WrappedInt<A>, WrappedInt<B>, Op, AndCombine>
{
  enum { val = A && B };
  typedef WrappedInt<val> Type_t;
  inline static
  Type_t combine(WrappedInt<A> , WrappedInt<B>, AndCombine)
  {
    return Type_t();
  }
};

//-----------------------------------------------------------------------------
// PETE functors to check the compression status of an expression.
//-----------------------------------------------------------------------------

/**
 * Compressible is used to query if an engine can be compressed, which
 * is a compile time trait.
 * ForEach<Expr,CompressibleTag,AndCombine>::Type_t is WrappedInt<true>
 * if Expr contains engines that are compressible.
 * (The apply() members here are never used.  Maybe they're not necessary.)
 */

/**
 * General Nodes are scalars which are compressible.
 */

template<class T>
struct EngineFunctorScalar<T, Compressible >
{
  typedef WrappedInt<true> Type_t;
  static Type_t apply(const T &, const Compressible &)
  {
    return Type_t();
  }
};

/**
 * Scalars are always compressed.
 */

template<class T>
struct EngineFunctorScalar<T, Compressed >
{
  typedef bool Type_t;
  static Type_t apply(const T &, const Compressed &)
  {
    return true;
  }
};

/**
 * PETE functors to perform compressed assignment.
 */

template<class T>
struct EngineFunctorScalar<T, CompressedRead >
{
  typedef T Type_t;
  static inline
  Type_t apply(const T& s, const CompressedRead &)
  {
    return s;
  }
};

/**
 * General engines are not compressible, so the general Compressible returns
 * false for compressed().  It is an error to attempt to perform compressed
 * reads of general engines, so your code should never attempt this.  (Make
 * a compile time switch first based on
 * ForEach<Expr,CompressibleTag,AndCombine>::Type_t::val, which is true
 * only if the expression contains engines that recognize compressibility.)
 */

template<class Engine>
struct EngineFunctorDefault<Engine,Compressible>
{
  typedef WrappedInt<false> Type_t;
};

template<class Engine>
struct EngineFunctorDefault<Engine,Compressed>
{
  typedef bool  Type_t;

  static inline
  Type_t apply(const Engine &, const Compressed &)
  {
    return false;
  }
};

/**
 * CompressibleBricks are the simplest case since all the functions of
 * Compressible just mirror the engine's functions.
 */

template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrick>,Compressible>
{
  typedef WrappedInt<true> Type_t;
};

template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrick>,Compressed>
{
  typedef Engine<Dim,T,CompressibleBrick> Engine_t;

  typedef bool Type_t;

  static inline
  Type_t apply(const Engine_t &e, const Compressed &)
  {
    return e.compressed();
  }
};

template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrick>,CompressedRead>
{
  typedef Engine<Dim,T,CompressibleBrick> Engine_t;

  typedef T Type_t;

  static inline
  Type_t apply(const Engine_t &e, const CompressedRead &)
  {
    return e.compressedRead();
  }
};

template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrick>,CompressedReadWrite>
{
  typedef Engine<Dim,T,CompressibleBrick> Engine_t;

  typedef T &Type_t;

  static inline
  Type_t apply(const Engine_t &e, const CompressedReadWrite &)
  {
    return e.compressedReadWrite();
  }
};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrick>,
  CompressedBrickIsWholeView>
{
  typedef Engine<Dim,T,CompressibleBrick> Engine_t;

  typedef bool Type_t;

  static inline
  bool apply(const Engine_t &e, const CompressedBrickIsWholeView &)
  {
    return e.compressedBrickIsWholeView();
  }

};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrick>,UnCompressedViewEngine>
{
  typedef Engine<Dim,T,CompressibleBrick> Engine_t;

  typedef Engine<Dim,T,BrickView> Type_t;

  static inline
  Type_t apply(const Engine_t &e, const UnCompressedViewEngine &)
  {
    return Type_t(e);
  }
};
    
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrickView>,Compressible>
{
  typedef WrappedInt<true> Type_t;
};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrickView >,Compressed>
{
  typedef Engine<Dim,T,CompressibleBrickView > Engine_t;

  typedef bool Type_t;

  static inline
  Type_t apply(const Engine_t &e, const Compressed &)
  {
    return e.compressed();
  }
};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrickView >,CompressedRead>
{
  typedef Engine<Dim,T,CompressibleBrickView > Engine_t;

  typedef T Type_t;

  static inline
  Type_t apply(const Engine_t &e, const CompressedRead &)
  {
    return e.compressedRead();
  }
};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrickView >,
  CompressedReadWrite>
{
  typedef Engine<Dim,T,CompressibleBrickView > Engine_t;

  typedef T &Type_t;

  static inline
  Type_t apply(const Engine_t &e, const CompressedReadWrite &)
  {
    return e.compressedReadWrite();
  }
};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrickView >,
  CompressedBrickIsWholeView>
{
  typedef Engine<Dim,T,CompressibleBrickView > Engine_t;

  typedef bool Type_t;

  static inline
  bool apply(const Engine_t &e, const CompressedBrickIsWholeView &)
  {
    return e.compressedBrickIsWholeView();
  }

};

template<int Dim, class T>
struct EngineFunctor<Engine<Dim,T,CompressibleBrickView >,
  UnCompressedViewEngine>
{
  typedef Engine<Dim,T,CompressibleBrickView > Engine_t;

  typedef Engine<Dim,T,BrickView> Type_t;

  static inline
  Type_t apply(const Engine_t &e, const UnCompressedViewEngine &)
  {
    return Type_t(e);
  }
};
    
/**
 * Constant function engine is definitely compressed. (You can't write to it,
 * though.)
 */

template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,ConstantFunction>,Compressible>
{
  typedef WrappedInt<true> Type_t;

  static inline
  Type_t apply(const Engine<Dim,T,ConstantFunction> &, const Compressible &)
  {
    return WrappedInt<true>();
  }
};
  
template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,ConstantFunction>,Compressed>
{
  typedef bool Type_t;

  static inline
  Type_t apply(const Engine<Dim,T,ConstantFunction> &, const Compressed &)
  {
    return true;
  }
};
  
template<int Dim,class T>
struct EngineFunctor<Engine<Dim,T,ConstantFunction>,CompressedRead>
{
  typedef T Type_t;

  static inline
  Type_t apply(const Engine<Dim,T,ConstantFunction> &e,
	       const CompressedRead &)
  {
    return e.constant();
  }
};
  
/**
 * Component forward engine can be compressed if it contains compressed
 * engines.
 */

template<int Dim, class T, class Eng, class Components>
struct EngineFunctor<Engine<Dim,T,CompFwd<Eng, Components> >,Compressible>
{
  typedef typename EngineFunctor<Eng,Compressible>::Type_t Comp_t;
  enum { compressible = Comp_t::val };
  typedef WrappedInt<compressible> Type_t;
};

template<int Dim, class T, class Eng, class Components>
struct EngineFunctor<Engine<Dim,T,CompFwd<Eng, Components> >,CompressedRead>
{
  typedef Engine<Dim,T,CompFwd<Eng, Components> > Engine_t;
  typedef typename Engine_t::CompAccess_t CompAccess_t;
  typedef typename CompAccess_t::Element_t Type_t;

  static inline
  Type_t apply(const Engine_t &e, const CompressedRead &tag)
  {
    return CompAccess_t::index(engineFunctor(e.elemEngine(),tag),
			       e.components());
  }
};

template<int Dim, class T, class Eng, class Components>
struct EngineFunctor<Engine<Dim,T,CompFwd<Eng, Components> >,CompressedReadWrite>
{
  typedef Engine<Dim,T,CompFwd<Eng, Components> > Engine_t;
  typedef typename Engine_t::CompAccess_t CompAccess_t;
  typedef typename CompAccess_t::ElementRef_t Type_t;

  static inline
  Type_t apply(const Engine_t &e, const CompressedReadWrite &tag)
  {
    return CompAccess_t::indexRef(engineFunctor(e.elemEngine(),tag),
				  e.components());
  }
};

template<int Dim, class T, class Eng, class Components>
struct EngineFunctor<Engine<Dim,T,CompFwd<Eng, Components> >,
  UnCompressedViewEngine>
{
  typedef Engine<Dim,T,CompFwd<Eng, Components> > Engine_t;

  typedef typename EngineFunctor<Eng,
    UnCompressedViewEngine>::Type_t CompEngine_t;
  typedef Engine<Dim,T,CompFwd<CompEngine_t, Components> > Type_t;

  static inline
  Type_t apply(const Engine_t &e,const UnCompressedViewEngine &tag)
  {
    return Type_t(engineFunctor(e.elemEngine(), tag), e.components());
  }
};
  
#endif     // POOMA_EVALUATOR_COMPRESSIBLEENGINES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CompressibleEngines.h,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
