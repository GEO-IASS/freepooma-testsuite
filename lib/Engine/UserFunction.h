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
//   UserFunction          - Base class for defining a user function.
//   UserFunctionEngine    - An tag for an engine for representing a function
//   Engine                - Specialization for UserFunctionEngine
//   NewEngine             - Specialization for UserFunctionEngine
//   IntersectEngine       - Specialization for UserFunctionEngine
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_USERFUNCTION_H
#define POOMA_ENGINE_USERFUNCTION_H

// namespace Pooma {

/** @file
 * @ingroup Engine
 * @brief
 * UserFunction objects are a way to build an object which applies a
 * function to an Array, and returns a new Array for the
 * expression.
 * 
 * This is the recommended way for users to make elementwise functions
 * apply to Arrays.
 *
 * UserFunction
 *
 *   A base class from which users would inherit to produce a specific
 *   functor.  This mainly implements operator()(expr), which
 *   constructs the expression with the function applied to the
 *   expression.
 *
 * UserFunctionEngine<D,T2,Expression>
 *
 *   An engine for Arrays which applies a user function.  This takes
 *   another engine as a template argument and applies the function to
 *   that engine.
 *
 * NewEngine
 *
 *   Defines the type of UserFunctionEngine you get when you subset it.
 *   It just subsets the engine inside of it.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include "Evaluator/EngineTraits.h"
#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Engine/IntersectEngine.h"
#include "Engine/DataObject.h"
#include "PETE/ErrorType.h"
#include "Pooma/View.h"
#include "Pooma/FunctorResult.h"

template<int D, class T, class E> class Array;


/**
 * This is just a tag class for the user function engine.
 * It is templated on:
 *  - UserFunction: The user function type.  This will be a class that
 *     inherits from UserFunction below.
 *  - Expression: The type of the expression to which the function
 *      is being applied.  This should be an Array.
 */

template<class UserFunction, class Expression> 
struct UserFunctionEngine 
{ };


/**
 * A specialization of Engine for UserFunctionEngine.
 *
 * This does all of the usual Engine things:
 *  - Typedefs for the tag, element types, domain and dimensions.
 *  - Operator() with integers to evaluate elements quickly.
 *  - Operator() with a domain to subset.
 *  - Accessor for the domain.
 */

template<int D, class T, class UserFunction, class Expression>
class Engine< D, T, UserFunctionEngine<UserFunction,Expression> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants
  //---------------------------------------------------------------------------

  typedef UserFunctionEngine<UserFunction,Expression>  Tag_t;
  typedef Engine<D,T,Tag_t>                  This_t;
  typedef This_t                             Engine_t;
  typedef typename Expression::Domain_t      Domain_t;
  typedef typename Expression::Layout_t      Layout_t;
  typedef T                                  Element_t;
  typedef ErrorType                          ElementRef_t;
  typedef typename Expression::Engine_t      ExprEngine_t;
  enum { dimensions = D };
  enum { hasDataObject = ExprEngine_t::hasDataObject };
  enum { dynamic = false };
  enum { zeroBased = ExprEngine_t::zeroBased };
  enum { multiPatch = ExprEngine_t::multiPatch };

  //---------------------------------------------------------------------------
  // Construct from a UserFunction object and an expression.
  //---------------------------------------------------------------------------

  Engine(const UserFunction& s, const Expression& e)
    : userFunction_m(s), expression_m(e) {}

  //---------------------------------------------------------------------------
  // Construct from a UserFunctionEngine and a domain.
  // Take a subset.
  //---------------------------------------------------------------------------

  template<class OtherE, class Domain>
  Engine(const Engine<D,T,UserFunctionEngine<UserFunction,OtherE> >& e,
	 const Domain& d)
    : userFunction_m(e.userFunction()), 
      expression_m( e.expression(), d ) {}

  //---------------------------------------------------------------------------
  // Element access via Loc and ints for speed.
  //---------------------------------------------------------------------------

  inline Element_t read(const Loc<D> &loc) const 
  {
    return userFunction_m(expression_m.read(loc));
  }
  inline Element_t read(int i) const 
  {
    return userFunction_m(expression_m.read(i));
  }
  inline Element_t read(int i, int j) const 
  {
    return userFunction_m(expression_m.read(i,j));
  }
  inline Element_t read(int i, int j, int k) const 
  {
    return userFunction_m(expression_m.read(i,j,k));
  }
  inline Element_t read(int i, int j, int k, int l) const 
  {
    return userFunction_m(expression_m.read(i,j,k,l));
  }
  inline Element_t read(int i, int j, int k, int l,int m) const 
  {
    return userFunction_m(expression_m.read(i,j,k,l,m));
  }
  inline Element_t read(int i, int j, int k, int l,int m,int n) const 
  {
    return userFunction_m(expression_m.read(i,j,k,l,m,n));
  }
  inline Element_t read(int i, int j, int k, int l,int m,int n,int o) const 
  {
    return userFunction_m(expression_m.read(i,j,k,l,m,n,o));
  }

  inline Element_t operator()(const Loc<D> &loc) const 
  {
    return userFunction_m(expression_m(loc));
  }
  inline Element_t operator()(int i) const 
  {
    return userFunction_m(expression_m(i));
  }
  inline Element_t operator()(int i, int j) const 
  {
    return userFunction_m(expression_m(i,j));
  }
  inline Element_t operator()(int i, int j, int k) const 
  {
    return userFunction_m(expression_m(i,j,k));
  }
  inline Element_t operator()(int i, int j, int k, int l) const 
  {
    return userFunction_m(expression_m(i,j,k,l));
  }
  inline Element_t operator()(int i, int j, int k, int l,int m) const 
  {
    return userFunction_m(expression_m(i,j,k,l,m));
  }
  inline Element_t operator()(int i, int j, int k, int l,int m,int n) const 
  {
    return userFunction_m(expression_m(i,j,k,l,m,n));
  }
  inline Element_t operator()(int i, int j, int k, int l,int m,int n,int o)
    const 
  {
    return userFunction_m(expression_m(i,j,k,l,m,n,o));
  }

  //---------------------------------------------------------------------------
  // Return the domain.
  //---------------------------------------------------------------------------

  inline const Domain_t& domain() const { return expression_m.domain(); }

  //---------------------------------------------------------------------------
  // Return first index in the specified direction.
  //---------------------------------------------------------------------------

  inline int first(int d) const
  {
    return expression_m.first(d);
  }
  
  //---------------------------------------------------------------------------
  // Accessors.
  //---------------------------------------------------------------------------

  const UserFunction& userFunction() const { return userFunction_m; }
  const Expression& expression() const { return expression_m; }

  //---------------------------------------------------------------------------
  // Need to pass lock requests to the contained engine.
  //---------------------------------------------------------------------------

  template<class RequestType>
  inline
  typename DataObjectRequest<RequestType>::Type_t
  dataObjectRequest(const DataObjectRequest<RequestType>& f) const
  {
    return engineFunctor(expression_m.engine(),f);
  }

private:
  UserFunction userFunction_m;
  Expression expression_m;
};


template<class Func> class UserFunction;

template<class Func,int D,class T,class E>
struct View1<UserFunction<Func>,Array<D,T,E> >
{
  typedef Array<D,T,E>                           Expr_t;
  typedef UserFunctionEngine<Func,Expr_t>        NewTag_t;
  typedef typename FunctorResult<Func,T>::Type_t NewT_t;
  typedef Engine<D,NewT_t,NewTag_t>              NewEngine_t;
  typedef Array<D,NewT_t,NewTag_t>               Type_t;
};


/**
 * To construct a user function class using UserFunction, define:
 *
 *   class MyUserFunction
 *
 * Give it the member function template:
 *
 *   template<class T> T operator()(const T&);
 * 
 * The input value is an element of an array, and the output is the
 * value from applying the user's function.
 *
 * Then UserFunction<MyUserFunction> can be applied to arrays.
 */

template<class Func>
class UserFunction
{
public:

  //---------------------------------------------------------------------------
  // Constructors
  //
  // UserFunction can be constructed using the default constructor, a function
  // objects, or from up to 7 arguments that are passed on to the function
  // object constructor.
  //---------------------------------------------------------------------------

  UserFunction()
  { }

  UserFunction(const Func &func)
    : function_m(func)
  { }

  template<class Init>
  UserFunction(const Init &init)
    : function_m(init)
  { }

  template<class I1, class I2>
  UserFunction(const I1 &i1, const I2 &i2)
    : function_m(i1,i2)
  { }

  template<class I1, class I2, class I3>
  UserFunction(const I1 &i1, const I2 &i2, const I3 &i3)
    : function_m(i1,i2,i3)
  { }

  template<class I1, class I2, class I3, class I4>
  UserFunction(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4)
    : function_m(i1,i2,i3,i4)
  { }

  template<class I1, class I2, class I3, class I4, class I5>
  UserFunction(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4,
	       const I5 &i5)
    : function_m(i1,i2,i3,i4,i5)
  { }

  template<class I1, class I2, class I3, class I4, class I5, class I6>
  UserFunction(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4,
	       const I5 &i5, const I6 &i6)
    : function_m(i1,i2,i3,i4,i5,i6)
  { }

  template<class I1, class I2, class I3, class I4, class I5, class I6,
    class I7>
  UserFunction(const I1 &i1, const I2 &i2, const I3 &i3, const I4 &i4,
	       const I5 &i5, const I6 &i6, const I7 &i7)
    : function_m(i1,i2,i3,i4,i5,i6,i7)
  { }

  //---------------------------------------------------------------------------
  // operator() 
  //---------------------------------------------------------------------------

  template<int D, class T, class E>
  typename View1<UserFunction<Func>,Array<D,T,E> >::Type_t
  operator()(const Array<D,T,E>& expr) const
  {
    typedef typename View1<UserFunction<Func>,Array<D,T,E> >::NewEngine_t
      NewEngine_t;
    typedef typename View1<UserFunction<Func>,Array<D,T,E> >::Type_t
      Type_t;
    return Type_t(NewEngine_t(function(),expr));
  }

  //---------------------------------------------------------------------------
  // accessors
  //---------------------------------------------------------------------------

  inline Func &function() { return function_m; }
  inline const Func &function() const { return function_m; }

private:

  // The function object.

  Func function_m;
};

/**
 * Specializations of NewEngine for subsetting a UserFunctionEngine with
 * an arbitrary domain.
 *
 * This just says that the subsetting operation is passed on to
 * the expression we're applying the function to.
 */

template <int Dim, class T, class S, class E, class Domain>
struct NewEngine< Engine<Dim,T,UserFunctionEngine<S,E> >, Domain >
{
  typedef typename View1<E,Domain>::Type_t               NewExpr_t;
  typedef UserFunctionEngine<S,NewExpr_t>                NewTag_t;
  typedef typename NewExpr_t::Element_t                  OldElement_t;
  typedef typename FunctorResult<S,OldElement_t>::Type_t NewElement_t;
  typedef Engine<Dim,NewElement_t,NewTag_t>              Type_t;
};

/**
 * Specializations for selecting the appropriate evaluator for the UserFunction
 * engine.  We just get the appropriate types from the Expression's engine.
 */

template<class UserFunction,class Expression>
struct EvaluatorEngineTraits<UserFunctionEngine<UserFunction,Expression> >
{
  typedef typename Expression::Engine_t Engine_t;
  typedef typename Engine_t::Tag_t Tag_t;
  typedef typename EvaluatorEngineTraits<Tag_t>::Evaluator_t Evaluator_t;
};

//---------------------------------------------------------------------------
// General version of engineFunctor to passes the request to
// the contained engine.
//---------------------------------------------------------------------------

template <int Dim, class T, class S, class E, class EFTag>
struct EngineFunctor<Engine<Dim, T, UserFunctionEngine<S, E> >, EFTag>
{
  typedef typename EngineFunctor<E, EFTag>::Type_t Type_t;

  static Type_t
  apply(const Engine<Dim, T, UserFunctionEngine<S, E> > &engine,
	const EFTag &tag)
  {
    return engineFunctor(engine.expression(), tag);
  }
};

template <int D, class T, class Func, class Expr, class Tag>
struct LeafFunctor<Engine<D, T, UserFunctionEngine<Func, Expr> >,
  EngineView<Tag> >
{
  typedef LeafFunctor<Expr, EngineView<Tag> > LeafFunctor_t;
  typedef typename LeafFunctor_t::Type_t NewViewed_t;
  typedef Engine<D, T, UserFunctionEngine<Func, NewViewed_t> > Type_t;

  static
  Type_t apply(const Engine<D, T, UserFunctionEngine<Func, Expr> > &engine,
	       const EngineView<Tag> &tag)
  {
    return Type_t(engine.userFunction(),
		  LeafFunctor_t::apply(engine.expression(), tag));
  }
};

template <int D, class T, class Func, class Expr, class Tag>
struct LeafFunctor<Engine<D, T, UserFunctionEngine<Func, Expr> >,
  ExpressionApply<Tag> >
{
  typedef LeafFunctor<Expr, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  static
  Type_t apply(const Engine<D, T, UserFunctionEngine<Func, Expr> > &engine,
	       const ExpressionApply<Tag> &tag)
  {
    return LeafFunctor_t::apply(engine.expression(), tag);
  }
};

#endif // POOMA_ENGINE_USERFUNCTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UserFunction.h,v $   $Author: richard $
// $Revision: 1.32 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
