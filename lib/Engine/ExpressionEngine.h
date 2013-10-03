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
//   DomainFunctorTag
//   LeafFunctor<DomainFunctorTag, T>
//   EvalLeaf
//   LeafFunctor<EvalLeaf<Domain>, Scalar<T> >
//   ViewFunctorTag
//   LeafFunctor< ViewFunctorTag<Domain>, Scalar<T> >
//   ExpressionTag<Expr>
//   Engine<Dim,T,ExpressionTag>
//   Combine2<Domain1, Domain2, Op, DomainFunctorTag>
//   NewEngine<Engine<Dim,T,ExpressionTag<Expr> >, Domain>
//   EngineFunctorTag
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_EXPRESSIONENGINE_H
#define POOMA_ENGINE_EXPRESSIONENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 * Engine<Dim,T,ExpressionTag<Expr> > (aka Expression-Engine) is the
 * engine that contains a PETE expression and provides all the Array
 * interfaces for it.
 */

#include "Domain/Loc.h"
#include "Domain/NullDomain.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include "Engine/Engine.h"
#include "Engine/DataObject.h"
#include "Engine/IntersectEngine.h"
#include "Engine/EngineFunctor.h"
#include "Engine/EnginePatch.h"
#include "PETE/PETE.h"
#include "Utilities/WrappedInt.h"

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template <int Dim>
class DomainLayout;

/**
 * EvalLeaf is used by Expression-Engine to evaluate itself at a point,
 * specified by a Loc<Dim>. We also supply a series of EvalLeaf<N> functors
 * for evaluating an engine from a set of indices.
 */

template<int Dim>
struct EvalLeaf { };

/**
 * We provide specializations for Scalars.
 */

template<class T, int Dim>
struct LeafFunctor<Scalar<T>, EvalLeaf<Dim> >
{
  typedef T Type_t;
  inline static
  Type_t apply(const Scalar<T> &s, const EvalLeaf<Dim> &) 
    {
      return s.value();
    }
};

/**
 * If an engine appears at the leaf, then we hand the engine to the EvalLeaf
 * object which determine whether or not to subtract off the first indices.
 */

template<int Dim, class T, class E>
struct LeafFunctor<Engine<Dim, T, E>, EvalLeaf<Dim> >
{
  typedef T Type_t;

  inline static
  Type_t apply(const Engine<Dim, T, E> &e, const EvalLeaf<Dim> &t) 
  {
    return t.eval(e);
  }
};

template<>
struct EvalLeaf<1>
{
  int i1_m;

  inline EvalLeaf(int i1) 
    : i1_m(i1)
  { }

  inline EvalLeaf(const Loc<1> &loc) 
    : i1_m(loc[0].first())
  { }

  inline int val1() const { return i1_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1());
  }
};

template<>
struct EvalLeaf<2>
{
  int i1_m, i2_m;

  inline EvalLeaf(int i1, int i2) 
    : i1_m(i1), i2_m(i2)
  { }

  inline EvalLeaf(const Loc<2> &loc) 
    : i1_m(loc[0].first()), i2_m(loc[1].first())
  { }

  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1(), val2());
  }
};

template<>
struct EvalLeaf<3>
{
  int i1_m, i2_m, i3_m;

  inline EvalLeaf(int i1, int i2, int i3) 
    : i1_m(i1), i2_m(i2), i3_m(i3)
  { }

  inline EvalLeaf(const Loc<3> &loc) 
    : i1_m(loc[0].first()), i2_m(loc[1].first()), i3_m(loc[2].first())
  { }

  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1(), val2(), val3());
  }
};

template<>
struct EvalLeaf<4>
{
  int i1_m, i2_m, i3_m, i4_m;

  inline EvalLeaf(int i1, int i2, int i3, int i4) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4)
  { }

  inline EvalLeaf(const Loc<4> &loc) 
    : i1_m(loc[0].first()), i2_m(loc[1].first()), i3_m(loc[2].first()),
      i4_m(loc[3].first())
  { }

  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1(), val2(), val3(), val4());
  }
};

template<>
struct EvalLeaf<5>
{
  int i1_m, i2_m, i3_m, i4_m, i5_m;

  inline EvalLeaf(int i1, int i2, int i3, int i4, int i5) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5)
  { }

  inline EvalLeaf(const Loc<5> &loc) 
    : i1_m(loc[0].first()), i2_m(loc[1].first()), i3_m(loc[2].first()),
      i4_m(loc[3].first()), i5_m(loc[4].first())
  { }

  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  inline int val5() const { return i5_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1(), val2(), val3(), val4(), val5());
  }
};

template<>
struct EvalLeaf<6>
{
  int i1_m, i2_m, i3_m, i4_m, i5_m, i6_m;

  inline EvalLeaf(int i1, int i2, int i3, int i4, int i5, int i6) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5), i6_m(i6)
  { }

  inline EvalLeaf(const Loc<6> &loc) 
    : i1_m(loc[0].first()), i2_m(loc[1].first()), i3_m(loc[2].first()),
      i4_m(loc[3].first()), i5_m(loc[4].first()), i6_m(loc[5].first())
  { }

  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  inline int val5() const { return i5_m; }
  inline int val6() const { return i6_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1(), val2(), val3(), val4(), val5(), val6());
  }
};

template<>
struct EvalLeaf<7>
{
  int i1_m, i2_m, i3_m, i4_m, i5_m, i6_m, i7_m;

  inline EvalLeaf(int i1, int i2, int i3, int i4, int i5, int i6, int i7) 
    : i1_m(i1), i2_m(i2), i3_m(i3), i4_m(i4), i5_m(i5), i6_m(i6), i7_m(i7)
  { }

  inline EvalLeaf(const Loc<7> &loc) 
    : i1_m(loc[0].first()), i2_m(loc[1].first()), i3_m(loc[2].first()),
      i4_m(loc[3].first()), i5_m(loc[4].first()), i6_m(loc[5].first()),
      i7_m(loc[6].first())
  { }

  inline int val1() const { return i1_m; }
  inline int val2() const { return i2_m; }
  inline int val3() const { return i3_m; }
  inline int val4() const { return i4_m; }
  inline int val5() const { return i5_m; }
  inline int val6() const { return i6_m; }
  inline int val7() const { return i7_m; }
  
  template<class Engine>
  inline typename Engine::Element_t eval(const Engine &e) const
  {
    return e.read(val1(), val2(), val3(), val4(), val5(), val6(), val7());
  }
};

/**
 * NewEngine<Engine<Dim,T,ExpressionTag<Expr> >, Domain >::Type_t is
 * supposed to give the type of ExpressionEngine that you would get by
 * taking views of all of the expressions leaves based on the domain
 * Domain.  To accomplish this, we use the ViewFunctorTag, which
 * contains the domain.  All classes that can appear as leaves in
 * Expressions that can go into an Expression-Engine should specialize
 * a version of LeafFunctor for this tag that takes a view the leaf
 * using the domain provided. We provide a specialization of this sort
 * here for Scalars.
 */

template<class Domain>
struct ViewFunctorTag
{
  const Domain &domain_m;
  inline ViewFunctorTag(const Domain &domain) : domain_m(domain) { }
};

template<class T, class Domain>
struct LeafFunctor<Scalar<T>, ViewFunctorTag<Domain> >
{
  typedef Scalar<T> Type_t;
  inline static
  Type_t apply(const Scalar<T> &s, const ViewFunctorTag<Domain> &) 
  {
    return s;
  }
};


/**
 * DomainFunctorTag is a functor tag used to divine domains. We ask the leaves
 * for their types and use domain() to get the domain.  For scalars, there
 * is no domain, so return NullDomain.
 */

struct DomainFunctorTag { };

template<class T>
struct LeafFunctor<Scalar<T>, DomainFunctorTag>
{
  typedef NullDomain Type_t;
  inline static
  Type_t apply(const Scalar<T> &, const DomainFunctorTag &) 
  {
    return NullDomain();
  }
};

template<class T>
struct LeafFunctor<T, DomainFunctorTag>
{
  typedef typename T::Domain_t Type_t;
  inline static
  Type_t apply(const T &leaf, const DomainFunctorTag &) 
  {
    return leaf.domain();
  }
};

//-----------------------------------------------------------------------------
// Here we supply combiners that work with DomainFunctorTag to combine
// domains. 
//
// These combiners are incomplete.  Right now, they just return the leftmost
// domain in an expression.  Eventually we want to add some kind of runtime
// check here that checks to see if the domains are compatible.
//-----------------------------------------------------------------------------

/**
 * In general adding two things on different domain types
 * is an error.
 */

template<class Domain1, class Domain2, class Op>
struct Combine2<Domain1, Domain2, Op, DomainFunctorTag>
{
  typedef ErrorDomain Type_t;
  inline static
  Type_t combine(const Domain1 &, const Domain2 &, const DomainFunctorTag &) 
  {
    return ErrorDomain();
  }
};

/**
 * Anything on a domain can be added to the something on the
 * same domain.
 */

template<class Domain, class Op>
struct Combine2<Domain, Domain, Op, DomainFunctorTag>
{
  typedef Domain Type_t;
  inline static
  Type_t combine(const Domain &a, const Domain &, const DomainFunctorTag &) 
  { 
    return a; 
  }
};

/**
 * Something plus a scalar should return the domain type
 * of the something.
 */

template<class Domain, class Op>
struct Combine2<Domain, NullDomain, Op, DomainFunctorTag>
{
  typedef Domain Type_t;
  inline static
  Type_t combine(const Domain &a, const NullDomain &,
		 const DomainFunctorTag &) 
  { 
    return a; 
  }
};

template<class Domain,class Op>
struct Combine2<NullDomain, Domain, Op, DomainFunctorTag>
{
  typedef Domain Type_t;
  inline static
  Type_t combine(const NullDomain &, const Domain &b, 
		 const DomainFunctorTag &) 
    { 
      return b; 
    }
};

/**
 * EngineFunctorTag is used to apply an EngineFunctor to an expression
 * containing several engines.
 */

template<class Tag>
struct EngineFunctorTag
{
  EngineFunctorTag() : tag_m() { }
  EngineFunctorTag(const Tag &tag) : tag_m(tag) { }
  inline const Tag &tag() const
  {
    return tag_m;
  }
  Tag tag_m;
};

/**
 * Tag class used to encode the type of an Expression for Engine.
 */

template<class Expr>
struct ExpressionTag
{
  typedef Expr Expression_t;
};

/**
 * This engine stores the expression tree Expr and acts like an
 * engine that lets you look at the values of the expression as if it
 * were an ordinary Brick type engine.
 */

template<int Dim, class T, class Expr>
class Engine<Dim, T, ExpressionTag<Expr> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants. Engine_t is my type.  Element_t is the
  // return type of the PETE expression.
  // Domain_t is the type of domain for the entire 
  // expression. Tag_t is my tag. Finally, dimensions is the number of
  // dimensions (comes from Domain_t)

  typedef Engine<Dim, T, ExpressionTag<Expr> > Engine_t;
  typedef ExpressionTag<Expr> Tag_t;
  typedef T Element_t;
  typedef ErrorType ElementRef_t;
  typedef typename ForEach<Expr, DomainFunctorTag, DomainFunctorTag>::Type_t
    Domain_t;
  typedef Expr Expression_t;
  typedef DomainLayout<Dim>   Layout_t;

  //---------------------------------------------------------------------------
  /// Was enum { dimensions = Domain_t::dimensions };
  /// It's possible for the dimension of an expression to be different from
  /// that of the domain.  For example, you can wrap a scalar in an array
  /// that has arbitrary dimension, but the domain of the scalar is
  /// NullDomain.

  enum { dimensions = Dim };

  //---------------------------------------------------------------------------
  /// Expressions might be multi-patch so we say they are to force code to check

  enum { multiPatch = true };

  //---------------------------------------------------------------------------
  /// We say that expression have a block so that functions that access data
  /// objects will call our message functor and we'll traverse the tree for them.

  enum { hasDataObject = true };
  enum { dynamic = false };

  //---------------------------------------------------------------------------
  /// Expression-engines are zero-based.
  
  enum { zeroBased = true };
  
  //---------------------------------------------------------------------------
  /// Expression constructor. Just stick the expression in local storage.

  inline Engine(const Expr &expr) : expr_m(expr),
    domain_m(forEach(expr_m, DomainFunctorTag(), DomainFunctorTag())) { }

  //---------------------------------------------------------------------------
  /// Copy constructor.

  inline Engine(const Engine_t &engine) : expr_m(engine.expression()),
    domain_m(engine.domain()) { }

  //---------------------------------------------------------------------------
  /// Subsetting Constructor. We build this expression engine, in place, from
  /// another expression engine and a domain. We pass a ViewFunctorTag since
  /// we will need to do some fiddling with the domain. Expression-engines are
  /// zero based, but can contain objects at their leaves that are not 
  /// zero-based. This means that when we get to the leaves, we must
  /// adjust the domain based on where the indices start for the leaf's engine.

  template<int Dim2, class T2, class Expr2, class Initializer>
  inline Engine(const Engine<Dim2, T2, ExpressionTag<Expr2> > &e, 
    const Initializer &i)
  : expr_m(e.expression(), i),
    domain_m(forEach(expr_m, DomainFunctorTag(), DomainFunctorTag()))
    { }

  template<int Dim2, class T2, class Expr2, class I1, class I2>
  inline Engine(const Engine<Dim2, T2, ExpressionTag<Expr2> > &e, 
                const I1 &i1, const I2 &i2)
    : expr_m(e.expression(), i1, i2),
      domain_m(forEach(expr_m, DomainFunctorTag(), DomainFunctorTag()))
  { }

  //---------------------------------------------------------------------------
  /// Construct from another expression without a domain.

  template<class Expr2>
  explicit inline Engine(const Engine<Dim,T,ExpressionTag<Expr2> > &e)
    : expr_m(e.expression()),
      domain_m(forEach(expr_m, DomainFunctorTag(), DomainFunctorTag()))
  { }

  //---------------------------------------------------------------------------
  /// Accessor functions that return the expression.

  inline const Expression_t &expression() const { return expr_m; }
  inline Expression_t &expression() { return expr_m; }

  //---------------------------------------------------------------------------
  /// Get a private copy of the expression

  Engine_t &makeOwnCopy();

  //---------------------------------------------------------------------------
  /// Accessor functions for a specific element. We recursively go through the
  /// expression tree asking each node and leaf to evaluate itself and combine
  /// results based on the OpCombine. 

  inline Element_t read(int i0) const 
    {
      return forEach(expr_m, EvalLeaf<1>(i0), OpCombine());
    }

  inline Element_t read(int i0, int i1) const 
    {
      return forEach(expr_m, EvalLeaf<2>(i0, i1), OpCombine());
    }

  inline Element_t read(int i0, int i1, int i2) const 
    {
      return forEach(expr_m, EvalLeaf<3>(i0, i1, i2), OpCombine());
    }

  inline Element_t read(int i0, int i1, int i2, int i3) const 
    {
      return forEach(expr_m, EvalLeaf<4>(i0, i1, i2, i3), 
        OpCombine());
    }

  inline Element_t read(int i0, int i1, int i2, int i3, int i4) const 
    {
      return forEach(expr_m, EvalLeaf<5>(i0, i1, i2, i3, i4), 
        OpCombine());
    }

  inline Element_t read(int i0, int i1, int i2, int i3, int i4, int i5) const 
    {
      return forEach(expr_m, EvalLeaf<6>(i0, i1, i2, i3, i4, i5), 
        OpCombine());
    }

  inline Element_t read(int i0, int i1, int i2, int i3, int i4, int i5,
    int i6) const 
    {
      return forEach(expr_m, EvalLeaf<7>(i0, i1, i2, i3, i4, i5, i6), 
        OpCombine());
    }

  inline Element_t read(const Loc<Dim> &loc) const 
  {
    return forEach(expr_m, EvalLeaf<Dim>(loc), OpCombine());
  }

  //---------------------------------------------------------------------------
  /// Function to return the common domain. We recursively go through the
  /// expression tree asking each node and leaf to return their domain and 
  /// combine the results based on the DomainFunctorTag. The DomainFunctorTag 
  /// combiners are above.

  inline const Domain_t& domain() const 
    {
      return domain_m;
    }

  //---------------------------------------------------------------------------
  /// Return a layout.

  inline Layout_t layout() const { return Layout_t(domain()); }

  //---------------------------------------------------------------------------
  /// Return the first value for the specified direction (always zero since this
  /// engine is zero-based).
  
  inline int first(int) const
  {
    return 0;
  }
  
  //---------------------------------------------------------------------------
  /// Need to pass lock requests to the leaves.

  template<class RequestType>
  inline
  typename DataObjectRequest<RequestType>::Type_t
  dataObjectRequest(const DataObjectRequest<RequestType>& f) const
  {
    typedef DataObjectRequest<RequestType> Tag_t;
    typedef EngineFunctorTag<Tag_t> Functor_t;
    typedef typename Tag_t::Combine_t Combine_t;

    return forEach(expr_m,Functor_t(f),Combine_t());
  }
 
private:

  /// The expression is stored here.

  Expr expr_m;

  // The domain of the expression.

  Domain_t domain_m;

};

/**
 * Here we supply the NewEngine traits class for Expression-Engine. We
 * go through the engine's expression recursively using the
 * ViewFunctorTag to divine the types that would result from taking
 * views of the leaves.  We use a TreeCombine to put these together
 * into an expression tree.
 */

template<int Dim, class T, class Expr, class Domain>
struct NewEngine<Engine<Dim, T, ExpressionTag<Expr> >, Domain>
{
  typedef ViewFunctorTag<Domain> FTag_t;
  typedef typename ForEach<Expr, FTag_t, TreeCombine>::Type_t ExprView_t;
  typedef Engine<Dim, T, ExpressionTag<ExprView_t> > Type_t;
};

template <int Dim, class T, class Expr, int SliceDim>
struct NewEngine<Engine<Dim, T, ExpressionTag<Expr> >, 
  SliceInterval<Dim,SliceDim> >
{
  typedef ViewFunctorTag<SliceInterval<Dim,SliceDim> > FTag_t;
  typedef typename ForEach<Expr, FTag_t, TreeCombine>::Type_t ExprView_t;
  typedef Engine<SliceDim, T, ExpressionTag<ExprView_t> > Type_t;
};

template <int Dim, class T, class Expr, int SliceDim>
struct NewEngine<Engine<Dim,T,ExpressionTag<Expr> >, 
  SliceRange<Dim,SliceDim> >
{
  typedef ViewFunctorTag<SliceRange<Dim,SliceDim> > FTag_t;
  typedef typename ForEach<Expr, FTag_t, TreeCombine>::Type_t ExprView_t;
  typedef Engine<SliceDim, T, ExpressionTag<ExprView_t> > Type_t;
};

//-----------------------------------------------------------------------------
// LeafFunctor specializations for EngineFunctorTag:

/**
 * EngineFunctorTag applied to a general node is an error.
 */

template<class Node, class Tag>
struct LeafFunctor<Node, EngineFunctorTag<Tag> >
{
};

/**
 * Scalar<T> specialization for LeafFunctor uses EngineFunctorScalar
 * to evaluate scalar leaf nodes.
 */

template<class T, class Tag>
struct LeafFunctor<Scalar<T>, EngineFunctorTag<Tag> >
{
  typedef typename EngineFunctorScalar<T,Tag>::Type_t Type_t;

  inline static
  Type_t apply(const Scalar<T> &scalar,  const EngineFunctorTag<Tag> &tag)
  {
    return EngineFunctorScalar<T,Tag>::apply(scalar.value(), tag.tag());
  }
};

/**
 * Engines could appear at expression leaves.
 */

template<int Dim, class T, class E, class Tag>
struct LeafFunctor<Engine<Dim, T, E>, EngineFunctorTag<Tag> >
{
  typedef Engine<Dim, T, E> Engine_t;
  typedef typename EngineFunctor<Engine_t, Tag>::Type_t Type_t;

  inline static
  Type_t apply(const Engine_t &engine, 
	       const EngineFunctorTag<Tag> &tag)
  {
    return EngineFunctor<Engine_t, Tag>::apply(engine, tag.tag());
  }
};

/**
 * EngineFunctors get applied to expressions using ForEach, with the
 * combiner supplied by the Tag.
 */

template<int Dim, class T, class Expr, class Tag>
struct EngineFunctor<Engine<Dim,T,ExpressionTag<Expr> >,Tag>
{
  typedef EngineFunctorTag<Tag> Functor_t;
  typedef typename Tag::Combine_t Combine_t;
  typedef typename ForEach<Expr,Functor_t,Combine_t>::Type_t Type_t;

  inline static 
  Type_t apply(const Engine<Dim, T, ExpressionTag<Expr> > &engine,
	       const Tag &tag)
  {
    return forEach(engine.expression(), Functor_t(tag), Combine_t());
  }
};

/**
 * EnginePatch is a special case because we want to package up the patched
 * expression in an expression engine.
 */

template<int Dim, class T, class Expr>
struct EngineFunctor<Engine<Dim, T, ExpressionTag<Expr> >, EnginePatch>
{
  typedef typename EnginePatch::Combine_t Combine_t;
  typedef typename ForEach<Expr, EnginePatch, Combine_t>::Type_t NewExpr_t;

  typedef Engine<Dim, T, ExpressionTag<NewExpr_t> > Type_t;

  inline static 
  Type_t apply(const Engine<Dim, T, ExpressionTag<Expr> > &engine,
	       const EnginePatch &tag)
  {
    return Type_t(forEach(engine.expression(), tag, Combine_t()));
  }
};

/**
 * EngineView acting on an expression returns an expression engine containing
 * the result of applying that functor to the leaves.
 */

template<int Dim, class T, class Expr, class Tag>
struct LeafFunctor<Engine<Dim, T, ExpressionTag<Expr> >, EngineView<Tag> >
{
  typedef EngineView<Tag> Functor_t;
  typedef typename Functor_t::Combine_t Combine_t;
  typedef typename ForEach<Expr, Functor_t, Combine_t>::Type_t NewExpr_t;

  typedef Engine<Dim, T, ExpressionTag<NewExpr_t> > Type_t;

  inline static 
  Type_t apply(const Engine<Dim, T, ExpressionTag<Expr> > &engine,
	       const EngineView<Tag> &tag)
  {
    return Type_t(forEach(engine.expression(), tag, Combine_t()));
  }
};

/**
 * ExpressionApply acting on an expression gets applied to the leaves.
 */

template<int Dim, class T, class Expr, class Tag>
struct LeafFunctor<Engine<Dim, T, ExpressionTag<Expr> >, ExpressionApply<Tag> >
{
  typedef int Type_t;

  inline static 
  Type_t apply(const Engine<Dim, T, ExpressionTag<Expr> > &engine,
	       const ExpressionApply<Tag> &tag)
  {
    return forEach(engine.expression(), tag, NullCombine());
  }
};

#endif // POOMA_ENGINE_EXPRESSIONENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ExpressionEngine.h,v $   $Author: richard $
// $Revision: 1.79 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
