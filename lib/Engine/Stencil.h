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
//   Stencil               - Base class for defining stencils
//   StencilEngine         - An tag for an engine for representing a stencil
//   View1                 - Specialization for Stencil
//   Engine                - Specialization for StencilEngine
//   NewEngine             - Specialization for StencilEngine
//   EvaluatorEngineTraits - Specialization for StencilEngine
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_STENCIL_H
#define POOMA_ENGINE_STENCIL_H

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

/** @file
 * @ingroup Engine
 * @brief
 * Stencil objects are a way to build an object which applies a stencil
 * to an array, and returns a new Array for the expression.
 * 
 * There are several reasons one might want to do this:
 * -# Abstraction.  Once a stencil like Laplace is constructed, you can say
 *    things like Laplace(A) to take the laplacian of A.  That way the 
 *    definition of the laplacian can be abstracted out and put in one
 *    place.  
 * -# Polymorphism. The Laplace object above can take different actions 
 *    depending on the type of A, giving it compile-time polymorphism.
 * -# Run-time efficiency.  Because the stencil object directly represents 
 *    what happens in the inner loop, more optimizations are available.
 *    Two particular ones are of greatest importance.  When an array appears
 *    several times in a stencil, it can recognize that the pointers
 *    are the same, saving registers, and the values of the integer offsets 
 *    from those pointers are visible and can be put in the instruction
 *    stream instead of registers.  Together, these two optimizations 
 *    allow a third, reusing values from the stencil from one loop
 *    iteration to the next.
 * -# Compile-time efficiency.  Stencil objects are much easier to compile
 *    than expression templates, so compilation goes much faster.
 *
 * Stencil
 *
 *   A base class from which users would inherit to produce a specific
 *   stencil.  This mainly implements operator()(expr), which constructs
 *   the expression with the stencil applied to the expression.
 *
 * StencilEngine<D,T2,Expression>
 *
 *   An engine for Arrays which applies a stencil.  This takes another
 *   engine as a template argument and applies the stencil to that 
 *   engine.
 *
 * NewEngine
 *
 *   Defines the type of StencilEngine you get when you subset it.
 *   It just subsets the engine inside of it.
 * 
 *
 * Stencil Concepts:
 * 
 * A stencil is a pattern repeatedly applied to elements in an input domain to
 * yield elements in the output domain.  For example, the simplest
 * stencil copies each element in the input domain to exactly the same
 * element in the output domain.  A second-order difference stencil can
 * be represented by the formula
 * 
 * 		 out(i) = 2 in(i-1) + in(i) + in(i+1)
 * 
 * where in(i) and out(i) indicate the ith input and output elements,
 * respectively.  This stencil illustrates that a stencil can use more
 * than one input element, but that all input elements must be
 * contiguous.
 * 
 * A StencilEngine is an engine applying a stencil to an input array.
 * When invoked, the result is an array filled with values from applying
 * the stencil to the input array.  We explain the engine's data members
 * and assumptions.  Even though a StencilEngine stores the data for its
 * computation, actually performing the computation only when requested,
 * we use the slang of its "output" to avoid writing "its output when the
 * computation is invoked."  Also, in the explanation below, we use
 * one-dimensional terminology.  The only supported domains and ranges
 * are Cartesian products so the one-dimensional terminology is easily
 * generalized.
 * 
 * When created, engines frequently are given the desired array output
 * range indices, e.g., -3, ..., 5.  Any such range can be shifted so the
 * leftmost element's index is zero, i.e., zero-based.  For example, 0,
 * ..., 8 with an offset of -3.  To return to the "original," desired
 * range, add the offset to each index.  The `domain_m' variable records
 * the number of output elements.
 * 
 * Assume the engine's stencil uses input array elements with indices
 * lowerExtent, lowerExtent+1, ..., 0, ..., upperExtent.  Thus, to
 * produce out(0) requires knowing in(lowerExtent), ..., in(upperExtent).
 * The input domain to consisting of the values used to compute the
 * zero-based output range is in(lowerExtent), ..., in(domain_m +
 * upperExtent).
 * 
 * The StencilEngine's data members are
 *  -# function_m representing the stencil
 *  -# expression_m which is approximately the input
 *  -# domain_m representing the indices for the output
 *  -# offset_m representing the 'shift' to yield zero-based output indices
 *
 * Note all members concern output, not input.
 * 
 * When reading the source code below, "domain" is used for both input
 * and output indices.  The reader must decide the meaning of each
 * occurrence.
 */ 

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Evaluator/EngineTraits.h"
#include "Engine/Engine.h"
#include "Engine/Intersector.h"
#include "Engine/IntersectEngine.h"
#include "Engine/EngineFunctor.h"
#include "PETE/ErrorType.h"
#include "Pooma/View.h"
#include "Pooma/FunctorResult.h"
#include "Engine/ViewEngine.h"
#include "Utilities/WrappedInt.h"

template<int D, class T, class E> class Array;
template<class M, class T, class E> class Field;
template<class ST> class Stencil;

/**
 * This is just a tag class for the stencil engine.
 * It is templated on:
 *  - Stencil:  The stencil type.  This will be a class that inherits
 *      from Stencil below.
 *  - Expression: The type of the expression to which the stencil
 *      is being applied.  This should be an Array<...>.
 *
 * This defines the type:
 *  - Element_t: The type of each element that is output from the stencil.
 *      This defaults to the same type as the expression, and should be
 *      specialized to something else if that is not the case.
 *  - ElementRef_t: A type to be used for referring to elements in a 
 *      stencil.  This will only have meaning if the stencil does something
 *      like selects a component from a vector.  This will not be a common
 *      case, so it is not defined by default.
 */

template<class Function, class Expression> 
struct StencilEngine 
{
  typedef typename Expression::Element_t                      ViewedElement_t;
  typedef typename FunctorResult<Function,ViewedElement_t>::Type_t  Element_t;
};

/// insetDomain() computes the insetDomain domain of the stencil for users
/// (it's not zero-based).  If you just got a random stencil from who knows
/// where and want to apply it to another array, you could say:
///
/// b(st.insetDomain(a.domain())) = st(a);
///
/// Note that you can always just say:
///
/// b(range) = st(a,range);
///
/// because that version doesn't inset.
///
/// In other words, given a stencil and an input domain, return the
/// resulting output indices.

template<class Function, int D>
inline
Interval<D> insetDomain(const Function &f, const Interval<D> &domain)
{
  Interval<D> ret;
  int d;
  for (d = 0; d < D; ++d)
  {
    ret[d] = Interval<1>(domain[d].first() + f.lowerExtent(d),
			 domain[d].last() - f.upperExtent(d));
  }
  return ret;    
}


/**
 * A specialization of Engine for StencilEngine.
 *
 * This does all of the usual Engine things:
 *  - Typedefs for the tag, element types, domain and dimensions.
 *  - Operator() with integers to evaluate elements quickly.
 *  - Operator() with a domain to subset.
 *  - Accessor for the domain.
 */

template<int D, class T, class Function, class Expression>
class Engine<D, T, StencilEngine<Function, Expression> >
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef StencilEngine<Function, Expression>   Tag_t;
  typedef Engine<D, T, Tag_t>                   This_t;
  typedef This_t                                Engine_t;
  typedef Interval<D>                           Domain_t;
  typedef DomainLayout<D>		        Layout_t;
  typedef T                                     Element_t;
  typedef ErrorType                             ElementRef_t;
  typedef typename Expression::Engine_t         ExprEngine_t;

  enum { dimensions = D };
  enum { hasDataObject = ExprEngine_t::hasDataObject };
  enum { dynamic = false };
  enum { multiPatch = ExprEngine_t::multiPatch };
  enum { zeroBased = true };

  // FIXME: using any of the two below disables using of
  // expression engines as Expression type, because these
  // are not default-constructible.
  // Only FieldEngine ever default-constructs these, though.

  Engine()
    : function_m(), expression_m(), domain_m(Pooma::NoInit())
  {
  }

  template <class Layout2>
  explicit Engine(const Layout2 &layout)
    : function_m(), expression_m(), domain_m(layout.domain())
  {
  }

  //============================================================
  // Construct from a Function object (effectively a stencil)
  // and an expression (effectively the input array), and
  // sometimes an output (not input) domain.
  //============================================================

  Engine(const Function &f, const Expression &e)
    : function_m(f), expression_m(e), domain_m(Pooma::NoInit())
  {
    // inset is the indices for the stencil's output.
    Interval<D> inset = insetDomain(f, e.domain());
    int d;
    for (d = 0; d < D; ++d)
    {
      domain_m[d] = Interval<1>(inset[d].length());
      offset_m[d] = function().lowerExtent(d);
    }
  }

  Engine(const Function &f, const Expression &e, const Interval<D> &domain)
    : function_m(f), expression_m(e), domain_m(Pooma::NoInit())
  {
    int d;
    for (d = 0; d < D; ++d)
    {
      domain_m[d] = Interval<1>(domain[d].length());
      offset_m[d] = domain[d].first();
    }
  }

  // Construct an engine for composing stencils, e.g.,
  // stencil1(stencil2(array)).
  template<class OtherE>
  Engine(const Engine<D, T, StencilEngine<Function, OtherE> > &e,
	 const INode<D> &node)
    : function_m(e.function()),
      expression_m(e.expression()(e.viewDomain(node))),
      domain_m(Pooma::NoInit())
  {
    int d;
    for (d = 0; d < D; ++d)
    {
      domain_m[d] = Interval<1>(node.domain()[d].length());
      offset_m[d] = function().lowerExtent(d);
    }
  }

  Engine(const Engine_t &e,
	 const Interval<D> &domain)
    : function_m(e.function()),
      expression_m(e.expression()),
      domain_m(Pooma::NoInit())
  {
    int d;
    for (d = 0; d < D; ++d)
    {
      domain_m[d] = Interval<1>(domain[d].length());
      offset_m[d] = e.offset_m[d] + domain[d].first();
    }
  }

  template <int Dim, class Tx, class EngineTag>
  void initExpressionFromModel(const Array<Dim, Tx, EngineTag>& model)
  {
    expression_m.engine() = model.engine();
  }

  template <class Mesh, class Tx, class EngineTag>
  void initExpressionFromModel(const Field<Mesh, Tx, EngineTag>& model)
  {
    expression_m.fieldEngine() = model.fieldEngine();
  }

  This_t &operator=(const This_t &model)
  {
    domain_m = model.domain();
    function_m = model.function();
    initExpressionFromModel(model.expression());
    for (int d = 0; d < D; ++d)
    {
      domain_m[d] = model.domain()[d];
      offset_m[d] = model.offset(d);
    }
    return *this;
  }

  //============================================================
  // Element access via ints for speed.  The arguments correspond to
  // output elements, not input elements.
  //============================================================

  inline Element_t read(int i) const 
  {
    // Input index `i + offset_m[0]' corresponds to output index `i'.
    return function()(expression_m,
		      i + offset_m[0]);
  }
  inline Element_t read(int i, int j) const 
  {
    return function()(expression_m,
		      i + offset_m[0],
		      j + offset_m[1]);
  }
  inline Element_t read(int i, int j, int k) const 
  {
    return function()(expression_m,
		      i + offset_m[0],
		      j + offset_m[1],
		      k + offset_m[2]);
  }

  inline Element_t read(const Loc<1> &loc) const 
  {
    return function()(expression_m,
		      loc[0].first() + offset_m[0]);
  }
  inline Element_t read(const Loc<2> &loc) const 
  {
    return function()(expression_m,
		      loc[0].first() + offset_m[0],
		      loc[1].first() + offset_m[1]);
  }
  inline Element_t read(const Loc<3> &loc) const 
  {
    return function()(expression_m,
		      loc[0].first() + offset_m[0],
		      loc[1].first() + offset_m[1],
		      loc[2].first() + offset_m[2]);
  }

  //============================================================
  // operator() are provided since users typically write stencils
  // as x(i, j) + x(i, j - 1), so for stencils of stencils to work
  // the engine needs this interface.

  inline Element_t operator()(int i) const
  {
    return read(i);
  }

  inline Element_t operator()(int i, int j) const
  {
    return read(i, j);
  }

  inline Element_t operator()(int i, int j, int k) const
  {
    return read(i, j, k);
  }

  //============================================================
  // Return the output domain.
  //============================================================

  inline const Domain_t &domain() const { return domain_m; }

  //============================================================
  // Return the output layout.
  //============================================================

  inline Layout_t layout() const { return Layout_t(domain_m); }

  //============================================================
  // Return the first output index value for the specified direction
  // (always zero since this engine is zero-based).
  //============================================================
  
  inline int first(int i) const
  {
#if POOMA_BOUNDS_CHECK
    PAssert(i >= 0 && i < D);
#endif
    return 0;
  }

  //---------------------------------------------------------------------------
  // viewDomain() gives the region of the expression needed to compute a given
  // region of the stencil.  That is, viewDomain(outputDomain) yields
  // the corresponding input domain.
  //---------------------------------------------------------------------------

  inline
  Interval<D> viewDomain(const Interval<D> &domain) const
  {
    Interval<D> ret;
    int d;
    for (d = 0; d < D; ++d)
    {
      // The computation subtracts and adds the stencil's extent from
      // the "original", unshifted output domain.
      ret[d] =
	Interval<1>(
		    domain[d].first() + offset_m[d]
		    - function().lowerExtent(d),
		    domain[d].last() + offset_m[d] + function().upperExtent(d)
		    );
    }
    return ret;
  }

  inline
  INode<D> viewDomain(const INode<D> &inode) const
  {
    return INode<D>(inode, viewDomain(inode.domain()));
  }

  //---------------------------------------------------------------------------
  // intersectDomain() gives the "original", unshifted output domain.
  //---------------------------------------------------------------------------

  inline
  Interval<D> intersectDomain() const
  {
    Interval<D> ret;
    int d;
    for (d = 0; d < D; ++d)
    {
      ret[d] =
	Interval<1>(
		    domain_m[d].first() + offset_m[d],
		    domain_m[d].last() + offset_m[d]
		    );
    }

    return ret;
  }

  //============================================================
  // Accessors.
  //============================================================

  inline const Function   &function() const   { return function_m; }
  inline const Expression &expression() const { return expression_m; }
  int offset(int d) const { return offset_m[d]; }

private:

  Function    function_m;
  Expression  expression_m;
  Interval<D> domain_m;
  int         offset_m[D];
};

//-----------------------------------------------------------------------------
// View types for stencil objects.  Stencils define operator() to return a
// stencil engine object, which, when invoked, yields the result of
// applying the stencil to the given array.
//
// If you wanted to store that object, you could write:
//
// A a;
// Laplace laplace;
// typename View1<Stencil<Laplace>,A>::Type_t b = laplace(a);
//-----------------------------------------------------------------------------

template<class Function, int D, class T, class E>
struct View1<Stencil<Function>, Array<D, T, E> >
{
  typedef Array<D,T,E> ArrayIn_t;
  typedef StencilEngine<Function, ArrayIn_t> NewTag_t;
  typedef typename NewTag_t::Element_t NewT_t;
  typedef Engine<D, NewT_t, NewTag_t> NewEngine_t;
  typedef Array<D, NewT_t, NewTag_t> Type_t;

  static inline
  Type_t make(const Stencil<Function> &s, const ArrayIn_t &a)
  {
    return Type_t(NewEngine_t(s.function(), a,
			      insetDomain(s.function(), a.domain())
			      ));
  }
};

//-----------------------------------------------------------------------------
// View2 is used to construct the return type for stencils where the
// output domain is given as well.
//-----------------------------------------------------------------------------

template<class Function, class ArrayIn, int Dim>
struct View2<Stencil<Function>,ArrayIn,Interval<Dim> >
{
  enum { dim = ArrayIn::dimensions };
  typedef Interval<Dim> ViewDom_t;
  typedef typename View1<ArrayIn,ViewDom_t>::Type_t Expression_t;
  //  typedef StencilEngine<Function, ArrayIn> NewTag_t;
  typedef StencilEngine<Function, Expression_t> NewTag_t;
  typedef typename NewTag_t::Element_t NewT_t;
  typedef Engine<dim,NewT_t,NewTag_t> NewEngine_t;
  typedef Array<dim,NewT_t,NewTag_t> Type_t;

  static inline
  Type_t make(const Stencil<Function> &s, const ArrayIn &a,
	      const Interval<Dim> &d)
  {
    return Type_t(NewEngine_t(s.function(), a(s.inputDomain(d))));
    //    return Type_t(NewEngine_t(s.function(), a, d));
  }
};

template<class Function, class ArrayIn, class Dom>
struct View2<Stencil<Function>, ArrayIn, Dom>
{
  enum { dim2 = ArrayIn::dimensions };
  enum { dim = Dom::dimensions };
  typedef Interval<dim2> ViewDom_t;
  typedef typename View1<ArrayIn, ViewDom_t>::Type_t Expression_t;
  typedef StencilEngine<Function, Expression_t> StencilTag_t;
  //  typedef StencilEngine<Function, ArrayIn> StencilTag_t;
  typedef typename StencilTag_t::Element_t NewT_t;
  typedef ViewEngine<dim2, StencilTag_t> NewTag_t;
  typedef Engine<dim,NewT_t,NewTag_t> NewEngine_t;
  typedef Array<dim,NewT_t,NewTag_t> Type_t;

  static inline
  Type_t make(const Stencil<Function> &s, const ArrayIn &a, const Dom &dom)
  {
    ViewDom_t viewDom = s.inputDomain(dom);
    ViewDom_t insetDom = insetDomain(s.function(), viewDom);
    ViewIndexer<dim,dim2> indexer(insetDom);
    Dom stView;
    indexer.baseToLocal(dom,stView);

    typedef typename NewEngine_t::ViewedEngine_t ViewedEngine_t;
    ViewedEngine_t viewed(s.function(),a(viewDom));
    //    ViewedEngine_t viewed(s.function(), a, insetDom);

    return Type_t(NewEngine_t(viewed,stView));
  }

};


/**
 * Stencil
 *
 * Base class for constructing stencil classes.  To construct 
 * a stencil class using Stencil, define:
 *
 *   class MyStencil
 *
 * Give it the member function template:
 *
 *   template<class A> T operator()(const A& expr, int, int, ...) const;
 * 
 * The argument 'expr' is the type of the expression the stencil 
 * is being applied to.  This will generally be some kind of Array.
 * The integer arguments have the location at which the stencil
 * is being applied.  (The const is important.  The stencil may be passed
 * to the evaluator as a const reference.)
 * 
 * The return type is whatever the stencil outputs.  If this is not
 * the same type as the elements of 'expr', you must specialize
 * the Pooma FunctorResult class (see Pooma/FunctorResult.h).
 *
 * To apply a stencil, create an instance of the Stencil<> class.
 *
 *   Stencil<MyStencil> myStencil;
 *
 * This class really only does one thing: defines operator()(expr),
 * and operator()(expr,domain).
 * When given an expression it wraps it in a StencilEngine and builds
 * an array with that engine, so that you can write:
 *
 * <PRE>
 *   b = myStencil(a);
 *   b(dom) = myStencil(a,dom);
 * </PRE>
 */

template<class Function>
class Stencil
{
public:

  Stencil()
  { }

  Stencil(const Stencil<Function> &model)
    : function_m(model.function_m)
  { }

  ~Stencil() {}

#if 0 // Unnecessary, but triggers CW 6 bug - SWH
  Stencil(const Function &function)
    : function_m(function)
  { }
#endif

  template<class Init>
  Stencil(const Init &init)
    : function_m(init)
  { }

  /// @name Array apply
  //@{

  template<int D, class T, class E>
  typename View1<Stencil<Function>,Array<D,T,E> >::Type_t
  operator()(const Array<D,T,E>& expr) const
  {
    typedef View1<Stencil<Function>,Array<D,T,E> > Ret_t;
    return Ret_t::make(*this,expr);
  }

  template<int D, class T, class E,class Dom>
  typename View2<Stencil<Function>,Array<D,T,E>,Dom>::Type_t
  operator()(const Array<D,T,E>& expr,const Dom &domain) const
  {
    CTAssert(D==Dom::dimensions); // maybe unnecessary
    typedef View2<Stencil<Function>,Array<D,T,E>,Dom> Ret_t;
    return Ret_t::make(*this,expr,domain);
  }

  //@}

  template<int D>
  inline
  Interval<D> insetDomain(const Interval<D> &domain)
  {
    return ::insetDomain(function(), domain);
  }

  //---------------------------------------------------------------------------
  // inputDomain() gives the region required to compute the stencil values on
  // a given subregion.
  //---------------------------------------------------------------------------

  template<int D, class DT>
  inline
  Interval<D> inputDomain(const Domain<D,DT> &domain) const
  {
    Interval<D> ret;

    int d;
    for (d = 0; d < D; ++d)
    {
      ret[d] =
	Interval<1>(
		    domain[d].first() - function().lowerExtent(d),
		    domain[d].last()  + function().upperExtent(d)
		    );
    }
    return ret;
  }

  inline Function &function() { return function_m; }
  inline const Function &function() const { return function_m; }

private:

  Function function_m;
};


/**
 * Specializations of NewEngine for subsetting a StencilEngine with
 * an arbitrary domain.
 *
 * This just says that the subsetting operation is passed on to
 * the expression we're applying the stencil to.
 */

template <int Dim, class T, class S, class E>
struct NewEngine<Engine<Dim, T, StencilEngine<S, E> >, Interval<Dim> >
{
  typedef Engine<Dim, T, StencilEngine<S, E> > Type_t;
};

template <int Dim, class T, class S, class E>
struct NewEngine<Engine<Dim, T, StencilEngine<S, E> >, INode<Dim> >
{
  typedef typename View1<E, INode<Dim> >::Type_t NewExpr_t;
  typedef StencilEngine<S, NewExpr_t> NewTag_t;
  typedef Engine<Dim, T, NewTag_t> Type_t;
};

template <int Dim, class T, class S, class E>
struct NewEngine<Engine<Dim, T, StencilEngine<S, E> >, Range<Dim> >
{
  typedef StencilEngine<S,E> SETag_t;
  typedef ViewEngine<Dim,SETag_t> NewTag_t;
  typedef Engine<Dim,T,NewTag_t> Type_t;
};

template <int Dim, class T, class S, class E, int SliceDim>
struct NewEngine< Engine<Dim,T,StencilEngine<S,E> >,
  SliceInterval<Dim,SliceDim> >
{
  typedef StencilEngine<S,E> SETag_t;
  typedef ViewEngine<Dim,SETag_t> NewTag_t;
  typedef Engine<SliceDim,T,NewTag_t> Type_t;
};

template <int Dim, class T, class S, class E, int SliceDim>
struct NewEngine< Engine<Dim,T,StencilEngine<S,E> >,
  SliceRange<Dim,SliceDim> >
{
  typedef StencilEngine<S,E> SETag_t;
  typedef ViewEngine<Dim,SETag_t> NewTag_t;
  typedef Engine<SliceDim,T,NewTag_t> Type_t;
};

/**
 * Specializations for selecting the appropriate evaluator for the Stencil
 * engine.  We just get the appropriate types from the Expression's engine.
 */

template<class UserFunction,class Expression>
struct EvaluatorEngineTraits<StencilEngine<UserFunction,Expression> >
{
  typedef typename Expression::Engine_t Engine_t;
  typedef typename Engine_t::Tag_t Tag_t;
  typedef typename EvaluatorEngineTraits<Tag_t>::Evaluator_t Evaluator_t;
};

//-----------------------------------------------------------------------------
// StencilIntersector is a special intersector that gets used when we come
// across a stencil object in an expression.
//-----------------------------------------------------------------------------

template<int Dim, class Intersect>
class StencilIntersector
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef typename Intersect::IntersectorData_t         IntersectorData_t;
  typedef StencilIntersector<Dim, Intersect>            This_t;
  typedef typename IntersectorData_t::const_iterator    const_iterator;
  typedef RefCountedPtr<IntersectorData_t>              DataPtr_t;
  
  enum { dimensions = Intersect::dimensions };
  
  StencilIntersector(const This_t &model)
    : domain_m(model.domain_m),
      stencilExtent_m(model.stencilExtent_m),
      intersector_m(model.intersector_m)
  { }

  StencilIntersector(const Interval<Dim> &domain, const Intersect &intersect,
		  const GuardLayers<Dim> &stencilExtent)
    : domain_m(domain),
      stencilExtent_m(stencilExtent),
      intersector_m(intersect)
  { }

  This_t &operator=(const This_t &model)
  {
    if (this != &model)
    {
      intersector_m = model.intersector_m;
      domain_m = model.domain_m;
      stencilExtent_m = model.stencilExtent_m;
    }
    return *this;
  }

  ~StencilIntersector() { }

  inline DataPtr_t &data() { return intersector_m.data(); }
  inline const DataPtr_t &data() const { return intersector_m.data(); }

  //===========================================================================
  // Accessors
  //===========================================================================

  // STL iterator support.
  
  inline const_iterator begin() const { return data()->inodes_m.begin(); }
  
  inline const_iterator end() const { return data()->inodes_m.end(); }
  

  //===========================================================================
  // Intersect routines
  //===========================================================================

  // All domains.
  
  template<class Engine>
  inline
  void intersect(const Engine &engine) 
  {
    typedef typename NewEngine<Engine,Interval<Dim> >::Type_t NewEngine_t;

    NewEngine_t newEngine(engine, domain_m);

    intersector_m.intersect(newEngine);

    data()->shared(engine.layout().ID(),newEngine.layout().ID());
  }

  template<class Engine, int Dim2>
  inline
  bool intersect(const Engine &engine, const GuardLayers<Dim2> &g,
		  GuardLayers<Dim> &usedGuards) 
  {
    intersect(engine);
    usedGuards = stencilExtent_m;
    return true;
  }

private:
  Interval<Dim> domain_m;
  GuardLayers<Dim> stencilExtent_m;
  Intersect     intersector_m;
};

//-----------------------------------------------------------------------------
// IntersectEngine specialization
//-----------------------------------------------------------------------------

template <int D, class T, class S, class E, class Intersect>
struct LeafFunctor<Engine<D,T,StencilEngine<S,E> >,
  ExpressionApply<IntersectorTag<Intersect> > >
{
  typedef int Type_t;

  static
  Type_t apply(const Engine<D,T,StencilEngine<S,E> > &engine,
	       const ExpressionApply<IntersectorTag<Intersect> > &tag)
  {
    typedef StencilIntersector<D, Intersect> NewIntersector_t;
    GuardLayers<D> stencilExtent;
    for (int i=0; i<D; ++i) {
      stencilExtent.lower(i) = engine.function().lowerExtent(i);
      stencilExtent.upper(i) = engine.function().upperExtent(i);
    }
    NewIntersector_t newIntersector(engine.intersectDomain(),
				    tag.tag().intersector_m,
				    stencilExtent);

    expressionApply(engine.expression(),
		    IntersectorTag<NewIntersector_t>(newIntersector));
    return 0;
  }
};

//---------------------------------------------------------------------------
// Specialization of DataObjectRequest engineFunctor to pass the request to
// the contained engine.
//---------------------------------------------------------------------------

template<class RequestType> class DataObjectRequest;

template <int D, class T, class S, class E, class RequestType>
struct EngineFunctor<Engine<D, T, StencilEngine<S, E> >,
  DataObjectRequest<RequestType> >
{
  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;

  static Type_t
  apply(const Engine<D, T, StencilEngine<S, E> > &engine,
	const DataObjectRequest<RequestType> &tag)
  {
    return engineFunctor(engine.expression(), tag);
  }
};

//-----------------------------------------------------------------------------
//
// The generic version of EngineView just accesses the contained engine and
// applies EngineView to it.
//
// The default version doesn't fiddle with the domain, since it is assumed
// that the typical view doesn't need to.  Specializations will be required
// for INode views etc...  Probably we should come up with a generic approach.
//
//-----------------------------------------------------------------------------

template <int D, class T, class S, class E, class Tag>
struct LeafFunctor<Engine<D, T, StencilEngine<S, E> >, EngineView<Tag> >
{
  typedef LeafFunctor<E, EngineView<Tag> > LeafFunctor_t;
  typedef typename LeafFunctor_t::Type_t NewViewed_t;
  typedef Engine<D, T, StencilEngine<S, NewViewed_t> > Type_t;

  static
  Type_t apply(const Engine<D, T, StencilEngine<S, E> > &engine,
	       const EngineView<Tag> &tag)
  {
    return Type_t(engine.function(), 
		  LeafFunctor_t::apply(engine.expression(), tag));
  }
};

//-----------------------------------------------------------------------------
//
// The generic version of ExpressionApply just accesses the contained engine
// and applies ExpressionApply to it.
//
//-----------------------------------------------------------------------------

template <int D, class T, class S, class E, class Tag>
struct LeafFunctor<Engine<D, T, StencilEngine<S, E> >, ExpressionApply<Tag> >
{
  typedef LeafFunctor<E, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  static
  Type_t apply(const Engine<D, T, StencilEngine<S, E> > &engine,
	       const ExpressionApply<Tag> &tag)
  {
    return LeafFunctor_t::apply(engine.expression(), tag);
  }
};

#endif // POOMA_ENGINE_STENCIL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Stencil.h,v $   $Author: richard $
// $Revision: 1.55 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
