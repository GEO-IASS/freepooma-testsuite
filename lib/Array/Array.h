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
//   Array
//   View0
//   View[1-7]<Array,various domains>
//   LeafFunctor<Array, DomainFunctorTag>
//   LeafFunctor<Array, ViewFunctorTag<Domain> >
//   LeafFunctor<Array, EvalLeaf[1-7] >
//   ElementProperties<Array>
// Functions:
//   assign()
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Array
 * @brief
 * Array classes.
 */

#ifndef POOMA_ARRAY_ARRAY_H
#define POOMA_ARRAY_ARRAY_H

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag> class Array;
template<int Dim, class T, class EngineTag> class Engine;
template<class Subject, class Sub1, bool SV>
struct View1Implementation;

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "Array/PrintArray.h"
#include "Domain/Loc.h"
#include "Domain/NewDomain.h"
#include "Domain/CombineDomainOpt.h"
#include "Engine/DataObject.h"
#include "Engine/Engine.h"
#include "Engine/BrickEngine.h"
#include "Engine/CompressedFraction.h"
#include "Engine/ConstantFunctionEngine.h"
#include "Engine/ExpressionEngine.h"
#include "Engine/ForwardingEngine.h"
#include "Engine/IndirectionEngine.h"
#include "Engine/EngineFunctor.h"
#include "Engine/EnginePatch.h"
#include "Evaluator/Evaluator.h"
#include "Evaluator/WhereProxy.h"
#include "Layout/INode.h"
#include "PETE/PETE.h"
#include "Pooma/PETEExtras.h"
#include "Pooma/PETE/ExpressionTraits.h"
#include "Pooma/View.h"
#include "Utilities/Conform.h"
#include "Utilities/PerformUpdate.h"
#include "Utilities/ElementProperties.h"
#include "Utilities/ModelElement.h"
#include "Utilities/NotifyPreRead.h"
#include "Utilities/PAssert.h"

#include "Array/ArrayOperators.h"
#include "Array/PoomaArrayOperators.h"
#include "Array/ArrayOperatorSpecializations.h"
#include "Array/VectorArrayOperators.h"
#include "Array/CreateLeaf.h"
#include "Functions/Reductions.h"

#include <iosfwd>

//-----------------------------------------------------------------------------
// Prototypes for the assign function used to assign an expression to an Array.
//
// Prototypes defined here:
//   Array = Array
//   Array = scalar
//
// If you wish to have Array work with other types of objects on the right-
// hand side (for example, other classes that derive from Array), define
// extra assign() functions that take the following arguments:
//
//   assign(Array<Dim,T,EngineTag>, yourclass, Operator)
//
// where "yourclass" is the class that you would like to work on the
// right-hand side in an expression with an Array on the left-hand side.
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag,
  int OtherDim, class OtherT, class OtherEngineTag, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs,
       const Array<OtherDim, OtherT, OtherEngineTag> &rhs,
       const Op &op);

template<int Dim, class T, class EngineTag, class T1, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs, const T1 &rhs, const Op &op);
    
//-----------------------------------------------------------------------------
// View specializations for Array.
//-----------------------------------------------------------------------------

// These are for views of the form a(b).

// Single-valued version. Handles scalars and Locs.

template<class Subject, class Sub1, bool SV>
struct View1Implementation;

template<int Dim, class T, class EngineTag, class Domain>
struct View1Implementation<Array<Dim, T, EngineTag>, Domain, true>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types are pretty simple here.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Read/write versions.
  
  template<class S1, class Combine>
  inline static 
  Type_t make(const Subject_t &a, const S1 &s1,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  template<class S1, class S2, class Combine>
  inline static 
  Type_t make(const Subject_t &a,
	      const S1 &s1, const S2 &s2,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1, s2));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  template<class S1, class S2, class S3,
    class Combine>
  inline static 
  Type_t make(const Subject_t &a,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  template<class S1, class S2, class S3, class S4,
    class Combine>
  inline static 
  Type_t make(const Subject_t &a,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const S4 &s4,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class Combine>
  inline static 
  Type_t make(const Subject_t &a,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const S4 &s4, const S5 &s5,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class Combine>
  inline static 
  Type_t make(const Subject_t &a,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const S4 &s4, const S5 &s5, const S6 &s6,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5, s6));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class S7, class Combine>
  inline static 
  Type_t make(const Subject_t &a,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const S4 &s4, const S5 &s5, const S6 &s6,
	      const S7 &s7,
	      const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5, s6, s7));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine()(s);
    }

  // Read only versions.
  
  template<class S1, class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a, const S1 &s1,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }

  template<class S1, class S2, class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a,
	                  const S1 &s1, const S2 &s2,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1, s2));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }

  template<class S1, class S2, class S3,
    class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a,
	                  const S1 &s1, const S2 &s2, const S3 &s3,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }

  template<class S1, class S2, class S3, class S4,
    class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a,
	                  const S1 &s1, const S2 &s2, const S3 &s3,
	                  const S4 &s4,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a,
	                  const S1 &s1, const S2 &s2, const S3 &s3,
	                  const S4 &s4, const S5 &s5,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a,
	                  const S1 &s1, const S2 &s2, const S3 &s3,
	                  const S4 &s4, const S5 &s5, const S6 &s6,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5, s6));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class S7, class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &a,
	                  const S1 &s1, const S2 &s2, const S3 &s3,
	                  const S4 &s4, const S5 &s5, const S6 &s6,
	                  const S7 &s7,
	                  const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5, s6, s7));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
        "Array view bounds error.");
#endif
      return a.engine().read(s);
    }
};

// Non-single-valued implementation. Works for general domains
// including Nodes and INodes.

template<int Dim, class T, class EngineTag, class Domain>
struct View1Implementation<Array<Dim, T, EngineTag>, Domain, false>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce the template parameters for the output type.
  
  typedef typename Subject_t::Engine_t Engine_t;
  typedef typename NewEngine<Engine_t, Domain>::Type_t NewEngine_t;
  enum { newDim = NewEngine_t::dimensions };
  typedef typename NewEngine_t::Tag_t NewEngineTag_t;
  
  // The output type.
  
  typedef Array<newDim, T, NewEngineTag_t> Type_t;
  typedef Type_t ReadType_t;
  
  typedef NewEngineEngine<Engine_t, Domain> NewEE_t;
  typedef NewEngineDomain<Engine_t, Domain> NewED_t;

  // Functions that create the view.

  // Read-write versions.
  
  template<class S1, class Combine>
  static 
  Type_t make(const Subject_t &a, const S1 &s1,
	          const Combine &)
    {
      Domain s(Combine::make(a, s1));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
		NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  template<class S1, class S2, class Combine>
  static 
  Type_t make(const Subject_t &a, const S1 &s1,
	          const S2 &s2, const Combine &)
    {
      Domain s(Combine::make(a, s1, s2));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
	    NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  template<class S1, class S2, class S3,
    class Combine>
  static 
  Type_t make(const Subject_t &a,
	          const S1 &s1, const S2 &s2, const S3 &s3,
	          const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
		NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  template<class S1, class S2, class S3, class S4,
    class Combine>
  static 
  Type_t make(const Subject_t &a,
	          const S1 &s1, const S2 &s2, const S3 &s3,
	          const S4 &s4,
	          const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
		NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class Combine>
  static 
  Type_t make(const Subject_t &a,
	          const S1 &s1, const S2 &s2, const S3 &s3,
	          const S4 &s4, const S5 &s5,
	          const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
		NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class Combine>
  static 
  Type_t make(const Subject_t &a,
	          const S1 &s1, const S2 &s2, const S3 &s3,
	          const S4 &s4, const S5 &s5, const S6 &s6,
	          const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5, s6));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
		NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class S7, class Combine>
  static 
  Type_t make(const Subject_t &a,
	          const S1 &s1, const S2 &s2, const S3 &s3,
	          const S4 &s4, const S5 &s5, const S6 &s6,
	          const S7 &s7,
	          const Combine &)
    {
      Domain s(Combine::make(a, s1, s2, s3, s4, s5, s6, s7));
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), s),
	    "Array view bounds error.");
#endif

      return Type_t(
		NewEE_t::apply(a.engine(), s),
		NewED_t::apply(a.engine(), s));
    }

  // Read-only versions.
  
  template<class S1, class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a, const S1 &s1,
	              const Combine &c)
    {
      return make(a, s1, c);
    }

  template<class S1, class S2, class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a, const S1 &s1,
	              const S2 &s2, const Combine &c)
    {
      return make(a, s1, s2, c);
    }

  template<class S1, class S2, class S3,
    class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a,
	              const S1 &s1, const S2 &s2, const S3 &s3,
	              const Combine &c)
    {
      return make(a, s1, s2, s3, c);
    }

  template<class S1, class S2, class S3, class S4,
    class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a,
	              const S1 &s1, const S2 &s2, const S3 &s3,
	              const S4 &s4, const Combine &c)
    {
      return make(a, s1, s2, s3, s4, c);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a,
	              const S1 &s1, const S2 &s2, const S3 &s3,
	              const S4 &s4, const S5 &s5, const Combine &c)
    {
      return make(a, s1, s2, s3, s4, s5, c);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6,
    class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a,
	              const S1 &s1, const S2 &s2, const S3 &s3,
	              const S4 &s4, const S5 &s5, const S6 &s6,
	              const Combine &c)
    {
      return make(a, s1, s2, s3, s4, s5, s6, c);
    }

  template<class S1, class S2, class S3, class S4, class S5,
    class S6, class S7,
    class Combine>
  inline static 
  Type_t makeRead(const Subject_t &a,
	              const S1 &s1, const S2 &s2, const S3 &s3,
	              const S4 &s4, const S5 &s5, const S6 &s6,
	              const S7 &s7, const Combine &c)
    {
      return make(a, s1, s2, s3, s4, s5, s6, s7, c);
    }
};

// General version.

template<int Dim, class T, class EngineTag, class Domain>
struct View1<Array<Dim, T, EngineTag>, Domain>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  // At some point, we need to fix NewDomain1; until then, use
  // the temporary version from Array.h.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef TemporaryNewDomain1<Domain_t, Domain> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;
  
  // The optimized domain combiner:
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Domain &s1)
    {
      return Dispatch_t::make(a, s1, Combine_t());
    }

  inline static
  ReadType_t makeRead(const Subject_t &a, const Domain &s1)
    {
      return Dispatch_t::makeRead(a, s1, Combine_t());
    }
};

// View0 deals with the special case of read() and
// operator().
// We don't use TemporaryNewDomain1 because the domain of
// an existing engine cannot be any kind of slice domain.
// Also, bounds checking would make no sense because it would
// reduce to contains(a.domain(), a.domain());

template<int Dim, class T, class EngineTag>
struct View0<Array<Dim, T, EngineTag> >
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  // At some point, we need to fix NewDomain1; until then, use
  // the temporary version from Array.h.
  
  typedef typename Subject_t::Engine_t Engine_t;
  typedef typename Subject_t::Domain_t Domain_t;

  // Deduce the template parameters for the output type.
  
  typedef typename NewEngine<Engine_t, Domain_t>::Type_t NewEngine_t;
  enum { newDim = NewEngine_t::dimensions };
  typedef typename NewEngine_t::Tag_t NewEngineTag_t;
  
  // The output types.
  
  typedef Array<newDim, T, NewEngineTag_t> Type_t;
  typedef Type_t ReadType_t;

  // Functions that create the view.
  
  static Type_t make(const Subject_t &a)
    {
      typedef NewEngineEngine<Engine_t, Domain_t> NewEE_t;
      typedef NewEngineDomain<Engine_t, Domain_t> NewED_t;

      return Type_t(
		NewEE_t::apply(a.engine(), a.engine().domain()),
		NewED_t::apply(a.engine(), a.engine().domain()));
    }
    
  inline static ReadType_t makeRead(const Subject_t &a)
    {
      return make(a);
    }
};

template<int Dim, class T, class EngineTag>
struct View1<Array<Dim, T, EngineTag>, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  Type_t make(const Subject_t &a, int s1)
    { 
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<1>(s1)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1);
    }

  inline static
  ReadType_t makeRead(const Subject_t &a, int s1)
    { 
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<1>(s1)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1);
    }
};

template<int D1, class T1, class E1, int D2, class T2, class E2>
struct View1<Array<D1, T1, E1>, Array<D2, T2, E2> >
{
  typedef Array<D1, T1, E1> Array1_t;
  typedef Array<D2, T2, E2> Array2_t;

  typedef IndirectionTag<Array1_t, Array2_t> Tag_t;

  // The return types.
  
  typedef Array<D2, T1, Tag_t> Type_t;
  typedef Type_t ReadType_t;

  // Functions that make the view.
  
  static 
  Type_t make(const Array1_t &a, const Array2_t &s)
    {
      return Type_t(a, s);
    }

  inline static 
  Type_t makeRead(const Array1_t &a, const Array2_t &s)
    {
      return make(a, s);
    }
};

// These are for views of the form a(b, c).

template<int Dim, class T, class EngineTag, 
  class Sub1, class Sub2>
struct View2<Array<Dim, T, EngineTag>, Sub1, Sub2>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain2<Sub1, Sub2> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;

  // Create an optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1, const Sub2 &s2)
    {
      return Dispatch_t::make(a, s1, s2, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, const Sub1 &s1, const Sub2 &s2)
    {
      return Dispatch_t::makeRead(a, s1, s2, Combine_t());
    }
};

template<int Dim, class T, class EngineTag>
struct View2<Array<Dim, T, EngineTag>, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  Type_t make(const Subject_t &a, int s1, int s2)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<2>(s1, s2)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1, s2);
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, int s1, int s2)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<2>(s1, s2)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1, s2);
    }
};

// These are for views of the form a(b, c, d).

template<int Dim, class T, class EngineTag, 
  class Sub1, class Sub2, class Sub3>
struct View3<Array<Dim, T, EngineTag>, Sub1, Sub2, Sub3>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain3<Sub1, Sub2, Sub3> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;
  
  // Create an optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3)
    {
      return Dispatch_t::make(a, s1, s2, s3, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3)
    {
      return Dispatch_t::makeRead(a, s1, s2, s3, Combine_t());
    }
};

template<int Dim, class T, class EngineTag>
struct View3<Array<Dim, T, EngineTag>, int, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  Type_t make(const Subject_t &a, int s1, int s2, int s3)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<3>(s1, s2, s3)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1, s2, s3);
    }

  inline static
  ReadType_t makeRead(const Subject_t &a, int s1, int s2, int s3)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<3>(s1, s2, s3)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1, s2, s3);
    }
};

// These are for views of the form a(b, c, d, e).

template<int Dim, class T, class EngineTag, 
  class Sub1, class Sub2, class Sub3, class Sub4>
struct View4<Array<Dim, T, EngineTag>, 
  Sub1, Sub2, Sub3, Sub4>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain4<Sub1, Sub2, Sub3, Sub4> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;
  
  // Create an optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4)
    {
      return Dispatch_t::make(a, s1, s2, s3, s4, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4)
    {
      return Dispatch_t::makeRead(a, s1, s2, s3, s4, Combine_t());
    }
};

template<int Dim, class T, class EngineTag>
struct View4<Array<Dim, T, EngineTag>, int, int, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  Type_t make(const Subject_t &a, int s1, int s2, int s3, int s4)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<4>(s1, s2, s3, s4)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1, s2, s3, s4);
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, int s1, int s2, int s3, int s4)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<4>(s1, s2, s3, s4)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1, s2, s3, s4);
    }
};

// These are for views of the form a(b, c, d, e, f).

template<int Dim, class T, class EngineTag, 
  class Sub1, class Sub2, class Sub3, class Sub4, class Sub5>
struct View5<Array<Dim, T, EngineTag>, 
  Sub1, Sub2, Sub3, Sub4, Sub5>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain5<Sub1, Sub2, Sub3, Sub4, Sub5> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;
  
  // Create an optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4, const Sub5 &s5)
    {
      return Dispatch_t::make(a, s1, s2, s3, s4, s5, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4, const Sub5 &s5)
    {
      return Dispatch_t::makeRead(a, s1, s2, s3, s4, s5, Combine_t());
    }
};

template<int Dim, class T, class EngineTag>
struct View5<Array<Dim, T, EngineTag>, int, int, int, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  ReadType_t makeRead(const Subject_t &a, int s1, int s2, int s3, 
    int s4, int s5)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<5>(s1, s2, s3, s4, s5)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1, s2, s3, s4, s5);
    }
  
  inline static
  Type_t make(const Subject_t &a, int s1, int s2, int s3, int s4, int s5)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<5>(s1, s2, s3, s4, s5)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1, s2, s3, s4, s5);
    }
};

// These are for views of the form a(b, c, d, e, f, g).

template<int Dim, class T, class EngineTag, 
  class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
  class Sub6>
struct View6<Array<Dim, T, EngineTag>, 
  Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain6<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;
  
  // Create an optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4, const Sub5 &s5, const Sub6 &s6)
    {
      return Dispatch_t::make(a, s1, s2, s3, s4, s5, s6, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4, const Sub5 &s5, const Sub6 &s6)
    {
      return Dispatch_t::makeRead(a, s1, s2, s3, s4, s5, s6, Combine_t());
    }
};

template<int Dim, class T, class EngineTag>
struct View6<Array<Dim, T, EngineTag>, int, int, int, int, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  Type_t make(const Subject_t &a, int s1, int s2, int s3, int s4, int s5,
	      int s6)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<6>(s1, s2, s3, s4, s5, s6)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1, s2, s3, s4, s5, s6);
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, int s1, int s2, int s3, 
    int s4, int s5, int s6)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<6>(s1, s2, s3, s4, s5, s6)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1, s2, s3, s4, s5, s6);
    }
};

// These are for views of the form a(b, c, d, e, f, g, h).

template<int Dim, class T, class EngineTag, 
  class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
  class Sub6, class Sub7>
struct View7<Array<Dim, T, EngineTag>, 
  Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain7<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7> 
    NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;
  
  // Create an optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;

  // The return types.
  
  typedef typename Dispatch_t::Type_t Type_t;
  typedef typename Dispatch_t::ReadType_t ReadType_t;

  // Functions that create the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4, const Sub5 &s5, const Sub6 &s6, 
    const Sub7 &s7)
    {
      return Dispatch_t::make(a, s1, s2, s3, s4, s5, s6, s7, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3, const Sub4 &s4, const Sub5 &s5, const Sub6 &s6, 
    const Sub7 &s7)
    {
      return Dispatch_t::makeRead(a, s1, s2, s3, s4, s5, s6, s7, 
        Combine_t());
    }
};

template<int Dim, class T, class EngineTag>
struct View7<Array<Dim, T, EngineTag>, int, int, int, int, int, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // Functions that perform the indexing.
  
  inline static
  Type_t make(const Subject_t &a, int s1, int s2, int s3, int s4, int s5,
	      int s6, int s7)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<7>(s1, s2, s3, s4, s5, s6, s7)),
	    "Array view bounds error.");
#endif
      return a.engine()(s1, s2, s3, s4, s5, s6, s7);
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &a, int s1, int s2, int s3, 
    int s4, int s5, int s6, int s7)
    {
#if POOMA_BOUNDS_CHECK
      PInsist(contains(a.domain(), Loc<7>(s1, s2, s3, s4, s5, s6, s7)),
	    "Array view bounds error.");
#endif
      return a.engine().read(s1, s2, s3, s4, s5, s6, s7);
    }
};

/**
 * Patch specialization for Array.
 */

template<int Dim, class T, class EngineTag>
struct Patch<Array<Dim, T, EngineTag> >
{
  typedef Array<Dim, T, EngineTag> Subject_t;
  typedef typename Subject_t::Engine_t OldEngine_t;
  typedef typename EngineFunctor<OldEngine_t, EnginePatch>::Type_t Engine_t;

  // We've assumed that Dim and T are the same for the patch engine.

  typedef Array<Dim, T, typename Engine_t::Tag_t> Type_t;

  inline static
  Type_t make(const Subject_t &subject, int i)
    {
      return Type_t(engineFunctor(subject.engine(), EnginePatch(i)));
    }
};


/**
 * ComponentView specialization for Array.
 */

template<class Components, int Dim, class T, class EngineTag>
struct ComponentView<Components, Array<Dim, T, EngineTag> >
{
  // Convenience typedef for the thing we're taking a component view of.
  
  typedef Array<Dim, T, EngineTag> Subject_t;

  // Deduce the template parameters for the output type.
  
  typedef Engine<Dim, T, EngineTag> Engine_t;
  typedef typename Engine_t::Element_t Element_t;
  typedef typename ComponentAccess<Element_t, Components>::Element_t NewT_t;
  typedef CompFwd<Engine_t, Components> NewEngineTag_t;
  
  // The output type.
  
  typedef Array<Dim, NewT_t, NewEngineTag_t> Type_t;

  // A function that creates the view.
  
  inline static
  Type_t make(const Subject_t &a, const Components &c)
    {
      return Type_t(a, ComponentWrapper<Components>(c));
    }
};


//-----------------------------------------------------------------------------
// Array
//-----------------------------------------------------------------------------

/** Arrays are used to apply operations to N-dimensional (N <= 7) logically 
 * rectangular, logically dense sets of elements. Array provides several 
 * services:
 *   - General subsetting operations. Using operator(), the user can ask the
 *     array to construct new arrays, which give views onto part 
 *     of an Array's domain.
 *   - Representation independence. Arrays work with arbitrary engines, 
 *     which fill requests for views with actual data.
 */

template<int Dim, class T = POOMA_DEFAULT_ELEMENT_TYPE,
  class EngineTag = POOMA_DEFAULT_ENGINE_TYPE>
class Array
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  //---------------------------------------------------------------------------
  /// This_t is a convenience typedef referring to this class.
  typedef Array<Dim, T, EngineTag> This_t;

  //---------------------------------------------------------------------------
  /// Engine_t is a typedef that maps the EngineTag to the Engine type.
  typedef Engine<Dim, T, EngineTag> Engine_t;
  /// Engine Tag_t is this tag type, exported to other classes.
  typedef EngineTag EngineTag_t;

  //---------------------------------------------------------------------------
  // This might be a reference or a proxy object. 

  /// Element_t is the type of elements managed by this array's engine.
  typedef typename Engine_t::Element_t Element_t;
  /// ElementRef_t is the type that is used to write to a single element. 
  typedef typename Engine_t::ElementRef_t ElementRef_t;
  /// The typedef Domain_t gives the type of domain this array is defined on.
  typedef typename Engine_t::Domain_t Domain_t;
  typedef typename Engine_t::Layout_t Layout_t;
  /// The enum 'rank' gives the rank of this array. 'dimensions' is the same
  /// thing.
  enum { dimensions = Engine_t::dimensions };
  enum { rank = Engine_t::dimensions };

  // Arrays dont support relations attached to them.

  enum { hasRelations = false };
  
  //===========================================================================
  // Constructors
  //===========================================================================
   
  //---------------------------------------------------------------------------
  /// Default constructor for Array. Exists so array can be resized.

  Array() { }
  
  //---------------------------------------------------------------------------
  /// @name Engine Array constructors.
  //@{

  inline explicit Array(const Engine_t &modelEngine)
  : engine_m(modelEngine)
    { }

  template<int Dim2, class T2, class EngineTag2, class Initializer>
  inline Array(const Engine<Dim2, T2, EngineTag2> &engine, 
    const Initializer &init)
  : engine_m(engine, init)
    { }
  //@}

  //---------------------------------------------------------------------------
  /// Indirection Array constructors. 

  template<int D1, class T1, class E1, int D2, class T2, class E2>
  inline Array(const Array<D1, T1, E1> &a1, const Array<D2, T2, E2> &a2)
  : engine_m(a1, a2)
    { }

  //---------------------------------------------------------------------------
  /// One way for a user to construct an Array is to use another Array as a 
  /// model. The non-templated version is the copy constructor.
  /// For the templated versions to work, Engine_t must possess a constructor
  /// that takes an OtherEngine and, perhaps, an OtherDomain.

  inline Array(const This_t &model)
  : engine_m(model.engine())
    { }
  
  /// This ctor is called, for example, to initialize a brick-view with a
  /// compressible brick.
  
  template<int OtherDim, class OtherT, class OtherEngineTag>
  inline explicit Array(const Array<OtherDim, OtherT, OtherEngineTag> &model)
  : engine_m(model.engine())
    { }
  
  template<int OtherDim, class OtherT, class OtherEngineTag, class OtherDomain>
  inline Array(const Array<OtherDim, OtherT, OtherEngineTag> &model, 
        const OtherDomain &domain)
  : engine_m(NewEngineEngine<Engine<OtherDim,OtherT,OtherEngineTag>, 
             OtherDomain>::apply(model.engine(),domain),
             NewEngineDomain<Engine<OtherDim,OtherT,OtherEngineTag>,
             OtherDomain>::apply(model.engine(),domain))
    { }
 
  template <class OtherT, class OtherEngineTag, class Components>
  Array(const Array<Dim, OtherT, OtherEngineTag> &a,
	const ComponentWrapper<Components>& c)
  : engine_m(a.engine(), c.components())
    { }

  //---------------------------------------------------------------------------
  /// @name constructors passing domain information to the engine
  /// All of these constructors pass domain information to the engine. 
  /// These domains are formed by combining up to 7 subdomains, which
  /// must yield an Interval<N>, where N must be the same as the
  /// rank of the array. These constructors call the default constructor
  /// for Element_t. 
  //@{
  template<class Sub1>
  explicit Array(const Sub1 &s1)
  : engine_m(NewDomain1<Sub1>::combine(s1))
    { }

  template<class Sub1, class Sub2>
  Array(const Sub1 &s1, const Sub2 &s2)
  : engine_m(NewDomain2<Sub1, Sub2>::combine(s1, s2))
    { }

  template<class Sub1, class Sub2, class Sub3>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3)
  : engine_m(NewDomain3<Sub1, Sub2, Sub3>::combine(s1, s2, s3))
    { }
    
  template<class Sub1, class Sub2, class Sub3, class Sub4>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4)
  : engine_m(NewDomain4<Sub1, Sub2, Sub3, Sub4>::
      combine(s1, s2, s3, s4))
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
        const Sub5 &s5)
  : engine_m(NewDomain5<Sub1, Sub2, Sub3, Sub4, Sub5>::
      combine(s1, s2, s3, s4, s5))
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
        const Sub5 &s5, const Sub6 &s6)
  : engine_m(NewDomain6<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>::
      combine(s1, s2, s3, s4, s5, s6))
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6, class Sub7>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
        const Sub5 &s5, const Sub6 &s6, const Sub7 &s7)
  : engine_m(NewDomain7<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>::
      combine(s1, s2, s3, s4, s5, s6, s7))
    { }
  //@}

  //---------------------------------------------------------------------------
  /// @name constructors passing domain information to the engine along with a model element
  /// All of these constructors pass domain information to the engine along
  /// with a model element, which is used to initialize all array elements.
  /// Therefore, it is assumed that Element_t either has deep-copy semantics
  /// by default, or can be made to have it using a wrapper class.
  /// These domains are formed by combining up to 7 subdomains, which
  /// must yield an Interval<N>, where N must be the same as the
  /// rank of the array. 
  //@{
  template<class Sub1>
  Array(const Sub1 &s1, const ModelElement<Element_t> &model)
  : engine_m(NewDomain1<Sub1>::combine(s1), model.element())
    { }

  template<class Sub1, class Sub2>
  Array(const Sub1 &s1, const Sub2 &s2, 
        const ModelElement<Element_t> &model)
  : engine_m(NewDomain2<Sub1, Sub2>::combine(s1, s2), model.element())
    { }

  template<class Sub1, class Sub2, class Sub3>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
        const ModelElement<Element_t> &model)
  : engine_m(NewDomain3<Sub1, Sub2, Sub3>::combine(s1, s2, s3), 
      model.element())
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4, 
        const ModelElement<Element_t> &model)
  : engine_m(NewDomain4<Sub1, Sub2, Sub3, Sub4>::
      combine(s1, s2, s3, s4), model.element())
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
        const Sub5 &s5, const ModelElement<Element_t> &model)
  : engine_m(NewDomain5<Sub1, Sub2, Sub3, Sub4, Sub5>::
      combine(s1, s2, s3, s4, s5), model.element())
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
        const Sub5 &s5, const Sub6 &s6, const ModelElement<Element_t> &model)
  : engine_m(NewDomain6<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>::
      combine(s1, s2, s3, s4, s5, s6), model.element())
    { }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6, class Sub7>
  Array(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
        const Sub5 &s5, const Sub6 &s6, const Sub7 &s7, 
        const ModelElement<Element_t> &model)
  : engine_m(NewDomain7<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>::
      combine(s1, s2, s3, s4, s5, s6, s7), model.element())
    { }
  //@}

  //===========================================================================
  /// @name Array initializers
  //===========================================================================
  //@{
  void initialize(const Engine_t &modelEngine)
    {
      engine_m = modelEngine;
    }
  
  template<int Dim2, class T2, class EngineTag2, class Initializer>
  void initialize(const Engine<Dim2, T2, EngineTag2> &engine, 
    const Initializer &init)
    {
      engine_m = Engine_t(engine, init);
    }
  //@}

  //---------------------------------------------------------------------------
  /// @name Initializers from array
  /// One way for a user to initialize an Array is to use another Array as a 
  /// model. The non-templated version is the copy constructor.
  /// For the templated versions to work, Engine_t must possess a constructor
  /// that takes an OtherEngine and, perhaps, an OtherDomain.
  //@{
  void initialize(const This_t &model)
    {
      engine_m = model.engine();
    }
  
  template<int OtherDim, class OtherT, class OtherEngineTag>
  void initialize(const Array<OtherDim, OtherT, OtherEngineTag> &model)
    {
      engine_m = Engine_t(model.engine());
    }
  
  template<int OtherDim, class OtherT, class OtherEngineTag, class OtherDomain>
  void initialize(const Array<OtherDim, OtherT, OtherEngineTag> &model, 
    const OtherDomain &domain)
    {
      engine_m = Engine_t(
        NewEngineEngine<Engine<OtherDim,OtherT,OtherEngineTag>, OtherDomain>::
        apply(model.engine(),domain),
        NewEngineDomain<Engine<OtherDim,OtherT,OtherEngineTag>, OtherDomain>::
        apply(model.engine(),domain));
    }
  //@}
  //---------------------------------------------------------------------------
  /// @name Initializers from domain
  /// All of these initializers pass domain information to the engine. 
  /// These domains are formed by combining up to 7 subdomains, which
  /// must yield an Interval<N>, where N must be the same as the
  /// rank of the array. These initializers call the default constructor
  /// for Element_t. 
  //@{
  template<class Sub1>
  void initialize(const Sub1 &s1)
    {
      engine_m = Engine_t(NewDomain1<Sub1>::combine(s1));  
    }

  template<class Sub1, class Sub2>
  void initialize(const Sub1 &s1, const Sub2 &s2)
    {
      engine_m = Engine_t(NewDomain2<Sub1, Sub2>::combine(s1, s2));
    }

  template<class Sub1, class Sub2, class Sub3>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3)
    { 
      engine_m = Engine_t(NewDomain3<Sub1, Sub2, Sub3>::combine(s1, s2, s3));
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4)
    {
      engine_m = Engine_t(NewDomain4<Sub1, Sub2, Sub3, Sub4>::
                          combine(s1, s2, s3, s4));
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
                  const Sub5 &s5)
    {
      engine_m = Engine_t(NewDomain5<Sub1, Sub2, Sub3, Sub4, Sub5>::
                          combine(s1, s2, s3, s4, s5));
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
                  const Sub5 &s5, const Sub6 &s6)
    {
      engine_m = Engine_t(NewDomain6<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>::
                          combine(s1, s2, s3, s4, s5, s6));
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6, class Sub7>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5, const Sub6 &s6, const Sub7 &s7)
    {
      engine_m = 
        Engine_t(NewDomain7<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>::
        combine(s1, s2, s3, s4, s5, s6, s7));
    }
  //@}
  //---------------------------------------------------------------------------
  /// @name Initializers from domain and model
  /// All of these constructors pass domain information to the engine along
  /// with a model element, which is used to initialize all array elements.
  /// Therefore, it is assumed that Element_t either has deep-copy semantics
  /// by default, or can be made to have it using a wrapper class.
  /// These domains are formed by combining up to 7 subdomains, which
  /// must yield an Interval<N>, where N must be the same as the
  /// rank of the array. 
  //@{
  template<class Sub1>
  void initialize(const Sub1 &s1, const ModelElement<Element_t> &model)
    {
      engine_m = Engine_t(NewDomain1<Sub1>::combine(s1), model.element());  
    }

  template<class Sub1, class Sub2>
  void initialize(const Sub1 &s1, const Sub2 &s2, 
    const ModelElement<Element_t> &model)
    {
      engine_m = Engine_t(NewDomain2<Sub1, Sub2>::combine(s1, s2), 
        model.element());
    }

  template<class Sub1, class Sub2, class Sub3>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const ModelElement<Element_t> &model)
    {
      engine_m = Engine_t(NewDomain3<Sub1, Sub2, Sub3>::combine(s1, s2, s3), 
        model.element());
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4, 
    const ModelElement<Element_t> &model)
    {
      engine_m = Engine_t(NewDomain4<Sub1, Sub2, Sub3, Sub4>::
        combine(s1, s2, s3, s4), model.element());
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
    const Sub5 &s5, const ModelElement<Element_t> &model)
    {
      engine_m = Engine_t(NewDomain5<Sub1, Sub2, Sub3, Sub4, Sub5>::
        combine(s1, s2, s3, s4, s5), model.element());
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, const Sub4 &s4,
    const Sub5 &s5, const Sub6 &s6, const ModelElement<Element_t> &model)
    {
      engine_m = Engine_t(NewDomain6<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>::
        combine(s1, s2, s3, s4, s5, s6), model.element());
    }

  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6, class Sub7>
  void initialize(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5, const Sub6 &s6, const Sub7 &s7, 
    const ModelElement<Element_t> &model)
    {
      engine_m = 
        Engine_t(NewDomain7<Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>::
        combine(s1, s2, s3, s4, s5, s6, s7), model.element());
    }
  //@}
  //===========================================================================
  // Destructor
  //===========================================================================
  
  //---------------------------------------------------------------------------
  /// The destructor is trivial since the engine is implicitly destroyed.

  ~Array() 
    { }
  
  //===========================================================================
  /// @name Accessors and mutators
  //===========================================================================
  //@{
  //---------------------------------------------------------------------------
  /// Patch accessor function returns the i'th patch.
 
  inline typename Patch<This_t>::Type_t
  patchLocal(int i) const
    {
      return Patch<This_t>::make(*this, i);
    }

  inline int
  numPatchesLocal() const
    {
      return engineFunctor(engine_m, EngineNumPatches());
    }

  //---------------------------------------------------------------------------
  /// @name Domain accessors
  /// Accessor functions that return this array's domain, which is obtained 
  /// from the engine. Note that arrays dont distinguish between vertex and
  /// cell domains and such the array's domains are vertex domains in the
  /// field sense.
  //@{

  /// Shortcut for totalDomain().

  inline const Domain_t& domain() const 
    {
      return engine_m.domain();
    }

  /// Returns the physical domain, i.e. the domain without external guards.

  inline Domain_t physicalDomain() const 
    {
      return engine_m.layout().innerDomain();
    }

  /// Returns the total domain, i.e. the domain with external guards.

  inline const Domain_t& totalDomain() const 
    {
      return engine_m.domain();
    }
  //@}
  //---------------------------------------------------------------------------
  /// Returns the Array's layout.
  
  inline typename Engine_t::Layout_t layout() const
    {
      return engine_m.layout();
    }
    
  //---------------------------------------------------------------------------
  // @name View-creation operations
  /// View-creation operations. These operator() functions take one or more
  /// sub-domains, which combine to form a domain with dimensionality identical
  /// to the rank of the array. Views based on up to 7 subdomains are supported.
  /// Indirection views will end up calling these functions.
  ///
  /// A zero-argument version of operator(), which takes a view of 
  /// array's domain, is also supplied.
  //@{
  typename View0<This_t>::ReadType_t 
  read() const
    {
      typedef View0<This_t> Ret_t;
      return Ret_t::makeRead(*this);
    }

  template<class Sub1> 
  inline typename View1<This_t, Sub1>::ReadType_t 
  read(const Sub1 &s1) const
    {
      typedef View1<This_t, Sub1> Ret_t;
      return Ret_t::makeRead(*this, s1);
    }
  
  template<class Sub1, class Sub2> 
  inline typename View2<This_t, Sub1, Sub2>::ReadType_t
  read(const Sub1 &s1, const Sub2 &s2) const
    {
      typedef View2<This_t, Sub1, Sub2> Ret_t;
      return Ret_t::makeRead(*this, s1, s2);
    }
  
  template<class Sub1, class Sub2, class Sub3> 
  inline typename View3<This_t, Sub1, Sub2, Sub3>::ReadType_t
  read(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3) const
    {
      typedef View3<This_t, Sub1, Sub2, Sub3> Ret_t;
      return Ret_t::makeRead(*this, s1, s2, s3);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4> 
  inline typename View4<This_t, Sub1, Sub2, Sub3, Sub4>::ReadType_t
  read(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4) const
    {
      typedef View4<This_t, Sub1, Sub2, Sub3, Sub4> Ret_t;
      return Ret_t::makeRead(*this, s1, s2, s3, s4);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5> 
  inline typename View5<This_t, Sub1, Sub2, Sub3, Sub4, Sub5>::ReadType_t
  read(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5) const
    {
      typedef View5<This_t, Sub1, Sub2, Sub3, Sub4, Sub5> Ret_t;
      return Ret_t::makeRead(*this, s1, s2, s3, s4, s5);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6> 
  inline typename View6<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>::ReadType_t
  read(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5, const Sub6 &s6) const
    {
      typedef View6<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6> Ret_t;
      return Ret_t::makeRead(*this, s1, s2, s3, s4, s5, s6);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6, class Sub7> 
  inline typename 
    View7<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>::ReadType_t
  read(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5, const Sub6 &s6, const Sub7 &s7) const
    {
      typedef View7<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7> Ret_t;
      return Ret_t::makeRead(*this, s1, s2, s3, s4, s5, s6, s7);
    }

  typename View0<This_t>::Type_t 
  operator()() const
    {
      typedef View0<This_t> Ret_t;
      return Ret_t::make(*this);
    }

  template<class Sub1> 
  inline typename View1<This_t,Sub1>::Type_t 
  operator()(const Sub1 &s1) const
    {
      typedef View1<This_t, Sub1> Ret_t;
      return Ret_t::make(*this, s1);
    }
  
  template<class Sub1, class Sub2> 
  inline typename View2<This_t, Sub1, Sub2>::Type_t
  operator()(const Sub1 &s1, const Sub2 &s2) const
    {
      typedef View2<This_t, Sub1, Sub2> Ret_t;
      return Ret_t::make(*this, s1, s2);
    }
  
  template<class Sub1, class Sub2, class Sub3> 
  inline typename View3<This_t, Sub1, Sub2, Sub3>::Type_t
  operator()(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3) const
    {
      typedef View3<This_t, Sub1, Sub2, Sub3> Ret_t;
      return Ret_t::make(*this, s1, s2, s3);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4> 
  inline typename View4<This_t, Sub1, Sub2, Sub3, Sub4>::Type_t
  operator()(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4) const
    {
      typedef View4<This_t, Sub1, Sub2, Sub3, Sub4> Ret_t;
      return Ret_t::make(*this, s1, s2, s3, s4);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5> 
  inline typename View5<This_t, Sub1, Sub2, Sub3, Sub4, Sub5>::Type_t
  operator()(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5) const
    {
      typedef View5<This_t, Sub1, Sub2, Sub3, Sub4, Sub5> Ret_t;
      return Ret_t::make(*this, s1, s2, s3, s4, s5);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6> 
  inline typename 
    View6<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6>::Type_t
  operator()(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5, const Sub6 &s6) const
    {
      typedef View6<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6> Ret_t;
      return Ret_t::make(*this, s1, s2, s3, s4, s5, s6);
    }
  
  template<class Sub1, class Sub2, class Sub3, class Sub4, class Sub5,
    class Sub6, class Sub7> 
  inline typename 
    View7<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7>::Type_t
  operator()(const Sub1 &s1, const Sub2 &s2, const Sub3 &s3, 
    const Sub4 &s4, const Sub5 &s5, const Sub6 &s6, const Sub7 &s7) const
    {
      typedef View7<This_t, Sub1, Sub2, Sub3, Sub4, Sub5, Sub6, Sub7> Ret_t;
      return Ret_t::make(*this, s1, s2, s3, s4, s5, s6, s7);
    }
  //@}

  //---------------------------------------------------------------------------
  /// @name Component forwarding functions
  /// Component forwarding function. This is used when Element_t has components,
  /// as is the case for vectors, tensors, or arrays. A Forwarding engine
  /// is made containing a shallow copy of this array's engine and a Loc<M>,
  /// where M is the number of components that are being set. For this to
  /// ultimately work, Element_t must have a working const operator()(Loc<M>)
  /// that returns a reference or proxy to the component.
  //@{
  inline typename ComponentView<Loc<1>, This_t>::Type_t
  comp(int i1) const
    {
      return ComponentView<Loc<1>, This_t>::make(*this, Loc<1>(i1));
    }

  inline typename ComponentView<Loc<2>, This_t>::Type_t
  comp(int i1, int i2) const
    {
      return ComponentView<Loc<2>, This_t>::make(*this, Loc<2>(i1, i2));
    }

  inline typename ComponentView<Loc<3>, This_t>::Type_t
  comp(int i1, int i2, int i3) const
    {
      return ComponentView<Loc<3>, This_t>::make(*this, Loc<3>(i1, i2, i3));
    }

  inline typename ComponentView<Loc<4>, This_t>::Type_t
  comp(int i1, int i2, int i3, int i4) const
    {
      return ComponentView<Loc<4>, This_t>::make(*this, 
        Loc<4>(i1, i2, i3, i4));
    }

  inline typename ComponentView<Loc<5>, This_t>::Type_t
  comp(int i1, int i2, int i3, int i4, int i5) const
    {
      return ComponentView<Loc<5>, This_t>::make(*this, 
        Loc<5>(i1, i2, i3, i4, i5));
    }

  inline typename ComponentView<Loc<6>, This_t>::Type_t
  comp(int i1, int i2, int i3, int i4, int i5, int i6) const
    {
      return ComponentView<Loc<6>, This_t>::make(*this, 
        Loc<6>(i1, i2, i3, i4, i5, i6));
    }

  inline typename ComponentView<Loc<7>, This_t>::Type_t
  comp(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
    {
      return ComponentView<Loc<7>, This_t>::make(*this, 
        Loc<7>(i1, i2, i3, i4, i5, i6, i7));
    }

  template<class Components>
  inline typename ComponentView<Components, This_t>::Type_t
  comp(const Components &components) const
    {
      return ComponentView<Components, This_t>::make(*this, components);
    }
  //@}
  
  //---------------------------------------------------------------------------
  /// Instruct the array to make its own copy of its data.

  inline void makeOwnCopy()
    { engine_m.makeOwnCopy(); }

  //---------------------------------------------------------------------------
  /// @name Short-circuit functions that return information from domain.
  //@{
  inline int first(int d) const
  {
#if POOMA_BOUNDS_CHECK
      PInsist2(d >= 0 && d < Dim,
               "Array<%d,...>::first() bounds error, index = %d.", Dim, d);
#endif
      return engine_m.first(d);
    }
  
  inline int last(int d) const
    {
#if POOMA_BOUNDS_CHECK
      PInsist2(d >= 0 && d < Dim,
               "Array<%d,...>::last() bounds error, index = %d.", Dim, d);
#endif
      return engine_m.domain()[d].last();
    }
  
  inline int length(int d) const
    {
#if POOMA_BOUNDS_CHECK
      PInsist2(d >= 0 && d < Dim,
               "Array<%d,...>::length() bounds error, index = %d.", Dim, d);
#endif
      return engine_m.domain()[d].length();
    }
  
  inline Loc<Dim> firsts() const
    {
      return engine_m.domain().firsts();
    }
  
  inline Loc<Dim> lasts() const
    {
      return engine_m.domain().lasts();
    }
  
  inline Loc<Dim> lengths() const
    {
      return engine_m.domain().lengths();
    }

  inline long size() const
    {
      return engine_m.domain().size();
    }
  //@}
  
  //---------------------------------------------------------------------------
  /// @name Copy assignment operators
  /// Copy assignment operators. We pack this assignment expression into a
  /// PETE binary expression tree node and then use this to construct an
  /// array with an expression engine. We then pass this on to an evaluator,
  /// which handles the computation. The first two versions handle assigning
  /// Arrays to Arrays and the fourth one handles assigning
  /// scalars.
  //@{
  This_t &operator=(const Array<Dim, T, EngineTag> &rhs)
    {
      assign(*this, rhs, OpAssign());
      return *this;
    }

  const This_t &operator=(const Array<Dim, T, EngineTag> &rhs) const
    {
      return assign(*this, rhs, OpAssign());
    }

  template<class T1>
  const This_t &operator=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpAssign());
    }
  //@}

  //---------------------------------------------------------------------------
  /// @name Op-assignment operators
  //@{

  /// Addition.

  template<class T1>
  const This_t &operator+=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpAddAssign());
    }

  /// Subtraction.

  template<class T1>
  const This_t &operator-=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpSubtractAssign());
    }

  /// Multiplication.

  template<class T1>
  const This_t &operator*=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpMultiplyAssign());
    }

  /// Division.

  template<class T1>
  const This_t &operator/=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpDivideAssign());
    }

  /// Modulus.

  template<class T1>
  const This_t &operator%=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpModAssign());
    }

  /// Bitwise-Or.

  template<class T1>
  const This_t &operator|=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpBitwiseOrAssign());
    }

  /// Bitwise-And.

  template<class T1>
  const This_t &operator&=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpBitwiseAndAssign());
    }

  /// Bitwise-Xor.

  template<class T1>
  const This_t &operator^=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpBitwiseXorAssign());
    }

  /// Left shift.

  template<class T1>
  const This_t &operator<<=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpLeftShiftAssign());
    }

  /// Right shift.

  template<class T1>
  const This_t &operator>>=(const T1 &rhs) const
    {
      return assign(*this, rhs, OpRightShiftAssign());
    }
  //@}
  
  //---------------------------------------------------------------------------
  /// @name Assessor functions that return the engine
  //@{
  inline Engine_t &engine() 
    { return engine_m; }
  inline const Engine_t &engine() const 
    { return engine_m; }
  //@}

private:    

  //---------------------------------------------------------------------------
  /// This array's engine. Vrroom.

  Engine_t engine_m;
};  

/**
 * This specialization of LeafFunctor is used to get the domain type or the
 * domain itself from an Array. Used only by Expression-Engine.
 */

template<int Dim, class T, class EngineTag>
struct LeafFunctor<Array<Dim, T, EngineTag>, DomainFunctorTag>
{
  typedef typename Engine<Dim, T, EngineTag>::Domain_t Type_t;
  static Type_t apply(const Array<Dim, T, EngineTag> &a, 
		      const DomainFunctorTag &)
  {
    return a.domain();
  }
};

/**
 * This specialization of LeafFunctor is used to apply a view (subsetting) 
 * operation to an Array. The domain will always be zero-based since this
 * is used only by Expression-Engine. This is why we add the firsts() to
 * the domain.
 */

template<int Dim, class T, class EngineTag, class Domain>
struct LeafFunctor<Array<Dim, T, EngineTag>, ViewFunctorTag<Domain> >
{
  typedef typename View1<Array<Dim, T, EngineTag>, Domain>::Type_t Type_t;
  inline static Type_t apply(const Array<Dim, T, EngineTag> &a, 
    const ViewFunctorTag<Domain> &t) 
    {
      typedef View1<Array<Dim, T, EngineTag>, Domain> Ret_t;
      return Ret_t::make(a, t.domain_m + a.firsts());
    }
};

/**
 * This version of LeafFunctor is used by Expression-Engines to used to 
 * evaluate an Array using indices. 
 */

template<int Dim, class T, class EngineTag>
struct LeafFunctor<Array<Dim, T, EngineTag>, EvalLeaf<Dim> >
{
  typedef typename Array<Dim, T, EngineTag>::Element_t Type_t;
  inline static
  Type_t apply(const Array<Dim, T, EngineTag> &a, const EvalLeaf<Dim> &t) 
    {
      return t.eval(a.engine());
    }
};

/**
 * EngineView functor acting on array.  The functor is applied to the contained
 * engine and the result is packed back inside the array.
 */

template<int Dim, class T, class E, class Tag>
struct LeafFunctor<Array<Dim, T, E>, EngineView<Tag> >
{
  typedef LeafFunctor<Engine<Dim, T, E>, EngineView<Tag> > LeafFunctor_t;
  typedef typename LeafFunctor_t::Type_t NewEngine_t;

  typedef typename NewEngine_t::Tag_t NewTag_t;
  typedef Array<Dim, T, NewTag_t> Type_t;

  inline static
  Type_t apply(const Array<Dim, T, E> &array, 
	       const EngineView<Tag> &tag)
  {
    return Type_t(LeafFunctor_t::apply(array.engine(), tag));
  }
};

/**
 * ExpressionApply functor acting on Array.  The functor is applied to the
 * contained engine.
 */

template<int Dim, class T, class E, class Tag>
struct LeafFunctor<Array<Dim, T, E>, ExpressionApply<Tag> >
{
  typedef LeafFunctor<Engine<Dim, T, E>, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  inline static
  Type_t apply(const Array<Dim, T, E> &array, 
	       const ExpressionApply<Tag> &tag)
  {
    return LeafFunctor_t::apply(array.engine(), tag);
  }
};

/**
 * Apply the ConformTag to the leaves of the tree.
 * Check to see if a given Array conforms.
 *
 * First we have the case where the rank of the Array is the
 * same as the rank of the ConformTag.
 * Just loop over all of the dimensions, making sure that the
 * size of each fits with the size of the left hand side.
 */

template<int Dim, class T, class EngineTag>
struct LeafFunctor<Array<Dim, T, EngineTag>, ConformTag<Dim> >
{
  typedef bool Type_t;
  static Type_t apply(const Array<Dim, T, EngineTag> &array,
    const ConformTag<Dim> &ct)
    {
      return conforms(array.domain(), ct);
    }
};

/**
 * Now the case where the rank of the Array is not the same
 * as the ConformTag.  These cannot conform.
 */

template<int Dim1, int Dim2, class T, class EngineTag>
struct LeafFunctor<Array<Dim1, T, EngineTag>, ConformTag<Dim2> >
{
  typedef bool Type_t;
  static Type_t apply(const Array<Dim1, T, EngineTag> &,
    const ConformTag<Dim2> &)
    {
      return false;
    }
};

/**
 * Do what needs to be done prior to reading. For Arrays, this 
 * means doing nothing.
 */

template<int Dim, class T, class EngineTag>
struct LeafFunctor<Array<Dim, T, EngineTag>, NotifyPreReadTag>
{
  typedef bool Type_t;
  static Type_t apply(const Array<Dim, T, EngineTag> &a, 
    const NotifyPreReadTag &)
    {
      return true;
    }
};

/**
 * Generalized Engine Functors.
 */

template<int Dim, class T, class E, class Tag>
struct LeafFunctor<Array<Dim, T, E>, EngineFunctorTag<Tag> >
{
  typedef typename Array<Dim,T,E>::Engine_t Engine_t;
  typedef typename EngineFunctor<Engine_t,Tag>::Type_t Type_t;
  inline static
  Type_t apply(const Array<Dim, T, E> &array, const EngineFunctorTag<Tag> &tag)
  {
    return EngineFunctor<Engine_t,Tag>::apply(array.engine(), tag.tag());
  }
};

/**
 * A specialization of EngineFunctor for arrays.  This permits slightly more
 * generic programming, since we can just say engineFunctor(a, foo) instead
 * of saying engineFunctor(a.engine(), foo).  (Then the same function works
 * if a is an array or an engine.)
 */

template<int Dim, class T, class E, class Tag>
struct EngineFunctor<Array<Dim, T, E>, Tag>
{
  typedef typename EngineFunctor<Engine<Dim, T, E>, Tag>::Type_t Type_t;

  inline static 
  Type_t apply(const Array<Dim, T, E> &array,
	       const Tag &tag)
  {
    return engineFunctor(array.engine(), tag);
  }
};

/// Overload the << operator to print an Array to a stream.  This
/// uses the 'PrintArray' class to perform the formatting of the data.
/// It will create a default printer, print out the array with it, and
/// return the stream object.

template <int Dim, class T, class EngineTag>
std::ostream &operator<<(std::ostream &o, 
                         const Array<Dim, T, EngineTag> &ca)
{
  Pooma::blockAndEvaluate();
  PrintArray().print(o, ca);
  return o;
}

template <int Dim, class T, class EngineTag>
std::fstream &operator<<(std::fstream &f, 
                         const Array<Dim, T, EngineTag> &ca)
{
  Pooma::blockAndEvaluate();
  PrintArray().print(f, ca);
  return f;
}

/** Traits class for expressions containing arrays. */

struct ExpressionIsArray { };

template<int Dim, class T, class EngineTag>
struct ExpressionTraits<Array<Dim, T, EngineTag> >
{
  typedef ExpressionIsArray Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsArray, ExpressionIsArray>
{
  typedef ExpressionIsArray Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsArray, ExpressionIsScalar>
{
  typedef ExpressionIsArray Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsScalar, ExpressionIsArray>
{
  typedef ExpressionIsArray Type_t;
};

//---------------------------------------------------------------------------
// Implementation for assignment operators. We pack this assignment 
// expression into a PETE binary expression tree node and then use this 
// to construct an array with an expression engine. We then pass this on 
// to an evaluator, which handles the computation. The first version 
// handles assigning Arrays to Arrays and the second one handles 
// assigning scalars.
//
// If you wish to have Array work with other classes that derive from
// Array, define your own version of assign() that takes the following
// arguments:
//
//   assign(Array<Dim,T,EngineTag>, yourclass, Operator)
//
// where "yourclass" is the class that you would like to work on the
// right-hand side in an expression with an Array on the left-hand side.
//---------------------------------------------------------------------------

/// assign for Array = Array

template<int Dim, class T, class EngineTag,
  int OtherDim, class OtherT, class OtherEngineTag, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs,
       const Array<OtherDim, OtherT, OtherEngineTag> &rhs,
       const Op &op)
{
  // Optionally check conformance. 

  PAssert(forEach(rhs, ConformTag<Dim>(lhs.domain()), AndCombine()));

  // Just evaluate the sucker.

  Evaluator<MainEvaluatorTag>().evaluate(lhs, op, rhs);
  
  return lhs;
}

/// assign for Array = scalar

template<int Dim, class T, class EngineTag, class T1, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs, const T1 &rhs, const Op &op)
{
  // Make an expression out of the scalar.

  Array<Dim, T1, ConstantFunction> rhsExpr(lhs.domain());
  rhsExpr.engine().setConstant(rhs);
    
  // Evaluate.

  Evaluator<MainEvaluatorTag>().evaluate(lhs, op, rhsExpr);
    
  return lhs;
}

template<class Tree>
struct ConvertWhereProxy<ExpressionIsArray, Tree>
{
  typedef MakeReturn<Tree> Make_t;
};

template<int Dim, class T, class EngineTag, class F, class B, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs,
       const WhereProxy<F,B> &rhs,
       const Op &op)
{
  assign(lhs, rhs.whereMask(), rhs.opMask(op));

  return lhs;
}

/// Compute the number of elements that are currently compressed.

template<int Dim, class T, class EngineTag>
inline long
elementsCompressed(const Array<Dim, T, EngineTag> &a)
{
  return elementsCompressed(a.engine());
}

/// Return whether or not all of the elements are currently compressed.

template<int Dim, class T, class EngineTag>
inline bool
compressed(const Array<Dim, T, EngineTag> &a)
{
  return compressed(a.engine());
}

/// (Try to) compress the array.

template<int Dim, class T, class EngineTag>
inline void
compress(Array<Dim, T, EngineTag> &a)
{
  compress(a.engine());
}

/// Manually uncompress the array.
/// Partial specializations that could conceivably return nonzero fractions:

template<int Dim, class T, class EngineTag>
inline void
uncompress(Array<Dim, T, EngineTag> &a)
{
  uncompress(a.engine());
}

/**
 * Traits class telling RefCountedBlockPointer that this class has
 * shallow semantics and a makeOwnCopy method.
 */

template <int Dim, class T, class EngineTag>
struct ElementProperties< Array<Dim, T, EngineTag> > 
  : public MakeOwnCopyProperties< Array<Dim, T, EngineTag> >
{ };

#endif // POOMA_ARRAY_ARRAY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Array.h,v $   $Author: richi $
// $Revision: 1.154 $   $Date: 2004/11/04 14:15:58 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
