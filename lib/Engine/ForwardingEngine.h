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
//   ComponentAccess<T>
//   CompFwd<Eng, N>
//   Engine<Dim, T, CompFwd<Eng, N> >
//   NewEngine< Engine<Dim, T, CompFwd<Eng, N> >, Domain >
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_FORWARDINGENGINE_H
#define POOMA_ENGINE_FORWARDINGENGINE_H

#include "Domain/Loc.h" 
#include "Functions/ComponentAccess.h"
#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Engine/EnginePatch.h"
#include "Engine/NotifyEngineWrite.h"
#include "PETE/PETE.h"

/** @file
 * @ingroup Engine
 * @brief
 * A ForwardingEngine is used to forward indices to the elements of another
 * engine.
 */

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template <int Dim> class DomainLayout;

/**
 * The component forwarding tag class.
 */

template<class Eng, class Components>
struct CompFwd { };

/**
 * A ForwardingEngine is used to forward indices to the elements of another
 * engine.
 */

template<int Dim, class T, class Eng, class Components>
class Engine<Dim, T, CompFwd<Eng, Components> >
{
public:

  //---------------------------------------------------------------------------
  // This_t is a convenience typedef referring to this class.

  typedef Engine<Dim, T, Eng> This_t;

  //---------------------------------------------------------------------------
  // Engine_t is a typedef that makes the template parameter Engine accessible 
  // to other classes. ElemEngine_t makes the engine we're forwarding to
  // visible.
  
  typedef This_t Engine_t;
  typedef Eng ElemEngine_t;

  //---------------------------------------------------------------------------
  // Element_t is the type of elements managed by this engine. 
  // ElementRef_t is the type that is used to write to a single element. 
  // This might be a reference or a proxy object. 
  // The typedef Domain_t gives the type of domain this engine is defined on. 
  // The enum 'dimensions' gives the dimensionality of this engine.
  // Tag_t is my tag.
  // Layout_t is the type of layout I'd cough up, if asked.

  typedef typename Eng::Element_t FwdElement_t;
  typedef ComponentAccess<FwdElement_t, Components> CompAccess_t;
  typedef typename CompAccess_t::Element_t  Element_t;
  typedef typename CompAccess_t::ElementRef_t  ElementRef_t;
  typedef typename Eng::Domain_t Domain_t;
  typedef CompFwd<Eng, Components> Tag_t;
  typedef typename Eng::Layout_t Layout_t;

  //---------------------------------------------------------------------------
  // required constants

  enum { dimensions = Eng::dimensions };
  enum { hasDataObject = Eng::hasDataObject };
  enum { dynamic = false };
  enum { zeroBased = Eng::zeroBased };
  enum { multiPatch = Eng::multiPatch };

  //---------------------------------------------------------------------------
  // Empty constructor required for containers of engines.

  Engine()
    : engine_m(), components_m()
  { }

  //---------------------------------------------------------------------------
  // This is the most basic way to build a forwarding engine. We take an
  // engine and a Components giving the components we're supposed to forward.

  Engine(const Eng &e, const Components &l)
    : engine_m(e), components_m(l)
  { }

  //---------------------------------------------------------------------------
  // Copy constructor.

  Engine(const This_t &e)
    : engine_m(e.elemEngine()), components_m(e.components()) { }

  //---------------------------------------------------------------------------
  // View constructor. This is used to take a view of another forwarding
  // engine. OtherEngine should be the type computed by taking a view of
  // Engine_t using Domain_t.

  template<class OtherEng, class Domain>
  Engine(const Engine< Dim, T, CompFwd<OtherEng, Components> > &e,
		const Domain &domain)
    : engine_m(NewEngineEngine<OtherEng,Domain>::apply(e.elemEngine(),domain),
	       NewEngineDomain<OtherEng,Domain>::apply(e.elemEngine(),domain)),
      components_m(e.components()) { }

  //---------------------------------------------------------------------------
  // Destructor. Trivial since the engine_m and components_m objects have their
  // destructors called automatically.

  ~Engine() { }

  //---------------------------------------------------------------------------
  // Element accessor. Get engine()'s element by passing the domain to
  // its operator(). Then, forward the components to the element using
  // Element_t's operator().

  inline ElementRef_t operator()(const Loc<dimensions> &eloc) const
  {
    return CompAccess_t::indexRef(elemEngine()(eloc), components());
  }

  inline ElementRef_t operator()(int i1) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1), components());
  }

  inline ElementRef_t operator()(int i1, int i2) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1, i2), components());
  }

  inline ElementRef_t operator()(int i1, int i2, int i3) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1, i2, i3), 
				  components());
  }

  inline ElementRef_t operator()(int i1, int i2, int i3, int i4) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1, i2, i3, i4), 
				  components());
  }

  inline ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1, i2, i3, i4, i5), 
				  components());
  }

  inline ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5,
				 int i6) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1, i2, i3, i4, i5, i6), 
				  components());
  }

  inline ElementRef_t operator()(int i1, int i2, int i3, int i4, int i5,
				 int i6, int i7) const
  {
    return CompAccess_t::indexRef(elemEngine()(i1, i2, i3, i4, i5, i6, i7),
				  components());
  }
      
  //---------------------------------------------------------------------------
  // Read-only element accessor. Get engine()'s element by passing the domain
  // to its operator(). Then, forward the components to the element using
  // Element_t's operator().

  inline Element_t read(const Loc<dimensions> &eloc) const
  {
    return CompAccess_t::index(elemEngine().read(eloc), components());
  }

  inline Element_t read(int i1) const
  {
    return CompAccess_t::index(elemEngine().read(i1), components());
  }

  inline Element_t read(int i1, int i2) const
  {
    return CompAccess_t::index(elemEngine().read(i1, i2), components());
  }

  inline Element_t read(int i1, int i2, int i3) const
  {
    return CompAccess_t::index(elemEngine().read(i1, i2, i3), 
			       components());
  }

  inline Element_t read(int i1, int i2, int i3, int i4) const
  {
    return CompAccess_t::index(elemEngine().read(i1, i2, i3, i4), 
			       components());
  }

  inline Element_t read(int i1, int i2, int i3, int i4, int i5) const
  {
    return CompAccess_t::index(elemEngine().read(i1, i2, i3, i4, i5), 
			       components());
  }

  inline Element_t read(int i1, int i2, int i3, int i4, int i5,
			int i6) const
  {
    return CompAccess_t::index(elemEngine().read(i1, i2, i3, i4, i5, i6), 
			       components());
  }

  inline Element_t read(int i1, int i2, int i3, int i4, int i5,
			int i6, int i7) const
  {
    return CompAccess_t::index(elemEngine().read(i1, i2, i3, i4, i5, i6, i7),
			       components());
  }

  //---------------------------------------------------------------------------
  // Returns the layout, which is acquired from the contained engine.

  inline const Layout_t& layout() const
  {
    return elemEngine().layout();
  }

  inline Layout_t& layout()
  {
    return elemEngine().layout();
  }

  //---------------------------------------------------------------------------
  // Returns the domain, which is acquired from the contained engine.
  
  inline const Domain_t& domain() const { return elemEngine().domain(); }

  //---------------------------------------------------------------------------
  // Return the first value for the specified direction.
  
  inline int first(int i) const
  {
    return elemEngine().first(i);
  }

  //---------------------------------------------------------------------------
  // Get a private copy of this engine. Simply forwards the 'makeOwnCopy'
  // request to the contained engine.

  This_t &makeOwnCopy() 
  { 
    elemEngine().makeOwnCopy(); 
    return *this; 
  }
  
  //---------------------------------------------------------------------------
  // Assessor functions that return the engine and components.
  
  Eng &elemEngine() { return engine_m; }
  const Eng &elemEngine() const { return engine_m; }
  
  const Components &components() const { return components_m; }

private:
  
  Eng engine_m;
  Components components_m;

};


/**
 * Specialization allowing a view to be taken of a forwarding engine.
 */

template <int Dim, class T, class Eng, class Components, class Domain>
struct NewEngine<Engine<Dim, T, CompFwd<Eng, Components> >, Domain>
{
  typedef typename NewEngine<Eng, Domain>::Type_t NewEngine_t;
  typedef Engine<NewEngine_t::dimensions, T, CompFwd<NewEngine_t, Components> > Type_t;
};

/**
 * General version of engineFunctor to passes the request to
 * the contained engine.
 */

template<int Dim, class T, class Eng, class Components, class EFTag>
struct EngineFunctor<Engine<Dim, T, CompFwd<Eng, Components> >, EFTag>
{
  typedef typename EngineFunctor<Eng, EFTag>::Type_t Type_t;

  static Type_t
  apply(const Engine<Dim, T, CompFwd<Eng, Components> > &engine,
	const EFTag &tag)
  {
    return engineFunctor(engine.elemEngine(), tag);
  }
};

template <int D, class T, class E, class Comp, class Tag>
struct LeafFunctor<Engine<D, T, CompFwd<E, Comp> >, EngineView<Tag> >
{
  typedef LeafFunctor<E, EngineView<Tag> > LeafFunctor_t;
  typedef typename LeafFunctor_t::Type_t NewViewed_t;
  typedef Engine<D, T, CompFwd<NewViewed_t, Comp> > Type_t;

  static
  Type_t apply(const Engine<D, T, CompFwd<E, Comp> > &engine,
	       const EngineView<Tag> &tag)
  {
    return Type_t(LeafFunctor_t::apply(engine.elemEngine(), tag),
		  engine.components());
  }
};

template <int D, class T, class E, class Comp, class Tag>
struct LeafFunctor<Engine<D, T, CompFwd<E, Comp> >, ExpressionApply<Tag> >
{
  typedef LeafFunctor<E, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  static
  Type_t apply(const Engine<D, T, CompFwd<E, Comp> > &engine,
	       const ExpressionApply<Tag> &tag)
  {
    return LeafFunctor_t::apply(engine.elemEngine(), tag);
  }
};

/**
 * Tell contained engine that it's dirty.
 */

template<int Dim, class T, class Eng, class Components>
struct NotifyEngineWrite<Engine<Dim,T,CompFwd<Eng,Components> > >
{
  inline static void
  notify(const Engine<Dim,T,CompFwd<Eng,Components> > &engine)
  {
    typedef typename Engine<Dim, T, 
      CompFwd<Eng, Components> >::ElemEngine_t Engine_t;
    NotifyEngineWrite<Engine_t>::notify(engine.elemEngine());
  }
};

/**
 * Version of EnginePatch that gets the patch from the viewed engine.
 */

template <int D, class T, class E, class Comp>
struct EngineFunctor<Engine<D, T, CompFwd<E, Comp> >, EnginePatch>
{
  typedef typename EngineFunctor<E, EnginePatch>::Type_t NewViewed_t;
  typedef Engine<D, T, CompFwd<NewViewed_t, Comp> > Type_t;

  static
  Type_t apply(const Engine<D, T, CompFwd<E, Comp> > &engine,
	       const EnginePatch &tag)
  {
    return Type_t(engineFunctor(engine.elemEngine(), tag),
		  engine.components());
  }
};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ForwardingEngine.h,v $   $Author: richard $
// $Revision: 1.50 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
