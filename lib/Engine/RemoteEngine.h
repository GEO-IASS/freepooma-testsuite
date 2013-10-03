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
//   Remote<Tag>                  - RemoteThing-Engine specialization tag.
//   Engine<Dim, T, Remote<Tag> > - the "RemoteThing-Engine" specialization.
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_REMOTEENGINE_H
#define POOMA_ENGINE_REMOTEENGINE_H

/** @file
 * @ingroup Engine
 * @brief Remote engine support.
 */

//-----------------------------------------------------------------------------
// Overview:
//
//   Remote<Tag>
//    - tag class used to select specializations of Engine
//
//   Engine<Dim, T, Remote<Tag> >
//    - A wrapper engine that remotifies an Engine<Dim, T, Tag>.
//      The remote version belongs to a particular context
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/Loc.h"
#include "Domain/Touches.h"
#include "Domain/SliceRange.h"
#include "Engine/Engine.h"
#include "Evaluator/EngineTraits.h"
#include "Evaluator/Evaluator.h"
#include "Evaluator/Reduction.h"
#include "Engine/EngineFunctor.h"
#include "Engine/ExpressionEngine.h"
#include "Engine/BrickEngine.h"
#include "Engine/CompressibleBrick.h"
#include "Engine/ForwardingEngine.h"
#include "Engine/Stencil.h"
#include "Layout/INode.h"
#include "Layout/Node.h"
#include "Layout/DomainLayout.h"
#include "Utilities/algorithms.h"
#include "Utilities/PAssert.h"
#include "Utilities/RefCounted.h"
#include "Utilities/RefCountedPtr.h"
#include "Tulip/Messaging.h"
#include "Tulip/SendReceive.h"
#include "Tulip/ReduceOverContexts.h"
#include "Tulip/RemoteProxy.h"

#include <algorithm>

// namespace Pooma {

/**
 * Tag class used to select the "RemoteBrick" and
 * "RemoteBrickView" specializations of the Engine class template.
 */

template<class Tag>
struct Remote
{ };

//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
///  Engine<Dim, T, Remote<Tag> > is
//-----------------------------------------------------------------------------

template <int Dim, class T, class Tag>
class Engine<Dim, T, Remote<Tag> > 
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================
  typedef Engine<Dim, T, Remote<Tag> >         This_t;
  typedef Engine<Dim, T, Remote<Tag> >         Engine_t;
  typedef Engine<Dim, T, Tag>                  LocalEngine_t;
  typedef Interval<Dim>                        Domain_t;
  typedef T                                    Element_t;
  //  typedef typename LocalEngine_t::ReadReturn_t ReadReturn_t;
  typedef RemoteProxy<T>                       ElementRef_t;
  typedef Remote<Tag>                          Tag_t;
  typedef DomainLayout<Dim>                    Layout_t;

  enum { dimensions = Dim };
  enum { hasDataObject = true };
  enum { dynamic    = false };
  enum { zeroBased  = false };
  enum { multiPatch = false };

  //============================================================
  // Internal convenience typedefs
  //============================================================

  typedef Shared<LocalEngine_t>        LocalShared_t;
  typedef RefCountedPtr<LocalShared_t> LocalPtr_t;

  //============================================================
  // Constructors and Factory Methods
  //============================================================

  /// Default constructor.

  Engine();

  /// These constructors take an Interval<Dim> and set the owning
  /// context to 0.  On context 0 we create a new LocalEngine_t.
  
  explicit 
  Engine(const Domain_t &domain);

  Engine(int owningContext, const Domain_t &domain);

  Engine(const Domain_t &domain, const T &elementModel);

  /// This constructor takes a Node object, extracts the domain, and
  /// creates a new LocalEngine_t on the context given by the Node.
  
  explicit
  Engine(const Node<Domain_t> &node);

  explicit
  Engine(const Layout_t &layout);

  /// Copy constructor should perform a SHALLOW copy.

  Engine(const Engine_t &model);
  Engine(const Engine_t &, const EngineConstructTag &);

  /// Subsetting Constructors.  All the work of the subsetting is
  /// deferred to the LocalEngine_t.

  template<class OtherEngine, class Domain>
  Engine(const OtherEngine &otherEngine, const Domain &domain);

  template<class OtherEngine, int D2>
  Engine(const OtherEngine &otherEngine,
	 const SliceRange<D2, Dim> &domain);

  //============================================================
  // Destructor
  //============================================================

  // All pointer members are smart pointers, so this is trivial.

  ~Engine(); 


  //============================================================
  // Assignment operators
  //============================================================

  /// Assigment should be SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);


  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  /// Element access via Loc.

  Element_t read(const Loc<Dim> &) const;
  ElementRef_t operator()(const Loc<Dim> &) const;

  /// Element access via ints for speed.

  Element_t read(int) const;
  Element_t read(int, int) const;
  Element_t read(int, int, int) const;
  Element_t read(int, int, int, int) const;
  Element_t read(int, int, int, int, int) const;
  Element_t read(int, int, int, int, int, int) const;
  Element_t read(int, int, int, int, int, int, int) const;

  ElementRef_t operator()(int) const;
  ElementRef_t operator()(int, int) const;
  ElementRef_t operator()(int, int, int) const;
  ElementRef_t operator()(int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int, int) const;
  ElementRef_t operator()(int, int, int, int, int, int, int) const;

  /// Return the domain.

  inline const Domain_t &domain() const
  {
    return domain_m;
  }

  /// Accessor function that checks if the engine is local.
  /// (Really you can get this from owningContext(), but most of the
  /// code we write is of the form if (local) ... else ...)

  inline
  bool engineIsLocal() const
  {
    return (Pooma::context() == owningContext_m);
  }

  /// The owningContext() is the context that actually allocates
  /// a local engine where the data is stored.

  inline
  int owningContext() const
  {
    return owningContext_m;
  }

  /// Return a reference to the localEngine().  This operation only
  /// makes sense on the context that owns the data.

  inline
  const LocalEngine_t &localEngine() const
  {
    PAssert(engineIsLocal());
    PAssert(localEnginePtr_m.isValid());
    return (*localEnginePtr_m).data();
  }

  inline
  LocalEngine_t &localEngine()
  {
    PAssert(engineIsLocal());
    PAssert(localEnginePtr_m.isValid());
    return (*localEnginePtr_m).data();
  }

  /// first() interface

  inline int first(int i) const
  {
    return domain_m[i].first();
  }

  /// Get a private copy of data viewed by this Engine.

  inline
  Engine_t &makeOwnCopy()
  {
    if (engineIsLocal() && localEnginePtr_m.isValid()) 
    {
      // Ideally this would be localEnginePtr_m.makeOwnCopy();
      // but Shared<> doesn't implement ElementProperties correctly.
      LocalEngine_t engine(localEngine());
      engine.makeOwnCopy();
      localEnginePtr_m = LocalPtr_t(new LocalShared_t(engine));
    }
 
    return *this;
  }

protected:

  /// The domain.  We don't just pull the domain out of the localEngine
  /// because it doesn't exist on every context.
  /// The domain is protected because RemoteDynamic engine is derived from
  /// Remote engine and needs to update the domain when sync() is called.

  Domain_t domain_m;

private:

  //============================================================
  // Private data
  //============================================================

  /// The remote-engine on owningContext_m actually owns the data.

  int owningContext_m;

  /// Pointer to the local engine which only gets new'd on the owning
  /// context.  Eventually this needs to be changed to some form of
  /// shared object (or perhaps all the private data here will be
  /// collected in a shared object).

  LocalPtr_t localEnginePtr_m;
};


//////////////////////////////////////////////////////////////////////
//
// Inline implementation of the functions for Engine<D, T, Remote<Tag> >
//
//////////////////////////////////////////////////////////////////////

/// Return the element specified by loc.

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(const Loc<Dim> &loc) const
{
  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(loc);
  }
  return ElementRef_t(value, owningContext());
}

/// Return the element specified by list of ints.

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1) const
{
  PAssert(Dim == 1);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1);
  }
  return ElementRef_t(value, owningContext());
}

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1, int i2) const
{
  PAssert(Dim == 2);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1, i2);
  }
  return ElementRef_t(value, owningContext());
}

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1, int i2, int i3) const
{
  PAssert(Dim == 3);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1, i2, i3);
  }
  return ElementRef_t(value, owningContext());
}

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1, int i2, int i3, int i4) const
{
  PAssert(Dim == 4);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1, i2, i3, i4);
  }
  return ElementRef_t(value, owningContext());
}

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1, int i2, int i3, int i4, int i5) const
{
  PAssert(Dim == 5);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1, i2, i3, i4, i5);
  }
  return ElementRef_t(value, owningContext());
}

template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1, int i2, int i3, int i4, int i5, int i6) const
{
  PAssert(Dim == 6);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1, i2, i3, i4, i5, i6);
  }
  return ElementRef_t(value, owningContext());
}


template <int Dim, class T, class Tag>
inline T Engine<Dim, T, Remote<Tag> >::
read(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
{
  PAssert(Dim == 7);

  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1, i2, i3, i4, i5, i6, i7);
  }
  return ElementRef_t(value, owningContext());
}


/// Return a reference to the element specified by loc.

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(const Loc<Dim> &loc) const
{
  if (engineIsLocal())
  {
    T &value = localEngine()(loc);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}

/// Return a reference to the element specified by list of ints..

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1) const
{
  PAssert(Dim == 1);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1, int i2) const
{
  PAssert(Dim == 2);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1, i2);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1, int i2, int i3) const
{
  PAssert(Dim == 3);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1, i2, i3);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1, int i2, int i3, int i4) const
{
  PAssert(Dim == 4);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1, i2, i3, i4);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1, int i2, int i3, int i4, int i5) const
{
  PAssert(Dim == 5);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1, i2, i3, i4, i5);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}

template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1, int i2, int i3, int i4, int i5, int i6) const
{
  PAssert(Dim == 6);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1, i2, i3, i4, i5, i6);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}


template <int Dim, class T, class Tag>
inline RemoteProxy<T> Engine<Dim, T, Remote<Tag> >::
operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const
{
  PAssert(Dim == 7);

  if (engineIsLocal())
  {
    T &value = localEngine()(i1, i2, i3, i4, i5, i6, i7);
    return ElementRef_t(value, owningContext());
  }
  else
  {
    T val;
    return ElementRef_t(val, owningContext());
  }
}


//-----------------------------------------------------------------------------
//
// Engine<Dim, T, Remote<Tag> >()
//
// Constructs an empty Engine<Dim, T, Remote<Tag> >.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::Engine()
  : owningContext_m(0)
{
  PAssert(owningContext_m < Pooma::contexts());

  // In this case, we leave a null localEnginePtr_m...
  // Do we want to create an empty local engine?
}

//-----------------------------------------------------------------------------
//
// Engine<Dim, T, Remote<Tag> >(const Node<Interval<Dim> > &node)
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::Engine(const Node<Domain_t> &node)
  : domain_m(node.allocated()),
    owningContext_m(node.context())
{
  PAssert(owningContext_m < Pooma::contexts());

  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(node)));
  }
}

/// This constructor basically ignores the context given by the DomainLayout,
/// because that context is currently bogus.  (It should be -1 when used for
/// Bricks and set to a specific context for RemoteBricks, not to
/// Pooma::context() which implies that everyone thinks that they own the data
/// and no one else owns the data.)

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::Engine(const Layout_t &layout)
  : domain_m(layout.node().allocated()),
    owningContext_m(0)
{
  PAssert(owningContext_m < Pooma::contexts());

  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(layout.
								  node())));
  }
}


/// Constructs a Remote-Engine holding type T elements with the
/// multidimensional domain given by Interval<Dim>. Elements are
/// initialized with the default constructor.

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::Engine(const Domain_t &dom)
  : domain_m(dom), owningContext_m(0)
{
  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(domain_m)));
  }
}


//-----------------------------------------------------------------------------
//
// Engine<Dim, T, Remote<Tag> >(int owningContext, const Domain_t &domain)
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::Engine(int owningContext, const Domain_t &dom)
  : domain_m(dom),
    owningContext_m(owningContext)
{
  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(domain_m)));
  }
}


/// Constructs a Remote-Engine holding type T elements with the
/// multidimensional domain given by Interval<Dim>. Initializes these
/// with a model.

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::Engine(const Domain_t &dom, const T& model)
  : domain_m(dom), owningContext_m(0)
{
  if (engineIsLocal())
  {
    localEnginePtr_m =
      LocalPtr_t(new LocalShared_t(LocalEngine_t(domain_m, model)));
  }
}

/// Copy constructor for Remote-Engine.

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::
Engine(const Engine<Dim, T, Remote<Tag> > &modelEngine)
  : domain_m(modelEngine.domain()),
    owningContext_m(modelEngine.owningContext()),
    localEnginePtr_m(modelEngine.localEnginePtr_m)
{
}

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::
Engine(const Engine<Dim, T, Remote<Tag> > &modelEngine,
       const EngineConstructTag &)
  : domain_m(modelEngine.domain()),
    owningContext_m(modelEngine.owningContext()),
    localEnginePtr_m(modelEngine.localEnginePtr_m)
{
}

template <int Dim, class T, class Tag>
template<class OtherEngine, class Domain>
Engine<Dim, T, Remote<Tag> >::
Engine(const OtherEngine &otherEngine, const Domain &domain)
  : owningContext_m(otherEngine.owningContext())
{
  if (engineIsLocal())
  {
    localEnginePtr_m =
      LocalPtr_t(new LocalShared_t(LocalEngine_t(otherEngine.localEngine(),
						 domain)));
  }

  int i;
  for (i = 0; i < Dim; ++i)
  {
    domain_m[i] = Interval<1>(domain[i].length());
  }
}

template <int Dim, class T, class Tag>
template<class OtherEngine, int D2>
Engine<Dim, T, Remote<Tag> >::
Engine(const OtherEngine &otherEngine, const SliceRange<D2, Dim> &domain)
  : owningContext_m(otherEngine.owningContext())
{
  if (engineIsLocal())
  {
    localEnginePtr_m =
      LocalPtr_t(new LocalShared_t(LocalEngine_t(otherEngine.localEngine(),
						 domain)));
  }

  int i;
  for (i = 0; i < Dim; ++i)
  {
    domain_m[i] = Interval<1>(domain.totalDomain()[i].length());
  }
}

/// Assignment operator for Remote-Engines.

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> > & 
Engine<Dim, T, Remote<Tag> >::
operator=(const Engine<Dim, T, Remote<Tag> > &modelEngine)
{
  // Can skip the rest if we're trying to assign to ourselves

  if (this == &modelEngine)
    return *this;

  // Copy values over from provided engine

  owningContext_m  = modelEngine.owningContext_m;
  domain_m         = modelEngine.domain_m;
  localEnginePtr_m = modelEngine.localEnginePtr_m;

  return *this;
}


//-----------------------------------------------------------------------------
//
// ~Engine<Dim, T, Remote<Tag> >()
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class Tag>
Engine<Dim, T, Remote<Tag> >::~Engine()
{
}

//-----------------------------------------------------------------------------
//
// NewEngine specializations for taking views.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class Tag, class Domain>
struct NewEngine<Engine<Dim, T, Remote<Tag> >, Domain>
{
  typedef typename NewEngine<Engine<Dim, T, Tag>, Domain>::Type_t Local_t;
  typedef typename Local_t::Tag_t NewTag_t;
  enum { newDim = Local_t::dimensions };
  typedef Engine<newDim, T, Remote<NewTag_t> > Type_t;
};

template<int Dim, class T, class Tag>
struct NewEngineDomain<Engine<Dim, T, Remote<Tag> >, INode<Dim> >
{
  typedef const Interval<Dim> &Type_t;
  static inline const Interval<Dim> &
  apply(const Engine<Dim, T, Remote<Tag> > &,
	const INode<Dim> &i)
  {
    return i.domain();
  }
};

//-----------------------------------------------------------------------------
//
// RemoteView, RemoteSend
//
// These two functor tags are used with engineFunctor() to generate brick-views
// from expressions or engines containing remote-brick-views.  On the receiving
// side you say:
//
// Engine<2, double, Brick> a = engineFunctor(remotebrick, RemoteView());
//
// On the side that owns the data you say:
//
// engineFunctor(remotebrick, RemoteSend(toContext));
//
// The receive operation (RemoteView) generates an engine or expression that
// contains the incoming data.  The send operation just sends the data and has
// no return.
//
//-----------------------------------------------------------------------------

struct RemoteView { };

template<>
struct EngineView<RemoteView>
{
  typedef TreeCombine Combine_t;
};

struct RemoteSend
{
  RemoteSend(int n) : toContext_m(n) { }

  // We're sending the remote brick information to this context:

  inline int toContext() const { return toContext_m; }
  int toContext_m;
};

template<class Engine>
struct DefaultExpressionApply<Engine, RemoteSend >
{
  typedef Engine Subject_t;
  typedef int Type_t;

  static inline
  Type_t apply(const Subject_t &engine,
	       const ExpressionApply<RemoteSend> &)
  {
    return 0;
  }
};

template<class Engine>
struct DefaultEngineView<Engine, RemoteView>
{
  typedef Engine Subject_t;
  typedef Engine Type_t;

  static inline
  Type_t apply(const Subject_t &engine,
	       const EngineView<RemoteView> &)
  {
    return engine;
  }
};

template<int Dim, class T, class Tag>
struct LeafFunctor<Engine<Dim, T, Remote<Tag> >, ExpressionApply<RemoteSend> >
{
  typedef Engine<Dim, T, Remote<Tag> > Subject_t;
  typedef int Type_t;

  static inline
  Type_t apply(const Subject_t &engine,
	       const ExpressionApply<RemoteSend> &sendTag)
  {
    if (engine.engineIsLocal())
    {
      if (sendTag.tag().toContext() == -1)
      {
	for (int i = 0; i < Pooma::contexts(); i++)
	  if (i != Pooma::context())
	    SendReceive::send(engine.localEngine(), i);
      }
      else
      {
	if (Pooma::context() != sendTag.tag().toContext())
	  SendReceive::send(engine.localEngine(), sendTag.tag().toContext());
      }
    }
      
    return 0;
  }
};

template<int Dim, class T, class Tag>
struct LeafFunctor<Engine<Dim, T, Remote<Tag> >, EngineView<RemoteView> >
{
};

template<int Dim, class T>
struct LeafFunctor<Engine<Dim, T, Remote<Brick> >, EngineView<RemoteView> >
{
  typedef Engine<Dim, T, Remote<Brick> > Subject_t;

  typedef Engine<Dim, T, Brick> Type_t;

  static inline
  Type_t apply(const Subject_t &engine,
	       const EngineView<RemoteView> &)
  {
    if (engine.engineIsLocal())
    {
      return engine.localEngine();
    }
    else
    {
      Type_t local(engine.domain());

      Receive<Type_t>::receive(local, engine.owningContext());
	  
      return local;
    }
  }
};

template<int Dim, class T>
struct LeafFunctor<Engine<Dim, T, Remote<BrickView> >,
  EngineView<RemoteView> >
{
  typedef Engine<Dim, T, Remote<BrickView> > Subject_t;

  typedef Engine<Dim, T, Brick> IncomingView_t;
  typedef Engine<Dim, T, Brick> Brick_t;
  typedef Engine<Dim, T, BrickView> Type_t;
 
  static inline
  Type_t apply(const Subject_t &engine,
	       const EngineView<RemoteView> &)
  {
    if (engine.engineIsLocal())
    {
      return engine.localEngine();
    }
    else
    {
      Interval<Dim> dom = engine.domain();

      Brick_t local(dom);

      Type_t view(local, dom);
	  
      Receive<IncomingView_t>::receive(view, engine.owningContext());

      return view;
    }
  }
};

template<int Dim, class T>
struct LeafFunctor<Engine<Dim, T, Remote<CompressibleBrick> >,
  EngineView<RemoteView>  >
{
  typedef Engine<Dim, T, Remote<CompressibleBrick> > Subject_t;

  typedef Engine<Dim, T, CompressibleBrick> Type_t;

  static inline
  Type_t apply(const Subject_t &engine, const EngineView<RemoteView>  &)
  {
    if (engine.engineIsLocal())
    {
      return engine.localEngine();
    }
    else
    {
      Type_t local(engine.domain());
	  
      Receive<Type_t>::receive(local, engine.owningContext());
	  
      return local;
    }
  }
};

template<int Dim, class T>
struct LeafFunctor<Engine<Dim, T, Remote<CompressibleBrickView> >,
  EngineView<RemoteView>  >
{
  typedef Engine<Dim, T, Remote<CompressibleBrickView> > Subject_t;

  typedef Engine<Dim, T, CompressibleBrick> CompBrick_t;
  typedef Engine<Dim, T, CompressibleBrickView> Type_t;

  static inline
  Type_t apply(const Subject_t &engine, const EngineView<RemoteView>  &)
  {
    if (engine.engineIsLocal())
    {
      return engine.localEngine();
    }
    else
    {
      Interval<Dim> dom = engine.domain();
      CompBrick_t local(dom);
      Type_t view(local, dom);
	  
      Receive<Type_t>::receive(view, engine.owningContext());
	  
      return view;
    }
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

struct EngineBlockSerialize
{
  template<class Op, class Eng>
  inline static int apply(Op &op, const Eng &engine)
    {
      typedef typename Eng::Domain_t Domain_t;
      return apply(op, engine, engine.domain(),
		   WrappedInt<Domain_t::dimensions>());
    }
  
  template<class Op, class Eng, class Dom>
  inline static int apply(Op &op, const Eng &engine, const Dom &domain)
    {
      return apply(op, engine, domain,
		   WrappedInt<Dom::dimensions>());
    }

  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain, WrappedInt<1>)
    {
      CTAssert(Domain::unitStride == 1);
      int f0 = domain[0].first();
      int e0 = domain[0].last();
      for (int i0 = f0; i0<=e0; ++i0)
	op(engine(i0));
      return op.total_m;
    }
  
  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain,WrappedInt<2>)
    {
      CTAssert(Domain::unitStride == 1);
      int f0 = domain[0].first();
      int f1 = domain[1].first();
      int e0 = domain[0].last();
      int e1 = domain[1].last();
      for (int i1 = f1; i1<=e1; ++i1)
	for (int i0 = f0; i0<=e0; ++i0)
	  op(engine(i0,i1));
      return op.total_m;
    }
  
  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain,WrappedInt<3>)
    {
      CTAssert(Domain::unitStride == 1);
      int f0 = domain[0].first();
      int f1 = domain[1].first();
      int f2 = domain[2].first();
      int e0 = domain[0].last();
      int e1 = domain[1].last();
      int e2 = domain[2].last();
      for (int i2 = f2; i2<=e2; ++i2)
	for (int i1 = f1; i1<=e1; ++i1)
	  for (int i0 = f0; i0<=e0; ++i0)
	    op(engine(i0,i1,i2));
      return op.total_m;
    }
  
  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain,WrappedInt<4>)
    {
      CTAssert(Domain::unitStride == 1);
      int f0 = domain[0].first();
      int f1 = domain[1].first();
      int f2 = domain[2].first();
      int f3 = domain[3].first();
      int e0 = domain[0].last();
      int e1 = domain[1].last();
      int e2 = domain[2].last();
      int e3 = domain[3].last();
      for (int i3 = f3; i3<=e3; ++i3)
	for (int i2 = f2; i2<=e2; ++i2)
	  for (int i1 = f1; i1<=e1; ++i1)
	    for (int i0 = f0; i0<=e0; ++i0)
	      op(engine(i0,i1,i2,i3));
      return op.total_m;
    }

  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain,WrappedInt<5>)
    {
      CTAssert(Domain::unitStride == 1);
      int f0 = domain[0].first();
      int f1 = domain[1].first();
      int f2 = domain[2].first();
      int f3 = domain[3].first();
      int f4 = domain[4].first();
      int e0 = domain[0].last();
      int e1 = domain[1].last();
      int e2 = domain[2].last();
      int e3 = domain[3].last();
      int e4 = domain[4].last();
      for (int i4 = f4; i4<=e4; ++i4)
	for (int i3 = f3; i3<=e3; ++i3)
	  for (int i2 = f2; i2<=e2; ++i2)
	    for (int i1 = f1; i1<=e1; ++i1)
	      for (int i0 = f0; i0<=e0; ++i0)
		op(engine(i0,i1,i2,i3,i4));
      return op.total_m;
    }

  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain,WrappedInt<6>)
    {
      CTAssert(Domain::unitStride == 1);
      int f0 = domain[0].first();
      int f1 = domain[1].first();
      int f2 = domain[2].first();
      int f3 = domain[3].first();
      int f4 = domain[4].first();
      int f5 = domain[5].first();
      int e0 = domain[0].last();
      int e1 = domain[1].last();
      int e2 = domain[2].last();
      int e3 = domain[3].last();
      int e4 = domain[4].last();
      int e5 = domain[5].last();
      for (int i5 = f5; i5<=e5; ++i5)
	for (int i4 = f4; i4<=e4; ++i4)
	  for (int i3 = f3; i3<=e3; ++i3)
	    for (int i2 = f2; i2<=e2; ++i2)
	      for (int i1 = f1; i1<=e1; ++i1)
		for (int i0 = f0; i0<=e0; ++i0)
		  op(engine(i0,i1,i2,i3,i4,i5));
      return op.total_m;
    }

  template<class Op,class Eng,class Domain>
  inline static int apply(Op &op,const Eng &engine,
			  const Domain &domain,WrappedInt<7>)
{
  CTAssert(Domain::unitStride == 1);
  int f0 = domain[0].first();
  int f1 = domain[1].first();
  int f2 = domain[2].first();
  int f3 = domain[3].first();
  int f4 = domain[4].first();
  int f5 = domain[5].first();
  int f6 = domain[6].first();
  int e0 = domain[0].last();
  int e1 = domain[1].last();
  int e2 = domain[2].last();
  int e3 = domain[3].last();
  int e4 = domain[4].last();
  int e5 = domain[5].last();
  int e6 = domain[6].last();
  for (int i6 = f6; i6<=e6; ++i6)
    for (int i5 = f5; i5<=e5; ++i5)
      for (int i4 = f4; i4<=e4; ++i4)
	for (int i3 = f3; i3<=e3; ++i3)
	  for (int i2 = f2; i2<=e2; ++i2)
	    for (int i1 = f1; i1<=e1; ++i1)
	      for (int i0 = f0; i0<=e0; ++i0)
		op(engine(i0,i1,i2,i3,i4,i5,i6));
  return op.total_m;
}

private:

};

#if POOMA_MESSAGING

#include "Tulip/Messaging.h"

struct EngineElemSerialize
{
  EngineElemSerialize(char *buffer)
    : buffer_m(buffer),
      total_m(0)
  { }


  template<class T>
  inline void operator()(const T &t)
  {
    int change = Cheetah::Serialize<Cheetah::CHEETAH, T>::pack(t, buffer_m);
    buffer_m   += change;
    total_m    += change;
  }

  char *buffer_m;
  int total_m;
};

struct EngineElemDeSerialize
{
  EngineElemDeSerialize(char *buffer)
    : buffer_m(buffer),
      total_m(0)
  { }

  template<class T>
  inline void operator()(T &t)
  {
    T *a;
    int change = Cheetah::Serialize<Cheetah::CHEETAH, T>::unpack(a, buffer_m);
    t = *a;
    buffer_m   += change;
    total_m    += change;
    Cheetah::Serialize<Cheetah::CHEETAH, T>::cleanup(a);
  }

  char *buffer_m;
  int total_m;
};


namespace Cheetah {

// All these serializers/deserializers share a common header,
// namely domain and compressed flag.

template<int Dim, class T>
class Serialize<CHEETAH, Engine<Dim, T, BrickView> >
{
public:
  typedef Engine<Dim, T, BrickView> Engine_t;
  typedef Interval<Dim> Domain_t;

  static inline int
  size(const Engine_t &a)
  {
    int nBytes=0;
    
    nBytes += Serialize<CHEETAH, Domain_t>::size(a.domain());
    bool compressed = false;
    nBytes += Serialize<CHEETAH, bool>::size(compressed);
    nBytes += a.domain().size() * Serialize<CHEETAH, T>::size(T());

    return nBytes;
  }

  static inline int
  pack(const Engine_t &a, char *buffer)
  {
    Domain_t dom = a.domain();

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::pack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool compressed = false;
    change = Serialize<CHEETAH, bool>::pack(compressed, buffer);
    buffer += change;
    nBytes += change;

    EngineElemSerialize op(buffer);

    change = EngineBlockSerialize::apply(op, a, dom);

    nBytes += change;

    return nBytes;
  }

  // We support a special unpack to avoid an extra copy.

  static inline int
  unpack(Engine_t &a, char *buffer)
  {
    Interval<Dim> *dom;

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::unpack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool *compressed;
    change = Serialize<CHEETAH, bool>::unpack(compressed, buffer);
    buffer += change;
    nBytes += change;

    // domains dont match probably, but at least their sizes must
    for (int i=0; i<Dim; ++i)
      PAssert((*dom)[i].size() == a.domain()[i].size());

    if (*compressed)
    {
      T *value;
      change = Serialize<CHEETAH, T>::unpack(value, buffer);

      // we can't use usual array assignment here, because this would
      // irritate the scheduler and lead to bogous results
      Array<Engine_t::dimensions, T, typename Engine_t::Tag_t> lhs;
      lhs.engine() = a;
      Array<Engine_t::dimensions, T, ConstantFunction> rhs(*dom);
      rhs.engine().setConstant(*value);
      KernelEvaluator<InlineKernelTag>::evaluate(lhs, OpAssign(), rhs);
    } else {
      EngineElemDeSerialize op(buffer);

      change = EngineBlockSerialize::apply(op, a, a.domain());
    }
    nBytes += change;

    Serialize<CHEETAH, Domain_t>::cleanup(dom);
    Serialize<CHEETAH, bool>::cleanup(compressed);

    return nBytes;
  }

};

template<int Dim, class T>
class Serialize<CHEETAH, Engine<Dim, T, Brick> >
{
public:
  typedef Engine<Dim, T, Brick> Engine_t;
  typedef Interval<Dim> Domain_t;

  static inline int
  size(const Engine_t &a)
  {
    int nBytes=0;
    
    nBytes += Serialize<CHEETAH, Domain_t>::size(a.domain());
    bool compressed = false;
    nBytes += Serialize<CHEETAH, bool>::size(compressed);
    nBytes += a.domain().size() * Serialize<CHEETAH, T>::size(T());

    return nBytes;
  }

  static inline int
  pack(const Engine_t &a, char *buffer)
  {
    Domain_t dom = a.domain();

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::pack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool compressed = false;
    change = Serialize<CHEETAH, bool>::pack(compressed, buffer);
    buffer += change;
    nBytes += change;

    EngineElemSerialize op(buffer);

    change = EngineBlockSerialize::apply(op, a, dom);

    nBytes += change;

    return nBytes;
  }

  // Old-style unpack with extra copy.

  static inline int
  unpack(Engine_t* &a, char *buffer)
  {
    Interval<Dim> *dom;

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::unpack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool *compressed;
    change = Serialize<CHEETAH, bool>::unpack(compressed, buffer);
    buffer += change;
    nBytes += change;
    PAssert(!*compressed);

    a = new Engine<Dim, T, Brick>(*dom);

    EngineElemDeSerialize op(buffer);

    change = EngineBlockSerialize::apply(op, *a, *dom);

    nBytes += change;

    Serialize<CHEETAH, Domain_t>::cleanup(dom);
    Serialize<CHEETAH, bool>::cleanup(compressed);

    return nBytes;
  }

  static inline void
  cleanup(Engine_t* a)
  {
    delete a;
  }

};

template<int Dim, class T>
class Serialize<CHEETAH, Engine<Dim, T, CompressibleBrick> >
{
public:
  typedef Engine<Dim, T, CompressibleBrick> Engine_t;
  typedef Interval<Dim> Domain_t;

  static inline int
  size(const Engine_t &a)
  {
    int nBytes=0;
    
    nBytes += Serialize<CHEETAH, Domain_t>::size(a.domain());

    // we cannot use a.compressed() here, because we need to
    // set up a big enough receive buffer and the compressed
    // flag is not valid across contexts.
    bool compressed = false;
    nBytes += Serialize<CHEETAH, bool>::size(compressed);

    if (compressed)
    {
      nBytes += Serialize<CHEETAH, T>::size(T());
    }
    else
    {
      nBytes += a.domain().size() * Serialize<CHEETAH, T>::size(T());
    }
    return nBytes;
  }

  static inline int
  pack(const Engine_t &a, char *buffer)
  {
    Domain_t dom = a.domain();

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::pack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool compressed = a.compressed();

    change = Serialize<CHEETAH, bool>::pack(compressed, buffer);
    buffer += change;
    nBytes += change;

    if (compressed)
    {
      change = Serialize<CHEETAH, T>::pack(a.compressedRead(), buffer);
    }
    else
    {
      EngineElemSerialize op(buffer);

      change = EngineBlockSerialize::apply(op, a, dom);
    }
    nBytes += change;

    return nBytes;
  }

  // Old-style unpack with extra copy.

  static inline int
  unpack(Engine_t* &a, char *buffer)
  {
    Interval<Dim> *dom;

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::unpack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool *compressed;
    change = Serialize<CHEETAH, bool>::unpack(compressed, buffer);
    buffer += change;
    nBytes += change;

    if (*compressed)
    {
      T *value;

      change = Serialize<CHEETAH, T>::unpack(value, buffer);

      a = new Engine_t(*dom, *value);
    }
    else
    {
      a = new Engine_t(*dom);

      EngineElemDeSerialize op(buffer);

      change = EngineBlockSerialize::apply(op, *a, *dom);
    }
    nBytes += change;

    Serialize<CHEETAH, Domain_t>::cleanup(dom);
    Serialize<CHEETAH, bool>::cleanup(compressed);

    return nBytes;
  }

  static inline void
  cleanup(Engine_t* a)
  {
    delete a;
  }

};

template<int Dim, class T>
class Serialize<CHEETAH, Engine<Dim, T, CompressibleBrickView> >
{
public:
  typedef Engine<Dim, T, CompressibleBrickView> Engine_t;
  typedef Interval<Dim> Domain_t;

  static inline int
  size(const Engine_t &a)
  {
    int nBytes=0;
    
    nBytes += Serialize<CHEETAH, Domain_t>::size(a.domain());

    // we cannot use a.compressed() here, because we need to
    // set up a big enough receive buffer and the compressed
    // flag is not valid across contexts.
    bool compressed = false;
    nBytes += Serialize<CHEETAH, bool>::size(compressed);

    if (compressed)
    {
      nBytes += Serialize<CHEETAH, T>::size(T());
    }
    else
    {
      nBytes += a.domain().size() * Serialize<CHEETAH, T>::size(T());
    }

    return nBytes;
  }

  static inline int
  pack(const Engine_t &a, char *buffer)
  {
    Domain_t dom = a.domain();

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::pack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool compressed = a.compressed();

    change = Serialize<CHEETAH, bool>::pack(compressed, buffer);
    buffer += change;
    nBytes += change;

    if (compressed)
    {
      change = Serialize<CHEETAH, T>::pack(a.compressedRead(), buffer);
    }
    else
    {
      EngineElemSerialize op(buffer);

      change = EngineBlockSerialize::apply(op, a, dom);
    }
    nBytes += change;

    return nBytes;
  }

  static inline int
  unpack(Engine_t* &a, char *buffer)
  {
    Interval<Dim> *dom;

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::unpack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool *compressed;

    change = Serialize<CHEETAH, bool>::unpack(compressed, buffer);
    buffer += change;
    nBytes += change;

    if (*compressed)
    {
      T *value;

      change = Serialize<CHEETAH, T>::unpack(value, buffer);

      Engine<Dim, T, CompressibleBrick> foo(*dom, *value);

      a = new Engine_t(foo, *dom);
    }
    else
    {
      Engine<Dim, T, CompressibleBrick> foo(*dom);

      EngineElemDeSerialize op(buffer);

      change = EngineBlockSerialize::apply(op, foo, *dom);

      a = new Engine_t(foo, *dom);
    }
    nBytes += change;

    return nBytes;
  }

  static inline void
  cleanup(Engine_t* a)
  {
    delete a;
  }

  // We support a special unpack to avoid an extra copy.

  static inline int
  unpack(Engine_t &a, char *buffer)
  {
    Interval<Dim> *dom;

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::unpack(dom, buffer);
    buffer += change;
    nBytes += change;

    bool *compressed;
    change = Serialize<CHEETAH, bool>::unpack(compressed, buffer);
    buffer += change;
    nBytes += change;

    // domains dont match probably, but at least their sizes must
    for (int i=0; i<Dim; ++i)
      PAssert((*dom)[i].size() == a.domain()[i].size());

    if (*compressed)
    {
      T *value;

      change = Serialize<CHEETAH, T>::unpack(value, buffer);

      // we can't use usual array assignment here, because this would
      // irritate the scheduler and lead to bogous results
      a.compressedReadWrite() = *value;
    }
    else
    {
      EngineElemDeSerialize op(buffer);

      change = EngineBlockSerialize::apply(op, a, *dom);
    }
    nBytes += change;

    Serialize<CHEETAH, Domain_t>::cleanup(dom);
    Serialize<CHEETAH, bool>::cleanup(compressed);

    return nBytes;
  }
};

} // namespace Cheetah

#endif // POOMA_MESSAGING


//-----------------------------------------------------------------------------
//
// Compressible support.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class Tag>
long elementsCompressed(const Engine<Dim, T, Remote<Tag> > &engine)
{
  return elementsCompressed(engine.localEngine());
}

template <int Dim, class T, class Tag>
void compress(Engine<Dim, T, Remote<Tag> > &engine)
{
  compress(engine.localEngine());
}

template <int Dim, class T, class Tag>
void uncompress(Engine<Dim, T, Remote<Tag> > &engine)
{
  uncompress(engine.localEngine());
}

template <int Dim, class T, class Tag>
bool compressed(const Engine<Dim, T, Remote<Tag> > &engine)
{
  return compressed(engine.localEngine());
}


/**
 * EngineFunctor for gathering up the contexts in an expression and
 * returning the most common. We need to use the PIMPL pattern below because
 * we need to retain on-board data and this tag can be wrapped in an
 * EngineFunctorTag object, which would normally trigger a copy.
 */

class GatherContexts
{
private:

  //===========================================================================
  // Nested class GatherContextsData. Allows retain our context
  // list as multiple copies of the GatherContexts tag are made.
  //===========================================================================
  
  class GatherContextsData : public RefCounted
  {
  public:
  
    //-------------------------------------------------------------------------
    // Trivial constructors and destructors.
    
    inline GatherContextsData() {}
    inline GatherContextsData(const GatherContextsData &model)
    : contexts_m(model.contexts_m) {}
    inline ~GatherContextsData() {}

    //-------------------------------------------------------------------------
    // Used to add a context to our list. If it is a real context (>=0), we
    // push it at the end of our vector. The empty() check is used to keep
    // the standard lib from using up unreasonably large amounts of memory 
    // as it expands the vector (some implementations use a default size of
    // 1024). If we have encountered an object that lives everywhere, indicated
    // by c == -q, we don't add it.
    
    void addContext(int c) const
      {
        if (c != -1)
          {
            if (contexts_m.empty())
              contexts_m.reserve(4);

	    contexts_m.push_back(c);
	  }
      }

    //-------------------------------------------------------------------------
    // Sorts the contexts and finds the most common one unless the vector
    // contained no entries, in which case we return -1. It is up to the caller
    // to decide if this makes any sense.
    
    int mostCommonContext() const
      {
        if (contexts_m.size() != 0)
          {
            std::sort(contexts_m.begin(), contexts_m.end());
            return
	      *Pooma::Algorithms::find_most_common(contexts_m.begin(), 
		    			           contexts_m.end());
  	  }
        else
          {
            return -1;
          }
      }
    
  private:

    //-------------------------------------------------------------------------
    // Our container is mutable so we can add to it via const member functions.
    // Tags are logically const so this is required.
    
    mutable std::vector<int> contexts_m;
  };    
    
public:

  //-------------------------------------------------------------------------
  // Required EngineFunctor typedef.
  
  typedef NullCombine Combine_t;

  //-------------------------------------------------------------------------
  // Simple constructors implementing shallow copy semantics for the data.

  inline GatherContexts()
  : data_m(new GatherContextsData) {}
  
  inline GatherContexts(const GatherContexts &model)
  : data_m(model.data_m) { }
  
  //-------------------------------------------------------------------------
  // Shallow assignment operator.
  
  GatherContexts &operator=(const GatherContexts &rhs)
    {
      data_m = rhs.data_m;
      return *this;
    }
  
  //-------------------------------------------------------------------------
  // The RefCountedPtr will do memory management when it goes away.
  
  inline ~GatherContexts() {}
     
  //-------------------------------------------------------------------------
  // Accessors and modifiers defer to the data object.
  
  inline void addContext(int c) const
    { data_m->addContext(c); }
   
  inline int mostCommonContext() const
    { return data_m->mostCommonContext(); }
      
private:

  //-------------------------------------------------------------------------
  // Our data, stored as a ref-counted pointer to simplify memory management.
  
  RefCountedPtr<GatherContextsData> data_m;

};

template<class T>
struct EngineFunctorScalar<T, GatherContexts>
{
  typedef int Type_t;
  static inline
  Type_t apply(const T &, const GatherContexts &)
    {
      return 0;
    }
};

template<class Engine>
struct EngineFunctorDefault<Engine, GatherContexts>
{
  typedef Engine Subject_t;
  typedef int Type_t;

  static inline
  Type_t apply(const Subject_t &, const GatherContexts &)
    {
      return 0;
    }
};

template<int Dim, class T, class Tag>
struct EngineFunctor<Engine<Dim, T, Remote<Tag> >, GatherContexts>
{
  typedef Engine<Dim, T, Remote<Tag> > Subject_t;
  typedef int Type_t;

  static inline
  Type_t apply(const Subject_t &engine, const GatherContexts &tag)
    {
      tag.addContext(engine.owningContext());
      
      return 0;
    }
};

//-----------------------------------------------------------------------------
//
// Specializations of evaluator for remote brick engines.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Single-patch Evaluator involving remote engines:
//
// This evaluator handles a single patch involving engines that may be remote.
//-----------------------------------------------------------------------------

template <>
struct Evaluator<RemoteSinglePatchEvaluatorTag>
{
  // Default ctor.

  Evaluator() { }

  // Destructor.

  ~Evaluator() { }

  // evaluate(expression)
  // Input an expression and cause it to be evaluated.
  // We just pass the buck to a special evaluator.

  template <class LHS, class RHS, class Op>
  void evaluate(const LHS &lhs, const Op &op, const RHS &rhs) const
  {
    GatherContexts gtag;
    engineFunctor(lhs.engine(), gtag);
    int lhsContext = gtag.mostCommonContext();

    expressionApply(rhs, RemoteSend(lhsContext));

    EngineView<RemoteView> view;

    if (lhsContext == -1 || Pooma::context() == lhsContext)
      {
        Evaluator<SinglePatchEvaluatorTag> speval;
        speval.evaluate(
		        forEach(lhs, view, TreeCombine()), op,
		        forEach(rhs, view, TreeCombine())
			);
      }
  }
};


//-----------------------------------------------------------------------------
// Multiple-patch Evaluator involving remote engines:
//
// The remote multiple patch version makes patches and sends them out to
// the remote single patch evaluator.
//-----------------------------------------------------------------------------

template <>
struct Evaluator<RemoteMultiPatchEvaluatorTag>
{
  // Default ctor.

  Evaluator() { }

  // Destructor.

  ~Evaluator() { }

  // evaluate(expression)
  // Input an expression and cause it to be evaluated.
  // We just pass the buck to a special evaluator.

  template <class LHS, class RHS, class Op>
  void evaluate(const LHS &lhs, const Op &op, const RHS &rhs) const
  {
    typedef Intersector<LHS::dimensions> Inter_t;
    Inter_t inter;

    expressionApply(lhs, IntersectorTag<Inter_t>(inter));
    expressionApply(rhs, IntersectorTag<Inter_t>(inter));
  
    typename Inter_t::const_iterator i = inter.begin();
    while (i != inter.end())
      {
        Evaluator<RemoteSinglePatchEvaluatorTag>().
          evaluate(lhs(*i), op, rhs(*i));
        ++i;
      }
  }
};

//-----------------------------------------------------------------------------
// Single-patch Reductions involving remote engines:
//
// This reduction handles a single patch involving engines that may be remote.
//-----------------------------------------------------------------------------

template <>
struct Reduction<RemoteSinglePatchEvaluatorTag>
{
  //---------------------------------------------------------------------------
  // Default ctor.

  Reduction() { }

  //---------------------------------------------------------------------------
  // Destructor

  ~Reduction() { }

  //---------------------------------------------------------------------------
  // Input an expression and cause it to be reduced. The procedure is as 
  // follows:
  //
  //   1. Decide which context the reduction will be performed on.
  //   2. If the current context is not the calculation context... 
  //        o and the data resides on this context, send it to the 
  //          calculation context.
  //        o look for the result of the reduction in a message from 
  //          the calculation context .
  //   3. If the current context is the calculation context... 
  //        o get a local view of the thing we're reducing and perform the 
  //          reduction.
  //        o send the result to the other contexts.
  
  template<class T, class Op, class Expr>
  void evaluate(T &ret, const Op &op, const Expr &e) const
  {
    GatherContexts gtag;
    engineFunctor(e.engine(), gtag);
    int computationContext = gtag.mostCommonContext();

    Pooma::CountingSemaphore csem;
    csem.height(1);

    Pooma::scheduler().beginGeneration();

    if (Pooma::context() != computationContext)
    {
      expressionApply(e, RemoteSend(computationContext));
      csem.incr();
    }
    else
    {
      EngineView<RemoteView> view;
      Reduction<SinglePatchEvaluatorTag>().
	evaluate(ret, op,
		 forEach(e, view, TreeCombine()), csem);
    }

    Pooma::scheduler().endGeneration();

    csem.wait();
#if POOMA_MPI
    // The above single thread waiting has the same problem as with
    // the MultiPatch variant.  So fix it.
    Pooma::blockAndEvaluate();
#endif

    RemoteProxy<T> globalRet(ret, computationContext);
    ret = globalRet;  
  }
};


//-----------------------------------------------------------------------------
// Multiple-patch Reduction involving remote engines:
//
// The multiple patch handles the case when some of the engines are remote.
//-----------------------------------------------------------------------------

template <>
struct Reduction<RemoteMultiPatchEvaluatorTag>
{
  //---------------------------------------------------------------------------
  // Default ctor.

  Reduction() { }

  //---------------------------------------------------------------------------
  // Destructor

  ~Reduction() { }

  //---------------------------------------------------------------------------
  // Input an expression and cause it to be reduced according to the 
  // computational scheme:
  //   1. Perform the intersection calculation to deduce the patches that 
  //      computation will proceed on.
  //   2. Determine the number of patches that will be computed on the current
  //      context and allocate an array that big.
  //   3. For each patch that is associated with the current context... 
  //        o if the current context is not the calculation context, send the 
  //          data (if necessary) to the calculation context.
  //        o if the current context is the calculation context, get a local 
  //          view of the thing we're reducing and perform the reduction.
  //   4. Perform the reduction over the patches local to this context. 
  //      This follows the pattern of the multi-patch reduction.
  //   5. Do a reduction over contexts on, say, context 0 by... 
  //        o performing an all-to-one communication to context 0.
  //        o doing the reduction
  //        o performing a boradcast back from context 0.
  
  template<class T, class Op, class Expr>
  void evaluate(T &ret, const Op &op, const Expr &e) const
  {
    typedef Intersector<Expr::dimensions> Inter_t;
    Inter_t inter;
    
    expressionApply(e, IntersectorTag<Inter_t>(inter));

    std::vector<bool> present(inter.size());
    std::vector<int> computationalContext(inter.size());
    typename Inter_t::const_iterator i = inter.begin();
    int j, k, n = 0;
    for (j = 0; j < inter.size(); j++)
      {
        present[j] = i->contextParticipates(Pooma::context());
	if (present[j])
          {
            computationalContext[j] = i->context();
            if (computationalContext[j] == Pooma::context())
              ++n;
          }
	++i;
      }

#if 0
#define SWH_C 1

    if (Pooma::context() == SWH_C)
      {
	for (j = 0; j < inter.size(); j++)
	  std::cerr << "p[" << j << "] = " << present[j] << std::endl;
	for (j = 0; j < inter.size(); j++)
	  std::cerr << "cc[" << j << "] = " << computationalContext[j] << std::endl;
	std::cerr << "n = " << n << std::endl;
      }
#endif

    Pooma::CountingSemaphore csem;
    csem.height(n);
    T *vals = new T[n];

    Pooma::scheduler().beginGeneration();

    i = inter.begin();
    k = 0;
    for (j = 0; j < inter.size(); j++)
      {
#if 0
	if (Pooma::context() == SWH_C)
	  {
	    std::cerr << "p[" << j << "] = " << present[j] << std::endl;
	    i->globalIDDataBase()->print(std::cerr);
	    std::cerr << std::endl;
	    std::cerr << *i << std::endl;
	  }
#endif
	if (present[j])
	  {
	    if (computationalContext[j] == Pooma::context())
	      {
		EngineView<RemoteView> view;
		Reduction<SinglePatchEvaluatorTag>().
		  evaluate(vals[k++], op,
		    forEach(e(*i).engine(), view, TreeCombine()), csem);
	      }
	    else
	    {
	      expressionApply(e(*i), RemoteSend(computationalContext[j]));
	    }
	  }
	
	++i;
      }

    Pooma::scheduler().endGeneration();
    csem.wait();
#if POOMA_MPI
    // We need to wait for Reductions on _all_ contexts to complete
    // here, as we may else miss to issue a igc update send iterate that a
    // remote context waits for.  Consider the 2-patch setup
    //  a,b     |         g|  |          g|
    // with the expressions
    //  a(I) = b(I+1);
    //  bool res = all(a(I) == 0);
    // here we issue the following iterates:
    //  0: guard receive from 1 (write request b)
    //  1: guard send to 0      (read request b)
    //  0/1: expression iterate (read request b, write request a)
    //  0/1: reduction (read request a)
    //  0/1: blocking MPI_XXX
    // here the guard send from 1 to 0 can be skipped starting the
    // blocking MPI operation prematurely while context 0 needs to
    // wait for this send to complete in order to execute the expression.
    //
    // The easiest way (and the only available) is to blockAndEvaluate().
    Pooma::blockAndEvaluate();
#endif

    if (n > 0)
      {
	ret = vals[0];
	for (j = 1; j < n; j++)
	  op(ret, vals[j]);
      }

    delete [] vals;

    ReduceOverContexts<T, Op> finalReduction(ret, 0, n > 0);
    if (Pooma::context() == 0)
      ret = finalReduction;

    RemoteProxy<T> broadcast(ret, 0);
    ret = broadcast;
  }
};


//---------------------------------------------------------------------------
// Specialization of EngineFunctor for EnginePatch for
// MultiPatch<Remote<Tag>>.
// (needed since you really want the local engine.)
//---------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class Tag>
struct EngineFunctor<Engine<Dim, T, MultiPatch<LayoutTag, Remote<Tag> > >,
  EnginePatch >
{
  typedef Engine<Dim, T, MultiPatch<LayoutTag, Remote<Tag> > > Subject_t;

  typedef Engine<Dim, T, Tag>               Type_t;

  static inline
  Type_t apply(const Subject_t &engine, const EnginePatch &tag)
  {
    return engine.localPatch(tag.patch_m).localEngine();
  }
};

//-----------------------------------------------------------------------------
// Traits class telling RefCountedBlockPointer that this class has
// shallow semantics and a makeOwnCopy method.
//-----------------------------------------------------------------------------

template <int Dim, class T, class Eng>
struct ElementProperties<Engine<Dim, T, Remote<Eng> > > 
  : public MakeOwnCopyProperties<Engine<Dim, T, Remote<Eng> > >
{ };

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_ENGINE_REMOTEENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RemoteEngine.h,v $   $Author: richard $
// $Revision: 1.45 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
