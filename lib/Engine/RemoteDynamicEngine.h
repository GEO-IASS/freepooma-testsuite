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
//   Engine<1, T, Remote<Dynamic> > - a specialization of remote that supports
//                                    dynamic operations
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_REMOTEDYNAMICENGINE_H
#define POOMA_ENGINE_REMOTEDYNAMICENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 *    A wrapper engine that remotifies an Engine<1, T, Dynamic>.
 *    The remote version belongs to a particular context.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Domain.h"
#include "Domain/Loc.h"
#include "Domain/SliceRange.h"
#include "Domain/IndirectionList.h"
#include "Engine/Engine.h"
#include "Engine/EnginePatch.h"
#include "Evaluator/EngineTraits.h"
#include "Evaluator/Evaluator.h"
#include "Evaluator/Reduction.h"
#include "Engine/EngineFunctor.h"
#include "Engine/RemoteEngine.h"
#include "Engine/DynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "Layout/Node.h"
#include "Utilities/algorithms.h"
#include "Utilities/PAssert.h"
#include "Utilities/RefCountedPtr.h"
#include "Tulip/Messaging.h"
#include "Tulip/SendReceive.h"
#include "Tulip/ReduceOverContexts.h"
#include "Tulip/RemoteProxy.h"


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
// Engine<1, T, Remote<Dynamic> > 
//
//  Engine<1, T, Remote<Dynamic> > is
//
//-----------------------------------------------------------------------------

template <class T>
class Engine<1, T, Remote<Dynamic> > 
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Engine<1, T, Remote<Dynamic> >       This_t;
  typedef Engine<1, T, Remote<Dynamic> >       Engine_t;
  typedef Engine<1, T, Dynamic>                LocalEngine_t;
  typedef DomainLayout<1>                      Layout_t;
  typedef Layout_t::PatchID_t                  PatchID_t;
  typedef Layout_t::CreateSize_t               CreateSize_t;
  typedef Interval<1>                          Domain_t;
  typedef T                                    Element_t;
  typedef T                                    ReadReturn_t;
  //  typedef typename LocalEngine_t::ReadReturn_t ReadReturn_t;
  typedef RemoteProxy<T>                       ElementRef_t;
  typedef Remote<Dynamic>                      Tag_t;

  enum { dimensions = 1 };
  enum { hasDataObject = true };
  enum { dynamic    = true };
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

  // Default constructor.

  Engine();

  // These constructors take an Interval<1> and set the owning
  // context to 0.  On context 0 we create a new LocalEngine_t.
  
  explicit 
  Engine(const Domain_t &domain);

  Engine(int owningContext, const Domain_t &domain);

  Engine(const Domain_t &domain, const T &elementModel);

  // This constructor takes a Node object, extracts the domain, and
  // creates a new LocalEngine_t on the context given by the Node.
  
  explicit
  Engine(const Node<Domain_t> &node);

  // Copy constructor should perform a SHALLOW copy.

  Engine(const Engine_t &model);
  Engine(const Engine_t &, const EngineConstructTag &);

  // Subsetting Constructors.  All the work of the subsetting is
  // deferred to the LocalEngine_t.

  template<class OtherEngine, class Domain>
  Engine(const OtherEngine &otherEngine, const Domain &domain);

  //============================================================
  // Destructor
  //============================================================

  // All pointer members are smart pointers, so this is trivial.

  ~Engine(); 


  //============================================================
  // Assignment operators
  //============================================================

  // Assigment should be SHALLOW, to be consistent with copy.

  Engine_t &operator=(const Engine_t &model);


  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  // Element access via Loc.

  ReadReturn_t read(const Loc<1> &) const;
  ElementRef_t operator()(const Loc<1> &) const;

  // Element access via ints for speed.

  ReadReturn_t read(int) const;
  ElementRef_t operator()(int) const;

  // Return the domain.

  inline const Domain_t &domain() const
  {
    return domain_m;
  }

  // Accessor function that checks if the engine is local.
  // (Really you can get this from owningContext(), but most of the
  // code we write is of the form if (local) ... else ...)

  inline
  bool engineIsLocal() const
  {
    return (Pooma::context() == owningContext_m);
  }

  // The owningContext() is the context that actually allocates
  // a local engine where the data is stored.

  inline
  int owningContext() const
  {
    return owningContext_m;
  }

  // Return a reference to the localEngine().  This operation only
  // makes sense on the context that owns the data.

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

  // Get a private copy of data viewed by this Engine.

  inline
  Engine_t &makeOwnCopy()
  {
    if (engineIsLocal() && localEnginePtr_m != NULL) 
    {
      // Ideally this would be localEnginePtr_m.makeOwnCopy();
      // but Shared<> doesn't implement ElementProperties correctly.
      LocalEngine_t engine(localEngine());
      engine.makeOwnCopy();
      localEnginePtr_m = LocalPtr_t(new LocalShared_t(engine));
    }
 
    return *this;
  }

  //============================================================
  // Dynamic interface methods.
  //============================================================

  // Create new elements by extending the current domain
  // on the local context by the requested number of elements.
  // Returns an Interval giving the domain of the newly created
  // elements.

  Interval<1> create(CreateSize_t num)
  {
    PAssert(engineIsLocal());
    Interval<1> newElems = localEngine().create(num);
    domain_m = localEngine().domain();
    return newElems;
  }

  // Delete the elements specified by the given domain. 
  // This backfills the deleted elements with elements from
  // the end of the list.
  
  template <class Dom>
  void destroy(const Dom &killList)
  {
    PAssert(engineIsLocal());
    localEngine().destroy(killList);
    domain_m = localEngine().domain();
  }
  
  template <class Iter>
  void destroy(Iter begin, Iter end)
  {
    PAssert(engineIsLocal());
    localEngine().destroy(begin, end);
    domain_m = localEngine().domain();
  }

  // Delete the elements specified by the given domain, or by a
  // pair of iterators into some sort of collection, and the
  // appropriate fill method. If offsetFlag is set to true, the 
  // domain is interpreted as a set of offsets rather than a set
  // of points in our domain.

  // Available fill mechanisms are backfill and shift-up, selected by
  // passing either a BackFill or ShiftUp object. BackFill will move 
  // elements from the end up to fill the holes. ShiftUp will shift
  // elements up to fill in holes. The latter is considerably slower,
  // but maintains the relative ordering of the elements, which may
  // be important for some applications.

  template <class Dom, class DeleteMethod>
  void destroy(const Dom &killList, 
               const DeleteMethod &method,
               bool offsetFlag = false)
  {
    PAssert(engineIsLocal());
    localEngine().destroy(killList, method, offsetFlag);
    domain_m = localEngine().domain();
  }
                      
  template <class Iter, class DeleteMethod>
  void destroy(Iter begin, 
               Iter end, 
               const DeleteMethod &method,
               bool offsetFlag = false)
  {
    PAssert(engineIsLocal());
    localEngine().destroy(begin, end, method, offsetFlag);
    domain_m = localEngine().domain();
  }

  // sync() function is a no-op for a single-patch engine.
  // This version of sync() may be called via the DynamicArray interface.

  void sync() { }

  // Modify the domain (but not the size) of this engine.
  // This version of sync() may be called by MultiPatchEngine on its patches.

  void sync(const Domain_t & d)
  {
    if (engineIsLocal())
    {
      localEngine().sync(d);
    }
    domain_m = d;
  }

#if POOMA_MESSAGING

  template <class Dom>
  int packSize(const Dom &packList) const
  {
    typedef Cheetah::Serialize<Cheetah::CHEETAH, T> Serializer_t;

    return packList.length() * Serializer_t::size(T());
  }

  int pack(const IndirectionList<int> &packList, char *buffer,
           bool zeroBasedDomain = true) const
  {
    typedef Cheetah::Serialize<Cheetah::CHEETAH, T> Serializer_t;

    const LocalEngine_t& lengine = localEngine();
    int length = packList.length();
    int i;
    int nBytes = 0;
    int change;
    int offset;

    // if the given domain is zero-based, add an offset
    offset = (zeroBasedDomain) ? lengine.domain().first() : 0;

    for (i = 0; i < length; ++i)
    {
      change = Serializer_t::pack(lengine(packList(i)+offset), buffer);
      buffer += change;
      nBytes += change;
    }

    return nBytes;
  }

  int unpack(const Interval<1> &unpackDomain, char *buffer,
             bool zeroBasedDomain = true)
  {
    typedef Cheetah::Serialize<Cheetah::CHEETAH, T> Serializer_t;

    LocalEngine_t& lengine = localEngine();
    int last = unpackDomain.last();
    int i;
    int nBytes = 0;
    int change;
    T *value;
    int offset;

    // if the given domain is zero-based, add an offset
    offset = (zeroBasedDomain) ? lengine.domain().first() : 0;

    for (i = unpackDomain.first(); i <= last; ++i)
    {
      change = Serializer_t::unpack(value, buffer);
      lengine(i+offset) = *value;
      buffer += change;
      nBytes += change;
    }

    return nBytes;
  }

#endif

private:

  //============================================================
  // Private data
  //============================================================

  // The remote-engine on owningContext_m actually owns the data.

  int owningContext_m;

  // Pointer to the local engine which only gets new'd on the owning
  // context.  Eventually this needs to be changed to some form of
  // shared object (or perhaps all the private data here will be
  // collected in a shared object).

  LocalPtr_t localEnginePtr_m;

  // The domain.  We don't just pull the domain out of the localEngine
  // because it doesn't exist on every context.

  Domain_t domain_m;
};


//////////////////////////////////////////////////////////////////////
//
// Inline implementation of the functions for Engine<D, T, Remote<Dynamic> >
//
//////////////////////////////////////////////////////////////////////

//
// Return the element specified by loc.
//

template <class T>
inline typename Engine<1, T, Remote<Dynamic> >::ReadReturn_t
Engine<1, T, Remote<Dynamic> >::read(const Loc<1> &loc) const
{
  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(loc);
  }
  return ElementRef_t(value, owningContext());
}

//
// Return the element specified by list of ints.
//

template <class T>
inline typename Engine<1, T, Remote<Dynamic> >::ReadReturn_t
Engine<1, T, Remote<Dynamic> >::read(int i1) const
{
  T value;
  if (engineIsLocal())
  {
    value = localEngine().read(i1);
  }
  return ElementRef_t(value, owningContext());
}

//
// Return a reference to the element specified by loc.
//

template <class T>
inline RemoteProxy<T> Engine<1, T, Remote<Dynamic> >::
operator()(const Loc<1> &loc) const
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

//
// Return a reference to the element specified by list of ints.
//

template <class T>
inline RemoteProxy<T> Engine<1, T, Remote<Dynamic> >::
operator()(int i1) const
{
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

//-----------------------------------------------------------------------------
//
// Engine<1, T, Remote<Dynamic> >()
//
// Constructs an empty Engine<1, T, Remote<Dynamic> >.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::Engine()
  : owningContext_m(0), localEnginePtr_m(NULL)
{
  PAssert(owningContext_m < Pooma::contexts());

  // In this case, we leave a null localEnginePtr_m...
  // Do we want to create an empty local engine?
}

//-----------------------------------------------------------------------------
//
// Engine<1, T, Remote<Dynamic> >(const Node<Interval<Dim> > &node)
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::Engine(const Node<Domain_t> &node)
  : owningContext_m(node.context()),
    domain_m(node.allocated())
{
  PAssert(owningContext_m < Pooma::contexts());

  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(node)));
  }
}


//-----------------------------------------------------------------------------
//
// Engine<1, T, Remote<Dynamic> >(const Domain_t &domain)
//
// Constructs a Remote-Engine holding type T elements with the
// multidimensional domain given by Interval<Dim>. Elements are
// initialized with the default constructor.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::Engine(const Domain_t &dom)
  : owningContext_m(0), domain_m(dom)
{
  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(domain_m)));
  }
}


//-----------------------------------------------------------------------------
//
// Engine<1, T, Remote<Dynamic> >(int owningContext, const Domain_t &domain)
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::Engine(int owningContext, const Domain_t &dom)
  : owningContext_m(owningContext),
    domain_m(dom)
{
  if (engineIsLocal())
  {
    localEnginePtr_m = LocalPtr_t(new LocalShared_t(LocalEngine_t(domain_m)));
  }
}


//-----------------------------------------------------------------------------
//
// Engine<1, T, Remote<Dynamic> >(const Interval<Dim> &domain, const T & model)
//
// Constructs a Remote-Engine holding type T elements with the
// multidimensional domain given by Interval<Dim>. Initializes these
// with a model.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::Engine(const Domain_t &dom, const T& model)
  : owningContext_m(0), domain_m(dom)
{
  if (engineIsLocal())
  {
    localEnginePtr_m =
      LocalPtr_t(new LocalShared_t(LocalEngine_t(domain_m, model)));
  }
}

//-----------------------------------------------------------------------------
//
// Engine<1, T, Remote<Dynamic> >(const Engine<D, T, Remote<Tag> > &)
//
// Copy constructor for Remote-Engine.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::
Engine(const Engine<1, T, Remote<Dynamic> > &modelEngine)
  : owningContext_m(modelEngine.owningContext()),
    localEnginePtr_m(modelEngine.localEnginePtr_m),
    domain_m(modelEngine.domain())
{
}

template <class T>
Engine<1, T, Remote<Dynamic> >::
Engine(const Engine<1, T, Remote<Dynamic> > &modelEngine,
       const EngineConstructTag &)
  : owningContext_m(modelEngine.owningContext()),
    localEnginePtr_m(modelEngine.localEnginePtr_m),
    domain_m(modelEngine.domain())
{
}

template <class T>
template<class OtherEngine, class Domain>
Engine<1, T, Remote<Dynamic> >::
Engine(const OtherEngine &otherEngine, const Domain &domain)
  : owningContext_m(otherEngine.owningContext())
{
  if (engineIsLocal())
  {
    localEnginePtr_m =
      LocalPtr_t(new LocalShared_t(LocalEngine_t(otherEngine.localEngine(),
						 domain)));
  }

  domain_m[1] = Interval<1>(domain[1].length());
}

//-----------------------------------------------------------------------------
//
// Engine<D,T,Remote> & operator=(const Engine<D,T,Remote> &)
//
// Assignment operator for Remote-Engines.
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> > & 
Engine<1, T, Remote<Dynamic> >::
operator=(const Engine<1, T, Remote<Dynamic> > &modelEngine)
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
// ~Engine<1, T, Remote<Dynamic> >()
//
//-----------------------------------------------------------------------------

template <class T>
Engine<1, T, Remote<Dynamic> >::~Engine()
{
}

//-----------------------------------------------------------------------------
//
// RemoteView, RemoteSend
//
// These two functor tags are used with engineFunctor() to generate brick-views
// from expressions or engines containing remote-brick-views.  On the receiving
// side you say:
//
// Engine<1, double, dynamic> a = engineFunctor(remotedynbrick, RemoteView());
//
// On the side that owns the data you say:
//
// expressionApply(remotedynbrick, RemoteSend(toContext));
//
// The receive operation (RemoteView) generates an engine or expression that
// contains the incoming data.  The send operation just sends the data and has
// no return.
//
//-----------------------------------------------------------------------------

template<class T>
struct LeafFunctor<Engine<1, T, Remote<Dynamic> >, EngineView<RemoteView> >
{
  typedef Engine<1, T, Remote<Dynamic> > Subject_t;

  typedef Engine<1, T, Dynamic> Type_t;

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

template<class T>
struct LeafFunctor<Engine<1, T, Remote<DynamicView> >,
  EngineView<RemoteView> >
{
  typedef Engine<1, T, Remote<DynamicView> > Subject_t;

  typedef Engine<1, T, DynamicView> Type_t;

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
      // Particle expressions should be aligned.  If there is a need
      // to do cross context assignment, then fill this space with
      // something that receives a remote engine.  (See RemoteEngine.h)

      PAssert(false);
      return engine.localEngine();  // (This will fail at run-time)
    }
  }
};

#if POOMA_MESSAGING

#include "Tulip/Messaging.h"

namespace Cheetah {

template<class T>
class Serialize<CHEETAH, Engine<1, T, Dynamic> >
{
public:
  typedef Engine<1, T, Dynamic> Engine_t;
  typedef Interval<1> Domain_t;

  static inline int
  size(const Engine_t &a)
  {
    int nBytes=0;
    
    nBytes += Serialize<CHEETAH, Domain_t>::size(a.domain());
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

    EngineElemSerialize op(buffer);

    change = EngineBlockSerialize::apply(op, a, dom);

    nBytes += change;

    return nBytes;
  }

  static inline int
  unpack(Engine_t* &a, char *buffer)
  {
    Interval<1> *dom;

    int change;
    int nBytes=0;

    change = Serialize<CHEETAH, Domain_t>::unpack(dom, buffer);
    buffer += change;
    nBytes += change;

    a = new Engine<1, T, Dynamic>(*dom);

    EngineElemDeSerialize op(buffer);

    change = EngineBlockSerialize::apply(op, *a, *dom);

    nBytes += change;

    return nBytes;
  }

  static inline void
  cleanup(Engine_t* a)
  {
    delete a;
  }
};

} // namespace Cheetah

#endif // POOMA_MESSAGING

/// checkDynamicID(obj, ID) is a specializable function that is used
/// by some classes to check the dynamic ID value stored in the first
/// argument by some means.  If it is the same as the given ID, this
/// returns false.  If it is not the same, it should return true and
/// changethe state of obj to indicate that it has "seen" the given ID.
///
/// That this function is required is very disturbing.

template<class T>
inline bool
checkDynamicID(Engine<1, T, Remote<Dynamic> > &be, ObserverEvent::ID_t did)
{
  PAssert(be.engineIsLocal());

  return checkDynamicID(be.localEngine(), did);
}

/// localPatchEngine() is a utility function used to perform operations on
/// multipatch engines where the patch engine could be a remote-engine.
/// Currently this function is used by the Multi-patch engine copy()
/// functions.

template<class T>
inline
Engine<1, T, Dynamic> &localPatchEngine(Engine<1, T, Remote<Dynamic> > &e)
{
  PAssert(e.engineIsLocal());
  return e.localEngine();
}

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_ENGINE_REMOTEDYNAMICENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RemoteDynamicEngine.h,v $   $Author: richard $
// $Revision: 1.23 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
