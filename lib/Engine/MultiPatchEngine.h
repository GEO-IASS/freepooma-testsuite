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
//   MultiPatch           
//     - MultiPatch-Engine specialization tag.
//   MultiPatchView       
//     - View specialization tag.
//   Engine<Dim,T,MultiPatch>  
//     - the "MultiPatch-Engine" specialization.
//   Engine<Dim,T,MultiPatchView>
//     - the "MultiPatchView-Engine" specialization.
//   NewEngine<Engine,SubDomain> 
//     - specializations for MultiPatchView-Engine.
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_MULTIPATCHENGINE_H
#define POOMA_ENGINE_MULTIPATCHENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 * MultiPatch engine:
 *  - MultiPatch & MultiPatchView
 *    tag classes used to select specializations of Engine
 *  - Engine<Dim,T,MultiPatch>
 *    an Engine that manages data stored in an array of PatchEngines.
 *  - Engine<Dim,T,MultiPatchView>
 *    an Engine that manages a view into a MultiPatch-Engine.
 *  - NewEngine
 *    Specialized traits class for creating MultiPatchView-Engines. 
 *    See Engine.h for the general version.
 *  - NewEngineEngine, NewEngineDomain
 *    fuctors that tell arrays how to get the right patch engine
 *    when making a view (so that Brick doesn't have to know about
 *    MultiPatchEngine.h).
 */

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Grid.h"
#include "Layout/MultiPatchLayoutTraits.h"
#include "Engine/Engine.h"
#include "Engine/Intersector.h"
#include "Engine/IntersectEngine.h"
#include "Engine/NotifyEngineWrite.h"
#include "Engine/EnginePatch.h"
#include "Utilities/Observer.h"
#include "Utilities/ObserverEvent.h"
#include "Utilities/WrappedInt.h"
#include "Layout/DynamicEvents.h"

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template <class DT> class SliceDomain;
template<class Tag> struct Remote;

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//-----------------------------------------------------------------------------
/// These are tag classes used to select the "MultiPatch" and
/// "MultiPatchView" specializations of the Engine class template.
//-----------------------------------------------------------------------------

template <class LayoutTag, class PatchTag>
struct MultiPatch
{
  MultiPatch(){}
  ~MultiPatch(){}
};

template <class LayoutTag, class PatchTag, int Dim2>
struct MultiPatchView 
{ 
  MultiPatchView(){}
  ~MultiPatchView(){}
};

//-----------------------------------------------------------------------------
// Forward declarations of the Engines.
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
class Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >;

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
class Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag,Dim2> >;

//-----------------------------------------------------------------------------
// NewEngine SPECIALIZATIONS
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Several specializations of NewEngine for combinations of 
// Engines and Domains that produce BrickView-Engines.
// Also specializations of NewEngineEngine and NewEngineDomain
// where we want to construct a patch engine from the MultiPatch
// engine.
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
struct NewEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  Interval<Dim> >
{
  typedef Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim> > Type_t;
};

template <int Dim, class T, class LayoutTag, class PatchTag>
struct NewEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  Range<Dim> >
{
  typedef Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim> > Type_t;
};

template <int Dim, class T, class LayoutTag, class PatchTag, class Domain>
struct NewEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  Node<Domain> >
{
  typedef typename 
    NewEngine<Engine<Dim, T, PatchTag>, Node<Domain> >::Type_t Type_t;
};

template <int Dim, class T, class LayoutTag, class PatchTag, class Domain>
struct NewEngineEngine<
  Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> >, 
  Node<Domain> >
{
  typedef Engine<Dim,T,PatchTag> &Type_t;
  static inline Engine<Dim,T,PatchTag> &
  apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &e,
	const Node<Domain> &i)
  {
    return e.globalPatch(i);
  }
};

template <int Dim, class T, class LayoutTag, class PatchTag, class Domain>
struct NewEngineDomain<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  Node<Domain> >
{
  typedef const Domain &Type_t;
  static inline const Domain &
  apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &,
	const Node<Domain> &i)
  {
    return i.domain();
  }
};

template <int Dim, class T, class LayoutTag, class PatchTag>
struct NewEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  INode<Dim> >
{
  typedef typename 
    NewEngine<Engine<Dim, T, PatchTag>, Interval<Dim> >::Type_t Type_t;
};

template <int Dim,  class T, class LayoutTag, class PatchTag>
struct NewEngineEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  INode<Dim> >
{
  typedef Engine<Dim,T,PatchTag> &Type_t;
  static inline Engine<Dim,T,PatchTag> &
  apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &e,
	const INode<Dim> &i)
  {
    return e.globalPatch(i);
  }
};

template <int Dim,  class T, class LayoutTag, class PatchTag>
struct NewEngineDomain<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  INode<Dim> >
{
  typedef const Interval<Dim> &Type_t;
  static inline const Interval<Dim> &
  apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &,
	const INode<Dim> &i)
  {
    return i.domain();
  }
};

template <int Dim, class T, class LayoutTag, class PatchTag, int SliceDim>
struct NewEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  SliceInterval<Dim,SliceDim> >
{
  typedef 
    Engine<SliceDim, T, MultiPatchView<LayoutTag, PatchTag, Dim> > Type_t;
};

template <int Dim, class T, class LayoutTag, class PatchTag, int SliceDim>
struct NewEngine<
  Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >, 
  SliceRange<Dim,SliceDim> >
{
  typedef 
    Engine<SliceDim, T, MultiPatchView<LayoutTag, PatchTag, Dim> > Type_t;
};

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
struct NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >, 
  Interval<Dim> >
{
  typedef 
    Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > Type_t;
};

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
struct NewEngine< 
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >, 
  Range<Dim> >
{
  typedef 
    Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > Type_t;
};

template <int Dim, class T, 
          class LayoutTag, class PatchTag, int Dim2, class Domain>
struct NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >, 
  Node<Domain> >
{
  typedef typename 
    NewEngine<Engine<Dim2, T, PatchTag>, SliceRange<Dim2, Dim> >::Type_t
      Type_t;
};

// Specialization of above for Dim == Dim2...
// Note that we always lose the unit stride information when doing this!!!

template <int Dim, class T, 
          class LayoutTag, class PatchTag, class Domain>
struct NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim> >, 
  Node<Domain> >
{
  typedef typename 
    NewEngine<Engine<Dim, T, PatchTag>, Range<Dim> >::Type_t
      Type_t;
};

template <int Dim,  class T, 
          class LayoutTag, class PatchTag, int Dim2, class Domain>
struct NewEngineEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag, Dim2> >, 
  Node<Domain> >
{
  typedef typename NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag, Dim2> >, 
    Node<Domain> >::Type_t Type_t;
  static inline Type_t
  apply(const Engine<Dim,T,MultiPatchView<LayoutTag,PatchTag,Dim2> > &e,
	const Node<Domain> &i)
  {
    return e.globalPatch(i);
  }
};

template <int Dim,  class T, 
          class LayoutTag, class PatchTag, int Dim2, class Domain>
struct NewEngineDomain<
  Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag, Dim2> >, 
  Node<Domain> >
{
  typedef EngineConstructTag Type_t;
  static inline EngineConstructTag
  apply(const Engine<Dim,T,MultiPatchView<LayoutTag,PatchTag,Dim2> > &,
	const Node<Domain> &)
  {
    return EngineConstructTag();
  }
};

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
struct NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >, 
  INode<Dim> >
{
  typedef typename 
    NewEngine<Engine<Dim2, T, PatchTag>, SliceRange<Dim2, Dim> >::Type_t
      Type_t;
};

// Specialization of above for Dim == Dim2.
// Note that we lose the unit stride info here.

template <int Dim, class T, class LayoutTag, class PatchTag>
struct NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim> >, 
  INode<Dim> >
{
  typedef typename 
    NewEngine<Engine<Dim, T, PatchTag>, Range<Dim> >::Type_t
      Type_t;
};

template <int Dim,  class T, class LayoutTag, class PatchTag, int Dim2>
struct NewEngineEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag, Dim2> >, 
  INode<Dim> >
{
  typedef typename NewEngine<
  Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag, Dim2> >, 
    INode<Dim> >::Type_t Type_t;
  static inline Type_t
  apply(const Engine<Dim,T,MultiPatchView<LayoutTag,PatchTag,Dim2> > &e,
	const INode<Dim> &i)
  {
    return e.globalPatch(i);
  }
};

template <int Dim,  class T, class LayoutTag, class PatchTag, int Dim2>
struct NewEngineDomain<
  Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag, Dim2> >, 
  INode<Dim> >
{
  typedef EngineConstructTag Type_t;
  static inline EngineConstructTag
  apply(const Engine<Dim,T,MultiPatchView<LayoutTag,PatchTag,Dim2> > &,
	const INode<Dim> &)
  {
    return EngineConstructTag();
  }
};

template <int Dim, class T, 
          class LayoutTag, class PatchTag, int Dim2, int SliceDim>
struct NewEngine< 
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >, 
  SliceInterval<Dim, SliceDim> >
{
  typedef 
    Engine<SliceDim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > Type_t;
};

template <int Dim, class T, 
          class LayoutTag, class PatchTag, int Dim2, int SliceDim>
struct NewEngine< 
  Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >, 
  SliceRange<Dim, SliceDim> >
{
  typedef 
    Engine<SliceDim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > Type_t;
};


/**
 * The multi-patch engine manages a Dim-dimensional logical "brick" of data 
 * of type T. The data is distributed amongst a list of patches of type 
 * Engine<Dim, T, PatchTag>. The union of the "owned" domains of the patches 
 * (the patch-engine domains ignoring the portions that are used for storing 
 * internal guard data) is equal to the domain of the multi-patch engine. 
 *
 * The Engine makes no assumptions about the type of T beyond that it have a
 * copy constructor and assignment operator. (Aren't these deferred to the
 * underlying RCBPtr???)
 *
 * The Domain of this engine is an Interval<Dim>, which is a tensor product 
 * of Dim 1-D intervals.
 *
 * The organization of the patches, and all services related to looking up 
 * which patches correspond to which domains, is deferred to the layout, whose
 * type is MultiPatchLayoutTraits<LayoutTag,Dim>::Layout_t. The type of the 
 * view layout is 
 *   MultiPatchLayoutTraits<LayoutTag,ViewDim>::View1<Dim>::Layout_t.
 *
 * Subsetting Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> > returns an
 * Engine<NewDim, T, MultiPatchView<LayoutTag,PatchTag,Dim> >. See below.
 */

template <int Dim, class T, class LayoutTag, class PatchTag>
class Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> > :
  public Observer<typename MultiPatchLayoutTraits<LayoutTag,Dim>::Layout_t>
{
public:
  
  //===========================================================================
  // Exported and convenience typedefs and constants
  //===========================================================================

  typedef MultiPatch<LayoutTag,PatchTag>                  Tag_t;
  typedef Engine<Dim,T,Tag_t>                             This_t;
  typedef Engine<Dim,T,Tag_t>                             Engine_t;
  typedef Interval<Dim>                                   Domain_t;
  typedef T                                               Element_t;
  typedef PatchTag                                        PatchTag_t;
  typedef Engine<Dim, T, PatchTag>                        PatchEngine_t;
  typedef typename PatchEngine_t::ElementRef_t            ElementRef_t;
  typedef RefCountedBlockPtr<PatchEngine_t>               PatchContainer_t;
  typedef MultiPatchLayoutTraits<LayoutTag,Dim>           LayoutTraits_t;
  typedef typename LayoutTraits_t::Layout_t               Layout_t;
  typedef Layout_t                                        Observable_t;
  typedef DynamicEvents::PatchID_t                        PatchID_t;
  typedef DynamicEvents::CreateSize_t                     CreateSize_t;
  typedef ObserverEvent::ID_t                             DynamicID_t;

  typedef typename NewEngine<PatchEngine_t, Domain_t>::Type_t PatchView_t;

  enum { brick = false };
  enum { dimensions = Dim };
  enum { hasDataObject = false };
  enum { dynamic = PatchEngine_t::dynamic };
  enum { zeroBased = false };
  enum { multiPatch = true };


  //===========================================================================
  // Constructors and factory methods
  //===========================================================================

  //---------------------------------------------------------------------------
  /// The default constructor is available, but must be followed by use
  /// of operator= to set this object up in a usable state.
  
  Engine();

  //---------------------------------------------------------------------------
  /// The constructor takes a Layout object.
  
  explicit Engine(const Layout_t &layout);

  //---------------------------------------------------------------------------
  /// Copy constructor.

  Engine(const Engine_t &model);


  //===========================================================================
  // Destructor
  //===========================================================================

  /// We need to detach from the layout that we're viewing.

  ~Engine();


  //===========================================================================
  /// Assignment operators
  //===========================================================================

  Engine_t &operator=(const Engine_t &model);


  //===========================================================================
  // Accessor and mutator functions
  //===========================================================================

  //---------------------------------------------------------------------------
  ///@name Element access via Loc.
  //@{

  Element_t read(const Loc<Dim> &) const;
  ElementRef_t operator()(const Loc<Dim> &) const;

  //@}

  //---------------------------------------------------------------------------
  ///@name Element access via ints for speed.
  //@{

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

  //@}

  //---------------------------------------------------------------------------
  /// Return a patch given a Node or INode.

  inline PatchEngine_t &globalPatch(const INode<Dim> &inode) const
  {
    // get the global id from the inode.

    int gid = inode.globalID(layout_m.ID());

    PatchEngine_t &pengine = data()[gid];
    
    // Make sure our initial assumption that we're fully contained
    // in the patch is correct.

    PAssert(contains(pengine.domain(), inode.domain()));
    
    return pengine;
  }

  template<class Domain>
  inline PatchEngine_t &globalPatch(const Node<Domain> &node) const
  {
    int gid = node.globalID();  

    PatchEngine_t &pengine = data()[gid];
    
    // Make sure that we're fully contained in the patch.

    PAssert(contains(pengine.domain(), node.domain()));
    
    return pengine;
  }

  inline PatchEngine_t &globalPatch(PatchID_t gpatchID) const
  {
    return data()[gpatchID];
  }
  
  inline PatchEngine_t &localPatch(PatchID_t lpatchID) const
  {
    int gID = (layout_m.nodeListLocal()[lpatchID])->globalID();
    return data()[gID];
  }

  inline bool patchEmpty(PatchID_t patch) const
  {
    return data()[patch].domain().empty();
  }

  //---------------------------------------------------------------------------
  /// Return the layout.

  inline const Layout_t &layout() const
  {
    return layout_m;
  }

  inline Layout_t &layout()
  {
    return layout_m;
  }

  //---------------------------------------------------------------------------
  /// Return the domain and base domain.

  inline const Domain_t &domain() const
  {
    return layout().domain();
  }

  inline const Domain_t &innerDomain() const
  {
    return layout().innerDomain();
  }

  //---------------------------------------------------------------------------
  /// Return the first index in the given direction.

  inline int first(int d) const
  {
    return layout().first(d);
  }

  //------------------------------------------------------------
  /// Return whether we have been initialized yet or not.

  inline bool initialized() const
  {
    return layout().initialized();
  }

  //---------------------------------------------------------------------------
  /// Return the data.

  inline const PatchContainer_t &data() const
  {
    return data_m;
  }

  //---------------------------------------------------------------------------
  /// Get a private copy of data viewed by this Engine.

  Engine_t &makeOwnCopy();
  
  //---------------------------------------------------------------------------
  /// Fill the internal guard cells.
  
  inline void fillGuards(const GuardLayers<Dim>& g) const 
  { 
    fillGuardsHandler(g, WrappedInt<Layout_t::supportsGuards>());
  }

  inline void fillGuards() const
  {
    fillGuards(layout().internalGuards());
  }
  
  inline void fillGuardsHandler(const GuardLayers<Dim>&, const WrappedInt<false>&) const { };
  void fillGuardsHandler(const GuardLayers<Dim>&, const WrappedInt<true>&) const ;
  
  //---------------------------------------------------------------------------
  /// Set the internal guard cells to a particular value.
  
  void setGuards(const T &val) const;
  
  //---------------------------------------------------------------------------
  /// Accumulate from the internal guards into 

  void accumulateFromGuards() const;

  //---------------------------------------------------------------------------
  /// Set and get the dirty flag (fillGuards is a no-op unless the 
  /// dirty flag is true).
    
  inline int dirty() const { return *pDirty_m; }

  inline void setDirty() const
  {
    *pDirty_m = (1<<(Dim*2))-1;
  }

  inline void clearDirty(int face = -1) const
  {
    if (face == -1)
      *pDirty_m = 0;
    else {
      PAssert(face >= 0 && face <= Dim*2-1);
      *pDirty_m &= ~(1<<face);
    }
  }

  inline bool isDirty(int face = -1) const
  {
    if (face == -1)
      return *pDirty_m != 0;
    else {
      PAssert(face >= 0 && face <= Dim*2-1);
      return *pDirty_m & (1<<face);
    }
  }

  //============================================================
  // Observer methods
  //============================================================

  /// Be notified of various events from the layout, including
  /// when the layout is deleted or when dynamic operations occur.

  virtual void notify(Observable_t &observed, const ObserverEvent &event);

  /// Handler for dynamic events. 
  /// Notify delegates dynamic event handling to one of these functions
  /// by calling dynamicHandler with WrappedInt<PatchEngine_t::dynamic>.
  /// The "true" version handles dynamic events, and requires that the
  /// patch engine have the dynamic interface.
  /// The "false" version should never be called, but has to be compilable.
  
  void dynamicHandler(Observable_t &, const ObserverEvent &, 
                      const WrappedInt<false> &);
  void dynamicHandler(Observable_t &, const ObserverEvent &, 
                      const WrappedInt<true> &);
  
  
  //============================================================
  // Dynamic interface methods.
  //============================================================

  /// Create new elements by extending the current domain
  /// on the local context by the requested number of elements, or the
  /// specified local patch.
  /// 'local' means on this same context.  The patch is refered to
  /// by local index, from 0 ... # local patches - 1.
  /// The default of -1 puts elements on the last local patch.

  inline void create(CreateSize_t num, PatchID_t localPatchID = (-1))
    {
      layout().create(num, localPatchID);
    }

  /// Delete the elements within our domain specified by the given
  /// Range, using either a backfill mechanism or a shift-up mechanism
  /// to delete the data.  The second argument should be one of the
  /// following types:
  ///  - BackFill() will move elements from the bottom up to fill the holes.
  ///  - ShiftUp() will shift elements up to fill in holes.
  /// The domain must be within the total domain of the DynamicArray.

  template<class Dom, class DeleteMethod>
  inline void destroy(const Dom &killlist, const DeleteMethod &method)
  {
    layout().destroy(killlist, method);
  }

  template<class Dom>
  inline void destroy(const Dom &killlist)
  {
    layout().destroy(killlist, BackFill());
  }

  /// Delete the elements within the specific local domain for the given
  /// patch.  The domain values in this case should all be zero-based,
  /// so that they are relative to the first element of the specified
  /// local patch.  The domain values should all be contained within the
  /// specified local patch as well.  To perform cross-patch destroy's,
  /// use the form where you do not specify a local patch number.

  template<class Dom, class DeleteMethod>
  inline void destroy(const Dom &killlist, PatchID_t frompatch,
		      const DeleteMethod &method)
    {
      layout().destroy(killlist, frompatch, method);
    }

  /// Copy all elements of domain n to the end of the last patch (default) or
  /// to the end of the specified patch.

  template<class Dom>
  inline void copy(const Dom &domain, PatchID_t topatch = (-1))
    {
      layout().copy(domain, topatch);
    }

  /// Copy all elements from the specified local patch to the end of
  /// the local to-patch.  The domain values in this case should all be zero-
  /// based, so that they are relative to the first element of the specified
  /// local patch.  The domain values should all be contained within the
  /// specified local patch as well.  To perform cross-patch copies,
  /// use the form where you do not specify a local from-patch number.

  template<class Dom>
  inline void copy(const Dom &killlist, 
                   PatchID_t frompatch, PatchID_t topatch)
    {
      layout().copy(killlist, frompatch, topatch);
    }

  /// Perform a "multiple patch" copy, using a list of IndirectionList's
  /// for a set of source patches, and an IndirectionList giving the
  /// patch ID for the source patches.  Copy data into the destination
  /// patch.  The source and desination patches must be specified, this
  /// is only for "zero-based" index lists.  If the last argument is
  /// true, storage is created at the end, otherwise elements are
  /// just copied to the end of the existing storage.

  inline void copy(const IndirectionList< IndirectionList<int> > &domlists,
		   const IndirectionList< int > &fromlist,
		   PatchID_t topatch,
		   bool docreate)
    {
      layout().copy(domlists, fromlist, topatch, docreate);
    }

  /// Synchronize all the contexts to update their domain information.
  /// This should be used after create/destroy operations have modified
  /// the domain of local context's data, and all contexts must be told
  /// of the new situation.  This should be an SPMD call.

  void sync()
    {
      layout().sync();
    }

private:

  //===========================================================================
  // Private dynamic routines
  //===========================================================================

  /// Carry out a request to perform a create operation in a particular
  /// patch.  The layout is responsible for figuring out what patch to do
  /// this in, so the patch number must be a valid index into our local
  /// patch list.

  void performCreate(CreateSize_t num, PatchID_t patch,
		     DynamicID_t did);

  /// Carry out the work to perform a destroy operation on a particular
  /// patch.  The layout is responsible for figuring out what patch to do
  /// this in, so the patch number must be a valid index into our local
  /// patch list.  Also, the domain must be a "relative" domain, with zero-
  /// based values.

  template<class Dom, class DeleteMethod>
  void performDestroy(const Dom &killlist, PatchID_t patch,
		      const DeleteMethod &method,
		      DynamicID_t did);

  /// Carry out the work to perform a copy of values from one patch
  /// to another (or perhaps to the same patch).
  /// The layout is responsible for figuring out what patch to do
  /// this in, so the patch number must be a valid index into our local
  /// patch list.  Also, the domain must be a "relative" domain, with zero-
  /// based values.

  template<class Dom>
  void performCopy(const Dom &domain, PatchID_t frompatch, PatchID_t topatch,
		   DynamicID_t did);

  /// Do the actual work of a multiple-list copy.

  void performPatchCopy(const IndirectionList< IndirectionList<int> > &dlists,
			const IndirectionList< int > &fromlist,
			PatchID_t topatch,
			bool docreate,
			DynamicID_t did);

  //===========================================================================
  /// Special Runnable class for allocating patches
  //===========================================================================

  template <class Node, class Counter>
  class PatchAllocator : public Pooma::Runnable_t
  {
  public:
  
    PatchAllocator(PatchEngine_t &dest, const Node &node, Counter &c)
      : Pooma::Runnable_t(node.affinity()), dest_m(dest),
	node_m(node), counter_m(c)
    {
      // std::cout<< " aff = " << node.affinity() <<std::endl;
    }
    
    ~PatchAllocator() { } // Trivial. 
    
    void run()
    {
      dest_m = PatchEngine_t(node_m);
      ++counter_m;
    }
    
  private:
  
    PatchEngine_t &dest_m;
    const Node    &node_m;
    Counter       &counter_m;
    
  };

  //===========================================================================
  // Data
  //===========================================================================

  /// Layout for this engine.

  Layout_t layout_m;

  /// The actual data is stored in a ref-counted block of Engines.

  PatchContainer_t data_m;
  
  /// Flag indicating whether internal guard cells need to be filled.
  /// We store a pointer to the flag since all copies and views
  /// must share the same flag. We use the reference count in
  /// data_m to decide whether to clean this up.

  int *pDirty_m;
};


/**
 * Helper traits class for MultiPatchView - should never need 
 * SliceRange if Dim == Dim2.
 */

template <int D1, int D2>
struct SubDomainTraits
{
  typedef SliceRange<D1,D2> LocalToBase_t;
  typedef Range<D1>         TotalDomain_t;
  inline static 
  TotalDomain_t totalDomain(const LocalToBase_t &d)
  {
    return d.totalDomain();
  }
};

template <int D1>
struct SubDomainTraits<D1, D1>
{
  typedef Range<D1>         LocalToBase_t;
  typedef Range<D1>         TotalDomain_t;
  inline static
  TotalDomain_t totalDomain(const LocalToBase_t &dom)
  {
    return dom;
  }
};

      
/**
 * Multi-patch-view-engine manages a view of a multi-patch engine.
 *
 * The Domain of this engine is an Interval<Dim>, which is a tensor
 * product of Dim 1-D intervals. For view-engines, these intervals
 * are 0-based (i.e. [0..N0]x[0..N1] etc.).
 *
 * The Dim2 parameter gives the dimension of the MP that ultimately
 * spawned this MPView. This does not necessarily need to equal Dim because
 * of the possibility of a slice.
 *  
 * Note that this is NOT the domain of the underlying data storage,
 * but rather it is the domain as presented to the outside world.
 */

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
class Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag,Dim2> >
{
public:
  
  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef MultiPatchView<LayoutTag, PatchTag, Dim2>       Tag_t;
  typedef Engine<Dim,T,Tag_t>                             This_t;
  typedef Engine<Dim,T,Tag_t>                             Engine_t;
  typedef Engine<Dim2,T,MultiPatch<LayoutTag,PatchTag> >  ViewedEngine_t;
  typedef Interval<Dim>                                   Domain_t;
  typedef T                                               Element_t;
  typedef PatchTag                                        PatchTag_t;
  typedef Engine<Dim2, T, PatchTag_t>                     PatchEngine_t;
  typedef typename PatchEngine_t::ElementRef_t            ElementRef_t;
  typedef RefCountedBlockPtr<PatchEngine_t>               PatchContainer_t;
  
  // The traits class is for the underlying layout as we need
  // to get the type of the appropriate view of that layout.
  
  typedef MultiPatchLayoutTraits<LayoutTag,Dim2>    ViewedLayoutTraits_t;
  typedef typename ViewedLayoutTraits_t::Layout_t   ViewedLayout_t;

  // For EGCS (doesn't like foo::bar::splat)

#ifdef __MWERKS__
  typedef typename ViewedLayoutTraits_t::View<Dim>
                                                    NestedViewTraits_t;
#else
  typedef typename ViewedLayoutTraits_t::template View<Dim>
                                                    NestedViewTraits_t;
#endif

  typedef typename NestedViewTraits_t::Layout_t     Layout_t;

  enum { dimensions = Dim };
  enum { hasDataObject = false };
  enum { dynamic = PatchEngine_t::dynamic };
  enum { zeroBased = true };
  enum { multiPatch = true };

  //===========================================================================
  // Constructors and factory methods
  //===========================================================================

  //---------------------------------------------------------------------------
  /// The default constructor is available, but must be followed by use
  /// of operator= to set this object up in a usable state.
  
  Engine() { }

  //---------------------------------------------------------------------------
  /// The constructors take a MultiPatch and a non-slice domain like an
  /// Interval<Dim2> or a Range<Dim2>.
  
  template<class DT>
  Engine(const ViewedEngine_t &engine, const Domain<Dim2, DT> &domain)
  : layout_m(engine.layout(), domain),
    baseEngine_m(engine)
  { }

  //---------------------------------------------------------------------------
  /// The constructors take a MultiPatch and a slice domain like a
  /// SliceInterval<Dim2,Dim> or a SliceRange<Dim2,Dim>.
  
  template<class DT>
  Engine(const ViewedEngine_t &engine, const SliceDomain<DT> &domain)
  : layout_m(engine.layout(), domain),
    baseEngine_m(engine)
  { }

  //---------------------------------------------------------------------------
  /// The constructors take a MultiPatchView and a non-slice domain like
  /// an Interval<Dim2> or a Range<Dim2>.

  template<class DT>
  Engine(const This_t &engine, const Domain<Dim, DT> &domain)
  : layout_m(engine.layout(), domain),
    baseEngine_m(engine.baseEngine_m)
  { }

  //---------------------------------------------------------------------------
  /// The constructors take a MultiPatchView and a slice domain like
  /// a SliceInterval<OrigDim,Dim> or a SliceRange<OrigDim,Dim>.

  template<int OrigDim, class DT>
  Engine(const Engine<OrigDim, T, Tag_t> &engine, 
    const SliceDomain<DT> &domain)
  : layout_m(engine.layout(), domain),
    baseEngine_m(engine.baseEngine())
  { }

  //---------------------------------------------------------------------------
  /// Copy constructor.

  Engine(const Engine_t &model);


  //===========================================================================
  // Destructor
  //===========================================================================

  //---------------------------------------------------------------------------
  /// non - Trivial.

  ~Engine();


  //===========================================================================
  /// Assignment operators
  //===========================================================================

  Engine_t &operator=(const Engine_t &model);


  //===========================================================================
  // Accessor and mutator functions
  //===========================================================================

  //---------------------------------------------------------------------------
  ///@name Element access via Loc.
  //@{

  Element_t read(const Loc<Dim> &) const;
  ElementRef_t operator()(const Loc<Dim> &) const;

  //@}

  //---------------------------------------------------------------------------
  ///@name Element access via ints for speed.
  //@{

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

  //@}

  //---------------------------------------------------------------------------
  /// Return a patch given a Node or INode. This is different than in
  /// MultiPatch because we really need to take a view here.

  typename NewEngine<This_t, INode<Dim> >::Type_t
  globalPatch(const INode<Dim> &inode) const
  {
    typedef SubDomainTraits<Dim2, Dim>                  SubDomainTraits_t;
    typedef typename SubDomainTraits_t::LocalToBase_t   LocalToBase_t; 

    int gid = inode.globalID(layout_m.ID());

    PatchEngine_t &pengine = data()[gid];

    // Now that we've got the patch's engine, we need to transform to base
    // coordinates, because that is the coordinate system for indexing the
    // engine.
    
    LocalToBase_t bdom = Pooma::NoInit();
    layout_m.localToBase(inode.domain(), bdom);
    
    // Make sure our initial assumption that we're fully contained
    // in the patch is correct.

    PAssert(contains(pengine.domain(), SubDomainTraits_t::totalDomain(bdom)));

    // OK, we're going to be making a new engine here.
    
    typedef typename NewEngine<This_t,INode<Dim> >::Type_t Ret_t;
    
    return Ret_t(pengine, bdom);
  }

  template<class Domain>
  typename NewEngine<This_t, Node<Domain> >::Type_t
  globalPatch(const Node<Domain> &node) const
  {
    typedef SubDomainTraits<Dim2, Dim>                  SubDomainTraits_t;
    typedef typename SubDomainTraits_t::LocalToBase_t   LocalToBase_t; 

    // First we need to transform the node's domain to the base coordinates.
    
    LocalToBase_t bdom = Pooma::NoInit();
    layout_m.localToBase(node.domain(), bdom);

    int gid = node.globalID();  

    PatchEngine_t &pengine = data()[gid];
    
    // Make sure our initial assumption that we're fully contained
    // in the patch is correct.

    PAssert(contains(pengine.domain(), SubDomainTraits_t::totalDomain(bdom)));
    
    // OK, we're going to be making a new engine here.
    
    typedef typename NewEngine<This_t,Node<Domain> >::Type_t Ret_t;
    
    return Ret_t(pengine, bdom);
  }

  //---------------------------------------------------------------------------
  /// Return the layout.
  
  inline Layout_t &layout()
  {
    return layout_m;
  }

  inline const Layout_t &layout() const
  {
    return layout_m;
  }

  //---------------------------------------------------------------------------
  /// Return the domain and base domain

  inline const Domain_t &domain() const
  {
    return layout().domain();
  }

  inline const Domain_t &innerDomain() const
  {
    return layout().innerDomain();
  }

  //---------------------------------------------------------------------------
  /// Return the first index in the given direction (always zero).

  inline int first(int) const
  {
    return 0;
  }

  //------------------------------------------------------------
  /// Return whether we have been initialized yet or not.

  inline bool initialized() const
  {
    return layout().initialized();
  }

  //===========================================================================
  // Guard layer support
  //===========================================================================
  
  // These operations are simply forwarded to the underlying viewed engine.
  
  //---------------------------------------------------------------------------
  /// Fill the internal guard cells.
  
  inline void fillGuards() const
  {
    baseEngine_m.fillGuards();
  }
  
  inline void fillGuards(const GuardLayers<Dim2>& g) const
  {
    baseEngine_m.fillGuards(g);
  }
  
  //---------------------------------------------------------------------------
  /// Set the internal guard cells to a particular value (default zero)
  
  inline void setGuards(const T &val) const
  {
    baseEngine_m.setGuards(val);
  }

  //---------------------------------------------------------------------------
  /// Accumulate from the internal guards into 

  inline void accumulateFromGuards() const
  {
    baseEngine_m.accumulateFromGuards();
  }

  //---------------------------------------------------------------------------
  /// Set and get the dirty flag (fillGuard is a no-op unless the 
  /// dirty flag is true).
  
  inline void setDirty() const
  {
    baseEngine_m.setDirty();
  }

  inline void clearDirty(int face=-1) const
  {
    baseEngine_m.clearDirty(face);
  }
  
  inline bool isDirty(int face=-1) const
  {
    return baseEngine_m.isDirty(face);
  }

  //---------------------------------------------------------------------------
  /// Return the data.

  inline const PatchContainer_t &data() const
  {
    return baseEngine_m.data();
  }

  //---------------------------------------------------------------------------
  /// Return whether or not this engine has the same controller
  /// as another.  We find this by finding if the beginning
  /// of the block for each of them is the same.

  template<int D1, class T1>
  inline bool sameController(const Engine<D1, T1, Tag_t> &e) const
  {
    return block() == e.block();
  }

  //---------------------------------------------------------------------------
  /// Return a (shallow) copy of the data block.
  ///
  /// WARNING: If a copy of the data block exists at the time that
  /// all views of the engine go away, the dirty flag will not be 
  /// deleted. So using this should be avoided, if possible.
  
  PatchContainer_t block() const
  {
    return baseEngine_m.data();
  }
          
  //---------------------------------------------------------------------------
  /// Return the base engine to be used by constructors.
  
  inline const ViewedEngine_t &baseEngine() const
  {
    return baseEngine_m;
  }

private:

  //===========================================================================
  // Data
  //===========================================================================

  /// Layout for this engine.

  Layout_t layout_m;

  /// Shallow copy of the underyling engine.
  /// We have to have this to support filling guard cells. There
  /// is only a single dirty flag, so guard cell fill requests must
  /// fill the guards for the entire engine, and we need a reference
  /// to that engine in order to carry this out.
  
  ViewedEngine_t baseEngine_m;

};


///////////////////////////////////////////////////////////////////////////////
//
// Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> > MEMBER FUNCTIONS
//
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
// Indexing operators.
//-----------------------------------------------------------------------------

// Read methods:

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(const Loc<Dim> &loc) const
{
  return data()[layout_m.globalID(loc)].read(loc);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0) const
{
  return data()[layout_m.globalID(i0)].read(i0);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0, int i1) const
{
  return data()[layout_m.globalID(i0, i1)].read(i0, i1);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0, int i1, int i2)
const
{
  return data()[layout_m.globalID(i0, i1, i2)].read(i0, i1, i2);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0, int i1, int i2,
  int i3) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3)].read(i0, i1, i2, i3);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0, int i1, int i2,
  int i3, int i4) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3, i4)].read
    (i0, i1, i2, i3, i4);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0, int i1, int i2,
  int i3, int i4, int i5) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3, i4, i5)].read
    (i0, i1, i2, i3, i4, i5);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::Element_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::read(int i0, int i1, int i2,
  int i3, int i4, int i5, int i6) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3, i4, i5, i6)].read
    (i0, i1, i2, i3, i4, i5, i6);
}

// Operator()'s:

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(const Loc<Dim> &loc)
  const
{
  return data()[layout_m.globalID(loc)](loc);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0) const
{
  return data()[layout_m.globalID(i0)](i0);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0, 
  int i1) const
{
  return data()[layout_m.globalID(i0, i1)](i0, i1);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0, int i1, 
  int i2) const
{
  return data()[layout_m.globalID(i0, i1, i2)](i0, i1, i2);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0, int i1, 
  int i2, int i3) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3)](i0, i1, i2, i3);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0, int i1,
  int i2, int i3, int i4) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3, i4)]
    (i0, i1, i2, i3, i4);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0, int i1, 
  int i2, int i3, int i4, int i5) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3, i4, i5)]
    (i0, i1, i2, i3, i4, i5);
}

template<int Dim, class T, class LayoutTag, class PatchTag>
inline typename Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::ElementRef_t
Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::operator()(int i0, int i1, 
  int i2, int i3, int i4, int i5, int i6) const
{
  return data()[layout_m.globalID(i0, i1, i2, i3, i4, i5, i6)]
    (i0, i1, i2, i3, i4, i5, i6);
}

// Dynamic event handler for non-dynamic patch engines:

template <int Dim, class T, class LayoutTag, class PatchTag>
inline void Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >::
dynamicHandler(Observable_t &, const ObserverEvent &, 
               const WrappedInt<false> &)
{ 
  PInsist(0,"This patch engine does not support dynamic events!");
}


///////////////////////////////////////////////////////////////////////////////
//
// Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag,Dim2> > MEMBER FUNCTIONS
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Copy constructor. This is pretty much a shallow copy since we use 
// our members' copy constructors to do all of the work. 
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
Engine(const Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > 
  &modelEngine)
: layout_m(modelEngine.layout_m), 
  baseEngine_m(modelEngine.baseEngine_m)
{ }


//-----------------------------------------------------------------------------
// Assignment operator. This is a shallow copy since we are just copying our
// smart pointer and layout.
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > &
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator=(const Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> > 
  &rhs)
{
  // Check for self-assignment.
  
  if (&rhs == this) return *this;
  
  // Assign our data and our layout.

  layout_m     = rhs.layout_m;
  baseEngine_m = rhs.baseEngine_m;
  return *this;  
}


//-----------------------------------------------------------------------------
//
// Destructor
// Trivial. Member destructors should handle everything.
//
//-----------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
Engine<Dim, T, MultiPatchView<LayoutTag,PatchTag,Dim2> >::
~Engine()
{ }

//-----------------------------------------------------------------------------
// Indexing operators.
//-----------------------------------------------------------------------------

// Read methods:

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline 
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(const Loc<Dim> &loc) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(loc, gloc);
  return data()[o].read(gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::read(int i0) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, gloc);
  return data()[o].read(gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(int i0, int i1) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, gloc);
  return data()[o].read(gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(int i0, int i1, int i2) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, gloc);
  return data()[o].read(gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(int i0, int i1, int i2, int i3) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, gloc);
  return data()[o].read(gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(int i0, int i1, int i2, int i3, int i4) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, i4, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(int i0, int i1, int i2, int i3, int i4, int i5) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, i4, i5, gloc);
  return data()[o].read(gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline
typename Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::Element_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
read(int i0, int i1, int i2, int i3, int i4, int i5, int i6) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, i4, i5, i6, gloc);
  return data()[o].read(gloc);
}

// Operator()'s:

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(const Loc<Dim> &loc) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(loc, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0, int i1) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0, int i1, int i2) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0, int i1, int i2, int i3) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0, int i1, int i2, int i3, int i4) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, i4, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0, int i1, int i2, int i3, int i4, int i5) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, i4, i5, gloc);
  return data()[o](gloc);
}

template<int Dim, class T, class LayoutTag, class PatchTag, int Dim2>
inline typename
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::ElementRef_t
Engine<Dim, T, MultiPatchView<LayoutTag, PatchTag, Dim2> >::
operator()(int i0, int i1, int i2, int i3, int i4, int i5, int i6) const
{
  Loc<Dim2> gloc = Pooma::NoInit();
  int o = layout_m.globalID(i0, i1, i2, i3, i4, i5, i6, gloc);
  return data()[o](gloc);
}

//---------------------------------------------------------------------------
// Specialization of IntersectEngine because these engines contain multiple
// patches.
// Respond to the IntersecEngineTag<Dim> message by intersecting our layout
// with the enclosed intersector.
//---------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag, class Intersect>
struct LeafFunctor<Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >,
  ExpressionApply<IntersectorTag<Intersect> > >
{
  typedef int Type_t;

  static Type_t
  apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &engine,
	const ExpressionApply<IntersectorTag<Intersect> > &tag)
  {
    GuardLayers<Dim> usedGuards;
    bool useGuards =
      tag.tag().intersector_m.intersect(engine,
				  engine.layout().internalGuards(), usedGuards);

    if (useGuards)
      engine.fillGuards(usedGuards);

    return 0;
  }
};

template <int Dim, class T, class LT, class PatchTag, int BD,
  class Intersect>
struct LeafFunctor<Engine<Dim, T, MultiPatchView<LT,PatchTag,BD> >,
  ExpressionApply<IntersectorTag<Intersect> > >
{
  typedef int Type_t;

  static Type_t
  apply(const Engine<Dim,T,MultiPatchView<LT,PatchTag,BD> > &engine,
	const ExpressionApply<IntersectorTag<Intersect> > &tag)
  {
    typedef typename MultiPatchLayoutTraits<LT,Dim>::Layout_t Layout_t;
    return applyHandler(engine, tag, WrappedInt<Layout_t::supportsGuards>());
  }

  inline static Type_t
  applyHandler(const Engine<Dim,T,MultiPatchView<LT,PatchTag,BD> > &engine,
	       const ExpressionApply<IntersectorTag<Intersect> > &tag,
	       const WrappedInt<true> &)
  {
    GuardLayers<BD> usedGuards;
    bool useGuards =
      tag.tag().intersector_m.
      intersect(engine,
		engine.layout().baseLayout().internalGuards(), usedGuards);

    if (useGuards)
      engine.fillGuards(usedGuards);

    return 0;
  }

  inline static Type_t
  applyHandler(const Engine<Dim,T,MultiPatchView<LT,PatchTag,BD> > &engine,
	       const ExpressionApply<IntersectorTag<Intersect> > &tag,
	       const WrappedInt<false> &)
  {
    tag.tag().intersector_m.intersect(engine); 
    return 0;
  }
};

//---------------------------------------------------------------------------
// Specialization of EngineFunctor for EnginePatch.
// EnginePatch(i) is a generic external method for obtaining the ith
// patch from an engine (even non-multipatch engines which are defined to
// have one patch).
//
// This functor is not specialized for MultiPatch engine views,
// because such an operation would be slightly strange.  (Do you want
// the view of the ith original patch or a view of the ith viewed patch
// or do you want the entire ith patch from the viewed engine?)
//---------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
struct EngineFunctor<Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >,
  EnginePatch >
{
  typedef Engine<Dim, T, MultiPatch<LayoutTag, PatchTag> > Subject_t;

  // typedef typename Subject_t::PatchView_t Type_t;
  
  typedef Engine<Dim,T,PatchTag>               Type_t;

  static inline
  Type_t apply(const Subject_t &engine, const EnginePatch &tag)
  {
    return engine.localPatch(tag.patch_m);
  }
};

//---------------------------------------------------------------------------
// Specialization of EngineFunctor for EngineNumPatches.
// Returns the number of patches in the engine.
//---------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
struct EngineFunctor<Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >,
  EngineNumPatches >
{
  typedef int Type_t;

  static inline
  Type_t apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &engine,
	       const EngineNumPatches &)
  {
    return engine.layout().sizeLocal();
  }
};

//---------------------------------------------------------------------------
// Tell multipatch engines that they're dirty.
//---------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
struct NotifyEngineWrite<Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> > >
{
  inline static void
  notify(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &engine)
  {
    engine.setDirty();
  }
};

template <int Dim, class T, class LT, class PatchTag, int BD>
struct NotifyEngineWrite<Engine<Dim, T, MultiPatchView<LT,PatchTag,BD> > >
{
  inline static void
  notify(const Engine<Dim,T,MultiPatchView<LT,PatchTag,BD> > &engine)
  {
    engine.setDirty();
  }
};

//---------------------------------------------------------------------------
// localPatchView() is a utility function used to perform operations on
// multipatch engines where the patch engine could be a remote-engine.
// Currently this function is used by the Multi-patch engine copy()
// functions.
//---------------------------------------------------------------------------

template<class PatchEngine>
inline
PatchEngine &localPatchEngine(PatchEngine &e)
{
  return e;
}

// Include .cpp file to get out-of-line functions.

#include "Engine/MultiPatchEngine.cpp"

// }; // namespace Pooma
///////////////////////////////////////////////////////////////////////////////


#endif // POOMA_ENGINE_MULTIPATCHENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MultiPatchEngine.h,v $   $Author: richi $
// $Revision: 1.122 $   $Date: 2004/11/10 21:58:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
