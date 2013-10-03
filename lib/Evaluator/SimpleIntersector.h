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
// SimpleIntersector
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Evaluator
 * @brief
 * Intersector that assumes matching layouts.
 */

#ifndef POOMA_EVALUATOR_SIMPLEINTERSECTOR_H
#define POOMA_EVALUATOR_SIMPLEINTERSECTOR_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/EngineFunctor.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<int Dim>
class SimpleIntersectorData
  : public RefCounted
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef SimpleIntersectorData<Dim>                    This_t;
  typedef INode<Dim>                                    INode_t;
  typedef std::vector<INode_t>                          INodeContainer_t;
  typedef typename INodeContainer_t::const_iterator     const_iterator;
  typedef Unique::Value_t                               LayoutID_t;
  
  enum { dimensions = Dim };
  
  
  //===========================================================================
  // Constructors
  //===========================================================================

  // Default constructor is trival.
  
  inline SimpleIntersectorData(const Interval<Dim> &domain, const GuardLayers<Dim> &extent)
    : seenFirst_m(false), domain_m(domain), extent_m(extent)
  {
  }

  //===========================================================================
  // Destructor
  //===========================================================================

  // Trival since members manage their own data.
  
  inline ~SimpleIntersectorData() { }
  
  template<class Engine>
  void intersect(const Engine &engine, bool useGuards) 
  {
    typedef typename Engine::Layout_t Layout_t;
    typedef typename NewEngine<Engine, Interval<Dim> >::Type_t NewEngine_t;
    const Layout_t &layout(engine.layout());

    // add an assertion that all layouts have the same base (probably
    // the same baseDomain)

    if (!seenFirst_m)
    {
      firstID_m = layout.ID();
      seenFirst_m = true;

      int key = GlobalIDDataBase::nullNodeKey();
      layout.touches(domain_m, std::back_inserter(inodes_m),
		     TouchesConstructINode<Dim>(firstID_m, key, &gidStore_m));
    }
    else
    {
      shared(layout.ID(), firstID_m);
    }
    // We need to process possible expression engines with different
    // guard needs here.  Modeled after StencilIntersector.
    if (useGuards) {
      expressionApply(NewEngine_t(engine, grow(domain_m, extent_m)),
		      IntersectorTag<Intersector<Dim> >(lhsi_m));
    } else {
      expressionApply(NewEngine_t(engine, domain_m),
		      IntersectorTag<Intersector<Dim> >(lhsi_m));
    }
  }

  inline
  void shared(LayoutID_t id1, LayoutID_t id2)
  {
    gidStore_m.shared(id1,id2);
  }

  //private:  

  // Don't ever want to copy one of these.
  SimpleIntersectorData(const This_t &);
  This_t &operator=(const This_t &);
  
  //===========================================================================
  // Data members
  //===========================================================================

  LayoutID_t firstID_m;
  bool seenFirst_m;
  INodeContainer_t inodes_m;
  GlobalIDDataBase gidStore_m;
  Interval<Dim> domain_m;
  GuardLayers<Dim> extent_m;
  Intersector<Dim> lhsi_m;
};

/**
 * This intersector handles matching layouts only.  It also assumes you
 * know in advance the amount of guards used.  But it allows differentiating
 * between engines that use or do not use guards.
 *
 * It doesnt intersect individual layouts but is done with creating INodes
 * from the first layout it sees by intersecting with the domain.
 * Use with care as no checks are present that assure layouts do really
 * match. So this is simple and cheap, but not safe.
 */

template<int Dim>
class SimpleIntersector
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef SimpleIntersectorData<Dim>                  SimpleIntersectorData_t;
  typedef SimpleIntersector<Dim>                              This_t;
  typedef typename SimpleIntersectorData_t::INode_t           INode_t;
  typedef typename SimpleIntersectorData_t::INodeContainer_t  INodeContainer_t;
  typedef typename SimpleIntersectorData_t::const_iterator    const_iterator;
  typedef RefCountedPtr<SimpleIntersectorData_t>              DataPtr_t;
  typedef NullCombine                                         Combine_t;
  
  enum { dimensions = Dim };
  
  SimpleIntersector(const Interval<Dim> &domain, const GuardLayers<Dim> &extent)
    : pdata_m(new SimpleIntersectorData_t(domain, extent)), useGuards_m(true)
  { }

  SimpleIntersector(const This_t &model)
    : pdata_m(model.pdata_m), useGuards_m(model.useGuards())
  { }

  This_t &operator=(const This_t &model)
  {
    if (this != &model) {
      pdata_m = model.pdata_m;
      useGuards_m = model.useGuards_m;
    }
    return *this;
  }

  ~SimpleIntersector() { }

  inline DataPtr_t &data() { return pdata_m; }
  inline const DataPtr_t &data() const { return pdata_m; }

  //===========================================================================
  // Accessors
  //===========================================================================

  // STL iterator support.
  
  inline const_iterator begin() const { return data()->inodes_m.begin(); }
  
  inline const_iterator end() const { return data()->inodes_m.end(); }

  inline int size() const { return data()->inodes_m.size(); }

  //===========================================================================
  // Intersect routines
  //===========================================================================

  // All domains.
  
  template<class Engine>
  inline
  void intersect(const Engine &l) const
  {
    data()->intersect(l, useGuards());

  }

  inline
  bool useGuards() const
  {
    return useGuards_m;
  }

  inline
  void useGuards(bool f) const
  {
    useGuards_m = f;
  }

  // Interface to be used by applyMultiArg()

  template<class A>
  void operator()(const A &a, bool f) const
  {
    useGuards(f);
    expressionApply(a, *this);
  }

private:

  DataPtr_t pdata_m;

  mutable bool      useGuards_m;
};


//-----------------------------------------------------------------------------
// The default behaviour for IntersectEngine is to simply return true.
// We assert that the engine is not multi-patch.
//-----------------------------------------------------------------------------

template<class Eng, int Dim>
struct DefaultExpressionApply<Eng, SimpleIntersector<Dim> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const Eng &,
	       const ExpressionApply<SimpleIntersector<Dim> > &)
  {
    // Engines that are multipatch must specialize this functor
    // to perform the correct intersection.
    CTAssert(!(Eng::multiPatch));
    return true;
  }
};

template<class LT, class PT>
struct MultiPatch;

//---------------------------------------------------------------------------
// Specialization of IntersectEngine because these engines contain multiple
// patches.
// Respond to the IntersecEngineTag<Dim> message by intersecting our layout
// with the enclosed intersector.
//---------------------------------------------------------------------------

template <int Dim, class T, class LayoutTag, class PatchTag>
struct LeafFunctor<Engine<Dim, T, MultiPatch<LayoutTag,PatchTag> >,
  ExpressionApply<SimpleIntersector<Dim> > >
{
  typedef int Type_t;

  static Type_t
  apply(const Engine<Dim,T,MultiPatch<LayoutTag,PatchTag> > &engine,
	const ExpressionApply<SimpleIntersector<Dim> > &apply)
  {
    apply.tag().intersect(engine);

    if (apply.tag().useGuards())
      engine.fillGuards(apply.tag().data()->extent_m);

    return 0;
  }
};

template <int Dim, class T, class LT, class PatchTag, int BD>
struct LeafFunctor<Engine<Dim, T, MultiPatchView<LT,PatchTag,BD> >,
  ExpressionApply<SimpleIntersector<Dim> > >
{
  typedef int Type_t;

  static Type_t
  apply(const Engine<Dim,T,MultiPatchView<LT,PatchTag,BD> > &engine,
	const ExpressionApply<SimpleIntersector<Dim> > &apply)
  {
    apply.tag().intersect(engine);

    if (apply.tag().useGuards())
      engine.fillGuards(apply.tag().data()->extent_m);

    return 0;
  }
};



//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_SIMPLEINTERSECTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SimpleIntersector.h,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
