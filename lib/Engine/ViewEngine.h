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

/** @file
 * @ingroup Engine
 * @brief
 * Generalized view engine that can handle intersections for contained
 * multi-patch engines.
 */

#ifndef POOMA_ENGINE_VIEWENGINE_H
#define POOMA_ENGINE_VIEWENGINE_H

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

#include "Utilities/ViewIndexer.h"
#include "Utilities/WrappedInt.h"
#include "Layout/Node.h"
#include "Layout/INode.h"
#include "Layout/DomainLayout.h"
#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Engine/DataObject.h"
#include "Engine/Intersector.h"
#include "Engine/IntersectEngine.h"
#include "Evaluator/EngineTraits.h"
#include "PETE/ErrorType.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<int OriginalDim, class ViewedEngineTag>
struct ViewEngine
{
  ViewEngine() {};
  ~ViewEngine() {};
 };

template<int Dim, class T, int OriginalDim, class ViewedEngineTag>
class Engine<Dim, T, ViewEngine<OriginalDim, ViewedEngineTag> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants

  typedef ViewEngine<OriginalDim, ViewedEngineTag> Tag_t;
  typedef Engine<Dim, T, Tag_t>                    This_t;
  typedef This_t                                   Engine_t;
  typedef Engine<OriginalDim, T, ViewedEngineTag>  ViewedEngine_t;
  typedef ViewIndexer<Dim, OriginalDim>            Indexer_t;
  typedef Interval<Dim>                            Domain_t;
  typedef DomainLayout<Dim>                        Layout_t;
  typedef T                                        Element_t;
  typedef ErrorType                                ElementRef_t;

  enum { dimensions = Dim };
  enum { hasDataObject = ViewedEngine_t::hasDataObject };
  enum { dynamic = false };
  enum { zeroBased = true };
  enum { multiPatch = ViewedEngine_t::multiPatch };

  //---------------------------------------------------------------------------
  // Default constructor allows engines to be used in containers.

  Engine()
    : eng_m()
  {
  }

  //---------------------------------------------------------------------------
  // Construct from an existing Engine and various sorts of domains 
  // (e.g., take a view).

  template<class DT>
  Engine(const Engine<Dim, T, ViewedEngineTag> &e, const Domain<Dim, DT> &dom)
    : eng_m(e), indexer_m(dom)
  {
    // We logically cannot be a slice.
    
    CTAssert(OriginalDim == Dim);
  }

  template<class DT>
  Engine(const Engine<OriginalDim, T, ViewedEngineTag> &e,
	 const SliceDomain<DT> &dom)
    : eng_m(e), indexer_m(dom)
  {
    // The domain's dimension should match ours.
    
    CTAssert(DT::sliceDimensions == Dim);
    CTAssert(DT::dimensions == OriginalDim);
  }

  template<class Domain>
  Engine(const Engine<Dim, T, ViewedEngineTag> &e, const Node<Domain> &node)
    : eng_m(e), indexer_m(node.domain())
  {
    // The nodes's dimension should match ours.
    
    CTAssert(Domain::dimensions == Dim);
  }

  Engine(const Engine<Dim, T, ViewedEngineTag> &e, const INode<Dim> &inode)
    : eng_m(e), indexer_m(inode.domain())
  { }

  //---------------------------------------------------------------------------
  // Construct from an existing ViewEngine and various sorts of domains 
  // (e.g., take a view).

  template<class DT>
  Engine(const Engine<Dim, T, ViewEngine<OriginalDim, ViewedEngineTag> > &e, 
	 const Domain<Dim, DT> &dom)
    : eng_m(e.viewedEngine()), indexer_m(e.indexer(), dom)
  {
  }

  template<int OrigDim, class DT>
  Engine(const Engine<OrigDim, T, ViewEngine<OriginalDim, ViewedEngineTag> > &e, 
	 const SliceDomain<DT> &dom)
    : eng_m(e.viewedEngine()), indexer_m(e,indexer(), dom)
  {
    // The domain's dimension should match ours.
    
    CTAssert(DT::sliceDimensions == Dim);
    CTAssert(DT::dimensions == OrigDim);
  }

  template<class Domain>
  Engine(const Engine<Dim, T, ViewEngine<OriginalDim, ViewedEngineTag> > &e, 
	 const Node<Domain> &node)
    : eng_m(e.viewedEngine()), indexer_m(e.indexer(), node.domain())
  {
    // The nodes's dimension should match ours.
    
    CTAssert(Domain::dimensions == Dim);
  }

  Engine(const Engine<Dim, T, ViewEngine<OriginalDim, ViewedEngineTag> > &e, 
	 const INode<Dim> &inode)
    : eng_m(e.viewedEngine()), indexer_m(e.indexer(), inode.domain())
  { }

  //---------------------------------------------------------------------------
  // Construct from another view-engine.

  Engine(const Engine<Dim, T, ViewEngine<OriginalDim, ViewedEngineTag> >
	 &model)
    : eng_m(model.viewedEngine()), indexer_m(model.indexer())
  { }

  //---------------------------------------------------------------------------
  // Assign one view-engine to another.

  This_t &operator=(const This_t &rhs)
  {
    eng_m = rhs.viewedEngine();
    indexer_m = rhs.indexer();
    
    return *this;
  }
    
  //---------------------------------------------------------------------------
  // Element access via ints for speed.

  inline Element_t read(int i0) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(int i0, int i1) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, i1, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(int i0, int i1, int i2) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, i1, i2, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(int i0, int i1, int i2, int i3) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, i1, i2, i3, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(int i0, int i1, int i2, int i3, int i4) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, i1, i2, i3, i4, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(int i0, int i1, int i2, int i3, int i4, int i5) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, i1, i2, i3, i4, i5, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(int i0, int i1, int i2, int i3, int i4, 
			int i5, int i6) const 
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(i0, i1, i2, i3, i4, i5, i6, oloc);
    return eng_m.read(oloc);
  }
  inline Element_t read(const Loc<Dim> &loc) const
  {
    Loc<OriginalDim> oloc;
    indexer_m.translate(loc, oloc);
    return eng_m.read(oloc);
  }

  //---------------------------------------------------------------------------
  // Return the domain.

  inline const Domain_t &domain() const { return indexer_m.domain(); }

  //---------------------------------------------------------------------------
  // Return the layout.

  inline Layout_t layout() const { return Layout_t(domain()); }

  //---------------------------------------------------------------------------
  // Return the first value for the specified direction (always zero since this
  // engine is zero-based).
  
  inline int first(int i) const
  {
    PAssert(i >= 0 && i < Dim);
    return 0;
  }

  //---------------------------------------------------------------------------
  // Accessors.

  inline const ViewedEngine_t &viewedEngine() const { return eng_m; }
  inline const Indexer_t &indexer() const { return indexer_m; }

  //---------------------------------------------------------------------------
  // Need to pass lock requests to the contained engine.

  template<class RequestType>
  inline
  typename DataObjectRequest<RequestType>::Type_t
  dataObjectRequest(const DataObjectRequest<RequestType>& f) const
  {
    return eng_m.dataObjectRequest(f);
  }

private:

  ViewedEngine_t eng_m;
  Indexer_t indexer_m;
};

/**
 * Specializations of NewEngine for subsetting a view-engines with
 * an arbitrary domain. 
 */

template <int Dim, class T, int D2, class ViewedTag>
struct NewEngine< Engine<Dim,T,ViewEngine<D2,ViewedTag> >, Interval<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<D2, ViewedTag> > Type_t;
};

template <int Dim, class T, int D2, class ViewedTag>
struct NewEngine< Engine<Dim,T,ViewEngine<D2,ViewedTag> >, Range<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<D2, ViewedTag> > Type_t;
};

template <int Dim, class T, int D2, class ViewedTag, int SliceDim>
struct NewEngine< Engine<Dim, T, ViewEngine<D2, ViewedTag> >,
  SliceInterval<Dim, SliceDim> >
{
  typedef Engine<SliceDim, T, ViewEngine<D2, ViewedTag> > Type_t;
};

template <int Dim, class T, int D2, class ViewedTag, int SliceDim>
struct NewEngine< Engine<Dim, T, ViewEngine<D2, ViewedTag> >,
  SliceRange<Dim, SliceDim> >
{
  typedef Engine<SliceDim, T, ViewEngine<D2, ViewedTag> > Type_t;
};

template <int Dim, class T, int D2, class ViewedTag>
struct NewEngine< Engine<Dim, T, ViewEngine<D2, ViewedTag> >, INode<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<D2, ViewedTag> > Engine_t;
  typedef typename Engine_t::ViewedEngine_t ViewedEngine_t;
  typedef typename NewEngine<ViewedEngine_t,
    INode<D2> >::Type_t NewViewedEngine_t;
  typedef typename NewEngine<NewViewedEngine_t,
    SliceRange<D2, Dim> >::Type_t Type_t;
};

template <int Dim, class T, int D2, class ViewedTag>
struct NewEngineEngine< Engine<Dim, T, ViewEngine<D2, ViewedTag> >,
  INode<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<D2, ViewedTag> > Engine_t;  
  typedef typename NewEngine< Engine_t, INode<Dim> >::NewViewedEngine_t Type_t;

  static inline Type_t
  apply(const Engine_t &e, const INode<Dim> &inode)
  {
    Range<D2> base;
    e.indexer().localToBase(inode.domain(), base);

    Interval<D2> baseInt;
    int i;
    for (i = 0; i < D2; ++i)
    {
      baseInt[i] = Interval<1>(base[i].first(), base[i].last());
    }

    INode<D2> viewNode(inode, baseInt);
    return Type_t(e.viewedEngine(), viewNode);
  }
};

template <int Dim, class T, int D2, class ViewedTag>
struct NewEngineDomain< Engine<Dim, T, ViewEngine<D2, ViewedTag> >,
  INode<Dim> >
{
  typedef SliceRange<D2, Dim> Type_t;
  typedef Engine<Dim, T, ViewEngine<D2, ViewedTag> > Engine_t;  

  static inline Type_t
  apply(const Engine_t &e, const INode<Dim> &inode)
  {
    SliceRange<D2, Dim> base;
    e.indexer().localToBase(inode.domain(), base);

    base.totalDomain() -= base.totalDomain().firsts();
    base.setSliceFromTotal();

    return base;
  }
};

template <int Dim, class T, class ViewedTag>
struct NewEngine<Engine<Dim, T, ViewEngine<Dim, ViewedTag> >, INode<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<Dim, ViewedTag> > Engine_t;
  typedef typename Engine_t::ViewedEngine_t ViewedEngine_t;
  typedef typename NewEngine<ViewedEngine_t,
    INode<Dim> >::Type_t NewViewedEngine_t;
  typedef typename NewEngine<NewViewedEngine_t,
    Range<Dim> >::Type_t Type_t;
};

template <int Dim, class T, class ViewedTag>
struct NewEngineEngine<Engine<Dim, T, ViewEngine<Dim, ViewedTag> >,
  INode<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<Dim, ViewedTag> > Engine_t;  
  typedef typename NewEngine<Engine_t, INode<Dim> >::NewViewedEngine_t Type_t;

  static inline Type_t
  apply(const Engine_t &e, const INode<Dim> &inode)
  {
    Range<Dim> base;
    e.indexer().localToBase(inode.domain(), base);

    Interval<Dim> baseInt;
    int i;
    for (i = 0; i < Dim; ++i)
    {
      baseInt[i] = Interval<1>(base[i].first(), base[i].last());
    }

    INode<Dim> viewNode(inode, baseInt);
    return Type_t(e.viewedEngine(), viewNode);
  }
};

template <int Dim, class T, class ViewedTag>
struct NewEngineDomain< Engine<Dim, T, ViewEngine<Dim, ViewedTag> >,
  INode<Dim> >
{
  typedef Range<Dim> Type_t;
  typedef Engine<Dim, T, ViewEngine<Dim, ViewedTag> > Engine_t;  

  static inline Type_t
  apply(const Engine_t &e, const INode<Dim> &inode)
  {
    Range<Dim> base;
    e.indexer().localToBase(inode.domain(), base);
    base -= base.firsts();

    return base;
  }
};

template<int OriginalDim, class ViewedEngineTag>
struct EvaluatorEngineTraits<ViewEngine<OriginalDim, ViewedEngineTag> >
{
  typedef typename EvaluatorEngineTraits<ViewedEngineTag>::Evaluator_t
  Evaluator_t;
};


template<int Dim, int ViewD1, int ViewD2>
class ViewIntersector
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef IntersectorData<Dim>                          IntersectorData_t;
  typedef ViewIntersector<Dim, ViewD1, ViewD2>          This_t;
  typedef typename IntersectorData_t::IDContainer_t     IDContainer_t;
  typedef typename IntersectorData_t::BaseDomain_t      BaseDomain_t;
  typedef typename IntersectorData_t::BaseDomainContainer_t
                                                        BaseDomainContainer_t;
  typedef typename IntersectorData_t::INode_t           INode_t;
  typedef typename IntersectorData_t::INodeContainer_t  INodeContainer_t;
  typedef typename IntersectorData_t::const_iterator    const_iterator;
  typedef RefCountedPtr<IntersectorData_t>              DataPtr_t;
  
  typedef ViewIndexer<ViewD1, ViewD2>                   Indexer_t;

  enum { dimensions = Dim };
  
  ViewIntersector(const ViewIndexer<ViewD1, ViewD2> &indexer,
		  const Intersector<Dim> &model)
    : pdata_m(model.data()), indexer_m(indexer)
  {
    // We haven't yet implemented the case where the view doesn't have
    // the same dimensions as the original expression
    CTAssert(Dim == ViewD1);
  }

  This_t &operator=(const This_t &model)
  {
    if (this != &model)
    {
      indexer_m = model.indexer_m;
      pdata_m = model.pdata_m;
    }
    return *this;
  }

  ~ViewIntersector() { }

  inline DataPtr_t &data() { return pdata_m; }
  inline const DataPtr_t &data() const { return pdata_m; }

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
  void intersect(const Engine &e) 
  {
    int n = data()->ids_m.size();

    typedef INode<ViewD2> INode2_t;
    typedef std::vector<INode2_t> INode2Container_t;

    INode2Container_t inodes2;

    int id = e.layout().ID();

    if (n == 0)
    {
      // If no intersections have been performed yet, then we need to
      // intersect the layout with the baseDomain of the indexer.
      // (We may only be viewing a portion of the engine with the
      // view engine.)

      Range<ViewD2> base = indexer_m.baseDomain();

      e.layout().touches(base,
			 std::back_inserter(inodes2),
			 TouchesConstructINode<ViewD2>(id, GlobalIDDataBase::
						       nullNodeKey(),
						       &(data()->gidStore_m))
			 );

      int i;
      int ni = inodes2.size();
      for (i = 0; i < ni; i++)
      {
	Interval<Dim> ival1;
	indexer_m.baseToLocalInterval(inodes2[i].domain(), ival1);
	INode<Dim> inode(inodes2[i], ival1);
	data()->inodes_m.push_back(inode);
      }
    }
    else
    {
      int ni = data()->inodes_m.size();
      int i;
      for (i = 0; i < ni; i++)
      {
	Range<ViewD2> range;
	indexer_m.localToBase(data()->inodes_m[i].domain(), range);

	e.layout().touches(range,
			   std::back_inserter(inodes2),
			   INode<ViewD2>::touchesConstructINode(id,
					    data()->inodes_m[i])
			   );
      }
        
      data()->inodes_m.erase(data()->inodes_m.begin(),
			     data()->inodes_m.begin() + ni);

      ni = inodes2.size();
      for (i = 0; i < ni; i++)
      {
	Interval<Dim> ival1;
	indexer_m.baseToLocalInterval(inodes2[i].domain(), ival1);
	INode<Dim> inode(inodes2[i], ival1);
	data()->inodes_m.push_back(inode);
      }
    }
  }

  template<class Engine, int Dim2>
  inline
  bool intersect(const Engine &l, const GuardLayers<Dim2> &guard) 
  {
    return (data()->intersect(l,guard));
  }

private:
  DataPtr_t pdata_m;
  Indexer_t indexer_m;

};

template <int Dim, class T, int D2, class ViewedTag, class Intersect>
struct LeafFunctor< Engine<Dim, T, ViewEngine<D2, ViewedTag> >,
  ExpressionApply<IntersectorTag<Intersect> > >
{
  typedef int Type_t;
  typedef Engine<Dim, T, ViewEngine<D2, ViewedTag> > Engine_t;

  static Type_t
  apply(const Engine_t &,
	const ExpressionApply<IntersectorTag<Intersect> > &,
	const WrappedInt<false> &)
  {
    return 0;
  }

  static Type_t
  apply(const Engine_t &engine,
	const ExpressionApply<IntersectorTag<Intersect> > &tag,
	const WrappedInt<true> &)
  {
    enum { d1 = Intersect::dimensions };
    ViewIntersector<d1, Dim, D2> newIntersector(engine.indexer(),
						tag.tag().intersector_m);
    ExpressionApply<IntersectorTag<ViewIntersector<d1, Dim, D2> > >
      newTag(newIntersector);
    
    forEach(engine.viewedEngine(), newTag, NullCombine());
    return 0;
  }

  static Type_t
  apply(const Engine_t &engine,
	const ExpressionApply<IntersectorTag<Intersect> > &tag)
  {
    enum { multiPatch =
	   Engine<Dim, T, ViewEngine<D2, ViewedTag> >::multiPatch };

    return apply(engine, tag, WrappedInt<multiPatch>());
  }
};

//---------------------------------------------------------------------------
// Specialization of  DataObjectRequest engineFunctor to pass the request to
// the contained engine.
//---------------------------------------------------------------------------

template<class RequestType> class DataObjectRequest;

template <int Dim, class T, int D2, class ViewedTag, class RequestType>
struct EngineFunctor< Engine<Dim, T, ViewEngine<D2, ViewedTag> >,
  DataObjectRequest<RequestType> >
{
  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;

  static Type_t
  apply(const Engine<Dim, T, ViewEngine<D2, ViewedTag> > &engine,
	const DataObjectRequest<RequestType> &tag)
  {
    return engineFunctor(engine.viewedEngine(), tag);
  }
};

template <int D, class T, int D2, class E, class Tag>
struct LeafFunctor<Engine<D, T, ViewEngine<D2, E> >, ExpressionApply<Tag> >
{
  typedef Engine<D, T, ViewEngine<D2, E> > Subject_t;
  typedef typename Subject_t::ViewedEngine_t Engine_t;
  typedef LeafFunctor<Engine_t, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  static
  Type_t apply(const Subject_t &engine, const ExpressionApply<Tag> &tag)
  {
    return LeafFunctor_t::apply(engine.viewedEngine(), tag);
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_VIEWENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ViewEngine.h,v $   $Author: richard $
// $Revision: 1.28 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
