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
 * @ingroup DiffOps
 * @brief
 * An engine capable of shifting the contents of a field.
 *   - FieldShift: Tag class for defining an engine capable of
 *                 shifting the contents of a field.
 *   - Engine: Specialization for FieldShift
 *   - NewEngine: Specializations for FieldShift
 */

#ifndef POOMA_FIELD_DIFFOPS_FIELDSHIFTENGINE_H
#define POOMA_FIELD_DIFFOPS_FIELDSHIFTENGINE_H

//-----------------------------------------------------------------------------
// Overview:
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Engine/Engine.h"
#include "Layout/INode.h"
#include "Layout/Node.h"
#include "PETE/ErrorType.h"

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template <int Dim>
class DomainLayout;

template<class Functor> class FieldShift;

/**
 * FieldShift is just a tag class for the field-stencil-application
 * engine.
 */

template <class Expression>
struct FieldShift;

/**
 * Engine<Dim, T, FieldShift<Expression> > is a specialization
 * of Engine for FieldShift<Expression>.
 * It makes an offset view of the input expression.
 */

template<int Dim, class T, class Expression>
class Engine<Dim, T, FieldShift<Expression> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants

  typedef FieldShift<Expression>   Tag_t;
  typedef Expression                               Expression_t;
  typedef Engine<Dim, T, Tag_t>                    This_t;
  typedef This_t                                   Engine_t;
  typedef Interval<Dim>                            Domain_t;
  typedef T                                        Element_t;
  typedef ErrorType                                ElementRef_t;
  typedef typename Expression_t::Engine_t          ExprEngine_t;
  typedef DomainLayout<Dim>                        Layout_t;

  enum { dimensions = Dim };
  enum { hasDataObject = ExprEngine_t::hasDataObject };
  enum { dynamic = false };
  enum { zeroBased = false };
  enum { multiPatch = ExprEngine_t::multiPatch };

  //---------------------------------------------------------------------------
  // Construct uninitialized Field stencil objects.  It's an error to use an
  // uninitialized engine, but we need to be able to create uninitialized
  // engines as placeholders to enable deferred initialization of fields.

  Engine()
    : domain_m(Pooma::NoInit()), exprEngine_m()
  {
  }

  //---------------------------------------------------------------------------
  // Construct an empty engine from a layout.

  template<class Layout>
  explicit Engine(const Layout &layout)
    : domain_m(layout.domain()), exprEngine_m()
  {
  }

  //---------------------------------------------------------------------------
  // Construct from a given field and an offset.

  Engine(const Expression_t &f, const Loc<Dim> &offset, Domain_t domain)
    : domain_m(domain),
      offset_m(offset),
      exprEngine_m(f)
  {
  }

  //---------------------------------------------------------------------------
  // Copy constructor.

  Engine(const This_t &model)
    : domain_m(model.domain()),
      offset_m(model.offset_m),
      exprEngine_m(model.exprEngine())
  {
  }    

  //---------------------------------------------------------------------------
  // Shallow assignment.
  
  This_t &operator=(const This_t &model)
  {
    domain_m = model.domain();
    offset_m = model.offset_m;
    exprEngine_m = model.exprEngine();

    return *this;
  }    

  //---------------------------------------------------------------------------
  // Element access via ints for speed.

  inline Element_t read(int i) const 
  {
    return exprEngine()(i + offset_m[0].first());
  }

  inline Element_t read(int i, int j) const 
  {
    return exprEngine()(i + offset_m[0].first(),
                        j + offset_m[1].first());
  }

  inline Element_t read(int i, int j, int k) const 
  {
    return exprEngine()(i + offset_m[0].first(),
                        j + offset_m[1].first(),
                        k + offset_m[2].first());
  }

  inline Element_t read(const Loc<Dim> &loc) const 
  {
    return exprEngine()(loc + offset_m);
  }

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

  inline Element_t operator()(const Loc<Dim> &loc) const 
  {
    return read(loc);
  }

  //---------------------------------------------------------------------------
  // Return the domain.

  inline const Domain_t &domain() const { return domain_m; }

  //---------------------------------------------------------------------------
  // 
  
  inline Loc<Dim> offset() const
  {
    return offset_m;
  }
  
  //---------------------------------------------------------------------------
  // Accessors.

  inline const Expression_t &exprEngine() const { return exprEngine_m; }

  //---------------------------------------------------------------------------
  // Need to pass lock requests to the contained engine.

  template<class RequestType>
  inline
  typename DataObjectRequest<RequestType>::Type_t
  dataObjectRequest(const DataObjectRequest<RequestType> &req) const
  {
    return exprEngine().engine().dataObjectRequest(req);
  }

  //---------------------------------------------------------------------------
  // viewDomain() gives the region of the expression needed to compute a given
  // region of the shift-engine.
  //---------------------------------------------------------------------------

  inline
  Interval<Dim> viewDomain(const Interval<Dim> &domain) const
  {
    Interval<Dim> ret;
    int d;
    for (d = 0; d < Dim; ++d)
    {
      ret[d] =
	Interval<1>(
		    domain[d].first() + offset_m[d].first(),
		    domain[d].last() + offset_m[d].first()
		    );
    }
    return ret;
  }

  inline
  INode<Dim> viewDomain(const INode<Dim> &inode) const
  {
    return INode<Dim>(inode, viewDomain(inode.domain()));
  }

  inline
  Interval<Dim> intersectDomain() const
  {
    Interval<Dim> ret;
    int d;
    for (d = 0; d < Dim; ++d)
    {
      ret[d] =
	Interval<1>(
		    domain_m[d].first() + offset_m[d].first(),
                    domain_m[d].last() + offset_m[d].first()
		    );
    }
    return ret;
  }

private:

  Interval<Dim> domain_m;
  Loc<Dim> offset_m;
  Expression_t exprEngine_m;
};

/**
 * NewEngine<Engine,SubDomain>
 *
 * Specializations of NewEngine for subsetting a constant-function-engine with
 * an arbitrary domain. 
 *
 * FIXME:We should consider making NewEngine a functor rather than using engine
 * constructors to initialize engine views. NewEngineEngine and NewEngineDomain
 * were introduced because a given engine doesn't know about all the engines
 * that can return a view implemented by that engine.  (For example, Brick
 * shouldn't know about MultiPatch, but MultiPatch<Brick> can return a view
 * of one of its patches via an INode<> view.)  Forwarding requests now
 * requires specializing 3 structs because of the constructor paradigm of
 * Engine(otherEngine, domain).
 */

template <int Dim, class T, class E>
struct NewEngine<Engine<Dim, T, FieldShift<E> >, Interval<Dim> >
{
  typedef typename NewEngine<E, Interval<Dim> >::Type_t Type_t;
};

template <int Dim, class T, class E>
struct NewEngineEngine<Engine<Dim, T, FieldShift<E> >, Interval<Dim> >
{
  typedef typename NewEngineEngine<E, Interval<Dim> >::Type_t Type_t;
  static inline
  Type_t apply(const Engine<Dim, T, FieldShift<E> > &e, const Interval<Dim> &d)
  {
    return NewEngineEngine<E, Interval<Dim> >::apply(e.exprEngine(),
                                                     e.viewDomain(d));
  }
};

template <int Dim, class T, class E>
struct NewEngineDomain<Engine<Dim, T, FieldShift<E> >, Interval<Dim> >
{
  typedef typename NewEngineDomain<E, Interval<Dim> >::Type_t Type_t;
  static inline
  Type_t apply(const Engine<Dim, T, FieldShift<E> > &e, const Interval<Dim> &d)
  {
    return NewEngineDomain<E, Interval<Dim> >::apply(e.exprEngine(),
                                                     e.viewDomain(d));
  }
};

template <int Dim, class T, class E>
struct NewEngine<Engine<Dim, T, FieldShift<E> >, INode<Dim> >
{
  typedef typename NewEngine<E, INode<Dim> >::Type_t Type_t;
};

template <int Dim, class T, class E>
struct NewEngineEngine<Engine<Dim, T, FieldShift<E> >, INode<Dim> >
{
  typedef typename NewEngineEngine<E, INode<Dim> >::Type_t Type_t;
  static inline
  Type_t apply(const Engine<Dim, T, FieldShift<E> > &e, const INode<Dim> &d)
  {
    return NewEngineEngine<E, INode<Dim> >::apply(e.exprEngine(),
                                                  e.viewDomain(d));
  }
};

template <int Dim, class T, class E>
struct NewEngineDomain<Engine<Dim, T, FieldShift<E> >, INode<Dim> >
{
  typedef typename NewEngineDomain<E, INode<Dim> >::Type_t Type_t;
  static inline
  Type_t apply(const Engine<Dim, T, FieldShift<E> > &e, const INode<Dim> &d)
  {
    return NewEngineDomain<E, INode<Dim> >::apply(e.exprEngine(),
                                                  e.viewDomain(d));
  }
};

/**
 * There are potentially many ways to construct field stencils.
 * FieldShiftSimple assumes that you just need to construct the output field
 * and stick ONE stencil engine into it.  Maybe this class can be generalized
 * for fields that contain multiple stencil engines.
 */

template<class Expression>
struct FieldShiftSimple
{
  typedef typename Expression::MeshTag_t MeshTag_t;
  typedef typename Expression::Element_t OutputElement_t;
  enum { outputDim = Expression::dimensions };

  typedef typename Expression::Engine_t InputEngine_t;
  typedef FieldShift<InputEngine_t> OutputEngineTag_t;

  typedef Field<MeshTag_t, OutputElement_t, OutputEngineTag_t> Type_t;

  typedef Engine<outputDim, OutputElement_t, OutputEngineTag_t> SEngine_t;

  static inline
  Type_t make(const Expression &f,
	      const FieldOffset<outputDim> &s1,
              const Centering<outputDim> &centering)
  {
    // Create a model field with the new centering.

    Type_t h(centering, f.layout(), f.mesh());
    h.fieldEngine().physicalCellDomain() = f.fieldEngine().physicalCellDomain();

    // Could change this to loop over centerings.

#   if POOMA_BOUNDS_CHECK
    if (f.numSubFields() > 0)
      {
	PInsist((s1.subFieldNumber() < f.numSubFields()) &&
		(s1.subFieldNumber() >= 0),
		"subField bounds error");
	PInsist(contains(f[s1.subFieldNumber()].totalDomain(),
			 f[s1.subFieldNumber()].domain() + s1.cellOffset()),
		"Field operator()(FieldOffset) bounds error");
      }
    else
      {
	PInsist(s1.subFieldNumber() == 0,
		"subField bounds error");
	PInsist(contains(f.totalDomain(), f.domain() + s1.cellOffset()),
		"Field operator()(FieldOffset) bounds error");
      }
#   endif   

    Expression fld = 
      (f.numSubFields() > 0) ? f[s1.subFieldNumber()] : f;
    const Loc<outputDim> &offset = s1.cellOffset();

    // FIXME: need to adjust guard layers based on centering???

    GuardLayers<outputDim> og(fld.fieldEngine().guardLayers());
    for (int d = 0; d < outputDim; d++)
      {
	og.lower(d) += offset[d].first();
	og.upper(d) -= offset[d].first();
      }

    h.fieldEngine().guardLayers() = og;
    h.fieldEngine().engine() = SEngine_t(fld.engine(), offset, fld.domain());

    return h;
  }

  static inline
  Type_t make(const Expression &f,
	      const std::vector<FieldOffset<outputDim> > &vs1,
              const Centering<outputDim> &centering)
  {
    typedef typename std::vector<FieldOffset<outputDim> >::size_type size_type;

    // Create a model field with the new centering.

    Type_t h(centering, f.layout(), f.mesh());
    h.fieldEngine().physicalCellDomain() = f.fieldEngine().physicalCellDomain();

    // Could change this to loop over centerings.

    PInsist(vs1.size() == centering.size(),
	    "The FieldOffset vector's length must match the centering's size.");

    // This code should simplify when unified access to fields with
    // one or more subfields is possible.

    for (size_type s1Index = 0; s1Index < vs1.size(); ++s1Index) {
      const FieldOffset<outputDim> s1 = vs1[s1Index];
      Type_t hField = (h.numSubFields() > 0) ? h[s1Index] : h;

#   if POOMA_BOUNDS_CHECK
      if (f.numSubFields() > 0)
	{
	  PInsist((s1.subFieldNumber() < f.numSubFields()) &&
		  (s1.subFieldNumber() >= 0),
		  "subField bounds error");
	  PInsist(contains(f[s1.subFieldNumber()].totalDomain(),
			   f[s1.subFieldNumber()].domain() + s1.cellOffset()),
		  "Field operator()(FieldOffset) bounds error");
	}
      else
	{
	  PInsist(s1.subFieldNumber() == 0,
		  "subField bounds error");
	  PInsist(contains(f.totalDomain(), f.domain() + s1.cellOffset()),
		  "Field operator()(FieldOffset) bounds error");
	}
#   endif   

      Expression fld = 
	(f.numSubFields() > 0) ? f[s1.subFieldNumber()] : f;
      const Loc<outputDim> &offset = s1.cellOffset();

      GuardLayers<outputDim> og(fld.fieldEngine().guardLayers());
      for (int d = 0; d < outputDim; d++)
	{
	  og.lower(d) += offset[d].first();
	  og.upper(d) -= offset[d].first();
	}

      hField.fieldEngine().guardLayers() = og;
      hField.fieldEngine().engine() =
	SEngine_t(fld.engine(), offset, fld.domain());
    }

    return h;
  }

};

/**
 * Specializations for selecting the appropriate evaluator for the Shift
 * engine.  We just get the appropriate types from the Expression's engine.
 */

template<class Expression>
struct EvaluatorEngineTraits<FieldShift<Expression> >
{
  typedef typename CreateLeaf<Expression>::Leaf_t Expr_t;
  typedef typename
    ForEach<Expr_t, EvaluatorTypeTag, EvaluatorCombineTag>::Type_t
      Evaluator_t;
};


/**
 * FieldShiftIntersector is a special intersector that gets used when we come
 * across a stencil object in an expression.
 */

template<int Dim, class Intersect>
class FieldShiftIntersector
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants

  typedef typename Intersect::IntersectorData_t         IntersectorData_t;
  typedef FieldShiftIntersector<Dim, Intersect>       This_t;
  typedef typename IntersectorData_t::const_iterator    const_iterator;
  typedef RefCountedPtr<IntersectorData_t>              DataPtr_t;
  typedef Interval<Dim>                                 Domain_t;
  
  enum { dimensions = Intersect::dimensions };
  
  //---------------------------------------------------------------------------
  // Constructors

  FieldShiftIntersector(const This_t &model)
    : domain_m(model.domain_m), intersector_m(model.intersector_m)
  { }

  FieldShiftIntersector(const Domain_t &dom, const Intersect &intersect)
    : domain_m(dom), intersector_m(intersect)
  { }

  This_t &operator=(const This_t &model)
  {
    if (this != &model)
    {
      domain_m = model.domain_m;
      intersector_m = model.intersector_m;
    }
    return *this;
  }

  ~FieldShiftIntersector() { }

  inline DataPtr_t &data() { return intersector_m.data(); }
  inline const DataPtr_t &data() const { return intersector_m.data(); }

  //---------------------------------------------------------------------------
  // Accessors

  // STL iterator support.
  
  inline const_iterator begin() const { return data()->inodes_m.begin(); }
  inline const_iterator end() const { return data()->inodes_m.end(); }

  //---------------------------------------------------------------------------
  // Intersect routines

  // All domains.
  
  template<class Engine>
  inline void intersect(const Engine &engine) 
  {
    typedef typename NewEngine<Engine, Interval<Dim> >::Type_t NewEngine_t;

    NewEngine_t newEngine(engine, domain_m);

    intersector_m.intersect(newEngine);

    data()->shared(engine.layout().ID(), newEngine.layout().ID());
  }

  template<class Engine, int Dim2>
  inline bool intersect(const Engine &engine, const GuardLayers<Dim2> &) 
  {
    intersect(engine);
    return true;
  }

private:

  
  Interval<Dim> domain_m;
  Intersect     intersector_m;
};


/**
 * IntersectEngine specialization
 */

template <int Dim, class T, class Expression, class Intersect>
struct LeafFunctor<Engine<Dim, T, FieldShift<Expression> >,
  ExpressionApply<IntersectorTag<Intersect> > >
{
  typedef int Type_t;

  static
  int apply(const Engine<Dim, T, FieldShift<Expression> > 
	    &engine, const ExpressionApply<IntersectorTag<Intersect> > &tag)
  {
    // We offset the domain to get a domain in the viewed engine that
    // the stencil looks at.  The intersection is performed with a view
    // of the contained engine over this domain.  The resulting answer works
    // even though the stencil looks beyond this domain, because the viewed
    // field guarantees enough guard layers for the stencil to work.
    // (Presently this assumption isn't checked anywhere, so a lack of guard
    // cells results in an error in the multipatch inode view.)

    typedef FieldShiftIntersector<Dim, Intersect> NewIntersector_t;
    NewIntersector_t newIntersector(engine.intersectDomain(),
				    tag.tag().intersector_m);

    expressionApply(engine.field(),
		    IntersectorTag<NewIntersector_t>(newIntersector));

    return 0;
  }
};

/**
 * Specialization of  DataObjectRequest engineFunctor to pass the request to
 * the contained engine.
 */

template<class RequestType> class DataObjectRequest;

template <int Dim, class T, class Expression, class RequestType>
struct EngineFunctor<Engine<Dim, T, FieldShift<Expression> >,
  DataObjectRequest<RequestType> >
{
  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;

  static Type_t
  apply(const Engine<Dim, T, FieldShift<Expression> > &engine,
	const DataObjectRequest<RequestType> &tag)
  {
    return engineFunctor(engine.field().engine(), tag);
  }
};

/**
 * The generic version of EngineView just accesses the contained engine and
 * applies EngineView to it.
 *
 * The default version doesn't fiddle with the domain, since it is assumed
 * that the typical view doesn't need to.  Specializations will be required
 * for INode views etc...  Probably we should come up with a generic approach.
 */

template <int Dim, class T, class Expression, class Tag>
struct LeafFunctor<Engine<Dim, T, FieldShift<Expression> >,
  EngineView<Tag> >
{
  typedef LeafFunctor<Expression, EngineView<Tag> > LeafFunctor_t;
  typedef typename LeafFunctor_t::Type_t NewViewed_t;
  typedef Engine<Dim, T, FieldShift<NewViewed_t> > Type_t;

  static
  Type_t apply(const Engine<Dim, T,
	       FieldShift<Expression> > &engine,
	       const EngineView<Tag> &tag)
  {
    return Type_t(LeafFunctor_t::apply(engine.field(), tag), engine);
  }
};

template <int Dim, class T, class Expression, class Tag>
struct LeafFunctor<Engine<Dim, T, FieldShift<Expression> >,
                   ExpressionApply<Tag> >
{
  typedef LeafFunctor<Expression, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  static
  Type_t apply(const Engine<Dim, T,
	       FieldShift<Expression> > &engine,
	       const ExpressionApply<Tag> &tag)
  {
    return LeafFunctor_t::apply(engine.field(), tag);
  }
};

/**
 * View2<Field, FieldOffset, Centering> specialization for indexing a
 * field with a FieldOffset.
 */

template<class MeshTag, class T, class EngineTag, int Dim>
struct View2<Field<MeshTag, T, EngineTag>, FieldOffset<Dim>,
             Centering<Dim> >
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<MeshTag, T, EngineTag> Subject_t;
  typedef typename Subject_t::Engine_t Engine_t;

  // The return types.

  typedef Field<MeshTag, T, FieldShift<Engine_t> > ReadType_t;
  typedef Field<MeshTag, T, FieldShift<Engine_t> > Type_t;

  // The functions that do the indexing.

  inline static
  Type_t make(const Subject_t &f, const FieldOffset<Dim> &s1,
              const Centering<Dim> &c)
  {
    CTAssert(Dim == Subject_t::dimensions);
    return FieldShiftSimple<Subject_t>::make(f, s1, c);
  }

  inline static
  ReadType_t makeRead(const Subject_t &f, const FieldOffset<Dim> &s1,
                      const Centering<Dim> &c)
  {
    CTAssert(Dim == Subject_t::dimensions);
    return FieldShiftSimple<Subject_t>::make(f, s1, c);
  }
};

/**
 * View2<Field, vector<FieldOffset>, Centering> specialization for indexing a
 * field with a vector<FieldOffset>.
 */

template<class MeshTag, class T, class EngineTag, int Dim>
struct View2<Field<MeshTag, T, EngineTag>, std::vector<FieldOffset<Dim> >,
             Centering<Dim> >
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<MeshTag, T, EngineTag> Subject_t;
  typedef typename Subject_t::Engine_t Engine_t;

  // The return types.

  typedef Field<MeshTag, T, FieldShift<Engine_t> > ReadType_t;
  typedef Field<MeshTag, T, FieldShift<Engine_t> > Type_t;

  // The functions that do the indexing.

  inline static
  Type_t make(const Subject_t &f,
	      const std::vector<FieldOffset<Dim> > &s1,
              const Centering<Dim> &c)
  {
    CTAssert(Dim == Subject_t::dimensions);
    return FieldShiftSimple<Subject_t>::make(f, s1, c);
  }

  inline static
  ReadType_t makeRead(const Subject_t &f,
		      const std::vector<FieldOffset<Dim> > &s1,
                      const Centering<Dim> &c)
  {
    CTAssert(Dim == Subject_t::dimensions);
    return FieldShiftSimple<Subject_t>::make(f, s1, c);
  }
};

#endif // POOMA_FIELD_DIFFOPS_FIELDSHIFTENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldShiftEngine.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:44 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
