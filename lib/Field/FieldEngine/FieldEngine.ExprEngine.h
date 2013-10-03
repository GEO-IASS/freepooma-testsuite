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
//   FieldEngine specialized for Expression-engines
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_FIELDENGINE_FIELDENGINEBASE_EXPRENGINE__H
#define POOMA_FIELD_FIELDENGINE_FIELDENGINEBASE_EXPRENGINE__H

/** @file
 * @ingroup Field
 * @brief
 * FieldEngine specialization for expression engine.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include "Engine/ExpressionEngine.h"
#include "Layout/GuardLayers.h"
#include "Utilities/PAssert.h"
#include "Field/FieldEngine/FieldEngine.h"


//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag> class Field;


//-----------------------------------------------------------------------------
// Handy combiners for getting the far-left field in an expression.
//-----------------------------------------------------------------------------

struct FarLeftTag
{
  POOMA_PURIFY_CONSTRUCTORS(FarLeftTag)
};

template<class G1, class T1, class E1, class Op>
struct Combine1<Field<G1, T1, E1>, Op, FarLeftTag>
{
  typedef Field<G1, T1, E1> Type_t;
  inline static
  const Type_t &combine(const Field<G1, T1, E1> &a,
			const FarLeftTag &) 
    { 
      return a; 
    }
};

template<class G1, class T1, class E1, class G2, class T2, class E2, 
  class Op>
struct Combine2<Field<G1, T1, E1>, Field<G2, T2, E2>, Op, FarLeftTag>
{
  typedef Field<G1, T1, E1> Type_t;
  inline static
  const Type_t &combine(const Field<G1, T1, E1> &a, 
                        const Field<G2, T2, E2> &, FarLeftTag)
    {
      return a;
    }
};

template<class T, class G2, class T2, class E2, class Op>
struct Combine2<T, Field<G2, T2, E2>, Op, FarLeftTag>
{
  typedef Field<G2, T2, E2> Type_t;
  inline static
  const Type_t &combine(const T &, 
                        const Field<G2, T2, E2> &b, FarLeftTag)
    {
      return b;
    }
};

template<class G1, class T1, class E1, class T, class Op>
struct Combine2<Field<G1, T1, E1>, T, Op, FarLeftTag>
{
  typedef Field<G1, T1, E1> Type_t;
  inline static
  const Type_t &combine(const Field<G1, T1, E1> &a, 
                        const T &, FarLeftTag)
    {
      return a;
    }
};

template<class A,class B,class C,class Op>
struct Combine3<A, B, C, Op, FarLeftTag>
{
  typedef typename Combine2<A, B, Op, FarLeftTag>::Type_t Type1_t;
  typedef typename Combine2<Type1_t, C, Op, FarLeftTag>::Type_t Type_t;
  inline static
  const Type_t &combine(const A& a, const B& b, const C& c,
			const FarLeftTag& t)
  {
    return
      Combine2<Type1_t, C,
      Op, FarLeftTag>::combine(Combine2<A, B, Op,
			       FarLeftTag>::combine(a, b, t), c, t);
  }
};

//-----------------------------------------------------------------------------
// This version of LeafFunctor is used by Expression-Engines to 
// get at the far-left field in an expression. 
//-----------------------------------------------------------------------------

template<class GeometryTag, class T, class EngineTag>
struct LeafFunctor<Field<GeometryTag, T, EngineTag>, FarLeftTag>
{
  typedef Field<GeometryTag, T, EngineTag> Type_t;
  inline static
  const Type_t &apply(const Field<GeometryTag, T, EngineTag> &f, 
    const FarLeftTag &) 
    {
      return f;
    }
};

template<class T>
struct LeafFunctor<Scalar<T>, FarLeftTag>
{
  typedef Scalar<T> Type_t;
  inline static
  const Type_t &apply(const Scalar<T> &s, const FarLeftTag &) 
    {
      return s;
    }
};


// ----------------------------------------------------------------------------
// FieldEngine<Dim, T, ExpressionTag<Expr> >
//
// This is a specialization of FieldEngine for expression-engines.
// ----------------------------------------------------------------------------

template<class Mesh, class T, class Expr>
class FieldEngine<Mesh, T, ExpressionTag<Expr> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.
    
  enum { dimensions = Mesh::dimensions };
  enum { Dim = dimensions };
  typedef ExpressionTag<Expr> EngineTag_t;
  typedef FieldEngine<Mesh, T, EngineTag_t> This_t;
  typedef Engine<Dim, T, EngineTag_t> Engine_t;
  typedef typename Engine_t::Domain_t Domain_t;
  typedef typename Engine_t::Layout_t Layout_t;
  typedef typename Engine_t::Element_t Element_t;
  typedef typename Engine_t::ElementRef_t ElementRef_t;
  typedef GuardLayers<Dim> GuardLayers_t;
  typedef typename ForEach<Expr, FarLeftTag, FarLeftTag>::Type_t 
    ReferenceField_t;

  //---------------------------------------------------------------------------
  // Constructors.

  // Expression constructor.  
  
  FieldEngine(const Engine_t &e)
  : engine_m(e),
    referenceField_m(
      forEachRef(engine_m.expression(), FarLeftTag(), FarLeftTag()))
    { }  
    
  // Domain view constructor. 

  template<class Expr2, class Domain>  
  FieldEngine(const FieldEngine<Mesh, T, ExpressionTag<Expr2> > &model, 
    const Domain &d)
  : engine_m(NewEngineEngine<Engine<Dim, T, ExpressionTag<Expr2> >, Domain>::apply(model.engine(), d),
	     NewEngineDomain<Engine<Dim, T, ExpressionTag<Expr2> >, Domain>::apply(model.engine(), d)
	     ),
    referenceField_m(
      forEachRef(engine_m.expression(), FarLeftTag(), FarLeftTag()))
    { }  

  ///@name Sub-field view constructors
  //@{

  /// This is when we want to construct a view of
  /// one of the subFields in our top-level list from material and centering.

  template<class Expr2>
  FieldEngine(const FieldEngine<Mesh, T, ExpressionTag<Expr2> > &model, int m, int c)
    : engine_m(Expr(model.engine().expression(), m, c)),
      referenceField_m(forEachRef(engine_m.expression(),
                                  FarLeftTag(), FarLeftTag()))
  { }

  /// Sub-field view for a material.

  template<class Expr2>
  FieldEngine(const FieldEngine<Mesh, T, ExpressionTag<Expr2> > &model, int m,
	      const Pooma::MaterialViewTag& tag)
    : engine_m(Expr(model.engine().expression(), m, tag)),
      referenceField_m(forEachRef(engine_m.expression(),
                                  FarLeftTag(), FarLeftTag()))
  { }  

  /// Sub-field view for a centering.

  template<class Expr2>
  FieldEngine(const FieldEngine<Mesh, T, ExpressionTag<Expr2> > &model, int c,
	      const Pooma::CenteringViewTag& tag)
    : engine_m(Expr(model.engine().expression(), c, tag)),
      referenceField_m(forEachRef(engine_m.expression(),
                                  FarLeftTag(), FarLeftTag()))
  { }  

  /// sub-material view. Deprecated.

  template<class Expr2>
  FieldEngine(const FieldEngine<Mesh, T, ExpressionTag<Expr2> > &model, int m)
    : engine_m(Expr(model.engine().expression, m)),
      referenceField_m(forEachRef(engine_m.expression(),
                                  FarLeftTag(), FarLeftTag()))
  { }  

  /// sub-center view. Deprecated.

  template<class Expr2>
  FieldEngine(int c, const FieldEngine<Mesh, T, ExpressionTag<Expr2> > &model)
    : engine_m(Expr(c, model.engine().expression())),
      referenceField_m(forEachRef(engine_m.expression(),
                                  FarLeftTag(), FarLeftTag()))
  { }

  //@}

  // Very important! Copy constructor is needed so that referenceField_m
  // doesn't refer to someone else's expression.      

  FieldEngine(const FieldEngine<Mesh, T, ExpressionTag<Expr> > &other)
  : engine_m(other.engine()),
    referenceField_m(forEachRef(engine_m.expression(),
				FarLeftTag(), FarLeftTag()))
  { }  
    
  //---------------------------------------------------------------------------
  // Accessors and modifiers.
    
  // FIXME: This function is deprecated.
  inline int numSubFields() const
    {
      return referenceField().numSubFields();
    }

  inline const Engine_t &engine() const
    {
      return engine_m;
    }

  inline Engine_t &engine()
    {
      return engine_m;
    }

  const ReferenceField_t &referenceField() const
  {
    return referenceField_m;
  }
  
  //---------------------------------------------------------------------------
  // Domain accessor functions. 
        
  inline const Domain_t physicalCellDomain() const
    {
      return referenceField().physicalCellDomain();
    }
        
  inline Domain_t totalCellDomain() const
    {
      return referenceField().totalCellDomain();
    }

  Domain_t physicalDomain() const
    {
      return referenceField().physicalDomain();
    }

  Domain_t totalDomain() const
    {
      return referenceField().totalDomain();
    }

  Domain_t physicalDomain(int iSubField) const
    {
      return referenceField().physicalDomain(iSubField);
    }

  Domain_t totalDomain(int iSubField) const
    {
      return referenceField().totalDomain(iSubField);
    }

  //---------------------------------------------------------------------------
  // Centering accessors.

  const Centering<Dim> &centering() const
  {
    return referenceField().centering();
  }        

  inline int centeringSize() const
    {
      return referenceField().centeringSize();
    }

  inline int numMaterials() const
    {
      return referenceField().numMaterials();
    }

  //---------------------------------------------------------------------------
  // Mesh accessors.

  Mesh &mesh()
  {
    return referenceField().mesh();
  }        

  const Mesh &mesh() const
  {
    return referenceField().mesh();
  }        

      
private:

  // There is no way to reinitialize the reference, so assignment
  // cannot be implemented.

  This_t &operator=(const This_t &other);
    
  Engine_t engine_m;
  const ReferenceField_t &referenceField_m;
};

template<class Mesh, class T, class Expr, class Tag>
struct LeafFunctor<FieldEngine<Mesh, T, ExpressionTag<Expr> >,
  ExpressionApply<Tag> >
{
  typedef FieldEngine<Mesh, T, ExpressionTag<Expr> > Subject_t;
  typedef typename Subject_t::Engine_t Engine_t;
  typedef LeafFunctor<Engine_t, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  inline static
  Type_t apply(const Subject_t &fieldEngineBase, 
	       const ExpressionApply<Tag> &tag)
  {
    LeafFunctor_t::apply(fieldEngineBase.engine(), tag);
    return 0;
  }
};

#endif // POOMA_FIELD_FIELDENGINE_FIELDENGINEBASE_EXPRENGINE__H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldEngine.ExprEngine.h,v $   $Author: richi $
// $Revision: 1.8 $   $Date: 2004/11/10 22:05:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
