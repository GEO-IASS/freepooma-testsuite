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
//   Field
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_FIELD_H
#define POOMA_FIELD_FIELD_H

/** @file
 * @ingroup Field
 * @brief
 * Ties together the notions of field-category and mesh.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Array/Array.h"
#include "Domain/CombineDomainOpt.h"
#include "Domain/NewDomain.h"
#include "Engine/ConstantFunctionEngine.h"
#include "Engine/Engine.h"
#include "Engine/EnginePatch.h"
#include "Engine/ExpressionEngine.h"
#include "Evaluator/Evaluator.h"
#include "PETE/PETE.h"
#include "Pooma/View.h"
#include "Utilities/PAssert.h"
#include "Utilities/RefCountedBlockPtr.h"
#include "Utilities/PerformUpdate.h"

// NOTE:  The current order of includes puts FieldCreateLeaf after the
// operators files to work around a bug with template instantiation in KCC.

#include "Field/FieldMakeReturn.h"
#include "Field/FieldOperators.h"
#include "Field/PoomaFieldOperators.h"
#include "Field/VectorFieldOperators.h"
#include "Field/FieldCreateLeaf.h"
#include "Field/FieldCentering.h"
#include "Field/NearestNeighbors.h"

#include "Field/PrintField.h"

#include "Field/FieldEngine/FieldEnginePatch.h"
#include "Field/FieldEngine/FieldEngine.h"

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

struct CompressibleBrick;

template<class Mesh, class T, class EngineTag>
class Field;

template<class LTag, class EngineTag>
struct MultiPatch;

template<int Dim> struct NoMesh;

struct POOMA_DEFAULT_ENGINE_TYPE;

template<class Subject> class SubFieldView;

template<class Subject, class Domain, bool SV>
struct View1Implementation;

class RelationListItem;

template <int Dim>
class FieldOffset;

template <int Dim>
class FieldOffsetList;

/// @name assign
/// Prototypes for the assign function used to assign an expression to a Field.
///
/// Prototypes defined here:
///  - Field = Field
///  - Field = Array
///  - Field = scalar
///  - Array = Field
///
/// If you wish to have Field work with other types of objects on the right-
/// hand side (for example, other classes that derive from Field), define
/// extra assign() functions that take the following arguments:
///
///   assign(Field<Mesh,T,EngineTag>, yourclass, Operator)
///
/// where "yourclass" is the class that you would like to work on the
/// right-hand side in an expression with a Field on the left-hand side.
//@{

template<class Mesh, class T, class EngineTag,
  class MeshTag2, class T2, class EngineTag2, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh,  T,  EngineTag> &lhs,
       const Field<MeshTag2, T2, EngineTag2> &rhs,
       const Op &op);

template<class Mesh, class T, class EngineTag, 
 int Dim2, class T2, class EngineTag2, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs, 
       const Array<Dim2, T2, EngineTag2> &rhs, const Op &op);

template<class Mesh, class T, class EngineTag, class T1, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs, 
       const T1 &rhs, const Op &op);

template<class Mesh, class T, class EngineTag, 
 int Dim2, class T2, class EngineTag2, class Op>
const Array<Dim2, T2, EngineTag2> &
assign(const Array<Dim2, T2, EngineTag2> &lhs,
       const Field<Mesh, T, EngineTag> &rhs, const Op &op);

template<class Mesh, class T, class EngineTag,
 class F, class B, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs,
       const WhereProxy<F, B> &rhs,
       const Op &op);

//@}

struct SubFieldViewFunctorTag;

/**
 * SubFieldView is used to implement the syntax f[i], which selects the
 * ith SubField for field f.
 */

template<class Mesh, class T, class EngineTag>
class SubFieldView<Field<Mesh, T, EngineTag> > {
  
public:
  
  // Use it to construct the output field type.

  typedef Field<Mesh, T, EngineTag> Type_t;

  // The function that actually creates the view.
  
  inline static Type_t make(const Type_t &s, int iSubField)
    {
      PBoundInsist(iSubField >= 0 && iSubField < s.numSubFields(),
		   "Field::operator[] indexing error.");
      return Type_t(s, iSubField);
    }

  inline static Type_t make(const Type_t &s, int m, int c)
  {
    PBoundInsist(m >= 0 && m < s.numMaterials()
		 && c >= 0 && c < s.centeringSize(),
		 "Field::subField(m, c) indexing error.");
    return Type_t(s, m, c);
  }

  inline static Type_t make(const Type_t &s, int c, const Pooma::CenteringViewTag &tag)
  {
    PBoundInsist(c >= 0 && c < s.centeringSize(),
		 "Field::center(c) indexing error.");
    return Type_t(s, c, tag);
  }

  inline static Type_t make(const Type_t &s, int m, const Pooma::MaterialViewTag &tag)
  {
    PBoundInsist(m >= 0 && m < s.numMaterials(),
		 "Field::material(m) indexing error.");
    return Type_t(s, m, tag);
  }
};

template<class Mesh, class T, class Expr>
class SubFieldView<Field<Mesh, T, ExpressionTag<Expr> > > {
  
public:
  
  // Use it to construct the output field type.

  typedef Field<Mesh, T, ExpressionTag<Expr> > Subject_t;
  typedef 
    typename ForEach<Expr, SubFieldViewFunctorTag, TreeCombine>::Type_t 
      Expr_t;
  typedef Field<Mesh, T, ExpressionTag<Expr_t> > Type_t;

  // The function that actually creates the view.
  
  inline static Type_t make(const Subject_t &s, int iSubField)
    {
      PBoundInsist(iSubField >= 0 && iSubField < s.numSubFields(),
		   "Field::operator[] indexing error.");
      return Type_t(s, iSubField);
    }

  inline static Type_t make(const Subject_t &s, int m, int c)
  {
    PBoundInsist(m >= 0 && m < s.numMaterials()
		 && c >= 0 && c < s.centeringSize(),
		 "Field::subField(m, c) indexing error.");
    return Type_t(s, m, c);
  }

  inline static Type_t make(const Subject_t &s, int c, const Pooma::CenteringViewTag &tag)
  {
    PBoundInsist(c >= 0 && c < s.centeringSize(),
		 "Field::center(c) indexing error.");
    return Type_t(s, c, tag);
  }

  inline static Type_t make(const Subject_t &s, int m, const Pooma::MaterialViewTag &tag)
  {
    PBoundInsist(m >= 0 && m < s.numMaterials(),
		 "Field::material(m) indexing error.");
    return Type_t(s, m, tag);
  }
};


/**
 * View1Implementation<Field, D, SV> specialization for indexing a field
 * with a single domain. There is a single-valued version (SV == true)
 * and a multi-valued version (SV == false).
 */

// Single-valued version. Handles scalars and Locs.

template<class Mesh, class T, class EngineTag, class Domain>
struct View1Implementation<Field<Mesh, T, EngineTag>, Domain, true>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // The return types are pretty simple here.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  template<class S1, class Combine>
  inline static 
  Type_t make(const Subject_t &f, const S1 &s1,
	      const Combine &)
    {
      PAssert(f.numSubFields() == 0);
      
      Domain s(Combine::make(f, s1));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");
      return f.engine()(s);
    }

  template<class S1, class S2, class Combine>
  inline static 
  Type_t make(const Subject_t &f,
	      const S1 &s1, const S2 &s2,
	      const Combine &)
    {
      PAssert(f.numSubFields() == 0);
      
      Domain s(Combine::make(f, s1, s2));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");
      return f.engine()(s);
    }

  template<class S1, class S2, class S3,
    class Combine>
  inline static 
  Type_t make(const Subject_t &f,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const Combine &)
    {
      PAssert(f.numSubFields() == 0);
      
      Domain s(Combine::make(f, s1, s2, s3));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");
      return f.engine()(s);
    }

  template<class S1, class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &f, const S1 &s1,
	      const Combine &)
    {
      PAssert(f.numSubFields() == 0);
      
      Domain s(Combine::make(f, s1));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");
      return f.engine().read(s);
    }

  template<class S1, class S2, class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &f,
	      const S1 &s1, const S2 &s2,
	      const Combine &)
    {
      PAssert(f.numSubFields() == 0);
      
      Domain s(Combine::make(f, s1, s2));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");
      return f.engine().read(s);
    }

  template<class S1, class S2, class S3,
    class Combine>
  inline static 
  ReadType_t makeRead(const Subject_t &f,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const Combine &)
    {
      PAssert(f.numSubFields() == 0);
      
      Domain s(Combine::make(f, s1, s2, s3));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");
      return f.engine().read(s);
    }
};

// Non-single-valued implementation. Works for general domains
// including Nodes and INodes.

// Use this little traits class to deduce the geometry tag for a view.
// It is always a NoGeometry unless the view is from an interval or
// an INode.

template<int Dim, class Mesh, class Domain>
struct NewMeshTag
{
  typedef NoMesh<Dim> Type_t;
};

template<int Dim, class Mesh>
struct NewMeshTag<Dim, Mesh, Interval<Dim> >
{
  typedef Mesh Type_t;
};

template<int Dim, class Mesh>
struct NewMeshTag<Dim, Mesh, INode<Dim> >
{
  typedef Mesh Type_t;
};

template<class Mesh, class T, class EngineTag, class Domain>
struct View1Implementation<Field<Mesh, T, EngineTag>, Domain, false>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Engine_t Engine_t;
  typedef typename NewEngine<Engine_t, Domain>::Type_t NewEngine_t;
  typedef typename NewEngine_t::Element_t NewT_t;
  typedef typename NewEngine_t::Tag_t NewEngineTag_t;
  
  // Deduce the new Mesh.
  
  typedef typename
    NewMeshTag<NewEngine_t::dimensions, Mesh, Domain>::Type_t 
      NewMeshTag_t;
  
  // The output types.
  
  typedef Field<NewMeshTag_t, NewT_t, NewEngineTag_t> ReadType_t;
  typedef Field<NewMeshTag_t, NewT_t, NewEngineTag_t> Type_t;

  template<class S1, class Combine>
  static 
  Type_t make(const Subject_t &f, const S1 &s1,
	      const Combine &)
    {
      Domain s(Combine::make(f, s1));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");

      return Type_t(f, s);
    }

  template<class S1, class S2, class Combine>
  static 
  Type_t make(const Subject_t &f, const S1 &s1,
	      const S2 &s2, const Combine &)
    {
      Domain s(Combine::make(f, s1, s2));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");

      return Type_t(f, s);
    }

  template<class S1, class S2, class S3,
    class Combine>
  static 
  Type_t make(const Subject_t &f,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const Combine &)
    {
      Domain s(Combine::make(f, s1, s2, s3));
      PBoundInsist(contains(f.totalDomain(), s),
		   "Field view bounds error.");

      return Type_t(f, s);
    }

  template<class S1, class Combine>
  inline static 
  Type_t makeRead(const Subject_t &f, const S1 &s1,
	      const Combine &c)
    {
      return make(f, s1, c);
    }

  template<class S1, class S2, class Combine>
  inline static 
  Type_t makeRead(const Subject_t &f, const S1 &s1,
	      const S2 &s2, const Combine &c)
    {
      return make(f, s1, s2, c);
    }

  template<class S1, class S2, class S3,
    class Combine>
  inline static 
  Type_t makeRead(const Subject_t &f,
	      const S1 &s1, const S2 &s2, const S3 &s3,
	      const Combine &c)
    {
      return make(f, s1, s2, s3, c);
    }
};


/**
 * View1<Field, S1> specialization for indexing a field with a single domain.
 */

template<class Mesh, class T, class EngineTag, class Sub1>
struct View1<Field<Mesh, T, EngineTag>, Sub1>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  // At some point, we need to fix NewDomain1; until then, use
  // the temporary version from NewDomain.h.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef TemporaryNewDomain1<Domain_t, Sub1> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;

  // The optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;
  
  // The return types.
  
  typedef typename Dispatch_t::ReadType_t ReadType_t;
  typedef typename Dispatch_t::Type_t Type_t;

  // The functions that create the view.
  
  inline static
  Type_t make(const Subject_t &f, const Sub1 &s1)
    {
      return Dispatch_t::make(f, s1, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &f, const Sub1 &s1)
    {
      return Dispatch_t::makeRead(f, s1, Combine_t());
    }
};


/**
 * View1<Field, int> specialization for indexing a field with an int.
 */

template<class Mesh, class T, class EngineTag>
struct View1<Field<Mesh, T, EngineTag>, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // The functions that do the indexing.

  inline static
  Type_t make(const Subject_t &f, int s1)
    {
      PAssert(f.numSubFields() == 0);

      PBoundInsist(contains(f.totalDomain(), Loc<1>(s1)),
		   "Field view bounds error.");
      return f.engine()(s1);
    }

  inline static
  ReadType_t makeRead(const Subject_t &f, int s1)
    {
      PAssert(f.numSubFields() == 0);

      PBoundInsist(contains(f.totalDomain(), Loc<1>(s1)),
		   "Field view bounds error.");
      return f.engine().read(s1);
    }
};


/**
 * View2<Field, S1, S2> specialization for indexing a field with two
 * domains.
 */

template<class Mesh, class T, class EngineTag, 
  class Sub1, class Sub2>
struct View2<Field<Mesh, T, EngineTag>, Sub1, Sub2>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain2<Sub1, Sub2> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;

  // The optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;
  
  // The return types.
  
  typedef typename Dispatch_t::ReadType_t ReadType_t;
  typedef typename Dispatch_t::Type_t Type_t;

  // The functions that create the view.
  
  inline static
  Type_t make(const Subject_t &f, const Sub1 &s1, const Sub2 &s2)
    {
      return Dispatch_t::make(f, s1, s2, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &f, const Sub1 &s1, const Sub2 &s2)
    {
      return Dispatch_t::makeRead(f, s1, s2, Combine_t());
    }
};


/**
 * View2<Field, int, int> specialization for indexing a field with two
 * integers.
 */

template<class Mesh, class T, class EngineTag>
struct View2<Field<Mesh, T, EngineTag>, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // The functions that do the indexing.

  inline static
  Type_t make(const Subject_t &f, int s1, int s2)
    {
      PAssert(f.numSubFields() == 0);
      
      PBoundInsist(contains(f.totalDomain(), Loc<2>(s1, s2)),
		   "Field view bounds error.");
      return f.engine()(s1, s2);
    }

  inline static
  ReadType_t makeRead(const Subject_t &f, int s1, int s2)
    {
      PAssert(f.numSubFields() == 0);
      
      PBoundInsist(contains(f.totalDomain(), Loc<2>(s1, s2)),
		   "Field view bounds error.");
      return f.engine().read(s1, s2);
    }
};


/**
 * View2<Field, FieldOffset<Dim>, Loc<Dim> > specialization for
 * indexing a field with a FieldOffset and a Loc.
 */

template<class Mesh, class T, class EngineTag, int Dim>
struct View2<Field<Mesh, T, EngineTag>,
             FieldOffset<Dim>,
             Loc<Dim> >
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // The field's dimension (i.e., the number of indices required to select a point).
  
  enum { dimensions = Subject_t::dimensions };

  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // The functions that do the indexing.

  inline static
  Type_t make(const Subject_t &f,
	      const FieldOffset<dimensions> &fo,
	      const Loc<dimensions> &loc)
    {
      CTAssert(dimensions == Dim);
      
      if (f.numSubFields() > 0) {
	PBoundInsist(contains(f[fo.subFieldNumber()].totalDomain(),
			      loc + fo.cellOffset()),
		     "Field view bounds error.");
	return f[fo.subFieldNumber()].engine()(loc + fo.cellOffset());
      }
      else {
	PBoundInsist(contains(f.totalDomain(), loc + fo.cellOffset()),
		     "Field view bounds error.");
	return f.engine()(loc + fo.cellOffset());
      }
    }

  inline static
  ReadType_t makeRead(const Subject_t &f,
		      const FieldOffset<dimensions> &fo,
		      const Loc<dimensions> &loc)
    {
      if (f.numSubFields() > 0) {
	PBoundInsist(contains(f[fo.subFieldNumber()].totalDomain(),
			      loc + fo.cellOffset()),
		     "Field view bounds error.");
	return f[fo.subFieldNumber()].engine().read(loc + fo.cellOffset());
      }
      else {
	PBoundInsist(contains(f.totalDomain(), loc + fo.cellOffset()),
		     "Field view bounds error.");
	return f.engine().read(loc + fo.cellOffset());
      }
    }
};


/**
 * View3<Field, S1, S2, S3> specialization for indexing a field with three
 * domains.
 */

template<class Mesh, class T, class EngineTag, 
  class Sub1, class Sub2, class Sub3>
struct View3<Field<Mesh, T, EngineTag>, Sub1, Sub2, Sub3>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef NewDomain3<Sub1, Sub2, Sub3> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Subject_t, SDomain_t, sv> Dispatch_t;

  // The optimized domain combiner.
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;
  
  // The return types.
  
  typedef typename Dispatch_t::ReadType_t ReadType_t;
  typedef typename Dispatch_t::Type_t Type_t;

  // The functions that create the view.
  
  inline static
  Type_t make(const Subject_t &f, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3)
    {
      return Dispatch_t::make(f, s1, s2, s3, Combine_t());
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &f, const Sub1 &s1, const Sub2 &s2, 
    const Sub3 &s3)
    {
      return Dispatch_t::makeRead(f, s1, s2, s3, Combine_t());
    }
};


/**
 * View3<Field, int, int, int> specialization for indexing a field with three
 * integers.
 */

template<class Mesh, class T, class EngineTag>
struct View3<Field<Mesh, T, EngineTag>, int, int, int>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;


  // The return types.
  
  typedef typename Subject_t::Element_t ReadType_t;
  typedef typename Subject_t::ElementRef_t Type_t;

  // The functions that do the indexing.
  
  inline static
  Type_t make(const Subject_t &f, int s1, int s2, int s3)
    {
      PAssert(f.numSubFields() == 0);
      
      PBoundInsist(contains(f.totalDomain(), Loc<3>(s1, s2, s3)),
		   "Field view bounds error.");
      return f.engine()(s1, s2, s3);
    }
  
  inline static
  ReadType_t makeRead(const Subject_t &f, int s1, int s2, int s3)
    {
      PAssert(f.numSubFields() == 0);
      
      PBoundInsist(contains(f.totalDomain(), Loc<3>(s1, s2, s3)),
		   "Field view bounds error.");
      return f.engine().read(s1, s2, s3);
    }
};


//-----------------------------------------------------------------------------
// Patch specialization for Field.
//-----------------------------------------------------------------------------

template<class Subject> struct Patch;

template<class Mesh, class T, class EngineTag>
struct Patch<Field<Mesh, T, EngineTag> >
{
  typedef Field<Mesh, T, EngineTag> Subject_t;
  typedef typename Subject_t::Engine_t OldEngine_t;
  typedef typename EngineFunctor<OldEngine_t, EnginePatch>::Type_t Engine_t;

  // We've assumed that Mesh and T are the same for the patch engine.

  typedef Field<Mesh, T, typename Engine_t::Tag_t> Type_t;

  enum { dim = OldEngine_t::dimensions };

  inline static
  Type_t make(const Subject_t &f, int i)
  {
    PAssert(f.numSubFields() == 0);

    return Type_t(f, FieldEnginePatch<dim>(i, f.physicalDomain()));
  }
};

template<class Mesh, class T, class LTag, class EngineTag>
struct Patch<Field<Mesh, T, MultiPatch<LTag, EngineTag> > >
{
  typedef Field<Mesh, T, MultiPatch<LTag, EngineTag> > Subject_t;
  typedef typename Subject_t::Engine_t OldEngine_t;
  typedef typename EngineFunctor<OldEngine_t, EnginePatch>::Type_t Engine_t;

  // We've assumed that Mesh and T are the same for the patch engine.

  typedef Field<Mesh, T, typename Engine_t::Tag_t> Type_t;

  enum { dim = OldEngine_t::dimensions };
  typedef typename OldEngine_t::Layout_t Layout_t;
  typedef typename Layout_t::Value_t Node_t;

  inline static
  Type_t make(const Subject_t &f, int i)
  {
    PAssert(f.numSubFields() == 0);

    Node_t *node = f.engine().layout().nodeListLocal()[i];
    
    return Type_t(f, FieldEnginePatch<dim>(i, intersect(f.physicalDomain(),
                                                        node->domain())));
  }
};


//-----------------------------------------------------------------------------
// ComponentView specialization for Field. Implements views of the form
// f.comp(loc).
//-----------------------------------------------------------------------------

template<class Components, class Mesh, class T, class EngineTag>
struct ComponentView<Components, Field<Mesh, T, EngineTag> >
{
  // Convenience typedef for the thing we're taking a component view of.
  
  typedef Field<Mesh, T, EngineTag> Subject_t;

  // Deduce the template parameters for the output type.
  
  typedef Engine<Mesh::dimensions, T, EngineTag> Engine_t;
  typedef typename Engine_t::Element_t Element_t;
  typedef typename ComponentAccess<Element_t, Components>::Element_t NewT_t;
  typedef CompFwd<Engine_t, Components> NewEngineTag_t;
  
  // The output type.
  
  typedef Field<Mesh, NewT_t, NewEngineTag_t> Type_t;

  // A function that creates the view.
  
  inline static
  Type_t make(const Subject_t &f, const Components &c)
    {
      return Type_t(f, ComponentWrapper<Components>(c));
    }
};


/**
 * Field.
 */

template<class Mesh, 
         class T = POOMA_DEFAULT_ELEMENT_TYPE,
         class EngineTag = POOMA_DEFAULT_ENGINE_TYPE>
class Field {
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.
    
  /// The specification type.
  
  typedef Mesh MeshTag_t;
  typedef Mesh Mesh_t;

  /// The type.
  
  typedef T T_t;
    
  /// The engine tag.
  
  typedef EngineTag EngineTag_t;
  
  /// This class.
  
  typedef Field<Mesh, T, EngineTag> This_t;
  
  /// The field engine type.
  
  typedef FieldEngine<Mesh, T, EngineTag> FieldEngine_t;
  
  /// The dimension (i.e., the number of indices required to select a point).
  
  enum { dimensions = FieldEngine_t::dimensions };
  
  /// The engine type.
  
  typedef Engine<dimensions, T, EngineTag> Engine_t;

  /// Element_t is the type of elements managed by this field's engine. 
  typedef typename Engine_t::Element_t Element_t;
  /// ElementRef_t is the writable version.
  typedef typename Engine_t::ElementRef_t ElementRef_t;
  
  /// Layout_t is the Engine's layout.
  
  typedef typename Engine_t::Layout_t Layout_t;
  
  /// The types of the our domains. 

  typedef typename Engine_t::Domain_t Domain_t;

  /// The types of the our centering.

  typedef Centering<dimensions> Centering_t;

  // Fields may have relations attached to them.

  enum { hasRelations = true };

  //---------------------------------------------------------------------------
  /// @name User-callable constructors
  /// These ctors are meant to be called by users.
  //@{

  /// Mesh/centering/layout constructors. We use the specified mesh 
  /// object to initialize our mesh and the layout to initialize 
  /// the engines. Clearly, these must be synchronized. This is appropriate 
  /// for multi-patch engines. We just store the centering.

  Field()
  : fieldEngine_m()
    { } 

  /// This version is used for expressions.

  template<class I1>  
  explicit Field(const I1 &i1)
    : fieldEngine_m(i1)
  { } 

  /// Layout is templated so you can use a compatible layout to construct the
  /// engine.

  template<class Layout2>
  Field(const Centering_t &centering, const Layout2 &layout, const Mesh_t &mesh)
    : fieldEngine_m(centering, layout, mesh)
  { } 

  template<class Layout2>
  Field(int materials, const Centering_t &centering, const Layout2 &layout, const Mesh_t &mesh)
    : fieldEngine_m(centering, layout, mesh, materials)
  { } 

  template<class I1, class I2>  
  Field(const Centering_t &centering, const Layout_t &layout, const I1 &i1, const I2 &i2)
    : fieldEngine_m(centering, layout, Mesh_t(layout, i1, i2))
  { } 

  Field(const Centering_t &centering, const Layout_t &layout)
    : fieldEngine_m(centering, layout, Mesh_t(layout))
  { } 

  template<class I1, class I2>  
  Field(int materials, const Centering_t &centering, const Layout_t &layout,
        const I1 &i1, const I2 &i2)
    : fieldEngine_m(centering, layout, Mesh_t(layout, i1, i2), materials)
  { } 

  /// Copy constructor.
  
  Field(const This_t &model)
  : fieldEngine_m(model.fieldEngine())
  { }

  /// Copy initializer.
  
  void initialize(const This_t &model)
  {
    fieldEngine_m = model.fieldEngine();
  }

  /// Initializers that are equivalent to the constructors.
  
  template<class Layout2>
  void
  initialize(const Centering_t &centering, const Layout2 &layout,
             const Mesh_t &mesh)
  {
    fieldEngine_m = FieldEngine_t(centering, layout, mesh);
  } 

  template<class Layout2>
  void
  initialize(int materials, const Centering_t &centering,
             const Layout2 &layout, const Mesh_t &mesh)
  {
    fieldEngine_m = FieldEngine_t(centering, layout, mesh, materials);
  } 

  void
  initialize(const Centering_t &centering, const Layout_t &layout)
  {
    fieldEngine_m = FieldEngine_t(centering, layout, Mesh_t(layout));
  } 

  //@}

  //---------------------------------------------------------------------------
  /// @name Internal POOMA constructors
  /// These ctors are used internally by POOMA.
  /// They are not really meant to be called by users.
  //@{

  /// Model-initializer constructor. Used by SubFieldView and 
  /// View1Implementation above and by MakeFieldReturn in FieldCreateLeaf.h.

  template<class GT2, class T2, class ET2, class Initializer>
  Field(const Field<GT2, T2, ET2> &model, const Initializer &i)
  : fieldEngine_m(model.fieldEngine(), i)
  { }

  template<class ET2>
  Field(const Field<Mesh, T, ET2> &model, int m, int c)
    : fieldEngine_m(model.fieldEngine(), m, c)
  { }

  template<class ET2>
  Field(int c, const Field<Mesh, T, ET2> &model)
    : fieldEngine_m(c, model.fieldEngine())
  { }

  template<class ET2>
  Field(const Field<Mesh, T, ET2> &model, int c, const Pooma::CenteringViewTag &tag)
    : fieldEngine_m(model.fieldEngine(), c, tag)
  { }

  template<class ET2>
  Field(const Field<Mesh, T, ET2> &model, int m, const Pooma::MaterialViewTag &tag)
    : fieldEngine_m(model.fieldEngine(), m, tag)
  { }

  //---------------------------------------------------------------------------
  /// Empty destructor is fine for us.
  
  ~Field() { }

  //@}

  //---------------------------------------------------------------------------
  /// @name Accessors
  //@{

  inline const Engine_t &engine() const
    {
      return fieldEngine_m.engine();
    }

  inline Engine_t &engine()
    {
      return fieldEngine_m.engine();
    }

  inline const FieldEngine_t &fieldEngine() const
    {
      return fieldEngine_m;
    }
    
  inline FieldEngine_t &fieldEngine()
    {
      return fieldEngine_m;
    }
    
  inline int numSubFields() const
    {
      return fieldEngine_m.numSubFields();
    }
        
  const Centering<dimensions> &centering() const
  {
    return fieldEngine().centering();
  }
        
  const Centering<dimensions> centering(int c) const
  {
    return fieldEngine().centering()[c];
  }

  inline int centeringSize() const
  {
    return fieldEngine().centeringSize();
  }

  inline int numMaterials() const
    {
      return fieldEngine().numMaterials();
    }

  /// Returns the physical cell domain (as opposed to the vertex or actual domain).

  inline const Domain_t physicalCellDomain() const
    {
      return fieldEngine_m.physicalCellDomain();
    }
        
  /// Returns the total cell domain (including external guards).

  inline Domain_t totalCellDomain() const
    {
      return fieldEngine_m.totalCellDomain();
    }

  /// Returns the actual physical domain of the specified subfield (which
  /// is a vertex or a cell domain dependend on the centering of the subfield).

  Domain_t physicalDomain(int iSubfield) const
    {
      return fieldEngine_m.physicalDomain(iSubfield);
    }

  /// Returns the actual total domain of the specified subfield (which
  /// is a vertex or a cell domain dependend on the centering of the subfield).

  Domain_t totalDomain(int iSubfield) const
    {
      return fieldEngine_m.totalDomain(iSubfield);
    }

  /// For centerings of size one this returns the actual physical domain of
  /// the field. For centerings of size greater than one this returns the
  /// physical cell domain (dont use in this case).

  Domain_t physicalDomain() const
    {
      return fieldEngine_m.physicalDomain();
    }

  /// For centerings of size one this returns the actual total domain of
  /// the field. For centerings of size greater than one this returns the
  /// total cell domain (dont use in this case).

  Domain_t totalDomain() const
    {
      return fieldEngine_m.totalDomain();
    }

  /// Alias for physicalDomain().

  Domain_t domain() const
    {
      return fieldEngine_m.physicalDomain();
    }

  inline
  const Mesh_t &mesh() const
  {
    return fieldEngine_m.mesh();
  }
        
  inline Layout_t layout() const
  {
    return fieldEngine_m.engine().layout();
  }

  //@}
        
  //---------------------------------------------------------------------------
  /// Instruct the field to make its own copy of its data.
  /// Recursively call ourself with subfield views of this field. When we're
  /// through, tell the fieldEngine to make a distinct copy of itself.

  void makeOwnCopy()
  {
    // Make a distinct copy of the fieldEngine.
          
    fieldEngine_m.makeOwnCopy(*this);
  }
      
  
  //---------------------------------------------------------------------------
  /// @name Sub-field view creation functions
  /// A field consists of (potentially) several sub-fields. This function
  /// returns a view of one of these.
  //@{
  inline typename SubFieldView<This_t>::Type_t
  operator[](int iSubfield) const
    {
      typedef SubFieldView<This_t> Ret_t;
      return Ret_t::make(*this, iSubfield);
    }

  inline typename SubFieldView<This_t>::Type_t
  subField(int m, int c) const
  {
    typedef SubFieldView<This_t> Ret_t;
    return Ret_t::make(*this, m, c);
  }

  inline typename SubFieldView<This_t>::Type_t
  center(int c) const
  {
    typedef SubFieldView<This_t> Ret_t;
    return Ret_t::make(*this, c, Pooma::CenteringViewTag());
  }

  inline typename SubFieldView<This_t>::Type_t
  material(int m) const
  {
    PAssert(numMaterials() > 1);
    typedef SubFieldView<This_t> Ret_t;
    return Ret_t::make(*this, m, Pooma::MaterialViewTag());
  }
  //@}

  //---------------------------------------------------------------------------
  /// @name View-creation operations
  /// These operator() and read() functions take 
  /// zero or more sub-domains, which combine to form a domain with 
  /// dimensionality identical to the rank of the field. The zero argument 
  /// version returns a view of the physical domain and the 'All'-suffixed
  /// versions return a view of the total domain.
  //@{

  inline typename View1<This_t, Domain_t>::ReadType_t 
  read() const
    {
      typedef View1<This_t, Domain_t> Ret_t;
      return Ret_t::makeRead(*this, physicalDomain());
    }

  inline typename View1<This_t, Domain_t>::ReadType_t 
  readAll() const
    {
      typedef View1<This_t, Domain_t> Ret_t;
      return Ret_t::makeRead(*this, totalDomain());
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

  inline typename View1<This_t, Domain_t>::Type_t 
  operator()() const
    {
      typedef View1<This_t, Domain_t> Ret_t;
      return Ret_t::make(*this, physicalDomain());
    }

  inline typename View1<This_t, Domain_t>::Type_t 
  all() const
    {
      typedef View1<This_t, Domain_t> Ret_t;
      return Ret_t::make(*this, totalDomain());
    }

  template<class Sub1> 
  inline typename View1<This_t, Sub1>::Type_t 
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
  //@}

  //---------------------------------------------------------------------------
  /// @name Component-forwarding functions
  /// These work quite similarly to the
  /// ones from Array except we produce a Field with the same Mesh.
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

  template<class Components>
  inline typename ComponentView<Components, This_t>::Type_t
  comp(const Components &loc) const
  {
    return ComponentView<Components, This_t>::make(*this, loc);
  }
  //@}

  //---------------------------------------------------------------------------
  /// @name Patch accessor functions returns the i'th patch.
  //@{
  inline typename Patch<This_t>::Type_t
  patchLocal(EnginePatch::PatchID_t i) const
    {
      return Patch<This_t>::make(*this, i);
    }

  inline int
  numPatchesLocal() const
  {
    return engineFunctor(engine(), EngineNumPatches());
  }
  //@}

  //---------------------------------------------------------------------------
  /// @name Copy assignment operators
  /// We pack this assignment expression into a
  /// PETE binary expression tree node and then use this to construct an
  /// array with an expression engine. We then pass this on to an evaluator,
  /// which handles the computation. The first three versions handle assigning
  /// Arrays and ConstArrays to Arrays and the fourth one handles assigning
  /// scalars.
  //@{

  This_t &operator=(const This_t &rhs)
    {
      assign(*this, rhs, OpAssign());
      return *this;
    }

  const This_t &operator=(const This_t &rhs) const
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
  /// @name Relation support
  //@{
 
  /// add a relation
  void addRelation(RelationListItem *item) const
  {
    PAssert(numSubFields() == 0);
    
    fieldEngine_m.relations().addRelation(item);
  }
 
  /// remove all relations
  void removeRelations()
  {
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centering().size(); ++ c)
      {
        fieldEngine_m.data(m, c).relations().erase();
      }
    }
  }
 
  /// trigger all relations dirty (or all, if makeDirty is set)
  void applyRelations(bool makeDirty = false) const
  {
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centering().size(); ++ c)
      {
        if (makeDirty)
          fieldEngine_m.data(m, c).relations().setDirty();
        fieldEngine_m.data(m, c).relations().notifyPreRead();
      }
    }
  }
 
  /// dirty field, dirtying all relations
  void setDirty() const
  {
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centering().size(); ++ c)
      {
        fieldEngine_m.data(m, c).relations().setDirty();
      }
    }
  }

  /// clear dirty flag of field, clearing all relations dirty flag
  void clearDirty() const
  {
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centering().size(); ++ c)
      {
        fieldEngine_m.data(m, c).relations().clearDirty();
      }
    }
  }

  //@}

private:

  FieldEngine_t fieldEngine_m;

};


//----------------------------------------------------------------------
// Set up a little traits class that distinguishes between OpAssign and
// other assignment operators that read the LHS.
//----------------------------------------------------------------------

template<class Op>
struct AssignOpReadWriteTraits
{
  enum { readLHS = true };
};

template<>
struct AssignOpReadWriteTraits<OpAssign>
{
  enum { readLHS = false };
};


//----------------------------------------------------------------------
// Apply the ConformTag to the leaves of the tree.
// Check to see if a given Field conforms.
//----------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, int Dim>
struct LeafFunctor<Field<Mesh, T, EngineTag>, ConformTag<Dim> >
{
  typedef bool Type_t;
  static Type_t apply1(const Interval<Dim> &d, 
    const ConformTag<Dim> &ct)
    {
      return conforms(d, ct);
    }
  template<int Dim2>
  static Type_t apply1(const Interval<Dim2> &d, 
    const ConformTag<Dim> &ct)
    {
      return false;
    }
  static Type_t apply(const Field<Mesh, T, EngineTag> &f,
    const ConformTag<Dim> &ct)
    {
      return apply1(f.physicalDomain(), ct);
    }
};


//----------------------------------------------------------------------
// This specialization of LeafFunctor is used to pass the 
// DataObjectRequest functor down into the FieldEngine. The default 
// behavior, given in the functor below, is to just pass it on to the 
// fieldEngine's engine.
//----------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, class RequestType>
struct LeafFunctor<Field<Mesh, T, EngineTag>,
  DataObjectRequest<RequestType> >
{
  typedef Field<Mesh, T, EngineTag> Subject_t;
  typedef typename Subject_t::FieldEngine_t FieldEngine_t;
  typedef LeafFunctor<FieldEngine_t, DataObjectRequest<RequestType> > 
    LeafFunctor_t;
  typedef typename LeafFunctor_t::Type_t Type_t;
  enum { dataObject = LeafFunctor_t::dataObject };
  
  inline static
  Type_t apply(const Subject_t &f,
	           const DataObjectRequest<RequestType> &functor)
    {
      return LeafFunctor_t::apply(f.fieldEngine(), functor);
    }
};

template<class Mesh, class T, class EngineTag, class RequestType>
struct LeafFunctor<FieldEngine<Mesh, T, EngineTag>,
  DataObjectRequest<RequestType> >
{
  typedef typename FieldEngine<Mesh, T, EngineTag>::Engine_t 
    Engine_t;
  enum { dataObject = Engine_t::dataObject };
  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;
  inline static
  Type_t apply(const FieldEngine<Mesh, T, EngineTag> &f,
	           const DataObjectRequest<RequestType> &functor)
    {
      return DataObjectApply<dataObject>::apply(f.engine(), functor);
    }
};


//-----------------------------------------------------------------------------
// This specialization of LeafFunctor is used to get the domain type or the
// (zero-based) domain itself from a Field. Used only by Expression-Engine.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
struct LeafFunctor<Field<Mesh, T, EngineTag>, DomainFunctorTag>
{
  typedef typename Field<Mesh, T, EngineTag>::Domain_t Type_t;

  inline static Type_t apply(const Field<Mesh, T, EngineTag> &f, 
    const DomainFunctorTag &)
    {
      // Return zero-based domain.
      
      return f.physicalDomain() - f.physicalDomain().firsts();
    }
};


//-----------------------------------------------------------------------------
// This specialization of LeafFunctor is used to pass the ExpressionApply
// functor
// down into the FieldEngine. The default behavior, given in the functor
// below, is to just pass it on to the fieldEngine's engine.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, class Tag>
struct LeafFunctor<Field<Mesh, T, EngineTag>, ExpressionApply<Tag> >
{
  typedef Field<Mesh, T, EngineTag> Subject_t;
  typedef typename Subject_t::FieldEngine_t FieldEngine_t;
  typedef LeafFunctor<FieldEngine_t, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  inline static
  Type_t apply(const Subject_t &field, 
	       const ExpressionApply<Tag> &tag)
    {
      return LeafFunctor_t::apply(field.fieldEngine(), tag);
    }
};

template<class Mesh, class T, class EngineTag, class Tag>
struct LeafFunctor<Field<Mesh, T, EngineTag>, EngineView<Tag> >
{
  typedef Field<Mesh, T, EngineTag> Subject_t;
  typedef typename Subject_t::Engine_t Engine_t;
  typedef typename LeafFunctor<Engine_t, EngineView<Tag> >::Type_t NewEngine_t;
  typedef typename NewEngine_t::Tag_t NewEngineTag_t;

  // Don't bother computing NewGeometry tag yet.
  // For now all EngineView operations are equivalent to Interval views.

  typedef Field<Mesh, T, NewEngineTag_t> Type_t;

  inline static
  Type_t apply(const Subject_t &field,
	       const EngineView<Tag> &tag)
  {
    return Type_t(field, tag);
  }
};


//-----------------------------------------------------------------------------
// Handle general engine functor requests.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, class Tag>
struct LeafFunctor<Field<Mesh, T, EngineTag>, EngineFunctorTag<Tag> >
{
  typedef typename Field<Mesh,T,EngineTag>::Engine_t Engine_t;
  typedef typename EngineFunctor<Engine_t,Tag>::Type_t Type_t;
  inline static
  Type_t apply(const Field<Mesh, T, EngineTag> &field,
	       const EngineFunctorTag<Tag> &tag)
  {
    return EngineFunctor<Engine_t,Tag>::apply(field.engine(), tag.tag());
  }
};


//---------------------------------------------------------------------------
// A specialization of EngineFunctor for field.
//---------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, class Tag>
struct EngineFunctor<Field<Mesh, T, EngineTag>, Tag>
{
  typedef typename Field<Mesh, T, EngineTag>::Engine_t Engine_t;
  typedef typename EngineFunctor<Engine_t, Tag>::Type_t Type_t;

  inline static 
  Type_t apply(const Field<Mesh, T, EngineTag> &field,
	           const Tag &tag)
    {
      return engineFunctor(field.engine(), tag);
    }
};


//-----------------------------------------------------------------------------
// This version of LeafFunctor is used by Expression-Engines to used to 
// evaluate a Field using indices. 
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, int Dim>
struct LeafFunctor<Field<Mesh, T, EngineTag>, EvalLeaf<Dim> >
{
  typedef typename Field<Mesh, T, EngineTag>::Element_t Type_t;
  inline static
  Type_t apply(const Field<Mesh, T, EngineTag> &f, 
    const EvalLeaf<Dim> &t) 
    {
      return t.eval(f.engine());
    }
};


//-----------------------------------------------------------------------------
// These leaf functor specializations are used to notify a field or expression
// that it is going to be read and, therefore, needs to update itself. 
//
// The second handles fields other than those with expression-engines by simply
// calling applyRelations(). The third passes the tag to the leaves.
//
// Fields with engines that store internal fields AND don't copy those
// fields' relations to its list must provide a specialization to get the 
// internal fields to update.
//
// NOTE: we don't use the ExpressionApply machinery here because this really
// operate on the engines.
//
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
struct LeafFunctor<Field<Mesh, T, EngineTag>, 
  PerformUpdateTag>
{
  typedef Field<Mesh, T, EngineTag> Subject_t;
  typedef int Type_t;

  inline static
  Type_t apply(const Subject_t &f, const PerformUpdateTag &)
    {
      f.applyRelations();
      return 0;
    }
};

template<class Mesh, class T, class Expr>
struct LeafFunctor<Field<Mesh, T, ExpressionTag<Expr> >, 
  PerformUpdateTag>
{
  typedef Field<Mesh, T, ExpressionTag<Expr> > Subject_t;
  typedef int Type_t;

  inline static
  Type_t apply(const Subject_t &f, const PerformUpdateTag &tag)
    {
      forEach(f.engine().expression(), tag, NullCombine());
      return 0;
    }
};


//-----------------------------------------------------------------------------
// This version of LeafFunctor is used to determine the type resulting from a
// sub-field view. 
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
struct LeafFunctor<Field<Mesh, T, EngineTag>, SubFieldViewFunctorTag>
{
  typedef Field<Mesh, T, EngineTag> Type_t;
};

template<class T>
struct LeafFunctor<Scalar<T>, SubFieldViewFunctorTag>
{
  typedef Scalar<T> Type_t;
};


//-----------------------------------------------------------------------------
// This specialization of LeafFunctor is used to apply a view (subsetting) 
// operation to a Field. The domain will always be zero-based since this
// is used only by Expression-Engine. This is why we add the firsts() to
// the domain.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, class Domain>
struct LeafFunctor<Field<Mesh, T, EngineTag>, ViewFunctorTag<Domain> >
{
  typedef typename View1<Field<Mesh, T, EngineTag>, Domain>::Type_t 
    Type_t;
};


//-----------------------------------------------------------------------------
// Overload the << operator to print a Field to a stream.  This
// uses the 'PrintField' class to perform the formatting of the data.
// It will create a default printer, print out the field with it, and
// return the stream object.
//-----------------------------------------------------------------------------

template <class Mesh, class T, class EngineTag>
std::ostream &operator<<(std::ostream &o, 
  const Field<Mesh, T, EngineTag> &cf)
{
  Pooma::blockAndEvaluate();
  PrintField().print(o, cf);
  return o;
}

template <class Mesh, class T, class EngineTag>
std::fstream &operator<<(std::fstream &f, 
  const Field<Mesh, T, EngineTag> &cf)
{
  Pooma::blockAndEvaluate();
  PrintField().print(f, cf);
  return f;
}


//-----------------------------------------------------------------------------
// Traits class for expressions containing fields.
//-----------------------------------------------------------------------------

struct ExpressionIsField { };

template<class Mesh, class T, class EngineTag>
struct ExpressionTraits<Field<Mesh, T, EngineTag> >
{
  typedef ExpressionIsField Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsField, ExpressionIsField>
{
  typedef ExpressionIsField Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsField, ExpressionIsScalar>
{
  typedef ExpressionIsField Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsScalar, ExpressionIsField>
{
  typedef ExpressionIsField Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsField, ExpressionIsArray>
{
  typedef ExpressionIsField Type_t;
};

template<>
struct CombineExpressionTraits<ExpressionIsArray, ExpressionIsField>
{
  typedef ExpressionIsField Type_t;
};


//-----------------------------------------------------------------------------
// assign() function for Field assign-op array.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, 
 int Dim2, class T2, class EngineTag2, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs, 
       const Array<Dim2, T2, EngineTag2> &rhs, const Op &op)
{
  for (int m = 0; m < lhs.numMaterials(); ++m)
    {
      for (int c = 0; c < lhs.centeringSize(); ++c)
        {
          const Field<Mesh, T, EngineTag> &l = lhs.subField(m, c);

          if (AssignOpReadWriteTraits<Op>::readLHS)
            l.applyRelations();
  
          // Evaluate.

          Evaluator<MainEvaluatorTag>().evaluate(l, op, rhs);

          // Having done the evaluation, we need to notify the LHS
          // that we've just written.
  
          l.setDirty();
        }
    }
        
  return lhs;
}


//-----------------------------------------------------------------------------
// assign() function for Field assign-op Field.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag,
  class Mesh2, class T2, class EngineTag2, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs,
       const Field<Mesh2, T2, EngineTag2> &rhs,
       const Op &op)
{
  PAssert(lhs.numMaterials() == rhs.numMaterials()
	  && lhs.centeringSize() == rhs.centeringSize());

  for (int m = 0; m < lhs.numMaterials(); ++m)
    {
      for (int c = 0; c < lhs.centeringSize(); ++c)
        {
          const Field<Mesh, T, EngineTag> &l = lhs.subField(m, c);
          const typename SubFieldView<Field<Mesh2, T2, EngineTag2> >::Type_t &r = 
            rhs.subField(m, c);

          forEach(r, PerformUpdateTag(), NullCombine());

          if (AssignOpReadWriteTraits<Op>::readLHS)
            l.applyRelations();
  
          // Evaluate.

          Evaluator<MainEvaluatorTag>().evaluate(l, op, r);

          // Having done the evaluation, we need to notify the LHS
          // that we've just written.
  
          l.setDirty();
        }
    }
        
  return lhs;
}


//-----------------------------------------------------------------------------
// assign() function for Field assign-op scalar.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, class T1, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs, const T1 &rhs,
       const Op &op)
{
  for (int m = 0; m < lhs.numMaterials(); ++m)
    {
      for (int c = 0; c < lhs.centeringSize(); ++c)
        {
          const Field<Mesh, T, EngineTag> &l = lhs.subField(m, c);
      
          if (AssignOpReadWriteTraits<Op>::readLHS)
            l.applyRelations();

          // Make an array out of the scalar.

          typedef Field<Mesh, T, EngineTag> LHS_t;
          Array<LHS_t::dimensions, T1, ConstantFunction>  rhsExpr(l.physicalDomain());
          rhsExpr.engine().setConstant(rhs);
     
          // Evaluate. 

          Evaluator<MainEvaluatorTag>().evaluate(l, op, rhsExpr);

          // Having done the evaluation, we need to notify the LHS
          // that we've just written.
  
          l.setDirty();
        }
    }
        
  return lhs;
}


//-----------------------------------------------------------------------------
// assign() function for Array assign-op Field.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag, 
 int Dim2, class T2, class EngineTag2, class Op>
const Array<Dim2, T2, EngineTag2> &
assign(const Array<Dim2, T2, EngineTag2> &lhs, 
       const Field<Mesh, T, EngineTag> &rhs, const Op &op)
{
  PAssert(rhs.numMaterials() == 1
          && rhs.centeringSize() == 1);

  forEach(rhs, PerformUpdateTag(), NullCombine());
  
  Evaluator<MainEvaluatorTag>().evaluate(lhs, op, rhs);
        
  return lhs;
}


//-----------------------------------------------------------------------------
// assign() function for Field assign-op WhereProxy.
//-----------------------------------------------------------------------------

template<class Tree>
struct ConvertWhereProxy<ExpressionIsField, Tree>
{
  typedef MakeFieldReturn<Tree> Make_t;
};

template<class Mesh, class T, class EngineTag,
 class F, class B, class Op>
const Field<Mesh, T, EngineTag> &
assign(const Field<Mesh, T, EngineTag> &lhs,
       const WhereProxy<F, B> &rhs, const Op &op)
{
  assign(lhs, rhs.whereMask(), rhs.opMask(op));

  return lhs;
}


//-----------------------------------------------------------------------------
// Compute whether or not a Field is currently compressed.
//
// This is only a sensible thing to do if there are no subfields, hence the
// assertion.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
inline bool compressed(const Field<Mesh, T, EngineTag> &f)
{
  PAssert(f.numSubFields() == 0);
  return compressed(f.engine());
}


//-----------------------------------------------------------------------------
// Compute the number of elements that are currently compressed.
//
// This is only a sensible thing to do if there are no subfields, hence the
// assertion.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
inline long elementsCompressed(const Field<Mesh, T, EngineTag> &f)
{
  PAssert(f.numSubFields() == 0);
  return elementsCompressed(f.engine());
}


//-----------------------------------------------------------------------------
// (Try to) compress all the patches of the Field. Only need to do work with
// multipatch engines.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class LTag>
void
compress(Field<Mesh, T, MultiPatch<LTag,CompressibleBrick> > &f)
{
  for (int m = 0; m < f.numMaterials(); ++m)
  {
    for (int c = 0; c < f.centeringSize(); ++c)
    {
      compress(f.fieldEngine().data(m, c).engine());
    }
  }
}


//-----------------------------------------------------------------------------
// Manually uncompress all the patches of the Field. Only need to do work with
// multipatch engines.
//-----------------------------------------------------------------------------

template<class Mesh, class T, class LTag>
void
uncompress(Field<Mesh, T, MultiPatch<LTag,CompressibleBrick> > &f)
{
  for (int m = 0; m < f.numMaterials(); ++m)
  {
    for (int c = 0; c < f.centeringSize(); ++c)
    {
      uncompress(f.fieldEngine().data(m, c).engine());
    }
  }
}


//-----------------------------------------------------------------------------
// Functions for getting the number of materials for Arrays and Fields.
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag>
inline int numMaterials(const Array<Dim, T, EngineTag> &a)
{
  return 1;
}

template<class Mesh, class T, class EngineTag>
inline int numMaterials(const Field<Mesh, T, EngineTag> &f)
{
  return f.numMaterials();
}


//-----------------------------------------------------------------------------
// Functions for getting the number of centering points for Arrays and Fields.
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag>
inline int centeringSize(const Array<Dim, T, EngineTag> &a)
{
  return 1;
}

template<class Mesh, class T, class EngineTag>
inline int centeringSize(const Field<Mesh, T, EngineTag> &f)
{
  return f.centeringSize();
}


//-----------------------------------------------------------------------------
// Functions for taking subfield views of Arrays and Fields.
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag>
inline Array<Dim, T, EngineTag> &subField(Array<Dim, T, EngineTag> &a, int, int)
{
  return a;
}

template<class Mesh, class T, class EngineTag>
inline typename SubFieldView<Field<Mesh, T, EngineTag> >::Type_t 
subField(Field<Mesh, T, EngineTag> &f, int m, int c)
{
  return f.subField(m, c);
}

#endif // POOMA_FIELD_FIELD_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Field.h,v $   $Author: richi $
// $Revision: 1.89 $   $Date: 2004/11/29 12:34:07 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
