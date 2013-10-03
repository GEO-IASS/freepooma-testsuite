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
//   FieldEngineBase
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_FIELDENGINE_FIELDENGINE_H
#define POOMA_FIELD_FIELDENGINE_FIELDENGINE_H

/** @file
 * @ingroup Field
 * @brief
 * FieldEngine and FieldEngineBaseData classes.
 *
 * POOMA supports a flexible form 
 * of "centering" that allows a hierarchy of multiple centering points per 
 * cell. The centering information, managed by the FieldEngine
 * class, is initialized using a flexible set of functors.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Domain/Shrink.h"
#include "Layout/INode.h"
#include "Layout/GuardLayers.h"
#include "Utilities/PAssert.h"
#include "Utilities/RefCountedBlockPtr.h"
#include "Engine/EnginePatch.h"
#include "Engine/EngineFunctor.h"
#include "Field/Relations/RelationList.h"
#include "Field/FieldCentering.h"
#include "Field/FieldEngine/FieldEnginePatch.h"

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag> class Engine;
template<class Components> class ComponentWrapper;

namespace Pooma {
  // Tag classes for taking sub-field views.
  struct MaterialViewTag {};
  struct CenteringViewTag {};
}


/**
 * FieldEngineBaseData holds an engine and the relations.
 */

template <int Dim, class T, class EngineTag>
class FieldEngineBaseData
{
public:

  FieldEngineBaseData()
    : engine_m()
  { }

  /// Initializer to be used with an engine compatible layout or
  /// similar initializer.

  template<class Initializer>
  FieldEngineBaseData(const Initializer &init)
    : engine_m(init)
  { }

  FieldEngineBaseData(const Pooma::NoInit &)
    : engine_m()
  { }

  template<class Initializer>
  FieldEngineBaseData(const Initializer &init, const RelationList &l)
    : engine_m(init),
      relations_m(l)
  { }

  template<class Engine, class Domain>
  FieldEngineBaseData(const Engine &e,
                      const Domain &d, const RelationList &l)
    : engine_m(NewEngineEngine<Engine, Domain>::apply(e, d), 
               NewEngineDomain<Engine, Domain>::apply(e, d)),
      relations_m(l)
  {
  }
    
  const Engine<Dim, T, EngineTag> &engine() const { return engine_m; }
  Engine<Dim, T, EngineTag> &engine() { return engine_m; }

  const RelationList &relations() const { return relations_m; }
  RelationList &relations() { return relations_m; }
  
private:

  Engine<Dim, T, EngineTag> engine_m;
  mutable RelationList relations_m;
};


/**
 * FieldEngineBase manages a hierarchy of engines, making it possible for
 * FieldEngine specializations to implement geometry-specific behavior only.
 */

template<class Mesh, class T, class EngineTag>
class FieldEngine
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.
    
  enum { dimensions = Mesh::dimensions };
  enum { Dim = dimensions };
  typedef FieldEngine<Mesh, T, EngineTag> This_t;
  typedef FieldEngineBaseData<Dim, T, EngineTag> Data_t;
  typedef Engine<Dim, T, EngineTag> Engine_t;
  typedef typename Engine_t::Domain_t Domain_t;
  typedef typename Engine_t::Layout_t Layout_t;
  typedef typename Engine_t::Element_t Element_t;
  typedef typename Engine_t::ElementRef_t ElementRef_t;
  typedef GuardLayers<Dim> GuardLayers_t;


  //---------------------------------------------------------------------------
  // Constructors.

  /// Default constructor.
  
  FieldEngine()
    : num_materials_m(0),
      physicalCellDomain_m(Pooma::NoInit()),
      guards_m(0)
  { }
  
  /// General version takes centering, layout, mesh, materials

  template<class Layout2>
  FieldEngine(const Centering<Dim> &centering, const Layout2 &layout,
              const Mesh &mesh, int materials = 1)
    : num_materials_m(materials),
      centering_m(centering),
      stride_m(centering.size()),
      physicalCellDomain_m(layout.domain()),
      guards_m(layout.externalGuards()),
      mesh_m(mesh)
  {
    shrinkInPlace(physicalCellDomain_m, guards_m);
    shrinkRightInPlace(physicalCellDomain_m, 1);
    addSubFields();
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centering.size(); ++ c)
      {
        data(m, c) = Data_t(layout);
      }
    }
  }

  /// Copy constructor.
  
  FieldEngine(const This_t &model)
    : num_materials_m(model.num_materials_m),
      centering_m(model.centering_m),
      stride_m(model.stride_m),
      data_m(model.data_m),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guards_m),
      mesh_m(model.mesh_m)
  {
  }

  ///@name Sub-field view constructors
  //@{

  /// Takes a view of
  ///  - the specified material including all centering points,
  ///    if there is more than one material
  ///  - the specified centering, if there is only one material
  /// These are weird semantics and thus this method is deprecated.
 
  FieldEngine(const This_t &model, int subField)
    : num_materials_m(1),
      stride_m(model.stride_m),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guards_m),
      mesh_m(model.mesh_m)
  {
    if (model.numMaterials() > 1)
    {
      centering_m = model.centering();
      data_m = model.data_m + model.stride_m * subField;
    }
    else
    {
      centering_m = model.centering()[subField];
      data_m = model.data_m + subField;
    }
  }

  /// Takes a view of the specified centering point of the specified material.

  FieldEngine(const This_t &model, int m, int c)
    : num_materials_m(1),
      centering_m(model.centering_m, c),
      stride_m(model.stride_m),
      data_m(model.data_m + model.stride_m * m + c),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guards_m),
      mesh_m(model.mesh_m)
  {
  }

  /// Takes a view of the specified centering point from all
  /// materials.

  FieldEngine(const This_t &model, int c, const Pooma::CenteringViewTag&)
    : num_materials_m(model.num_materials_m),
      centering_m(model.centering_m, c),
      stride_m(model.stride_m),
      data_m(model.data_m + c),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guards_m),
      mesh_m(model.mesh_m)
  {
  }

  /// Takes a view of the specified material retaining all centering points.

  FieldEngine(const This_t &model, int m, const Pooma::MaterialViewTag&)
    : num_materials_m(1),
      centering_m(model.centering_m),
      stride_m(model.stride_m),
      data_m(model.data_m + m * model.stride_m),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guards_m),
      mesh_m(model.mesh_m)
  {
  }

  /// Takes a view of the specified centering point of the first material.
  /// This is useless for fields with multiple materials and thus this
  /// method is deprecated. Use FieldEngine(field, 0, c).

  FieldEngine(int c, const This_t &model)
    : num_materials_m(1),
      centering_m(model.centering_m, c),
      stride_m(model.stride_m),
      data_m(model.data_m + c),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guards_m),
      mesh_m(model.mesh_m)
  {
  }

  //@}

  /// View constructors.  

  template<class T2, class EngineTag2>
  FieldEngine(const FieldEngine<Mesh, T2, EngineTag2> &model,
              const Domain_t &d)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      guards_m(0),
      mesh_m(model.mesh(),
             inputDomainToVertexDomain(d))
  {
    addSubFields();
    physicalCellDomain_m = d - d.firsts();
    if (centeringSize() == 1)
    {
      physicalCellDomain_m =
        centeringDomainToCellDomain(physicalCellDomain_m, centering_m, 0);
    }
    for (int m = 0; m < numMaterials(); ++m)
    {
      if (centeringSize() == 1)
      {
        data(m, 0) = Data_t(model.data(m, 0).engine(), d,
                            model.data(m, 0).relations());
      }
      else
      {
        for (int c = 0; c < centeringSize(); ++ c)
        {
          data(m, c) = Data_t(model.data(m, c).engine(),
                              cellDomainToCenteringDomain(d, centering_m, c),
                              model.data(m, c).relations());
        }
      }
    }
  }

  /// This constructor handle weird things like range views.

  template<class Mesh2, class T2, class EngineTag2, class Domain>
  FieldEngine(const FieldEngine<Mesh2, T2, EngineTag2> &model,
              const Domain &d)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      guards_m(0)
  {
    addSubFields();
    // FIXME: Does this ever happen to fields with multiple centering points?
    // (or event to fields with multiple materials???)
    PAssert(model.centeringSize() == 1);
    for (int m = 0; m < numMaterials(); ++m)
    {
      data(m, 0) = Data_t(model.data(m, 0).engine(), d,
                          model.data(m, 0).relations());
    }
    // FIXME: how do we construct the mesh?????
    mesh_m = Mesh(DomainLayout<Dim>(inputDomainToVertexDomain(data(0,0).engine().domain())));
    physicalCellDomain_m = mesh_m.physicalCellDomain();
  }

  template<class T2, class EngineTag2>
  FieldEngine(const FieldEngine<Mesh, T2, EngineTag2> &model,
              const INode<Dim> &i)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      guards_m(0),
      mesh_m(model.mesh(),
             inputDomainToVertexDomain(i.domain())) // FIXME: should hand INode to mesh?
  {
    addSubFields();
    physicalCellDomain_m = i.domain() - i.domain().firsts();
    if (centeringSize() == 1)
    {
      physicalCellDomain_m =
        centeringDomainToCellDomain(physicalCellDomain_m, centering_m, 0);
    }
    for (int m = 0; m < numMaterials(); ++m)
    {
      if (centeringSize() == 1)
      {
        data(m, 0) = Data_t(model.data(m, 0).engine(), i,
                            model.data(m, 0).relations());
      }
      else
      {
        for (int c = 0; c < centeringSize(); ++ c)
        { 
          data(m, c) =
            Data_t(model.data(m, c).engine(),
                   INode<Dim>(i, cellDomainToCenteringDomain(i.domain(),
                                                             centering_m, c)),
                   model.data(m, c).relations());
        }
      }
    }
  }

  template<class Mesh2, class T2, class EngineTag2, class Tag>
  FieldEngine(const FieldEngine<Mesh2, T2, EngineTag2> &model,
              const EngineView<Tag> &ev)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      physicalCellDomain_m(model.physicalCellDomain()),
      guards_m(model.guardLayers()),
      mesh_m(model.mesh())
  {
    typedef typename FieldEngine<Mesh2, T2, EngineTag2>::Engine_t EngIn_t;
    typedef LeafFunctor<EngIn_t, EngineView<Tag> > Functor_t;
    addSubFields();
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centeringSize(); ++ c)
      {
        data(m, c)
          = Data_t(Functor_t::apply(model.data(m, c).engine(), ev),
                   model.data(m, c).relations());
      }
    }
  }

  template<class EngineTag2>
  FieldEngine(const FieldEngine<Mesh, T, EngineTag2> &model,
              const FieldEnginePatch<Dim> &p)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      guards_m(model.guardLayers()),
      mesh_m(model.mesh()) // FIXME: should take a view of the mesh???
  {
    // FIXME: should we copy the relations for patch?  Do we want
    // to take patch views of composite fields?
    PAssert((model.numMaterials() == 1) && (model.centeringSize() == 1));
    addSubFields();
    data(0, 0) = Data_t(engineFunctor(model.engine(), EnginePatch(p.patch_m)));
    physicalCellDomain_m =
      centeringDomainToCellDomain(p.domain_m, centering_m, 0);
  }

  template<class Mesh2, class T2, class EngineTag2, class Components>
  FieldEngine(const FieldEngine<Mesh2, T2, EngineTag2> &model, 
              const ComponentWrapper<Components> &cw)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      physicalCellDomain_m(model.physicalCellDomain()),
      guards_m(model.guardLayers()),
      mesh_m(model.mesh())
  {
    addSubFields();
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centeringSize(); ++ c)
      {
        data(m, c) =
          Data_t(Engine_t(model.data(m, c).engine(), cw.components()),
                 model.data(m, c).relations());
      }
    }
  }

  FieldEngine(const This_t &model, 
              const Pooma::DontCopyRelations &d)
    : num_materials_m(model.numMaterials()),
      centering_m(model.centering()),
      stride_m(model.centeringSize()),
      physicalCellDomain_m(model.physicalCellDomain_m),
      guards_m(model.guardLayers()),
      mesh_m(model.mesh())
  {
    addSubFields();
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centeringSize(); ++ c)
      {
        data(m, c) = Data_t(model.data(m, c).engine());
      }
    }
  }
      
  //---------------------------------------------------------------------------
  // Initialize functions. 

  void initialize(const This_t &model)
  {
    num_materials_m = model.num_materials_m;
    stride_m = model.stride_m;
    centering_m = model.centering_m;
    data_m = model.data_m;
    physicalCellDomain_m = model.physicalCellDomain_m;
    guards_m = model.guards_m;
    mesh_m = model.mesh_m;
  }


  //---------------------------------------------------------------------------
  //@name Accessors and modifiers.
  //@{
    
  void addSubFields()
  {
    PAssert(data_m.size() == 0);

    int size = numMaterials() * centeringSize();

    data_m.reserve(size);
    data_m.resize(size);
  }

  /// FIXME: This function is deprecated.
  inline int numSubFields() const
  {
    if (numMaterials() > 1)
    {
      return numMaterials();
    }
    else
    {
      if (centering().size() > 1)
      {
        return centering().size();
      }
      else
      {
        return 0;
      }
    }
  }

  // FIXME: these should assert that there is 1 subfield.
  Engine_t &engine()
  {
    PAssert(data_m.isValid());
    return data_m->engine();
  }
    
  const Engine_t &engine() const
  {
    PAssert(data_m.isValid());
    return data_m->engine();
  }

  Engine_t &engine(int m, int c)
  {
    PAssert(data_m.isValid());
    return data(m,c).engine();
  }

  const Engine_t &engine(int m, int c) const
  {
    PAssert(data_m.isValid());
    return data(m,c).engine();
  }

  RelationList &relations() const
  {
    PAssert(data_m.isValid());
    return data_m->relations();
  }
    
  RelationList &relations(int m, int c) const
  {
    PAssert(data_m.isValid());
    return data(m, c).relations();
  }
    
  const GuardLayers_t &guardLayers() const
  {
    return guards_m;
  }

  GuardLayers_t &guardLayers()
  {
    return guards_m;
  }

  inline int numMaterials() const
  {
    return num_materials_m;
  }

  //@}

  //---------------------------------------------------------------------------
  //@name Domain accessor functions. 
  //@{
        
  /// The physical cell domain of all the sub-fields. Can be converted to
  /// the centering physical domain by means of cellDomainToCenteringDomain().

  Domain_t &physicalCellDomain()
  {
    return physicalCellDomain_m;
  }
        
  inline const Domain_t &physicalCellDomain() const
  {
    return physicalCellDomain_m;
  }
        
  inline Domain_t totalCellDomain() const
  {
    return grow(physicalCellDomain_m, guards_m);
  }

  /// Returns the physical domain suitable for viewing regardless of centering
  /// point count.

  Domain_t physicalDomain() const
  {
    if (centeringSize() == 1)
      return cellDomainToCenteringDomain(physicalCellDomain_m, centering_m, 0);
    else
      return physicalCellDomain_m;
  }

  /// Returns the physical domain of the specified centering.

  Domain_t physicalDomain(int i) const
  {
    return cellDomainToCenteringDomain(physicalCellDomain_m, centering_m, i);
  }

  /// Returns the total domain suitable for viewing regardless of centering
  /// point count.

  Domain_t totalDomain() const
  {
    if (centeringSize() == 1)
      return cellDomainToCenteringDomain(totalCellDomain(), centering_m, 0);
    else
      return totalCellDomain();
  }

  /// Returns the total domain of the specified centering.

  Domain_t totalDomain(int i) const
  {
    return cellDomainToCenteringDomain(totalCellDomain(), centering_m, i);
  }

  //@}

  //---------------------------------------------------------------------------
  //@name Centering accessors.
  //@{

  const Centering<Dim> &centering() const
  {
    return centering_m;
  }        

  inline int centeringSize() const
    {
      return centering_m.size();
    }

  //@}

  //---------------------------------------------------------------------------
  //@name Mesh accessors.
  //@{

  Mesh &mesh()
  {
    return mesh_m;
  }        

  const Mesh &mesh() const
  {
    return mesh_m;
  }        

  //@}

  //---------------------------------------------------------------------------
  /// Make a distinct copy of this fieldEngineBase.
 
  template<class Subject>
  void makeOwnCopy(const Subject &s)
  {
    PAssert(data_m.isValid());

    // Remember data_m as model

    RefCountedBlockPtr<Data_t> model = data_m;

    // Create a blank slate of engines:

    data_m = RefCountedBlockPtr<Data_t>();
    stride_m = centeringSize();
    addSubFields();

    // Copy then engines and relations and
    // Deepen the copies of the engine & relations list.
    
    for (int m = 0; m < numMaterials(); ++m)
    {
      for (int c = 0; c < centeringSize(); ++ c)
      {
        data(m, c) = model[m*stride_m + c];
        data(m, c).engine().makeOwnCopy();
        data(m, c).relations().makeOwnCopy(s.subField(m, c));
      }
    }
  }


  //---------------------------------------------------------------------------
  /// Domain translation function.
  /// FIXME: This function should go away.  Currently it's only used by
  /// the lagrangian field engine.

  inline Domain_t
  translateToVertexDomain(const Domain_t &d) const
  {
    if (centeringSize() == 1)
      return d;
    else
      return growRight(d, 1);
  }

  /// Converts an input domain (which is a cell domain for fields with
  /// multiple centering points and a centering domain for one centering point)
  /// to the corresponding vertex domain.

  Domain_t
  inputDomainToVertexDomain(const Domain_t &d) const
  {
    if (centeringSize() == 1)
      return growRight(centeringDomainToCellDomain(d, centering(), 0), 1);
    else
      return growRight(d, 1);
  }

  //---------------------------------------------------------------------------
  //@name Access material, centering subfield data.
  //@{
  
  inline Data_t &
  data(int material, int centering)
  {
    PAssert(data_m.isValid());
    return data_m[material * stride_m + centering];
  }

  inline const Data_t &
  data(int material, int centering) const
  {
    PAssert(data_m.isValid());
    return data_m[material * stride_m + centering];
  }

  //@}      
      
private:

  unsigned int num_materials_m;
  Centering<Dim> centering_m;
  int stride_m;
  RefCountedBlockPtr<Data_t> data_m;

  /// The physical cell domain of all the sub-fields. Can be converted to
  /// the centering physical domain by means of cellDomainToCenteringDomain().
  Domain_t physicalCellDomain_m;
  GuardLayers_t guards_m;

  Mesh mesh_m;
};

template<class Mesh, class T, class EngineTag, class Tag>
struct LeafFunctor<FieldEngine<Mesh, T, EngineTag>,
                   ExpressionApply<Tag> >
{
  typedef FieldEngine<Mesh, T, EngineTag> Subject_t;
  typedef typename Subject_t::Engine_t Engine_t;
  typedef LeafFunctor<Engine_t, ExpressionApply<Tag> > LeafFunctor_t;
  typedef int Type_t;

  inline static
  Type_t apply(const Subject_t &fieldEngineBase, 
	       const ExpressionApply<Tag> &tag)
  {
    for (int m = 0; m < fieldEngineBase.numMaterials(); ++m)
    {
      for (int c = 0; c < fieldEngineBase.centering().size(); ++ c)
      {
        LeafFunctor_t::apply(fieldEngineBase.data(m, c).engine(), tag);
      }
    }
    return 0;
  }
};

#endif // POOMA_FIELD_FIELDENGINE_FIELDENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldEngine.h,v $   $Author: richi $
// $Revision: 1.11 $   $Date: 2004/11/10 22:05:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
