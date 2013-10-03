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
//   FieldStencilSimple    - A wrapper class for a user-defined stencil.
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_DIFFOPS_FIELDSTENCIL_H
#define POOMA_FIELD_DIFFOPS_FIELDSTENCIL_H

/** @file
 * @ingroup DiffOps
 * @brief
 * This file contains the equipment required to write differential operators
 * that take the form of stencil objects using Fields.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Engine/Engine.h"
#include "Engine/Stencil.h"
#include "Layout/INode.h"
#include "Layout/Node.h"
#include "PETE/ErrorType.h"
#include "Field/FieldOffset.h"


/**
 * There are potentially many ways to construct field stencils.
 * FieldStencilSimple assumes that you just need to construct the output field
 * and stick ONE stencil engine into it.  Maybe this class can be generalized
 * for fields that contain multiple stencil engines.
 *
 * FieldStencil is used to wrap a user-defined field-based stencil class.
 * The idea is to encapsulate the majority of the crazy type manipulations
 * required to generate the output Field.
 *
 * To create a stencil, users must create a class similar to the one below,
 * which computes a central difference divergence of a vertex-centered Field
 * and maps it to a cell-centered Field:
 *
 * <PRE>
 * template<class T2, int Dim, class TM>
 * class DivVertToCell<Vector<Dim, T2>, UniformRectilinearMesh<Dim, TM> >
 * {
 * public:
 *  
 * typedef T2   OutputElement_t;
 *     
 * Centering<Dim> outputCentering() const 
 * {
 *   return canonicalCentering<Dim>(CellType, Continuous, AllDim);
 * }
 *
 * Centering<Dim> inputCentering() const 
 * {
 *   return canonicalCentering<Dim>(VertexType, Continuous, AllDim);
 * }
 *                           
 * // Constructors.
 *
 * // default version is required by default stencil engine constructor.
 *
 * DivVertToCell()
 * {
 *   for (int d = 0; d < Dim; ++d)
 *   {
 *      fact_m(d) = 1.0;
 *   }
 * }
 *
 * template<class FE>
 * DivVertToCell(const FE &fieldEngine)
 * {
 *   for (int d = 0; d < Dim; ++d)
 *   {
 *      fact_m(d) = 1 / fieldEngine.mesh().spacings()(d);
 *   }
 * }
 *
 * // Methods.
 *
 * int lowerExtent(int d) const { return 0; }
 * int upperExtent(int d) const { return 1; }
 *
 * template<class F>
 * inline OutputElement_t
 * operator()(const F &f, int i1) const
 * {
 *   return OutputElement_t
 *     (fact_m(0)*(f.read(i1+1)(0) - f.read(i1)(0)));
 * }
 *
 * // and versions for 2d and 3d
 * 
 * private:
 * Vector<Dim, TM> fact_m;
 * };
 * </PRE>
 *
 * There is one required typedefs: OutputElement_t. 
 * These export the type of the type resulting 
 * from applying the stencil at a point. 
 *
 * There are two required methods returning the input and
 * output centering.
 *
 * Then, there are two accessors: lowerExtent(int dir) and 
 * upperExtent(int dir). These return the extent of the stencil as a function 
 * of direction. As another example, a forward difference would have a lower
 * extent of 0 and an upper extent of 1. Finally, a series of inline apply()
 * functions, which take a Field of some sort and a set indices, must be
 * supplied. This is what actually computes the stencil.
 *
 * A Field that contains a StencilEngine that operates on
 * a Field f, is constructed by using make() from FieldStencilSimple:
 *
 * FieldStencilSimple<DivVertToCell<T, Mesh>, Field<Mesh, T, EngineTag> >
 *   ::make(DivVertToCell<T, Mesh>(f.fieldEngine()), f);
 */ 

template<class Functor, class Expression>
struct FieldStencilSimple
{
  typedef typename Expression::MeshTag_t MeshTag_t;
  enum { outputDim = Expression::dimensions };

  typedef typename Functor::OutputElement_t OutputElement_t;

  typedef StencilEngine<Functor, Expression> OutputEngineTag_t;
  typedef Field<MeshTag_t, OutputElement_t, OutputEngineTag_t> Type_t;

  typedef Engine<outputDim, OutputElement_t, OutputEngineTag_t> SEngine_t;

  static inline
  Type_t make(const Functor &stencil, const Expression &f)
  {
	// FIXME: need to add comparison for centerings.
	//    PAssert(f.centering() == stencil.inputCentering());

	// We need to use the centering, layout, mesh constructor.
	// The FieldEngine part initializes physicalCellDomain
	// and guards from the layout.

	Type_t h(stencil.outputCentering(), f.layout(), f.mesh());

	// Initialize engine with appropriate StencilEngine

	Interval<outputDim> domain = insetDomain(stencil, f.physicalDomain());
	h.fieldEngine().engine() = SEngine_t(stencil, f, domain);

	return h;
  }

  static inline
  Type_t make(const Functor &stencil, const Expression &f, const Interval<outputDim> &domain)
  {
	// FIXME: need to add comparison for centerings.
	//    PAssert(f.centering() == stencil.inputCentering());

	// We need to use the centering, layout, mesh constructor.
	// The FieldEngine part initializes physicalCellDomain
	// and guards from the layout.

	Type_t h(stencil.outputCentering(), f.layout(), f.mesh());

	// Initialize engine with appropriate StencilEngine

	h.fieldEngine().engine() = SEngine_t(stencil, f, domain);

	return h;
  }

  template<class Accumulate>
  static inline
  Type_t make(const Expression &f,
              const std::vector<FieldOffsetList<outputDim> > &nn,
              const Centering<outputDim> &outputCentering,
              Accumulate accumulate = Accumulate())
  {
    PAssert(nn.size() == outputCentering.size());

    Type_t h(outputCentering, f.layout(), f.mesh());
    h.fieldEngine().physicalCellDomain() = f.fieldEngine().physicalCellDomain();

    // FIXME: The guard layers are wrong; we need to find the maximum
    // offsets from all the functors below.  (Should the individual
    // sub-fields have their own guard layers???)

    h.fieldEngine().guardLayers() = f.fieldEngine().guardLayers();

    if (outputCentering.size() == 1)
    {
      h.fieldEngine().engine()
        = SEngine_t(Functor(nn[0], outputCentering, f.centering(),
                            accumulate),
                    f, h.physicalDomain());
    }
    else
    {
      int oc;

      for (oc = 0; oc < nn.size(); ++oc)
      {
        h[oc].fieldEngine().guardLayers() = f.fieldEngine().guardLayers();
        h[oc].fieldEngine().engine()
          = SEngine_t(Functor(nn[oc], outputCentering[oc], f.centering(),
                              accumulate),
                      f, h[oc].physicalDomain());
      }
    }

    return h;
  }
};


#endif // POOMA_FIELD_DIFFOPS_FIELDSTENCIL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldStencil.h,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:44 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
