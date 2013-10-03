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
// 
// FieldOffsetReduction
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_DIFFOPS_FIELDOFFSETREDUCTION_H
#define POOMA_FIELD_DIFFOPS_FIELDOFFSETREDUCTION_H

//-----------------------------------------------------------------------------
// Overview: 
//
// Support for reductions that take the result of a nearest neighbor
// calculation.  For example:
//
// sum(f, nearestNeighbor(f.centering(), centering2), centering2);
//
// results in a field with centering2, where values from f at the closest
// centering points are summed onto each centering point in the resulting
// field.
//
// Note that this function isn't really that general, because the nearest
// neighbor computation must use the input fields centering for the input
// centering and same output centering as the centering handed to sum.
// 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Field/DiffOps/FieldStencil.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
//
// FieldOffsetReduction<T, Dim, Accumulate>
//
// This class is a functor object used to create a field stencil engine that
// accumulates values from a FieldOffsetList.  The functor is constructed with
// a neighbor list, output centering, input centering and a binary function.
//
// Applying the functor to a field and location will accumulate the neighbors
// of the field at the location using the neighbor list.
//
// ----------------------------------------------------------------------------

template<class T, int Dim, class Accumulate>
class FieldOffsetReduction
{
public:

  //
  // Interface required to create a field stencil object:
  //

  typedef T   OutputElement_t;

  const Centering<Dim> &outputCentering() const
  {
    return outputCentering_m;
  }

  const Centering<Dim> &inputCentering() const
  {
    return inputCentering_m;
  }

  int lowerExtent(int d) const { return lower_m[d]; }
  int upperExtent(int d) const { return upper_m[d]; }
      
  // 
  // Constructors.
  // 

  FieldOffsetReduction()
  {
  }

  FieldOffsetReduction(const FieldOffsetList<Dim> &neighbors,
                       const Centering<Dim> &outputCentering,
                       const Centering<Dim> &inputCentering,
                       Accumulate accumulate = Accumulate())
    : neighbors_m(neighbors),
      outputCentering_m(outputCentering),
      inputCentering_m(inputCentering),
      accumulate_m(accumulate)
  {
    PInsist(neighbors.size() > 0, "no support for empty accumulation");
    PAssert(outputCentering.size() == 1);

    int d, i;

    for (d = 0; d < Dim; ++d)
    {
      upper_m[d] = 0;
      lower_m[d] = 0;
    }
    for (i = 0; i < neighbors_m.size(); ++i)
    {
      for (d = 0; d < Dim; ++d)
      {
        if (-neighbors_m[i].cellOffset()[d].first() > lower_m[d])
        {
          lower_m[d] = -neighbors_m[i].cellOffset()[d].first();
        }
        if (neighbors_m[i].cellOffset()[d].first() > upper_m[d])
        {
          upper_m[d] = neighbors_m[i].cellOffset()[d].first();
        }
      }
    }
  }

  //
  // Stencil operations.
  //

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1) const
  {
    T ret = f(neighbors_m[0], Loc<1>(i1));

    for (int i = 1; i < neighbors_m.size(); ++i)
    {
      ret = accumulate_m(ret, f(neighbors_m[i], Loc<1>(i1)));
    }
    return ret;
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2) const
  {
    T ret = f(neighbors_m[0], Loc<2>(i1, i2));

    for (int i = 1; i < neighbors_m.size(); ++i)
    {
      ret = accumulate_m(ret, f(neighbors_m[i], Loc<2>(i1, i2)));
    }
    return ret;
  }

  template<class F>
  inline OutputElement_t
  operator()(const F &f, int i1, int i2, int i3) const
  {
    T ret = f(neighbors_m[0], Loc<3>(i1, i2, i3));

    for (int i = 1; i < neighbors_m.size(); ++i)
    {
      ret = accumulate_m(ret, f(neighbors_m[i], Loc<3>(i1, i2, i3)));
    }
    return ret;
  }

private:

  FieldOffsetList<Dim> neighbors_m;
  Centering<Dim> outputCentering_m;
  Centering<Dim> inputCentering_m;
  Accumulate accumulate_m;
  int lower_m[Dim];
  int upper_m[Dim];
};


// ----------------------------------------------------------------------------
// Overloaded functions to perform reductions.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Sum up the results from a nearest neighbor computation
// ----------------------------------------------------------------------------

template<class GeometryTag, class T, class EngineTag, int Dim>
typename FieldStencilSimple<FieldOffsetReduction<T, Dim, OpAdd>,
                            Field<GeometryTag, T, EngineTag> >::Type_t
sum(const Field<GeometryTag, T, EngineTag> &f,
    const std::vector<FieldOffsetList<Dim> > &nn,
    const Centering<Dim> &outputCentering)
{
  typedef FieldOffsetReduction<T, Dim, OpAdd> Functor_t;
  typedef Field<GeometryTag, T, EngineTag> Field_t;

  return FieldStencilSimple<Functor_t, Field_t>::make(f, nn, outputCentering,
                                                      OpAdd());
}


#endif     // POOMA_FIELD_DIFFOPS_FIELDOFFSETREDUCTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldOffsetReduction.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:44 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
