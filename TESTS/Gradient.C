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
// An example computing the gradient using field-offsets and the neighbor
// function.
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/Fields.h"

#include <iostream>
#include <cmath>

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
#endif

template<class Geom, class T, class Eng, int Dim>
Field<Geom, Vector<Dim, T>, Eng>
gradient(const Field<Geom, T, Eng> &input,
         const Centering<Dim> &outputCentering,
	 Pooma::Tester &tester)
{
  // Need a centering, layout, mesh constructor:
  Field<Geom, Vector<Dim, T>, Eng>
    ret(outputCentering, input.engine().layout(),
	Vector<Dim>(0.0), Vector<Dim>(1.0));

  // Just do the 1 material - 1 centering point case for now:
  PAssert(input.numMaterials() <= 1);
  PAssert(input.centering().size() == 1);
  PAssert(outputCentering.size() == 1);

  std::vector<FieldOffsetList<Dim> > nn =
    nearestNeighbors(input.centering(), outputCentering);

  PAssert(nn.size() == 1);
  FieldOffsetList<Dim> offsets = nn[0];
  const int points = offsets.size();

  std::vector<Vector<Dim, double> > coeff;
  coeff.reserve(points);
  coeff.resize(points);

  const Vector<Dim, double> outputLocation = outputCentering.position(0);
  Vector<Dim, double> norm(0.0);

  int i, j;
  for (i = 0; i < points; ++i)
  {
    coeff[i] =
      inputPosition(input.centering(), offsets[i]) - outputLocation;

    for (j = 0; j < Dim; ++j)
    {
      norm(j) += coeff[i](j) * coeff[i](j);
    }
  }

  tester.out() << "Coefficients:" << std::endl;

  for (j = 0; j < Dim; ++j)
  {
    PInsist(norm(j) > 0.0,
	    "The gradient's neighborhood doesn't have enough points.");
  }
  for (i = 0; i < points; ++i)
  {
    for (j = 0; j < Dim; ++j)
    {
      coeff[i](j) /= norm(j);
    }
    tester.out() << i << ":" << coeff[i] << " for offset: " << offsets[i] << std::endl;
  }

  ret = Vector<Dim, double>(0.0);

  for (i = 0; i < points; ++i)
  {
    ret += coeff[i] * input(offsets[i], outputCentering);
  }

  return ret;
}


const int dim = 3;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<dim> physicalVertexDomain;

  for(int i = 0; i < dim; ++i)
    physicalVertexDomain[i] = Interval<1>(10);
  
  Loc<dim> blocks(2);

  UniformGridPartition<dim> partition(blocks, GuardLayers<dim>(1));   
  UniformGridLayout<dim> layout(physicalVertexDomain, partition,
                                LayoutTag_t());
  
  Centering<dim> cell = canonicalCentering<dim>(CellType, Continuous);
  Centering<dim> vertex = canonicalCentering<dim>(VertexType, Continuous);
  Centering<dim> discVertex = canonicalCentering<dim>(VertexType, Discontinuous);

  typedef UniformRectilinearMesh<dim> Geometry_t;

  typedef 
    Field<Geometry_t, double, MultiPatch<UniformTag, BrickTag_t> >
    Field_t;

  typedef 
    Field<Geometry_t, Vector<dim>, MultiPatch<UniformTag, BrickTag_t> >
    VField_t;

  Vector<dim> origin(0.0);
  Vector<dim> spacings(1.0);

  VField_t vfield(vertex, layout, origin, spacings);
  Field_t cfield(cell, layout, origin, spacings);

  cfield.all() = iota(cfield.all().domain()).comp(0);
  for (int i = 1; i < dim; ++i)
  {
    cfield.all() *= iota(cfield.all().domain()).comp(i);
  }

  vfield = gradient(cfield, vertex, tester);

  tester.out() << "input field" << std::endl
               << cfield.all() << std::endl;

  tester.out() << "output field" << std::endl
               << vfield.all() << std::endl;

  int ret = tester.results("Gradient");
  Pooma::finalize();
  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Gradient.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
