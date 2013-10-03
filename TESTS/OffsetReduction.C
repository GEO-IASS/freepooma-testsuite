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

#include "Field/DiffOps/FieldOffsetReduction.h"

#include <iostream>
#include <cmath>

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
#endif

const int dim = 2;

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
  
  Centering<dim> cell = canonicalCentering<dim>(CellType, Continuous, AllDim);
  Centering<dim> vertex = canonicalCentering<dim>(VertexType, Continuous, AllDim);
  Centering<dim> discVertex = canonicalCentering<dim>(VertexType, Discontinuous, AllDim);

  typedef UniformRectilinearMesh<dim> Geometry_t;

  typedef 
    Field<Geometry_t, double, MultiPatch<UniformTag, BrickTag_t> >
    Field_t;

  typedef 
    Field<Geometry_t, Vector<dim>, MultiPatch<UniformTag, BrickTag_t> >
    VField_t;

  Vector<dim> origin(0.0);
  Vector<dim> spacings(1.0);

  Field_t cfield(cell, layout, origin, spacings);
  Field_t r1(vertex, layout, origin, spacings);
  Field_t r2(discVertex, layout, origin, spacings);

  cfield.all() = iota(cfield.all().domain()).comp(0);
  for (int i = 1; i < dim; ++i)
  {
    cfield.all() *= iota(cfield.all().domain()).comp(i);
  }

  r1 = sum(cfield, nearestNeighbors(cfield.centering(), r1.centering()), r1.centering());
  r2 = sum(cfield, nearestNeighbors(cfield.centering(), r2.centering()), r2.centering());

  tester.out() << "input field" << std::endl
               << cfield.all() << std::endl;

  tester.out() << "r1" << std::endl
               << r1.all() << std::endl;

  tester.out() << "r2" << std::endl
               << r2.all() << std::endl;

  int ret = tester.results("OffsetReduction");
  Pooma::finalize();
  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: OffsetReduction.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
