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
// Test of the new Centerings class.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Field/FieldCentering.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  Centering<2> edges
    = canonicalCentering<2>(EdgeType, Continuous, XDim | YDim);

  std::cout << edges << std::endl;

  Interval<2> physicalVertexDomain(4, 4);
  DomainLayout<2> layout(physicalVertexDomain, GuardLayers<2>(1));
  typedef Field<UniformRectilinearMesh<2>, double, Brick> Field_t;

  // Create a field with edge-centered values for the x- and y-directions.
  Field_t f(edges, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));

  // Create a 3-material field with edge-centered values for the x-
  // and y-directions.
  Field_t g(3, edges, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));

  // Set some data in the field.
  
  f[0].all() = 2.0; f[0] = -1.0; 
  f[1].all() = 3.0; f[1] = -2.0; 
  
  std::cout << f.all() << std::endl;

  Pooma::finalize();
  return 0; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldTour3.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
