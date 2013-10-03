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
// A 2nd tour of the new Field class featuring multi-material fields.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  // To declare a field, you first need to set up a layout. This requires
  // knowing the physical vertex-domain and the number of external guard
  // cell layers. Vertex domains contain enough points to hold all of the
  // rectilinear centerings that POOMA is likely to support for quite
  // awhile. Also, it means that the same layout can be used for all
  // fields, regardless of centering.
  
  Interval<2> physicalVertexDomain(4, 4);
  DomainLayout<2> layout(physicalVertexDomain, GuardLayers<2>(1));
  
  // Now, we can declare a field. Let's make this a multi-material
  // field with cell centering.

  Centering<2> cell = canonicalCentering<2>(CellType, Continuous, AllDim);

  typedef Field<UniformRectilinearMesh<2>, double, Brick> Field_t;
  Field_t f(3, cell, layout, // 3 fields each Cell-centered
    Vector<2>(0.0), Vector<2>(1.0, 2.0));

  // Set some data in the field.
  
  f[0].all() = 2.0; f[0] = -1.0; 
  f[1].all() = 3.0; f[1] = -2.0; 
  f[2].all() = 4.0; f[2] = -3.0; 
  
  std::cout << f.all() << std::endl;

  // Try some reductions.

  std::cout << sum(f[0]) << std::endl;
  std::cout << sum(f[1] + f[2]) << std::endl;
  
  // Take a range-based view. Note: the only views allowed for
  // fields with sub-fields are those constructed using Intervals
  // and INodes. The reason is that a Range of cells can lead
  // to non-constant strides through the sub-field elements. One
  // can construct Range-based views of Fields with no subfields.
  // The result is a field with a NoGeometry, GeometryTag.
  
  Range<1> r(-1, 3, 2);
  Range<2> R(r, r);
  std::cout << f[2](R) << std::endl;

  Pooma::finalize();
  return 0; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldTour2.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
