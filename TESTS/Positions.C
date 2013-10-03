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
// Positions: illustrates the xField function which returns position values.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int NX = 5, NY = 5;
  Interval<1> I(NX), J(NY);
  Interval<2> physicalVertexDomain(I, J),
    td(Interval<1>(-1,NX-1), Interval<1>(-1,NY-1)), pd(NX-1, NY-1);
  
  // Create the mesh.
  
  Vector<2,double> origin;
  Vector<2,double> spacings;
  for (int d = 0; d < 2; d++) 
  {
    origin(d) = d;
    spacings(d) = d + 1;
  }
      
  // Make a Brick-Engine-based field.

  DomainLayout<2> layout1(physicalVertexDomain, GuardLayers<2>(1));

  Centering<2> cell = canonicalCentering<2>(CellType, Continuous);
  Centering<2> vert = canonicalCentering<2>(VertexType, Continuous);
  Centering<2> allFace = canonicalCentering<2>(FaceType, Continuous);

  typedef UniformRectilinearMesh<2> Mesh_t;
  typedef Field<Mesh_t, double, Brick> Field_t;

  typedef XField<Mesh_t>::Type_t XField_t;
  
  Field_t f(cell, layout1, origin, spacings);
  XField_t x(cell, layout1, origin, spacings);
  setXField(x);

  f = x.comp(0);

  tester.out() << f << std::endl;

  f = x.comp(1);

  tester.out() << f << std::endl;

  tester.out() << xField(f, vert) << std::endl;
  tester.out() << xField(f, allFace) << std::endl;

  // This statement fails because component view isn't implemented
  // for field expressions.

  //  tester.out() << (2 + x).comp(0) << std::endl;

  int ret = tester.results("Positions");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Positions.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
