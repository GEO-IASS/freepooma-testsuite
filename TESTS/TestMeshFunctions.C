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
// Test Mesh functors.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Arrays.h"
#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

//-----------------------------------------------------------------------------
// Globals
//-----------------------------------------------------------------------------

// Meshes are initialized with vertex-based PHYSICAL domains. The total domain
// should be the physical domain, expanded by the guard layers in each.
// The physical and total cell domains are shrunk by 1 on the right. When
// taking a view, the physical and total domains should be zero-based and 
// the same. Again, the physical and total cell domains are shrunk by 1 on 
// the right.

const int NX = 4, NY = 4;
GuardLayers<2> gl(1);
Interval<1> I(NX), J(NY);
Interval<2> physicalVertexDomain(I, J);
Vector<2> origin(0.0), spacings(1.0, 2.0);
  

//-----------------------------------------------------------------------------
// Test positions function.
//-----------------------------------------------------------------------------

template<class Field>
void testPositions(Pooma::Tester &tester, const Field &f)
{   
  tester.out() << positions(f) << std::endl;
}  

//-----------------------------------------------------------------------------
// Test normals functor.
//-----------------------------------------------------------------------------

template<class Field>
void testNormals(Pooma::Tester &tester, const Field &f)
{   
  tester.out() << outwardNormals(f) << std::endl;
  tester.out() << coordinateNormals(f).all() << std::endl;
}  

//-----------------------------------------------------------------------------
// Test cell volumes functor.
//-----------------------------------------------------------------------------

template<class Field>
void testCellVolumes(Pooma::Tester &tester, const Field &f)
{   
  tester.out() << cellVolumes(f) << std::endl;
}  

//-----------------------------------------------------------------------------
// Test face areas functor.
//-----------------------------------------------------------------------------

template<class Field>
void testFaceAreas(Pooma::Tester &tester, const Field &f)
{   
  tester.out() << faceAreas(f) << std::endl;
}  

//-----------------------------------------------------------------------------
// Test face areas functor.
//-----------------------------------------------------------------------------

template<class Field>
void testEdgeLengths(Pooma::Tester &tester, const Field &f)
{   
  tester.out() << edgeLengths(f) << std::endl;
}  

template<class Mesh>
void test(Pooma::Tester &tester)
{
  // Create a mesh using a DomainLayout.

  DomainLayout<2> layout1(physicalVertexDomain, gl);
  tester.out() << layout1 << std::endl;
  Mesh mesh1(layout1, origin, spacings);
  
  // Set up some centerings.
  
  Centering<2> cell = canonicalCentering<2>(CellType, Continuous);
  
  // Initialize a field.
  
  Field<Mesh> f1(cell, layout1, mesh1);

  // Do the tests.
  
  testPositions(tester, f1);
  testNormals(tester, f1);
  testCellVolumes(tester, f1);
  testFaceAreas(tester, f1);
  testEdgeLengths(tester, f1);
}

//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  // Test for uniform rectilinear mesh.
  test<UniformRectilinearMesh<2> >(tester);

  // Test for rectilinear mesh.
  test<RectilinearMesh<2> >(tester);
    
  int ret = tester.results("TestMeshFunctions");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestMeshFunctions.cpp,v $   $Author: richi $
// $Revision: 1.3 $   $Date: 2004/11/22 16:54:46 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
