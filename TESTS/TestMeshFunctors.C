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
GuardLayers<2> gl(Loc<2>(1, 2), Loc<2>(2, 1));

Interval<1> I(NX), J(NY), IV(NX - 2), JV(NY - 1);
Interval<2> 
  physicalVertexDomain(I, J),
  totalVertexDomain(Interval<1>(-gl.lower(0), NX+gl.upper(0)-1), 
    Interval<1>(-gl.lower(1),NY+gl.upper(1)-1)), 
  physicalCellDomain(shrinkRight(physicalVertexDomain, 1)),
  totalCellDomain(shrinkRight(totalVertexDomain, 1));
Interval<2>
  viewDomain(IV + 1, JV - 1),
  viewPhysVertexDomain(IV, JV), 
  viewPhysCellDomain(shrinkRight(viewPhysVertexDomain, 1));

Vector<2> origin(0.0), spacings(1.0, 2.0), viewOrigin(1.0, -2.0);
  

//-----------------------------------------------------------------------------
// Test positions functor.
//-----------------------------------------------------------------------------

template<class Mesh>
void testPositions(Pooma::Tester &tester, const Mesh &m)
{   
  Array<Mesh::dimensions, 
        typename Mesh::PointType_t, 
        typename Mesh::PositionsEngineTag_t> 
    a(m.physicalCellDomain());

  m.initializePositions(a.engine(), 
    canonicalCentering<Mesh::dimensions>(CellType, Continuous));
  tester.out() << a << std::endl;
}  

//-----------------------------------------------------------------------------
// Test normals functor.
//-----------------------------------------------------------------------------

template<class Mesh>
void testNormals(Pooma::Tester &tester, const Mesh &m, bool outward)
{   
  Array<Mesh::dimensions, 
        typename Mesh::VectorType_t, 
        typename Mesh::NormalsEngineTag_t> 
    a(m.physicalCellDomain());

  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(FaceType, Discontinuous);
  tester.out() << c << std::endl;
  
  m.initializeNormals(a.engine(), c[2], outward);
  tester.out() << a << std::endl;
}

//-----------------------------------------------------------------------------
// Test cell volumes functor.
//-----------------------------------------------------------------------------

template<class Mesh>
void testCellVolumes(Pooma::Tester &tester, const Mesh &m)
{   
  Array<Mesh::dimensions, 
        typename Mesh::T_t, 
        typename Mesh::CellVolumesEngineTag_t> 
    a(m.physicalCellDomain());

  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(CellType, Continuous);
  tester.out() << c << std::endl;

  m.initializeCellVolumes(a.engine(), c);
  tester.out() << a << std::endl;
}

//-----------------------------------------------------------------------------
// Test face areas functor.
//-----------------------------------------------------------------------------

template<class Mesh>
void testFaceAreas(Pooma::Tester &tester, const Mesh &m)
{   
  Array<Mesh::dimensions, 
        typename Mesh::T_t, 
        typename Mesh::FaceAreasEngineTag_t> 
    a(m.physicalCellDomain());

  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(FaceType, Continuous);
  tester.out() << c << std::endl;

  m.initializeFaceAreas(a.engine(), c[0]);
  tester.out() << a << std::endl;
}

//-----------------------------------------------------------------------------
// Test face areas functor.
//-----------------------------------------------------------------------------

template<class Mesh>
void testEdgeLengths(Pooma::Tester &tester, const Mesh &m)
{   
  Array<Mesh::dimensions, 
        typename Mesh::T_t, 
        typename Mesh::EdgeLengthsEngineTag_t> 
    a(m.physicalCellDomain());

  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(EdgeType, Continuous);
  tester.out() << c << std::endl;

  m.initializeEdgeLengths(a.engine(), c[0]);
  tester.out() << a << std::endl;
}


template <class Mesh>
void test(Pooma::Tester &tester)
{
  // Create a mesh using a DomainLayout and test.

  DomainLayout<2> layout1(physicalVertexDomain, gl);
  tester.out() << layout1 << std::endl;
  Mesh m1(layout1, origin, spacings);

  testPositions(tester, m1);
  testNormals(tester, m1, true);
  testNormals(tester, m1, false);
  testCellVolumes(tester, m1);
  testFaceAreas(tester, m1);
  testEdgeLengths(tester, m1);
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
    
  int ret = tester.results("TestMeshFunctors");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestMeshFunctors.cpp,v $   $Author: richi $
// $Revision: 1.3 $   $Date: 2004/11/22 16:54:46 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
