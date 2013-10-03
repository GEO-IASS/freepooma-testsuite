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
// Mesh Test 1: mesh constructors and accessors.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

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

const int NX = 8, NY = 12;
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
// Confirm that basic ctors and accessors are correct.
//-----------------------------------------------------------------------------

// Uniform rectilinear mesh.

void urmTest(Pooma::Tester &tester)
{
  // Create a uniform rectilinear mesh using a DomainLayout and test.

  DomainLayout<2> layout1(physicalVertexDomain, gl);
  tester.out() << layout1 << std::endl;
  UniformRectilinearMesh<2> m(layout1, origin, spacings);

  tester.check("URM.PVD", m.physicalVertexDomain(), physicalVertexDomain);
  tester.check("URM.TVD", m.totalVertexDomain(), totalVertexDomain);
  tester.check("URM.PCD", m.physicalCellDomain(), physicalCellDomain);
  tester.check("URM.TCD", m.totalCellDomain(), totalCellDomain);
  
  tester.check("URM.Origin", m.origin(), origin);
  tester.check("URM.Spacings", m.spacings(), spacings);

  // Test a view.
    
  UniformRectilinearMesh<2> m2(m, viewDomain);

  tester.check("V.URM.PVD", m2.physicalVertexDomain(), viewPhysVertexDomain);
  tester.check("V.URM.TVD", m2.totalVertexDomain(), viewPhysVertexDomain);
  tester.check("V.URM.PCD", m2.physicalCellDomain(), viewPhysCellDomain);
  tester.check("V.URM.TCD", m2.totalCellDomain(), viewPhysCellDomain);
  
  tester.check("V.URM.Origin", m2.origin(), viewOrigin);
  tester.check("V.URM.Spacings", m2.spacings(), spacings);
}

// No mesh

void nmTest(Pooma::Tester &tester)
{
  // Create a no-mesh using a DomainLayout and test.

  DomainLayout<2> layout1(physicalVertexDomain, gl);
  tester.out() << layout1 << std::endl;
  NoMesh<2> m(layout1);

  tester.check("NM.PVD", m.physicalVertexDomain(), physicalVertexDomain);
  tester.check("NM.TVD", m.totalVertexDomain(), totalVertexDomain);
  tester.check("NM.PCD", m.physicalCellDomain(), physicalCellDomain);
  tester.check("NM.TCD", m.totalCellDomain(), totalCellDomain);

  // Test a view.
    
  NoMesh<2> m2(m, viewDomain);

  tester.check("V.NM.PVD", m2.physicalVertexDomain(), viewPhysVertexDomain);
  tester.check("V.NM.TVD", m2.totalVertexDomain(), viewPhysVertexDomain);
  tester.check("V.NM.PCD", m2.physicalCellDomain(), viewPhysCellDomain);
  tester.check("V.NM.TCD", m2.totalCellDomain(), viewPhysCellDomain);
}

//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  urmTest(tester);
  nmTest(tester);
    
  int ret = tester.results("MeshTest1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MeshTest1.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
