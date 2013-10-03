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

// ----------------------------------------------------------------------------
// Various tests of det(Tensor<>) global function (determinant).
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

// Forward declarations:
template<int Dim>
void testDet(Pooma::Tester &tester);


//-----------------------------------------------------------------------------
// Main program:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  testDet<1>(tester);
  testDet<2>(tester);
  testDet<3>(tester);

  int ret = tester.results("TestDet");
  Pooma::finalize();
  return ret;
}

template<int Dim>
void testDet(Pooma::Tester &tester)
{
  // Create the physical Domains:
  const int nVerts = 6;
  const int nCells = nVerts - 1;
  int nCellsTot = 1;
  Interval<Dim> vertexDomain;
  for (int d = 0; d < Dim; d++) {
    vertexDomain[d] = Interval<1>(nVerts);
    nCellsTot *= nCells;
  }

  // Create the (uniform, logically rectilinear) mesh.
  Vector<Dim> origin(0.0), spacings(0.2);
  typedef UniformRectilinearMesh<Dim> Mesh_t;
  DomainLayout<Dim> layout(vertexDomain, GuardLayers<Dim>(0));

  // Create the Fields:
  Centering<Dim> cell = canonicalCentering<Dim>(CellType, Continuous);

  // Full, Antisymmetric, Symmetric, Diagonal Tensor Fields:
  Field<Mesh_t,Tensor<Dim,double,Full> > tff(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<Dim,double,Antisymmetric> >
    tfa(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<Dim,double,Symmetric> >
    tfs(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<Dim,double,Diagonal> >
    tfd(cell, layout, origin, spacings);

  // Assign values:
  Tensor<Dim,double,Full> tf(0.0);
  Tensor<Dim,double,Symmetric> ts(0.0);
  Tensor<Dim,double,Antisymmetric> ta(0.0);
  Tensor<Dim,double,Diagonal> td(0.0);
  double diagDet = 1.0;
  for (int i = 0; i < Dim; i++) {
    for (int j = 0; j < Dim; j++) {
      tf(i,j) = (i+1)*(i+1) + (j+1)*(j+1) + (i+4)*(j+4) + i;
      if (i == j) diagDet *= tf(i,j);
    }
  }
  ts = symmetrize<Symmetric>(tf);
  ta = symmetrize<Antisymmetric>(tf);
  td = symmetrize<Diagonal>(tf);
  tff = tf;
  tfs = ts;
  tfa = ta;
  tfd = td;

  // Test det of Full Tensor:
  double detValue;
  detValue = sum(det(tff));
  switch (Dim) {
  case 1:
    if (!tester.check("detValue", detValue, 18.0*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tff)) = " << detValue
                   << " != 18.0*nCellsTot = " << 18.0*nCellsTot << std::endl;
    }
    break;
  case 2:
    if (!tester.check("detValue", detValue, -38.0*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tff)) = " << detValue
                   << " != -38.0*nCellsTot = " << -38.0*nCellsTot << std::endl;
    }
    break;
  case 3:
    if (!tester.check("detValue", detValue, -4.0*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tff)) = " << detValue
                   << " != -4.0*nCellsTot = " << -4.0*nCellsTot << std::endl;
    }
    break;
  default:
    PInsist(Dim<4, "Tensor det() function not implemented for D>3!");
    break;
  }

  // Test det of Symmetric Tensor:
  detValue = sum(det(tfs));
  switch (Dim) {
  case 1:
    if (!tester.check("detValue", detValue, 18.0*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tfs)) = " << detValue
                   << " != 18.0*nCellsTot = " << 18.0*nCellsTot << std::endl;
    }
    break;
  case 2:
    if (!tester.check("detValue", detValue, -38.25*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tfs)) = " << detValue
                   << " != -38.25*nCellsTot = " 
                   << -38.25*nCellsTot << std::endl;
    }
    break;
  case 3:
    if (!tester.check("detValue", detValue, -4.0*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tfs)) = " << detValue
                   << " != -4.0*nCellsTot = " << -4.0*nCellsTot << std::endl;
    }
    break;
  default:
    PInsist(Dim<4, "Tensor det() function not implemented for D>3!");
    break;
  }

  // Test det of Antiymmetric Tensor:
  detValue = sum(det(tfa));
  switch (Dim) {
  case 1:
    if (!tester.check("detValue", detValue, 0.0)) {
      tester.out() << Dim << "D, sum(det(tfa)) = " << detValue
                   << " != 0.0 ." << std::endl;
    }
    break;
  case 2:
    if (!tester.check("detValue", detValue, ta(1,0)*ta(1,0)*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tfa)) = " << detValue
                   << " != ta(1,0)*ta(1,0)*nCellsTot = " 
                   << ta(1,0)*ta(1,0)*nCellsTot << std::endl;
    }
    break;
  case 3:
    if (!tester.check("detValue", detValue, 0.0)) {
      tester.out() << Dim << "D, sum(det(tfa)) = " << detValue
                   << " != 0.0 ." << std::endl;
    }
    break;
  default:
    PInsist(Dim<4, "Tensor det() function not implemented for D>3!");
    break;
  }

  // Test det of Diagonal Tensor:
  detValue = sum(det(tfd));
  if (Dim <= 3) {
    if (!tester.check("detValue", detValue, diagDet*nCellsTot)) {
      tester.out() << Dim << "D, sum(det(tfd)): " << detValue 
                   << " != diagDet*nCellsTot = " 
                   << diagDet*nCellsTot << std::endl;
    }
  } else {
    PInsist(Dim<4, "Tensor det() function not implemented for D>3!");
  }

}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestDet.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
