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
// Various tests of transpose(Tensor<>) global function (transposeerminant).
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

// Forward declarations:
template<int Dim>
void testTranspose(Pooma::Tester &tester);


//-----------------------------------------------------------------------------
// Main program:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  testTranspose<1>(tester);
  testTranspose<2>(tester);
  testTranspose<3>(tester);

  int ret = tester.results("TestTranspose");
  Pooma::finalize();
  return ret;
}

template<int Dim>
void testTranspose(Pooma::Tester &tester)
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
  Field<Mesh_t,Tensor<Dim,double,Symmetric> >
    tfs(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<Dim,double,Antisymmetric> >
    tfa(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<Dim,double,Diagonal> >
    tfd(cell, layout, origin, spacings);

  // Assign values:
  Tensor<Dim,double,Full> tf(0.0), tfTranspose;
  Tensor<Dim,double,Symmetric> ts(0.0), tsTranspose;
  Tensor<Dim,double,Antisymmetric> ta(0.0), taTranspose;
  Tensor<Dim,double,Diagonal> td(0.0), tdTranspose;
  for (int i = 0; i < Dim; i++) {
    for (int j = 0; j < Dim; j++) {
      tf(i,j) = (i+1)*(i+1) + (j+1)*(j+1) + (i+4)*(j+4) + i;
      tfTranspose(j,i) = tf(i,j);
    }
  }
  ts = symmetrize<Symmetric>(tf);
  ta = symmetrize<Antisymmetric>(tf);
  td = symmetrize<Diagonal>(tf);
  tff = tf;
  tfs = ts;
  tfa = ta;
  tfd = td;
  tsTranspose = ts;
  tdTranspose = td;
  for (int i = 1; i < Dim; i++) {
    for (int j = 0; j < i;  j++) {
      taTranspose(i,j) = -ta(i,j);
    }
  }
  tester.out() << "tf = " << tf << std::endl;
  tester.out() << "ts = " << ts << std::endl;
  tester.out() << "ta = " << ta << std::endl;
  tester.out() << "td = " << td << std::endl;

  // Test transpose of Full Tensor:
  Tensor<Dim,double,Full> transposeValF;
  transposeValF = sum(transpose(tff));
  if (!tester.check("transposeValF", transposeValF, tfTranspose*nCellsTot)) {
    tester.out() << Dim << "D, sum(transpose(tff)) = " << transposeValF
                 << " != tfTransposenCellsTot = " 
                 << tfTranspose*nCellsTot << std::endl;
  }

  // Test transpose of Symmetric Tensor:
  Tensor<Dim,double,Symmetric> transposeValS;
  transposeValS = sum(transpose(tfs));
  if (!tester.check("transposeValS", transposeValS, tsTranspose*nCellsTot)) {
    tester.out() << Dim << "D, sum(transpose(tfs)) = " << transposeValS
                 << " != tsTranspose*nCellsTot = " 
                 << tsTranspose*nCellsTot << std::endl;
  }

  // Test transpose of Antiymmetric Tensor:
  Tensor<Dim,double,Antisymmetric> transposeValA;
  transposeValA = sum(transpose(tfa));
  if (!tester.check("transposeValA", transposeValA, taTranspose*nCellsTot)) {
    tester.out() << Dim << "D, sum(transpose(tfa)) = " << transposeValA
                 << " != taTranspose*nCellsTot = " 
                 << taTranspose*nCellsTot << std::endl;
  }

  // Test transpose of Diagonal Tensor:
  Tensor<Dim,double,Diagonal> transposeValD;
  transposeValD = sum(transpose(tfd));
  if (!tester.check("transposeValD", transposeValD, tdTranspose*nCellsTot)) {
    tester.out() << Dim << "D, sum(transpose(tfd)) = " << transposeValD
                 << " != tdTranspose*nCellsTot = " 
                 << tdTranspose*nCellsTot << std::endl;
  }

}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestTranspose.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
