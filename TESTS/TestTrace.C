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
// Various tests of trace(Tensor<>) global function.
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

// Forward declarations:
template<int Dim>
void testTrace(Pooma::Tester &tester);


//-----------------------------------------------------------------------------
// Main program:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  testTrace<1>(tester);
  testTrace<2>(tester);
  testTrace<3>(tester);

  int ret = tester.results("TestTrace");
  Pooma::finalize();
  return ret;
}

template<int Dim>
void testTrace(Pooma::Tester &tester)
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
  Tensor<Dim,double,Antisymmetric> ta(0.0);
  Tensor<Dim,double,Symmetric> ts(0.0);
  Tensor<Dim,double,Diagonal> td(0.0);
  double fullSymDiagTrace = 0.0;
  for (int i = 0; i < Dim; i++) {
    for (int j = 0; j < Dim; j++) {
      tf(i,j) = (i+1)*(i+1) + (j+1)*(j+1) + (i+4)*(j+4) + i;
      if (i == j) fullSymDiagTrace += tf(i,j);
    }
  }
  //  tester.out() << "tf = " << tf << std::endl;
  ta = symmetrize<Antisymmetric>(tf);
  ts = symmetrize<Symmetric>(tf);
  td = symmetrize<Diagonal>(tf);
  tff = tf;
  tfa = ta;
  tfs = ts;
  tfd = td;

  // Test trace of Full Tensor:
  double traceValue;
  traceValue = sum(trace(tff));
  tester.out() << "traceValue = sum(trace(tff)): " << traceValue << std::endl;
  tester.check("traceValue", traceValue, fullSymDiagTrace*nCellsTot);

  // Test trace of Symmetric Tensor:
  traceValue = sum(trace(tfs));
  if (!tester.check("traceValue", traceValue, fullSymDiagTrace*nCellsTot)) {
    tester.out() << Dim << "D, sum(trace(tfs)): " << traceValue 
                 << " != fullSymDiagTrace*nCellsTot = " 
                 << fullSymDiagTrace*nCellsTot << std::endl;
  }

  // Test trace of Antisymmetric Tensor:
  traceValue = sum(trace(tfa));
  if (!tester.check("traceValue", traceValue, 0.0*nCellsTot)) {
    tester.out() << Dim << "D, sum(trace(tfa)): " << traceValue 
                 << " != 0.0" << std::endl;
  }

  // Test trace of Diagonal Tensor:
  traceValue = sum(trace(tfd));
  if (!tester.check("traceValue", traceValue, fullSymDiagTrace*nCellsTot)) {
    tester.out() << Dim << "D, sum(trace(tfd)): " << traceValue 
                 << " != fullSymDiagTrace*nCellsTot = " 
                 << fullSymDiagTrace*nCellsTot << std::endl;
  }

}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestTrace.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
