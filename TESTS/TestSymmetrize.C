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
// Various tests of the symmetrize<>() template function on Tensors
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

//-----------------------------------------------------------------------------
// Simplistic named functions that return a Full Tensor with specific symmetry.
// Use these for checking correctness of symmetrize<>() template functions.
//-----------------------------------------------------------------------------

template<int D, class T, class EngineTag>
Tensor<D,T,Full>
makeSymmetric(Tensor<D,T,EngineTag> &x)
{
  Tensor<D,T,Full> y(0.0);
  for (int i = 0; i < D; i++) {
    y(i,i) = x(i,i);
    for (int j = i + 1; j < D; j++) {
      y(i,j) = (x(i,j) + x(j,i))*0.5;
      y(j,i) = y(i,j);
    }
  }
  return y;
}

template<int D, class T, class EngineTag>
Tensor<D,T,Full>
makeAntisymmetric(Tensor<D,T,EngineTag> &x)
{
  Tensor<D,T,Full> y(0.0);
  for (int i = 1; i < D; i++) {
    for (int j = 0; j < i; j++) {
      y(i,j) = (x(i,j) - x(j,i))*0.5;
      y(j,i) = -y(i,j);
    }
  }
  return y;
}

template<int D, class T, class EngineTag>
Tensor<D,T,Full>
makeDiagonal(Tensor<D,T,EngineTag> &x)
{
  Tensor<D,T,Full> y(0.0);
  for (int i = 0; i < D; i++) {
    y(i,i) = x(i,i);
  }
  return y;
}

template<int D, class T, class EngineTag>
Tensor<D,T,Full>
makeFull(Tensor<D,T,EngineTag> &x)
{
  Tensor<D,T,Full> y(0.0);
  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++) {
      y(i,j) = x(i,j);
    }
  }
  return y;
}


//-----------------------------------------------------------------------------
// testSymmetrize function; main program calls this for various D values:
//-----------------------------------------------------------------------------

template<int D>
void testSymmetrize(Pooma::Tester &tester)
{

  tester.out() << std::endl << "========= " << D << "D =========" << std::endl;

  // Create Full, Antisymmetric, Symmetric, and Diagonal Tensors as inputs:
  Tensor<D,double,Full> tf;
  double value = 1.0;
  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++) {
      tf(i,j) = value;
      value++;
    }
  }
  tester.out() << "tf: " << tf << std::endl;
  Tensor<D,double,Antisymmetric> ta;
  value = 1.0;
  for (int i = 0; i < TensorStorageSize<D,Antisymmetric>::Size; i++) {
    ta(i) = i + 1.0;
  }
  tester.out() << "ta: " << ta << std::endl;
  Tensor<D,double,Symmetric> ts;
  for (int i = 0; i < TensorStorageSize<D,Symmetric>::Size; i++) {
    ts(i) = i + 1.0;
  }
  tester.out() << "ts: " << ts << std::endl;
  Tensor<D,double,Diagonal> td;
  for (int i = 0; i < TensorStorageSize<D,Diagonal>::Size; i++) {
    td(i) = i + 1.0;
  }
  tester.out() << "td: " << td << std::endl;

  //---------------------------------------------------------------------------
  // Make Fields of these types, to test forwarding of symmetrize<>() function:
  // Create the physical Domains:
  const int nVerts = 6;
  const int nCells = nVerts - 1;
  int nCellsTot = 1;
  Interval<D> vertexDomain;
  for (int d = 0; d < D; d++) {
    vertexDomain[d] = Interval<1>(nVerts);
    nCellsTot *= nCells;
  }

  // Create the (uniform, logically rectilinear) mesh.
  Vector<D> origin(0.0), spacings(0.2);
  typedef UniformRectilinearMesh<D> Mesh_t;
  DomainLayout<D> layout(vertexDomain, GuardLayers<D>(0));

  // Create the Fields:
  Centering<D> cell = canonicalCentering<D>(CellType, Continuous);

  // Full, Antisymmetric, Symmetric, Diagonal Tensor Fields:
  Field<Mesh_t,Tensor<D,double,Full> > tff(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<D,double,Symmetric> >
    tfs(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<D,double,Antisymmetric> >
    tfa(cell, layout, origin, spacings);
  Field<Mesh_t,Tensor<D,double,Diagonal> > tfd(cell, layout, origin, spacings);

  // Assign to the single-Tensor values:
  tff = tf;
  tfs = ts;
  tfa = ta;
  tfd = td;

#ifndef __MWERKS__ // This whole module is too much for CW5 to compile
  // --------------------------------------------------------------------------
  // Symmetrize from Full tensor to {Antisymmetric, Symmetric, Diagonal}:

  // To Antisymmetric:
  if (!tester.check("symmetrize<Antisymmetric>(tf): ", 
                    (symmetrize<Antisymmetric>(tf) == 
                     makeAntisymmetric(tf)))) {
    tester.out() << "symmetrize<Antisymmetric>(tf) = " 
                 << symmetrize<Antisymmetric>(tf) 
                 << " != makeAntisymmetric(tf) = " 
                 << makeAntisymmetric(tf) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Antisymmetric>(tff)): ", 
                    (sum(symmetrize<Antisymmetric>(tff)) == 
                     nCellsTot*makeAntisymmetric(tf)))) {
    tester.out() << "sum(symmetrize<Antisymmetric>(tff)) = " 
                 << sum(symmetrize<Antisymmetric>(tff)) 
                 << " != nCellsTot*makeAntisymmetric(tf) = "
                 << nCellsTot*makeAntisymmetric(tf) << std::endl;
  }


  // To Symmetric:
  if (!tester.check("symmetrize<Symmetric>(tf): ", 
                    (symmetrize<Symmetric>(tf) == 
                     makeSymmetric(tf)))) {
    tester.out() << "symmetrize<Symmetric>(tf) = " 
                 << symmetrize<Symmetric>(tf) 
                 << " != makeSymmetric(tf) = " 
                 << makeSymmetric(tf) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Symmetric>(tff)): ", 
                    (sum(symmetrize<Symmetric>(tff)) == 
                     nCellsTot*makeSymmetric(tf)))) {
    tester.out() << "sum(symmetrize<Symmetric>(tff)) = " 
                 << sum(symmetrize<Symmetric>(tff)) 
                 << " != nCellsTot*makeSymmetric(tf) = "
                 << nCellsTot*makeSymmetric(tf) << std::endl;
  }

  // To Diagonal::
  if (!tester.check("symmetrize<Diagonal>(tf): ", 
                    (symmetrize<Diagonal>(tf) == 
                     makeDiagonal(tf)))) {
    tester.out() << "symmetrize<Diagonal>(tf) = " 
                 << symmetrize<Diagonal>(tf) 
                 << " != makeDiagonal(tf) = " 
                 << makeDiagonal(tf) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Diagonal>(tff)): ", 
                    (sum(symmetrize<Diagonal>(tff)) == 
                     nCellsTot*makeDiagonal(tf)))) {
    tester.out() << "sum(symmetrize<Diagonal>(tff)) = " 
                 << sum(symmetrize<Diagonal>(tff)) 
                 << " != nCellsTot*makeDiagonal(tf) = "
                 << nCellsTot*makeDiagonal(tf) << std::endl;
  }

  // --------------------------------------------------------------------------
  // Symmetrize from Antisymmetric tensor to {Full, Symmetric, Diagonal}:

  // To Full:
  if (!tester.check("symmetrize<Full>(ta): ", 
                    (symmetrize<Full>(ta) == 
                     makeFull(ta)))) {
    tester.out() << "symmetrize<Full>(ta) = " 
                 << symmetrize<Full>(ta) 
                 << " != makeFull(ta) = " 
                 << makeFull(ta) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Full>(tfa)): ", 
                    (sum(symmetrize<Full>(tfa)) == 
                     nCellsTot*makeFull(ta)))) {
    tester.out() << "sum(symmetrize<Full>(tfa)) = " 
                 << sum(symmetrize<Full>(tfa)) 
                 << " != nCellsTot*makeFull(ta) = "
                 << nCellsTot*makeFull(ta) << std::endl;
  }

  // To Symmetric:
  if (!tester.check("symmetrize<Symmetric>(ta): ", 
                    (symmetrize<Symmetric>(ta) == 
                     makeSymmetric(ta)))) {
    tester.out() << "symmetrize<Symmetric>(ta) = " 
                 << symmetrize<Symmetric>(ta) 
                 << " != makeSymmetric(ta) = " 
                 << makeSymmetric(ta) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Symmetric>(tfa)): ", 
                    (sum(symmetrize<Symmetric>(tfa)) == 
                     nCellsTot*makeSymmetric(ta)))) {
    tester.out() << "sum(symmetrize<Symmetric>(tfa)) = " 
                 << sum(symmetrize<Symmetric>(tfa)) 
                 << " != nCellsTot*makeSymmetric(ta) = "
                 << nCellsTot*makeSymmetric(ta) << std::endl;
  }

  // To Diagonal::
  if (!tester.check("symmetrize<Diagonal>(ta): ", 
                    (symmetrize<Diagonal>(ta) == 
                     makeDiagonal(ta)))) {
    tester.out() << "symmetrize<Diagonal>(ta) = " 
                 << symmetrize<Diagonal>(ta) 
                 << " != makeDiagonal(ta) = " 
                 << makeDiagonal(ta) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Diagonal>(tfa)): ", 
                    (sum(symmetrize<Diagonal>(tfa)) == 
                     nCellsTot*makeDiagonal(ta)))) {
    tester.out() << "sum(symmetrize<Diagonal>(tfa)) = " 
                 << sum(symmetrize<Diagonal>(tfa)) 
                 << " != nCellsTot*makeDiagonal(ta) = "
                 << nCellsTot*makeDiagonal(ta) << std::endl;
  }
#endif // __MWERKS__ 

  // --------------------------------------------------------------------------
  // Symmetrize from Symmetric tensor to {Full, Antisymmetric, Diagonal}:

  // To Full:
  if (!tester.check("symmetrize<Full>(ts): ", 
                    (symmetrize<Full>(ts) == 
                     makeFull(ts)))) {
    tester.out() << "symmetrize<Full>(ts) = " 
                 << symmetrize<Full>(ts) 
                 << " != makeFull(ts) = " 
                 << makeFull(ts) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Full>(tfs)): ", 
                    (sum(symmetrize<Full>(tfs)) == 
                     nCellsTot*makeFull(ts)))) {
    tester.out() << "sum(symmetrize<Full>(tfs)) = " 
                 << sum(symmetrize<Full>(tfs)) 
                 << " != nCellsTot*makeFull(ts) = "
                 << nCellsTot*makeFull(ts) << std::endl;
  }

  // To Antisymmetric:
  if (!tester.check("symmetrize<Antisymmetric>(ts): ", 
                    (symmetrize<Antisymmetric>(ts) == 
                     makeAntisymmetric(ts)))) {
    tester.out() << "symmetrize<Antisymmetric>(ts) = " 
                 << symmetrize<Antisymmetric>(ts) 
                 << " != makeAntisymmetric(ts) = " 
                 << makeAntisymmetric(ts) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Antisymmetric>(tfs)): ", 
                    (sum(symmetrize<Antisymmetric>(tfs)) == 
                     nCellsTot*makeAntisymmetric(ts)))) {
    tester.out() << "sum(symmetrize<Antisymmetric>(tfs)) = " 
                 << sum(symmetrize<Antisymmetric>(tfs)) 
                 << " != nCellsTot*makeAntisymmetric(ts) = "
                 << nCellsTot*makeAntisymmetric(ts) << std::endl;
  }

  // To Diagonal::
  if (!tester.check("symmetrize<Diagonal>(ts): ", 
                    (symmetrize<Diagonal>(ts) == 
                     makeDiagonal(ts)))) {
    tester.out() << "symmetrize<Diagonal>(ts) = " 
                 << symmetrize<Diagonal>(ts) 
                 << " != makeDiagonal(ts) = " 
                 << makeDiagonal(ts) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Diagonal>(tfs)): ", 
                    (sum(symmetrize<Diagonal>(tfs)) == 
                     nCellsTot*makeDiagonal(ts)))) {
    tester.out() << "sum(symmetrize<Diagonal>(tfs)) = " 
                 << sum(symmetrize<Diagonal>(tfs)) 
                 << " != nCellsTot*makeDiagonal(ts) = "
                 << nCellsTot*makeDiagonal(ts) << std::endl;
  }

  // --------------------------------------------------------------------------
  // Symmetrize from Diagonal tensor to {Full, Antisymmetric, Symmetric}:

  // To Full:
  if (!tester.check("symmetrize<Full>(td): ", 
                    (symmetrize<Full>(td) == 
                     makeFull(td)))) {
    tester.out() << "symmetrize<Full>(td) = " 
                 << symmetrize<Full>(td) 
                 << " != makeFull(td) = " 
                 << makeFull(td) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Full>(tfd)): ", 
                    (sum(symmetrize<Full>(tfd)) == 
                     nCellsTot*makeFull(td)))) {
    tester.out() << "sum(symmetrize<Full>(tfd)) = " 
                 << sum(symmetrize<Full>(tfd)) 
                 << " != nCellsTot*makeFull(td) = "
                 << nCellsTot*makeFull(td) << std::endl;
  }

  // To Antisymmetric:
  if (!tester.check("symmetrize<Antisymmetric>(td): ", 
                    (symmetrize<Antisymmetric>(td) == 
                     makeAntisymmetric(td)))) {
    tester.out() << "symmetrize<Antisymmetric>(td) = " 
                 << symmetrize<Antisymmetric>(td) 
                 << " != makeAntisymmetric(td) = " 
                 << makeAntisymmetric(td) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Antisymmetric>(tfd)): ", 
                    (sum(symmetrize<Antisymmetric>(tfd)) == 
                     nCellsTot*makeAntisymmetric(td)))) {
    tester.out() << "sum(symmetrize<Antisymmetric>(tfd)) = " 
                 << sum(symmetrize<Antisymmetric>(tfd)) 
                 << " != nCellsTot*makeAntisymmetric(td) = "
                 << nCellsTot*makeAntisymmetric(td) << std::endl;
  }

  // To Symmetric:
  if (!tester.check("symmetrize<Symmetric>(td): ", 
                    (symmetrize<Symmetric>(td) == 
                     makeSymmetric(td)))) {
    tester.out() << "symmetrize<Symmetric>(td) = " 
                 << symmetrize<Symmetric>(td) 
                 << " != makeSymmetric(td) = " 
                 << makeSymmetric(td) << std::endl;
  }
  if (!tester.check("sum(symmetrize<Symmetric>(tfd)): ", 
                    (sum(symmetrize<Symmetric>(tfd)) == 
                     nCellsTot*makeSymmetric(td)))) {
    tester.out() << "sum(symmetrize<Symmetric>(tfd)) = " 
                 << sum(symmetrize<Symmetric>(tfd)) 
                 << " != nCellsTot*makeSymmetric(td) = "
                 << nCellsTot*makeSymmetric(td) << std::endl;
  }

}

//-----------------------------------------------------------------------------
// Main program:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  testSymmetrize<3>(tester);
  testSymmetrize<2>(tester);
  testSymmetrize<1>(tester);

  int ret = tester.results("TestSymmetrize");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestSymmetrize.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
