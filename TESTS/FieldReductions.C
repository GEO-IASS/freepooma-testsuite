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
// Test field reductions.
//-----------------------------------------------------------------------------


#include "Pooma/Fields.h"
#include "Utilities/Tester.h"
#include <iostream>


// Check the sum, average, minimum, and maximum functions for a
// specified position.

template <class Geometry, class T, class Engine, int Dim>
inline bool
checkFieldPosition(const Field<Geometry, T, Engine> &f,
		   const FieldOffsetList<Dim> &fol,
		   const Loc<Dim> &loc,
		   const T sumAnswer, const T avAnswer,
		   const T minAnswer, const T maxAnswer,
		   const double tolerance)
{
  return 
    std::abs(sum(f, fol, loc) - sumAnswer) < tolerance &&
    std::abs(av(f, fol, loc) - avAnswer) < tolerance &&
    std::abs(min(f, fol, loc) - minAnswer) < tolerance &&
    std::abs(max(f, fol, loc) - maxAnswer) < tolerance;
}


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const double eps = 1.0e-08;	// checking tolerance
  const int Dim = 2;
  std::vector<FieldOffsetList<Dim> > nn;
  std::vector<FieldOffsetList<3> > nn3;

  Centering<2> inputCenteringTwo, outputCenteringTwo;
  Centering<3> inputCenteringThree, outputCenteringThree;

  Interval<Dim> physicalVertexDomain(4, 4);
  DomainLayout<Dim> layout(physicalVertexDomain, GuardLayers<Dim>(1));
  typedef Field<UniformRectilinearMesh<Dim>, double, Brick> Field_t;


  // Test 2D Discontinuous Vertex -> Continuous Vertex.

  inputCenteringTwo = canonicalCentering<2>(VertexType, Discontinuous, AllDim);
  outputCenteringTwo = canonicalCentering<2>(VertexType, Continuous, AllDim);
  nn = nearestNeighbors(inputCenteringTwo, outputCenteringTwo);
  Field_t g(inputCenteringTwo, layout,
	    Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));

  g.all() = 2.0;
  g = -1.0;
  Pooma::blockAndEvaluate();
  g(FieldOffset<Dim>(Loc<Dim>(1,1), 0), Loc<Dim>(0,0)) = 17.0;
  tester.check("discontinuous vertex->continuous vertex",
	       checkFieldPosition(g, nn[0], Loc<Dim>(1),
				  14.0, 3.5, -1.0, 17.0, eps));


  // Test 2D Continuous Cell -> Continuous Cell.

  inputCenteringTwo = canonicalCentering<2>(CellType, Continuous, AllDim);
  outputCenteringTwo = canonicalCentering<2>(CellType, Continuous, AllDim);
  nn = nearestNeighbors(inputCenteringTwo, outputCenteringTwo);
  Field_t f(inputCenteringTwo, layout,
	    Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));

  f.all() = 2.0;
  f = -1.0;
  Pooma::blockAndEvaluate();
  f(FieldOffset<Dim>(Loc<Dim>(1,1), 0), Loc<Dim>(0,0)) = 17.0;
  tester.check("cell->cell",
	       checkFieldPosition(f, nn[0], Loc<Dim>(1,1),
				  17.0, 17.0, 17.0, 17.0, eps));


  // Test 2D Discontinuous Face -> Continuous Edge.

  inputCenteringTwo = canonicalCentering<2>(FaceType, Discontinuous, AllDim);
  outputCenteringTwo = canonicalCentering<2>(EdgeType, Continuous, AllDim);
  nn = nearestNeighbors(inputCenteringTwo, outputCenteringTwo);
  Field_t h(inputCenteringTwo, layout,
	    Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));

  h.all() = 2.0;
  h = -1.0;
  Pooma::blockAndEvaluate();
  h(FieldOffset<Dim>(Loc<Dim>(1), 0), Loc<Dim>(0)) = 17.0;
  tester.check("discontinuous face->edge",
	       checkFieldPosition(h, nn[0], Loc<Dim>(1),
				  -2.0, -1.0, -1.0, -1.0, eps));

  // Test 3D Discontinuous Vertex -> Continuous Cell.

  inputCenteringThree = canonicalCentering<3>(VertexType, Discontinuous, AllDim);
  outputCenteringThree = canonicalCentering<3>(CellType, Continuous, AllDim);
  nn3 = nearestNeighbors(inputCenteringThree, outputCenteringThree);

  Interval<3> physicalVertexDomain3(4, 4, 4);
  DomainLayout<3> layout3(physicalVertexDomain3, GuardLayers<3>(1));
  Field<UniformRectilinearMesh<3>, double, Brick>
    G(inputCenteringThree, layout3, Vector<3>(0.0), Vector<3>(1.0, 2.0, 0.0));

  G.all() = 2.0;
  G = -1.0;
  Pooma::blockAndEvaluate();
  G(FieldOffset<3>(Loc<3>(1), 0), Loc<3>(0)) = 17.0;
  tester.check("discontinuous vertex->cell",
	       checkFieldPosition(G, nn3[0], Loc<3>(1),
				  -46.0, -46.0/64.0, -1.0, 17.0, eps));

  int ret = tester.results("FieldReductions");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldReductions.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
