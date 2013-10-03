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
// Test replicating field values.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"
#include <vector>


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const double epsilon = 1.0e-08;
  const int Dim = 2;
  Centering<Dim> inputCenteringTwo, outputCenteringTwo;
  Interval<Dim> physicalVertexDomain(4, 4);
  DomainLayout<Dim> layout(physicalVertexDomain, GuardLayers<Dim>(1));
  typedef Field<UniformRectilinearMesh<Dim>, double, Brick> Field_t;

  // Test 2D Continuous Cell -> Discontinuous Edge.

  inputCenteringTwo = canonicalCentering<Dim>(CellType, Continuous, AllDim);
  outputCenteringTwo = canonicalCentering<Dim>(EdgeType, Discontinuous, AllDim);
  Field_t g(outputCenteringTwo, layout,
	    Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  Field_t f(inputCenteringTwo, layout,
	    Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  f.all() = 2.0;
  g.all() = 1.0;

  g = replicate(f, nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo, true),
		outputCenteringTwo);

  Pooma::blockAndEvaluate();
  tester.check("cell->discontinuous edge",
	       g(FieldOffset<Dim>(Loc<Dim>(0), 0), Loc<Dim>(0)),
	       2.0, epsilon);

  // Test 2D Continuous Vertex -> Discontinuous Vertex.

  inputCenteringTwo =
    canonicalCentering<Dim>(VertexType, Continuous, AllDim);
  outputCenteringTwo =
    canonicalCentering<Dim>(VertexType, Discontinuous, AllDim);

  Field_t g2(outputCenteringTwo, layout,
	     Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  Field_t f2(inputCenteringTwo, layout,
	     Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  f2.all() = 2.0;
  g2.all() = 1.0;

  g2 = replicate(f2, nearestNeighbors(inputCenteringTwo, outputCenteringTwo),
		 outputCenteringTwo);
  Pooma::blockAndEvaluate();
  tester.check("vertex->discontinuous vertex",
	       g2(FieldOffset<Dim>(Loc<Dim>(0), 0), Loc<Dim>(0)),
	       2.0, epsilon);

  int ret = tester.results("Replicate");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Replicate.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
