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
// Test of the new Centerings class.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const int Dim = 2;

  Centering<Dim> edges
    = canonicalCentering<Dim>(EdgeType, Continuous, XDim | YDim);

  Centering<Dim> cell
    = canonicalCentering<Dim>(CellType, Continuous);

  Interval<Dim> physicalVertexDomain(4, 4);
  DomainLayout<Dim> layout(physicalVertexDomain, GuardLayers<Dim>(1));
  typedef Field<UniformRectilinearMesh<Dim>, double, Brick> Field_t;
  Field_t f(edges, layout, Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  Field_t fS(cell, layout, Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  Field_t g(3, edges, layout, Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));

  // Set some data in the field.
  
  f[0].all() = 2.0; f[0] = -1.0;
  f[1].all() = 3.0; f[1] = -2.0;

  Pooma::blockAndEvaluate();
  
  // Test a field with subfields.

  tester.check("f[0](0,0)",
	       f(FieldOffset<Dim>(Loc<Dim>(0), 0), Loc<Dim>(0)),
	       -1.0, 1.0e-8);
  tester.check("f[0](0,0)",
	       f(FieldOffset<Dim>(Loc<Dim>(2,1), 0), Loc<Dim>(-2,-1)),
	       -1.0, 1.0e-8);
  tester.check("f[0](2,1)",
	       f(FieldOffset<Dim>(Loc<Dim>(2,1), 0), Loc<Dim>(0)),
	       -1.0, 1.0e-8);
  tester.check("f[1](0,0)",
	       f(FieldOffset<Dim>(Loc<Dim>(0), 1), Loc<Dim>(0)),
	       -2.0, 1.0e-8);
  tester.check("f[1](1,2)",
	       f(FieldOffset<Dim>(Loc<Dim>(1,2), 1), Loc<Dim>(0)),
	       -2.0, 1.0e-8);
  f(FieldOffset<Dim>(Loc<Dim>(3,2), 0), Loc<Dim>(-1,-1)) = 1.3;
  f(FieldOffset<Dim>(Loc<Dim>(3,2), 1), Loc<Dim>(-1,-1)) = 10.3;
  tester.check("f[0](2,1)",
	       f(FieldOffset<Dim>(Loc<Dim>(2,1), 0), Loc<Dim>(0)),
	       1.3, 1.0e-08);
  tester.check("f[1](2,1)",
	       f(FieldOffset<Dim>(Loc<Dim>(2,1), 1), Loc<Dim>(0)),
	       10.3, 1.0e-08);
  tester.check("f[0].read(2,1)",
	       f.read(FieldOffset<Dim>(Loc<Dim>(2,1), 0), Loc<Dim>(0)),
	       1.3, 1.0e-08);
  tester.check("f[1].read(2,1)",
	       f.read(FieldOffset<Dim>(Loc<Dim>(2,1), 1), Loc<Dim>(0)),
	       10.3, 1.0e-08);

  // Test a field with no subfields.

  Field_t h(canonicalCentering<Dim>(CellType, Continuous, AllDim),
	    layout, Vector<Dim>(0.0), Vector<Dim>(1.0, 2.0));
  h(FieldOffset<Dim>(Loc<Dim>(0,0)), Loc<Dim>(0,0)) = 1.3;
  h(FieldOffset<Dim>(Loc<Dim>(0,0)), Loc<Dim>(0,1)) = 2.3;
  h(FieldOffset<Dim>(Loc<Dim>(0,0)), Loc<Dim>(1,0)) = 2.8;
  h(FieldOffset<Dim>(Loc<Dim>(1,0)), Loc<Dim>(0,1)) = 3.3;

  Pooma::blockAndEvaluate();
  
  tester.check("h(0,0)",
	       h(FieldOffset<Dim>(Loc<Dim>(-1,-1)), Loc<Dim>(1,1)),
	       1.3, 1.0e-08);
  tester.check("h(0,1)",
	       h(FieldOffset<Dim>(Loc<Dim>(0,1)), Loc<Dim>(0,0)),
	       2.3, 1.0e-08);
  tester.check("h(1,0)",
	       h(FieldOffset<Dim>(Loc<Dim>(0,1)), Loc<Dim>(1,-1)),
	       2.8, 1.0e-08);
  tester.check("h(1,1)",
	       h(FieldOffset<Dim>(Loc<Dim>(0,0)), Loc<Dim>(1,1)),
	       3.3, 1.0e-08);
  tester.check("h.read(1,0)",
	       h.read(FieldOffset<Dim>(Loc<Dim>(0,1)), Loc<Dim>(1,-1)),
	       2.8, 1.0e-08);
  tester.check("h.read(1,1)",
	       h.read(FieldOffset<Dim>(Loc<Dim>(0,0)), Loc<Dim>(1,1)),
	       3.3, 1.0e-08);

  f[0] = iota(f[0].domain()).comp(1) * iota(f[0].domain()).comp(1);
  f[1] = iota(f[1].domain()).comp(0) * iota(f[1].domain()).comp(0);

  // Test the data-parallel uses.

  // Example data-parallel use:
  // result_field = f(FieldOffset, result_centering).  The FieldOffset
  //   should specify a location in the field f.  The second parameter
  //   specifies the desired output centering.

  FieldOffset<2> lowerXEdge(Loc<2>(0, 0), 0), upperXEdge(Loc<2>(0, 1), 0);
  FieldOffset<2>  leftYEdge(Loc<2>(0, 0), 1), rightYEdge(Loc<2>(1, 0), 1);

  tester.out() << f[0].fieldEngine().centering() << std::endl;
  tester.out() << f(upperXEdge, cell).physicalDomain() << std::endl;
  tester.out() << f(upperXEdge, cell) << std::endl;

  fS = (f(upperXEdge, cell) - f(lowerXEdge, cell))
    * (f(rightYEdge, cell) - f(leftYEdge, cell));

  tester.out() << "f" << std::endl << f << std::endl;
  tester.out() << "fS" << std::endl << fS << std::endl;

  f(upperXEdge, cell) = fS;

  // FIXME: Direct assignment to elements should probably work.
  // The following line fails.
  //  f(upperXEdge, cell)(1, 1) = 4.5;

  tester.out() << "f" << std::endl << f << std::endl;

  int ret = tester.results("FieldOffset");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldOffset.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
