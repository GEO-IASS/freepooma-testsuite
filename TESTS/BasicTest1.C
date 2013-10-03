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
// Basic Test 1: declaring, viewing, and indexing fields.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

//-----------------------------------------------------------------------------
// Globals
//-----------------------------------------------------------------------------

const int NX = 5, NY = 5;
Interval<1> I(NX), J(NY);
Interval<2> physicalVertexDomain(I, J),
  td(Interval<1>(-1,NX-1), Interval<1>(-1,NY-1)), pd(NX-1, NY-1);
  

//-----------------------------------------------------------------------------
// Test function
//-----------------------------------------------------------------------------

template<class Geom, class T, class EngineTag>
void doTest(Pooma::Tester &tester, Field<Geom, T, EngineTag> &f)
{
  tester.check("PD", f.physicalDomain(), pd);
  tester.check("TD", f.totalDomain(), td);

  Pooma::addAllConstantFaceBC(f, 0.0);    
    
  for (int i = 0; i <= f.physicalDomain()[0].last(); i++)
    for (int j = 0; j <= f.physicalDomain()[1].last(); j++)
      f(i, j) = i + j;
      
  tester.out() << f << std::endl;
  tester.out() << f.read() << std::endl;
  tester.out() << f() << std::endl;
  
  tester.check("f(4,4)", f(4,4), 0.0);
  tester.check("f.read(4,4)", f.read(4,4), 0.0);
  tester.check("f(0,3)", f(0,3), 3.0);
  
  Loc<2> loc(2, 3);
  f(loc) = 1.0;
  tester.out() << f << std::endl;
  
  Range<1> R(0,2,2);
  Range<2> RR(R, R);
  tester.out() << f(R, R) << std::endl;
  tester.check("sum(f(R,R))", sum(f(R,R)), 8.0);
}


//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  // Create the mesh.
  
  Vector<2,double> origin;
  Vector<2,double> spacings;
  for (int d = 0; d < 2; d++) 
    {
      origin(d) = d;
      spacings(d) = d + 1;
    }
      
  // Make a Brick-Engine-based field.

  DomainLayout<2> layout1(physicalVertexDomain, GuardLayers<2>(1));
  Centering<2> cell = canonicalCentering<2>(CellType, Continuous, AllDim);
  
  Field<UniformRectilinearMesh<2>, double, Brick> 
    f(cell, layout1, origin, spacings);

  doTest(tester, f);
  
  // Make a CompressibleBrick-Engine-based field.
  
  Field<UniformRectilinearMesh<2>, double, CompressibleBrick> 
    fc(cell, layout1, origin, spacings);

  doTest(tester, fc);
  
  // Make a Nonuniform Multipatch-Engine-based field.

  Loc<2> blocks(2,2);
  GridLayout<2> layout2(physicalVertexDomain, blocks, 
                        GuardLayers<2>(0), GuardLayers<2>(1),
                        ReplicatedTag());
  typedef MultiPatch<GridTag, Brick> MP2_t;
  Field<UniformRectilinearMesh<2>, double, MP2_t> 
    fg(cell, layout2, origin, spacings);

  doTest(tester, fg);
    
  int ret = tester.results("BasicTest1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BasicTest1.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
