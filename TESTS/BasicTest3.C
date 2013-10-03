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
// Basic Test3: Copying
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

#include <iostream>
#include <cmath>

//-----------------------------------------------------------------------------
// Main program:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Create the physical domain.
  
  const int NX = 5, NY = 5;
  Interval<1> I(NX), J(NY);
  
  Vector<2,double> origin;
  Vector<2,double> spacings;
  for (int d = 0; d < 2; d++) 
    {
      origin(d) = d;
      spacings(d) = d + 1;
    }

  DomainLayout<2> layout1(Interval<2>(I, J), GuardLayers<2>(1));
  Centering<2> vert = canonicalCentering<2>(VertexType, Continuous);
  
  Field<UniformRectilinearMesh<2>, double, Brick> f(vert, layout1, origin, spacings);
  f.all() = 2.0;

  // make a shallow copy
  Field<UniformRectilinearMesh<2>, double, Brick> g(f);
  tester.check("f == g", all(f == g));

  // write to shared engine
  f.all() = 5.0;
  tester.check("f == g", all(f == g));

  // break relation between f and g, write to now private g
  g.makeOwnCopy();
  g.all() = 1.0;
  tester.check("f != g", all(f == 5.0) && all(g == 1.0));

  int ret = tester.results("BasicTest3");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BasicTest3.cpp,v $   $Author: richi $
// $Revision: 1.3 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
