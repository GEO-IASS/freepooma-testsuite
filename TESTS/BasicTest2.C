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
// Basic Test2: Simple data parallel expressions.
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
  
  Field<UniformRectilinearMesh<2>, double, Brick> 
    f(vert, layout1, origin, spacings),
    g(vert, layout1, origin, spacings),
    h(vert, layout1, origin, spacings);

  Pooma::addAllConstantFaceBC(f, 0.0);
  Pooma::addAllConstantFaceBC(g, 0.0);
  Pooma::addAllConstantFaceBC(h, 0.0);
    
  for (int i = 0; i <= f.physicalDomain()[0].last(); i++)
    for (int j = 0; j <= f.physicalDomain()[1].last(); j++)
      {
        g(i, j) = 2.0 * i + j;
        h(i, j) = 4.0 - i - 3.0 * j;
      }

  tester.out() << "f = 1.0..." << std::endl;  
  f = 1.0;
  tester.out() << f << std::endl;
  tester.check("f = 1.0", sum(f), 25 * 1.0);

  tester.out() << "f -= g..." << std::endl;  
  f -= g;
  tester.out() << f << std::endl;
  tester.check("f -= g", sum(f), -125.0);

  tester.out() << sum(f) << std::endl;

  tester.out() << "f += sin(g) + 2.0 * h..." << std::endl;  
  f += sin(g) + 2.0 * h;
  tester.out() << f << std::endl;
  tester.check("f += sin(g) + 2.0 * h", sum(f), -324.60252, 1.0e-4);

  tester.out() << "TD f += sin(g) + 2.0 * h..." << std::endl;  
  f(f.totalDomain()) += sin(g.all()) + 2.0 * h(h.totalDomain());
  tester.out() << f << std::endl;
  tester.check("TD f += sin(g) + 2.0 * h", sum(f), -524.20503, 1.0e-4);

  tester.out() << "TD f += sin(g) + 2.0 * h..." << std::endl;  
  f.all() += sin(g(g.totalDomain())) + 2.0 * h(h.totalDomain());
  tester.out() << f << std::endl;
  tester.check("TD f += sin(g) + 2.0 * h", sum(f), -723.80755, 1.0e-4);

  int ret = tester.results("BasicTest2");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BasicTest2.cpp,v $   $Author: richi $
// $Revision: 1.3 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
