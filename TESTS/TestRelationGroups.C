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
// Test of relation groups.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

typedef Field<UniformRectilinearMesh<2>, double, Brick> Field_t;

const double g = 9.8;

void computeTotalEnergy(const Field_t &E, const Field_t &K, const Field_t &U)
{
  E = K + U;
}

class ComputeKineticEnergy {
public:

  ComputeKineticEnergy() { }
  
  ComputeKineticEnergy(const ComputeKineticEnergy &, const Field_t &) { }

  void operator()(const Field_t &K, const Field_t &m, const Field_t &v)
  {
    K = m * v * v / 2;
  }
};

void computePotentialEnergy(const Field_t &U, const Field_t &m, const Field_t &h)
{
  U = m * g * h;
}

struct ComputeVelocity
{
  void doit(const Field_t &v, const Field_t &p, const Field_t &m)
  {
    v = p /  m;
  }
};

void morePotentialEnergy(const Field_t &U)
{
  U += 3.0;
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<2> physicalVertexDomain(4, 4);
  DomainLayout<2> layout(physicalVertexDomain);
  
  // Now, we can declare a field.

  Centering<2> cell = canonicalCentering<2>(CellType, Continuous);

  Field_t E(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  Field_t K(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  Field_t U(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  Field_t v(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  Field_t p(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  Field_t m(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  Field_t h(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));

  p = 2.0;
  m = 3.0;
  h = 4.0;  

  ComputeVelocity obj;
  
  Pooma::newRelation(Pooma::functionPtr(computeTotalEnergy), E, K, U);
  Pooma::newRelation(ComputeKineticEnergy(), K, m, v);
  Pooma::newRelation(Pooma::functionPtr(computePotentialEnergy), U, m, h);
  Pooma::newRelation(Pooma::memberPtr(obj, &ComputeVelocity::doit), v, p, m);
  
  unsigned int g2 = Pooma::newRelationGroup();
  Pooma::deactivateRelationGroup(1);
  
  Pooma::newRelation(Pooma::functionPtr(morePotentialEnergy), U);
  Pooma::activateRelationGroup(1);
  Pooma::deactivateRelationGroup(g2);
  
  tester.out() << E << std::endl;
  
  Pooma::activateRelationGroup(g2);
  E.setDirty();
  
  tester.out() << E << std::endl;
    
  int ret = tester.results("TestRelationGroups");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestRelationGroups.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
