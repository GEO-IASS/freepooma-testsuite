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
// Tests of relations.
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

template<int Dim>
struct computePotentialEnergy {
  computePotentialEnergy() {}
  void scalarCodeInfo(ScalarCodeInfo& info) const
  {
    info.dimensions(Dim);
    info.arguments(3);
    info.write(0, true);
    info.write(1, false);
    info.write(2, false);
    info.useGuards(0, false);
    info.useGuards(1, false);
    info.useGuards(2, false);
    for (int i=0; i<Dim; ++i)
    {
      info.lowerExtent(i) = 0;
      info.upperExtent(i) = 0;
    }
  }
  template <class F1, class F2, class F3>
  void operator()(const F1& U, const F2& m, const F3& h, const Loc<Dim>& loc) const
  {
    U(loc) = m(loc) * g * h(loc);
  }
};

struct ComputeVelocity
{
  void doit(const Field_t &v, const Field_t &p, const Field_t &m)
  {
    v = p / m;
  }
};

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<2> physicalVertexDomain(4, 4);
  DomainLayout<2> layout(physicalVertexDomain);
  
  // Now, we can declare a field.

  Centering<2> cell = canonicalCentering<2>(CellType, Continuous);

  // total energy
  Field_t E(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  // kinetic energy
  Field_t K(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  // potential energy
  Field_t U(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  // velocity
  Field_t v(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  // momentum
  Field_t p(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  // mass
  Field_t m(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  // height
  Field_t h(cell, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));

  p = 2.0;
  m = 3.0;
  h = 4.0;

  ComputeVelocity obj;
  
  Pooma::newRelation(Pooma::functionPtr(computeTotalEnergy), E, K, U);
  Pooma::newRelation(ComputeKineticEnergy(), K, m, v);
  Pooma::newRelation(ScalarCode<computePotentialEnergy<2> >(), U, m, h);
  Pooma::newRelation(Pooma::memberPtr(obj, &ComputeVelocity::doit), v, p, m);
  
  tester.out() << E << std::endl;
  tester.check("Total energy at h=4.0", all(E == 3.0*g*4.0 + 0.5*3.0*pow(2.0/3.0, 2)));
  
  h = 0;
  
  tester.out() << E << std::endl;
  tester.check("Total energy at h=0.0", all(E == 3.0*g*0.0 + 0.5*3.0*pow(2.0/3.0, 2)));
    
  int ret = tester.results("TestBasicRelations");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestBasicRelations.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
