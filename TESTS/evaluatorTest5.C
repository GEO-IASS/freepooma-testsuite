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
// evaluatorTest5 - testing ScalarCode and boundary update
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Pooma/Fields.h" // for PerformUpdateTag() only!
#include "Evaluator/ScalarCode.h"
#include "Utilities/Tester.h"
#include <iostream>


// dummy operation

struct DirtyRelations
{
  DirtyRelations() {}

  template<class A>
  inline void operator()(const A &a, const Loc<1> &i) const
  {
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(1);
    i.dimensions(1);
    i.lowerExtent(0) = 0;
    i.upperExtent(0) = 0;
    i.write(0, true);
    i.useGuards(0, false);
  }
};
struct TriggerRelations
{
  TriggerRelations() {}

  template<class A>
  inline void operator()(const A &a, const Loc<1> &i) const
  {
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(1);
    i.dimensions(1);
    i.lowerExtent(0) = 1;
    i.upperExtent(0) = 1;
    i.write(0, false);
    i.useGuards(0, true);
  }
};
struct TriggerAndDirtyRelations
{
  TriggerAndDirtyRelations() {}

  template<class A>
  inline void operator()(const A &a, const Loc<1> &i) const
  {
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(1);
    i.dimensions(1);
    i.lowerExtent(0) = 1;
    i.upperExtent(0) = 1;
    i.write(0, true); // umm - _and_ read...
    i.useGuards(0, true);
  }
};

// boundary condition just incementing a global counter

static int bupd = 0;

class DummyBC
{
public:
  DummyBC() {}
  DummyBC(const DummyBC &) {}
  template <class Target>
  DummyBC(const DummyBC &, const Target &) {}
  DummyBC& operator=(const DummyBC&) {}
  template <class Target>
  void operator()(const Target&) const
  {
     bupd++;
  }
};


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Pooma::blockingExpressions(true);

  int size = 120;

  Interval<1> domain(size);
  DomainLayout<1> layout(domain, GuardLayers<1>(1));
  UniformRectilinearMesh<1> mesh(layout);
  Centering<1> cell = canonicalCentering<1>(CellType, Continuous);

  Field<UniformRectilinearMesh<1>, double, Brick>
     a(cell, layout, mesh), b(cell, layout, mesh);

  tester.out() << "Adding relation\n";
  Pooma::newRelation(DummyBC(), a);
  RelationListItem *rel = a.fieldEngine().data(0, 0).relations()(0);

  tester.check("a has dirty relation", rel->dirty());
  tester.check("a did not have relations applied", bupd == 0);

  bupd = 0;
  rel->setDirty();
  tester.out() << "Applying DirtyRelations()\n";
  ScalarCode<DirtyRelations>()(a);
  // not applying relations here is an optimization we're not able to do right now
  //tester.check("a did not have relations applied", bupd == 0);
  tester.check("a has dirty relation", rel->dirty());

  bupd = 0;
  rel->setDirty();
  tester.out() << "Applying TriggerRelations()\n";
  ScalarCode<TriggerRelations>()(a);
  tester.check("a did have relations applied", bupd == 1);
  tester.check("a has clean relation", !rel->dirty());

  bupd = 0;
  rel->clearDirty();
  tester.out() << "Applying TriggerAndDirtyRelations()\n";
  ScalarCode<TriggerAndDirtyRelations>()(a);
  tester.check("a did not have relations applied", bupd == 0);
  tester.check("a has dirty relation", rel->dirty());

  bupd = 0;
  rel->setDirty();
  tester.out() << "Reading from a.all()\n";
  b.all() = a.all();
  tester.check("a did have relations applied", bupd == 1);
  tester.check("a has clean relation", !rel->dirty());

  int retval = tester.results("evaluatorTest5 (ScalarCode)");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: evaluatorTest5.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
