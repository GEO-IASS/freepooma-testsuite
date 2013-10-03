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
// evaluatorTest9 - testing ScalarCode and custom evaluation domain
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Pooma/Fields.h" // for PerformUpdateTag() only!
#include "Evaluator/ScalarCode.h"
#include "Utilities/Tester.h"
#include <iostream>


// dummy operation

template <int Dim>
struct Copy
{
  Copy(int val) : val_m(val) {}

  template<class A>
  inline void operator()(const A &a, const Loc<Dim> &i) const
  {
     a(i) = val_m;
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(1);
    i.dimensions(Dim);
    i.write(1, true);
    i.useGuards(0, false);
  }

  const int val_m;
};


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Pooma::blockingExpressions(true);

  Interval<2> domain(16, 16);
  Loc<2> blocks(4, 4);
  UniformGridLayout<2> layout(domain, blocks, GuardLayers<2>(1), DistributedTag());
  UniformRectilinearMesh<2> mesh(layout);
  Centering<2> cell = canonicalCentering<2>(CellType, Continuous);

  Field<UniformRectilinearMesh<2>, int, MultiPatch<UniformTag, Remote<Brick> > >
    a(cell, layout, mesh),
    b(cell, layout, mesh);

  // initialize with zero
  a.all() = 0;
  b.all() = 0;

  // do assignments to various subdomains with both expression engine
  // and scalar code functor and compare the full results.
  Interval<2> I;

  (ScalarCode<Copy<2> >(1))(a);
  b = 1;
  tester.check("default (physical) domain", all(a.all() == b.all()));

  I = Interval<2>(Interval<1>(8, 14), Interval<1>(0, 14));
  (ScalarCode<Copy<2> >(2))(a, I);
  b(I) = 2;
  tester.check("partial set of physical patches", all(a.all() == b.all()));

  I = Interval<2>(Interval<1>(6, 9), Interval<1>(6, 9));
  (ScalarCode<Copy<2> >(3))(a, I);
  b(I) = 3;
  tester.check("arbitrary physical domain", all(a.all() == b.all()));

  I = Interval<2>(Interval<1>(0, 15), Interval<1>(-1, 2));
  (ScalarCode<Copy<2> >(4))(a, I);
  b(I) = 4;
  tester.check("arbitrary domain", all(a.all() == b.all()));

  int retval = tester.results("evaluatorTest9 (ScalarCode, evaluation domain)");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: evaluatorTest9.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
