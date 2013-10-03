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
// evaluatorTest5 - testing ScalarCode and expression arguments
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Evaluator/ScalarCode.h"
#include "Utilities/Tester.h"
#include <iostream>


// ScalarCode just evaluating/assigning an expression

struct EvaluateExpr
{
  EvaluateExpr() {}

  template<class LHS, class RHS>
  inline void operator()(const LHS &a, const RHS &b, const Loc<1> &i) const
  {
	  a(i) = b.read(i);
  }

  void scalarCodeInfo(ScalarCodeInfo& i) const
  {
    i.arguments(2);
    i.dimensions(1);
    i.write(0, true);
    i.write(1, false);
    i.useGuards(0, false);
    i.useGuards(1, false);
  }
};


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Pooma::blockingExpressions(true);

  Interval<1> domain(8);
  UniformGridLayout<1> layout(domain, Loc<1>(2), GuardLayers<1>(1), DistributedTag());

  Array<1, int, MultiPatch<UniformTag, Remote<Brick> > >
     a(layout), b(layout), c(layout);

  a = 0;
  b = 1;
  c = 2;
  ScalarCode<EvaluateExpr>()(a, c-b);
  tester.check("a = c - b", all(a(domain) == 1));
  tester.out() << a(domain) << std::endl;

  a = 0;
  ScalarCode<EvaluateExpr>()(a, b(domain-1)+c(domain+1));
  tester.check("a = b(i-1) + c(i+1)", all(a(domain) == 3));
  tester.out() << a(domain) << std::endl;

  tester.out() << "Manually triggering igc fill" << std::endl;
  b.engine().fillGuards();
  c.engine().fillGuards();
  a = 0;
  ScalarCode<EvaluateExpr>()(a, b(domain-1)+c(domain+1));
  tester.check("a = b(i-1) + c(i+1)", all(a(domain) == 3));
  tester.out() << a(domain) << std::endl;

  int retval = tester.results("evaluatorTest10 (ScalarCode with expressions)");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: evaluatorTest10.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
