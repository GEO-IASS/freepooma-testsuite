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
// Simple Reductions on Remote Multipatch arrays.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Arrays.h"
#include "Pooma/Indices.h"
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  Loc<1> blocks2(2), blocks5(5);
  UniformGridPartition<1> partition2(blocks2), partition5(blocks5);   
  UniformGridLayout<1> layout2(Interval<1>(10), partition2, DistributedTag()),
    layout5(Interval<1>(10), partition5, DistributedTag());
  Array<1, int, MultiPatch<UniformTag, Remote<Brick> > > a(layout2),
    b(layout5);
  Array<1, int> c(10);

  for (int i = 0; i < 10; i++)
    {
      a(i) = i + 1;
      b(i) = 2 * i;
      c(i) = i % 5;
    }

  Pooma::blockAndEvaluate();

  int ret;
  bool bret;

  // Test various sorts of reductions with a single array.

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), a);
  tester.check("sum", ret, 55);
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpMultiplyAssign(), 
    a(Interval<1>(9)));
  tester.check("prod", ret, 362880);
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, FnMinAssign(), a - 2);
  tester.check("min", ret, -1);
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(bret, FnAndAssign(), a - 1);
  tester.check("all", bret, false);
  tester.out() << bret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpBitwiseOrAssign(), a);
  tester.check("bitOr", ret, 15);
  tester.out() << ret << std::endl;

  // Test something with an expression engine (remote2 + remote5).

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), a + b);
  tester.check("sum(a + b)", ret, 55 + 90);
  tester.out() << ret << std::endl;

  // Test something with an expression engine (remote5 + remote2).

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), b + a);
  tester.check("sum(b + a)", ret, 90 + 55);
  tester.out() << ret << std::endl;

  // Test something with a brick (remote2 + remote5 + brick).

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), a + b + c);
  tester.check("sum(a + b + c)", ret, 90 + 55 + 20);
  tester.out() << ret << std::endl;

  // Test something with a brick (brick + remote5 + remote2).

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), c + b + a);
  tester.check("sum(c + b + a)", ret,  20 + 55 + 90);
  tester.out() << ret << std::endl;

  // Finish.

  int return_status = tester.results("ReductionTest4");

  Pooma::finalize();

  return return_status;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ReductionTest4.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
