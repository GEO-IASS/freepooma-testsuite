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
// Simple Reductions on Multipatch arrays.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Arrays.h"
#include "Pooma/Indices.h"
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  Loc<1> blocks(2);
  UniformGridPartition<1> partition(blocks);   
  UniformGridLayout<1> layout(Interval<1>(10), partition, ReplicatedTag());
  Array<1, int, MultiPatch<UniformTag, Brick> > a(layout);
  
  for (int i = 0; i < 10; i++)
    a(i) = i + 1;
  
  int ret;
  bool bret;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), a);
  tester.check("sum", ret, 55);
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpMultiplyAssign(), a(Interval<1>(9)));
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
  
  int return_status = tester.results("ReductionTest2");

  Pooma::finalize();

  return return_status;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ReductionTest2.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
