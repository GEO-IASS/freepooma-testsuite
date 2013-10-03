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
// Simple Reductions of compressible things.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Arrays.h"
#include "Pooma/Indices.h"
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  Array<2, int, CompressibleBrick> a(4, 4);
  a = 2;
  
  int ret;
  bool bret;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpAddAssign(), a);
  tester.check("sum", ret, int(2 * a.domain().size()));
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpMultiplyAssign(), a);
  tester.check("prod", ret, 65536);
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, FnMinAssign(), a);
  tester.check("min", ret, 2);
  tester.out() << ret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(bret, FnAndAssign(), a);
  tester.check("all", bret, true);
  tester.out() << bret << std::endl;

  Reduction<MainEvaluatorTag>().evaluate(ret, OpBitwiseOrAssign(), a);
  tester.check("bitOr", ret, 2);
  tester.out() << ret << std::endl;
  
  int return_status = tester.results("ReductionTest3");

  Pooma::finalize();

  return return_status;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ReductionTest3.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
