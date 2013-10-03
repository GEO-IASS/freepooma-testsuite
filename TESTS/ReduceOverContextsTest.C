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
// Test of ReduceOverContexts<T, ReductionOp>
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Tulip/ReduceOverContexts.h"
#include "Tulip/RemoteProxy.h"
#include "Utilities/Tester.h"

#include <vector>

#define BARRIER

#ifndef BARRIER
# if POOMA_CHEETAH
#  define BARRIER Pooma::controller()->barrier()
# else
#  define BARRIER
# endif
#endif

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const int numContexts = Pooma::contexts();
  const int myContext   = Pooma::context();

  tester.out() << "Sum test #1" << std::endl;

  typedef ReduceOverContexts<int, OpAddAssign> SumReduction_t;

  std::vector<int> foo(numContexts);

  for (int i = 0; i < numContexts; ++i)
    foo[i] = 3;

  int result1;

  SumReduction_t goo(foo[myContext]);
  goo.broadcast(result1);  

  BARRIER;

  tester.out().setOutputContext(0);
  tester.out() << "This should print three times the number of contexts,"
               << std::endl << "or " << 3*numContexts 
               << ", on all contexts." << std::endl;

  tester.check(result1 == 3*numContexts);

  BARRIER;

  tester.out().setOutputContext(Inform::allContexts);
  tester.out() << result1 << std::endl;

  BARRIER;

  tester.out().setOutputContext(0);

  tester.out() << "Sum test #2" << std::endl;

  std::vector<int> ans(numContexts);
  for (int i = 0; i < numContexts; i++)
    ans[i] = i + 1;

  tester.out() << "This should print context number plus one on each context."
               << std::endl;

  BARRIER;

  tester.out().setOutputContext(Inform::allContexts);
  tester.out() << ans[myContext] << std::endl;

  BARRIER;

  tester.out().setOutputContext(0);
  tester.out() << "Now reduce the values on contexts 0 and 1 only." 
               << std::endl;

  SumReduction_t reduce(ans[myContext], 0, myContext < 2);

  if (myContext == 0)
    ans[0] = reduce;

  BARRIER;

  RemoteProxy<int> broadcast(ans[myContext], 0);
  int final = broadcast;

  BARRIER;

  tester.out().setOutputContext(Inform::allContexts);
  tester.out() << final << std::endl;

  BARRIER;

  tester.check(final == (numContexts > 1 ? 3 : 1));

  int ret = tester.results("ReduceOverContexts Test");

  Pooma::finalize();

  return ret;
}

