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
// Test of PatchSizeSyncer
//-----------------------------------------------------------------------------

// Include files

#include "Tulip/PatchSizeSyncer.h"
#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Grid.h"

#define BARRIER

#ifndef BARRIER
#if POOMA_CHEETAH
# define BARRIER Pooma::controller()->barrier()
#else
# define BARRIER
#endif
#endif

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const int numContexts = Pooma::contexts();
  const int myContext   = Pooma::context();

  tester.out() << "Testing PatchSizeSyncer class . . .\n";
  tester.out() << "Running with " << numContexts << " contexts." << std::endl;

  int start = myContext * 10;
  int end   = (myContext + 1) * 10;

  Range<1> r(start, end, 2);

  Grid<1> foo(r);

  BARRIER;

  tester.out().setOutputContext(-1);
  tester.out() << "Initializing PatchSizeSyncer with grid " << foo 
	       << std::endl;

  Pooma::PatchSizeSyncer dls(myContext,foo);

  Grid<1> bar;

  dls.calcGlobalGrid(bar);

  BARRIER;

  tester.out().setOutputContext(0);
  tester.out() << "Here's the global grid (everyone's should be the same): " 
	       << std::endl;
  tester.out().setOutputContext(-1);
  tester.out() << bar << std::endl;

  BARRIER;

  tester.check(bar == Grid<1>(Range<1>(0,numContexts*10,2)));

  // Now we test one that mimics adding some particles:

  if (myContext == 1)
    {
      // Add 20 "elements" to the last domain on this context. 
      IndirectionList<int> tmp(foo.storage());
      tmp(foo.size()-1) += 20;
      foo = Grid<1>(tmp);
    }

  if (myContext == 2)
    {
      // Remove 1 "element" from all domains on this context. 
      foo = Grid<1>(Interval<1>(20,25));
    }

  tester.out().setOutputContext(0);
  tester.out() << "This test actually involves some changes..." << std::endl;
  tester.out().setOutputContext(-1);
  tester.out() << "Initializing PatchSizeSyncer with grid " << foo 
	       << std::endl;

  Pooma::PatchSizeSyncer dls2(myContext,foo);

  dls2.calcGlobalGrid(bar);

  BARRIER;

  tester.out().setOutputContext(0);
  tester.out() << "Here's the global grid: " << std::endl;
  tester.out().setOutputContext(-1);
  tester.out() << bar << std::endl;

  BARRIER;

  IndirectionList<int> idata(bar.size());

  int val = 0;
  for (int i = 0; i < bar.size(); ++i)
    {
      idata(i) = val;
      if (i < 9)  val+=2;
      if (i == 9) val+=22;
      if (i > 9 && i < 15) ++val;
      if (i >= 15) val+=2;
    }

  Grid<1> ans2(idata);

  tester.check(bar == ans2);

  int ret = tester.results("PatchSizeSyncer Test");
  Pooma::finalize();

  return ret;
}

