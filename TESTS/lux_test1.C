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
// Lux test 1: Display a 2D and a 3D fixed-size Array with Lux.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Lux.h"
#include "Pooma/Arrays.h"
#include "Pooma/DynamicArrays.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Lux Array display test" << std::endl;
  tester.out() << "-----------------------------------" << std::endl;

#if POOMA_LUX

  // Arrays to display

  Loc<3> blocks3D(1,2,2);
  Interval<3> domain3D(32, 32, 64);
  Interval<2> domain2D(100, 100);
  Interval<1> domain1D(20);
  GridLayout<3> layout3D(domain3D, blocks3D);
  Array<3, double, MultiPatch<GridTag,Brick> > a3D(layout3D);
  Array<2, int, Brick> a2D(domain2D);
  DynamicArray<float, SharedBrick> a1D(domain1D);

  // Initialize the arrays

  a3D = 100 * (iota(domain3D).comp(2) + 1) +
    10 * (iota(domain3D).comp(1) + 1) +
    iota(domain3D).comp(0) + 1;
  a2D = 1 + iota(domain2D).comp(1);
  a1D = 1 + 10 * iota(domain1D).comp(0);
  Pooma::blockAndEvaluate();

  // Create a Lux connection

  tester.out() << "Creating LuxConnection object ..." << std::endl;
  Connection<Lux> lux("test1");
  tester.out() << "Finished creating LuxConnection object." << std::endl;

  // Establish connections for the two arrays; also connect up a view of
  // the first array

  tester.out() << "Connecting a3D for display ..." << std::endl;
  lux.connect("a3D", a3D);
  tester.out() << "Connecting a2D for display ..." << std::endl;
  lux.connect("a2D", a2D);
  tester.out() << "Connecting a1D for display ..." << std::endl;
  lux.connect("a1D", a1D);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  lux.ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Do, in a loop, updates of the arrays, and redisplay/interact.

  int myiters = 20;
  while (myiters-- > 0)
    {
      tester.out() << "Incrementing for iters = " << myiters << std::endl;
      a3D -= 1;
      a2D += 1;
      Pooma::blockAndEvaluate();

      tester.out() << "Resizing dynamic for iters = " << myiters << std::endl;
      a1D.create(5);
      a1D = 1 + 10 * iota(a1D.domain()).comp(0);
      Pooma::blockAndEvaluate();

      tester.out() << "Updating for iters = " << myiters << std::endl;
      lux.update();

      tester.out() << "Interacting for iters = " << myiters << std::endl;
      lux.interact();
    }

  // Delete LUX connection, closing the window.

  tester.out() << "Closing LUX connection ..." << std::endl;
  lux.close();

#else // POOMA_LUX

  tester.out() << "Please configure with --lux to use this test code!"
	       << std::endl;

#endif // POOMA_LUX

  // Finish up and report results

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Lux Array display test");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: lux_test1.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
