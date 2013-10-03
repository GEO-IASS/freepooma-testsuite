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
// UniformGridPartition test
//-----------------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Layout/UniformGridLayout.h"
#include "Layout/GridLayout.h"
#include "Partition/UniformGridPartition.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[]) {
  int i;

  // Initialize POOMA and output stream, using Tester class
	
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Partition operations." << std::endl;
  tester.out() << "---------------------------------------------" << std::endl;

  // Create a UniformGridLayout with a non-empty domain and a set of
  // blocks.

  Loc<2> blocks(2,3);
  Interval<2> domain(12, 12);
  tester.out() << "Creating UniformGridLayout with blocks=" << blocks;
  tester.out() << ", domain=" << domain << std::endl;
  UniformGridLayout<2> ugrid1(domain, blocks, ReplicatedTag());
  tester.out() << "Layout = " << ugrid1 << std::endl;
  tester.check(ugrid1.sizeLocal() == 6);

  // Create a GridLayout with an empty domain and the same set of blocks

  Interval<2> domain2;
  tester.out() << "Creating GridLayout with blocks=" << blocks;
  tester.out() << ", domain=" << domain2 << std::endl;
  GridLayout<2> grid2(domain2, blocks,ReplicatedTag());
  tester.out() << "Layout = " << grid2 << std::endl;
  tester.check(grid2.sizeLocal() == 6);
  for (i=0; i < grid2.sizeLocal(); ++i)
    tester.check(grid2.domain(i).empty());

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Partition operations");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ugp_test.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:17:02 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
