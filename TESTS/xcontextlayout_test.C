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
// Cross context layout test
//-----------------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Pooma/GMPArrays.h"
#include "Layout/UniformGridLayout.h"
//#include "Layout/GridLayout.h"
//#include "Layout/SparseTileLayout.h"
#include "Partition/UniformGridPartition.h"
//#include "Partition/GridPartition.h"
//#include "Partition/TilePartition.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
 {
  // Initialize POOMA and output stream, using Tester class
	
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": UniformGridPartition operations." << std::endl;
  tester.out() << "---------------------------------------------" << std::endl;

  // Create a UniformGridLayout with a non-empty domain and a set of
  // blocks.

  Loc<2> blocks(20,30);
  Interval<2> domain(120, 120);
  tester.out() << "Creating UniformGridLayout with blocks=" << blocks;
  tester.out() << ", domain=" << domain << std::endl;
  UniformGridLayout<2> ugrid1(domain, blocks, ReplicatedTag());
  UniformGridPartition<2> par(Loc<2>(20,30));

  ugrid1.repartition(par);
  ugrid1.repartition(par,DistributedMapper<2>(par));
 
  tester.out() << "Layout = " << ugrid1 << std::endl;
  
  int b = blocks[0].last()*blocks[1].last();
  for ( int j = 0; j < Pooma::contexts(); ++j)
    {
      if(j == Pooma::context()  )
	{
	  tester.out() << " # local domains " << ugrid1.sizeLocal() <<std::endl;
	  tester.check(ugrid1.sizeLocal() == b/Pooma::contexts());
	  tester.out() << " # remote domains " << ugrid1.sizeRemote() <<std::endl;
	  tester.check(ugrid1.sizeRemote() == (Pooma::contexts()-1)*b/Pooma::contexts() );
	}
    }
  tester.out()<<std::endl;

  ugrid1.repartition(par, LocalMapper<2>(par));

 tester.out() << "Layout with LocalMapper = " << ugrid1 << std::endl;
 


  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Cross Context Layout tests");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: xcontextlayout_test.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
