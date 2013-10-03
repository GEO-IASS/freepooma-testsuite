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
// GridLayout test: Create and use GridLayout objects
//-----------------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Pooma/GMPArrays.h"
#include "Partition/ContextMapper.h"
#include "Partition/SpatialPartition.h"
#include "Partition/DistributedMapper.h"
#include "Utilities/Tester.h"
#include <iostream>
#include <iterator>


int main(int argc, char *argv[]) {
  int i;

  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": GridLayout operations." << std::endl;
  tester.out() << "----------------------------------------" << std::endl;


  Interval<1> I1(0,9);
  Interval<2> I2(I1,I1);
  
  GridLayout<2> tgl(I2,GridPartition<2>(Loc<2>(2,2),		
					GuardLayers<2>(2),
					GuardLayers<2>(2)),
		    ReplicatedTag() );
  
  tester.out()<<tgl<<std::endl;



  Interval<1> dom(1, 20);
  Range<1> back(19,2);
  Loc<1> mlocks(2);
  GridPartition<1> martition(mlocks);
  GridLayout<1> mayout2(dom, martition, ReplicatedTag() );
  Array<1, double, MultiPatch<GridTag, Brick> > aa(mayout2);

  aa = 9;

  aa(back) = 3.0;

  tester.out() << " testing reversed range view of GridLayout " <<std::endl;
  tester.out() << aa<<std::endl;


  // Create a grid, and a simple set of blocks

  IndirectionList<int> g1(5), g2(4);
  for (i=0; i < 5; ++i)
    {
      g1(i) = (i < 2 ? i + 1 : g1(i-1) + g1(i-2));
      if (i < 4)
	g2(i) = 2 * (i + 1);
    }
  Grid<2> grid(g1, g2);
  tester.out() << "Creating Grid<2> from indirection lists: Grid<2> = \n";
  tester.out() << grid << std::endl;

  // Create a GridLayout from the Grid

  tester.out() << "Creating empty GridLayout<2>:" << std::endl;
  GridLayout<2> gridlayout;
  tester.out() << gridlayout << std::endl;

  // Initialize the GridLayout

  tester.out() << "Initializing GridLayout<2>:" << std::endl;
  gridlayout.initialize(grid, GuardLayers<2>(2), GuardLayers<2>(1),ReplicatedTag());
  tester.out() << "Initialized; GridLayout<2>:" << std::endl;
  tester.out() << gridlayout << std::endl;

  // Find global ID of nodes at some points

  Loc<2> a1(4, 3);
  Loc<2> a2(2, 4);
  tester.out() << "Global ID of Node at pos " << a1 << ": ";
  tester.out() << gridlayout.globalID(a1) << std::endl;
  tester.check(gridlayout.globalID(a1) == 2);
  tester.out() << "Global ID of Node at pos " << a2 << ": ";
  tester.out() << gridlayout.globalID(a2) << std::endl;
  tester.check(gridlayout.globalID(a2) == 5);

  // Find the nodes touching a given domain

  Interval<2> test(5, 5);

#if POOMA_NO_OSTREAM_ITERATOR_1ARG
  std::ostream_iterator<Node< Interval<2> >, char> op(tester.out().stream(),
						      "\n");
#else
  std::ostream_iterator<Node< Interval<2> > > op(tester.out().stream(),
						 "\n");
#endif

  tester.out() << "Finding touching nodes for " << test << std::endl;
  int c = gridlayout.touches(test, op);
  tester.out() << "Result of touches: " << c << std::endl;
  tester.check(c == 6);

  tester.out() << "Finding touchingAlloc nodes for " << test << std::endl;
  c = gridlayout.touchesAlloc(test, op);
  tester.out() << "Result of touchesAlloc: " << c << std::endl;
  tester.check(c == 2);

  GridPartition<2> GPartition3(grid,GuardLayers<2>(1),GuardLayers<2>(0));
  DistributedMapper<2> dgpm(GPartition3);

  tester.out() << grid << std::endl;
  
  tester.out() << grid[0]<<std::endl;
  tester.out() << grid[1]<<std::endl;
  
  tester.out() << grid[0].first()<<std::endl;
  tester.out() << grid[1].first()<<std::endl;
  
  tester.out() << grid[0].last()<<std::endl;
  tester.out() << grid[1].last()<<std::endl;
  
  Interval<1> ii1(grid[0].first(),grid[0].last()-1);
  
  tester.out()<< ii1<<std::endl;
  

  Interval<1> ii2(grid[1].first(),grid[1].last()-1);
  
  tester.out()<< ii2<<std::endl;
  
  Interval<1> gridI0(grid[0].first(),grid[0].last()-1);
  Interval<1> gridI1(grid[1].first(),grid[1].last()-1);
  Interval<2> gridInterval(gridI0,gridI1);
			   

  GridLayout<2> GL3(gridInterval,GPartition3,ReplicatedTag());
  
  GridLayout<2> GL3m(gridInterval,GPartition3,dgpm);
 
  // test all of the constructors:

  tester.out()<<std::endl<<std::endl<<std::endl;

  Interval<5> dom5(20,20,20,20,20); 
  tester.out()<< std::endl<<"   Interval is "<<dom5<<std::endl;


  {
    GridLayout<5> GL5R(dom5,ReplicatedTag());
    GridLayout<5> GL5D(dom5,DistributedTag());
  }
  
  {
    GridLayout<5> GL5R(dom5,GuardLayers<5>(2),ReplicatedTag());
    GridLayout<5> GL5D(dom5,GuardLayers<5>(2),DistributedTag());
  }
    
  {
    GridLayout<5> GL5R(dom5,
		       GuardLayers<5>(2),
		       ReplicatedTag());

    GridLayout<5> GL5D(dom5,
		       GuardLayers<5>(2),
		       DistributedTag());
  }
  // with loc divisor
  {
    GridLayout<5> GL5R(dom5,Loc<5>(2),ReplicatedTag());
    GridLayout<5> GL5D(dom5,Loc<5>(2),DistributedTag());
  }
  
  {
    GridLayout<5> GL5R(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),ReplicatedTag());
    GridLayout<5> GL5D(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),DistributedTag());
  }
    
  {
    GridLayout<5> GL5R(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),
		       GuardLayers<5>(2),
		       ReplicatedTag());
    GridLayout<5> GL5D(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),
		       GuardLayers<5>(2),
		       DistributedTag());
  }
 

 
  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("GridLayout operations");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: gridlayout_test.cpp,v $   $Author: richard $
// $Revision: 1.18 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
