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
// UniformGridLayout test: Create and use UniformGridLayout objects
//-----------------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Pooma/UMPArrays.h"
#include "Partition/ContextMapper.h"
#include "Partition/SpatialPartition.h"
#include "Partition/DistributedMapper.h"
#include "Utilities/Tester.h"
#include <iostream>
#include <iterator>


int main(int argc, char *argv[]) {

  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": UniformGridLayout operations." << std::endl;
  tester.out() << "----------------------------------------" << std::endl;


  Interval<1> I1(0,9);
  Interval<2> I2(I1,I1);
  
  UniformGridLayout<2> tgl(I2,UniformGridPartition<2>(Loc<2>(2,2),		
					GuardLayers<2>(2),
					GuardLayers<2>(2)),
		    ReplicatedTag() );
  
  tester.out()<<tgl<<std::endl;



  Interval<1> dom(1, 20);
  Range<1> back(19,2);
  Loc<1> mlocks(2);
  UniformGridPartition<1> martition(mlocks);
  UniformGridLayout<1> mayout2(dom, martition, ReplicatedTag() );
  Array<1, double, MultiPatch<UniformTag, Brick> > aa(mayout2);

  aa = 9;

  aa(back) = 3.0;

  tester.out() << " testing reversed range view of UniformGridLayout " <<std::endl;
  tester.out() << aa<<std::endl;


  // Create a UniformGridLayout from the Grid

  tester.out() << "Creating empty UniformGridLayout<2>:" << std::endl;
  UniformGridLayout<2> ugridlayout;
  tester.out() << ugridlayout << std::endl;

  // Initialize the UniformGridLayout

  Interval<2> udom(Interval<1>(0,19),Interval<1>(0,19));

  tester.out() << "Initializing UniformGridLayout<2>:" << std::endl;
  ugridlayout.initialize(udom,
			 Loc<2>(2), 
			 GuardLayers<2>(2), 
			 GuardLayers<2>(1),
			 ReplicatedTag());
  tester.out() << "Initialized; UniformGridLayout<2>:" << std::endl;
  tester.out() << ugridlayout << std::endl;

  // Find global ID of nodes at some points

  Loc<2> a1(4, 3);
  Loc<2> a2(11, 14);
  tester.out() << "Global ID of Node at pos " << a1 << ": ";
  tester.out() << ugridlayout.globalID(a1) << std::endl;
  tester.check(ugridlayout.globalID(a1) == 0);
  tester.out() << "Global ID of Node at pos " << a2 << ": ";
  tester.out() << ugridlayout.globalID(a2) << std::endl;
  tester.check(ugridlayout.globalID(a2) == 3);

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
  int c = ugridlayout.touches(test, op);
  tester.out() << "Result of touches: " << c << std::endl;
  tester.check(c == 1);

  tester.out() << "Finding touchingAlloc nodes for " << test << std::endl;
  c = ugridlayout.touchesAlloc(test, op);
  tester.out() << "Result of touchesAlloc: " << c << std::endl;
  tester.check(c == 1);

  UniformGridPartition<2> UPartition3(Loc<2>(2),GuardLayers<2>(1),GuardLayers<2>(0));
  DistributedMapper<2> dgpm(UPartition3);

  
  Interval<1> ii1(0,19);
  
  tester.out()<< ii1<<std::endl;
  

  Interval<1> ii2(0,19);
  
  tester.out()<< ii2<<std::endl;
  
  //  Interval<1> gridI0(grid[0].first(),grid[0].last()-1);
  // Interval<1> gridI1(grid[1].first(),grid[1].last()-1);
  //  Interval<2> gridInterval(gridI0,gridI1);
	
  Interval<1> gridI0(0,19),gridI1(0,19);
  Interval<2> gridInterval(gridI0,gridI1);
		   

  UniformGridLayout<2> GL3(gridInterval,UPartition3,ReplicatedTag());
  
  UniformGridLayout<2> GL3m(gridInterval,UPartition3,dgpm);
 
  // test all of the constructors:

  tester.out()<<std::endl<<std::endl<<std::endl;

  Interval<5> dom5(20,20,20,20,20); 
  tester.out()<< std::endl<<"   Interval is "<<dom5<<std::endl;


  {
    UniformGridLayout<5> GL5R(dom5,ReplicatedTag());
    UniformGridLayout<5> GL5D(dom5,DistributedTag());
  }
  
  {
    UniformGridLayout<5> GL5R(dom5,GuardLayers<5>(2),ReplicatedTag());
    UniformGridLayout<5> GL5D(dom5,GuardLayers<5>(2),DistributedTag());
  }
    
  {
    UniformGridLayout<5> GL5R(dom5,
		       GuardLayers<5>(2),
		       ReplicatedTag());

    UniformGridLayout<5> GL5D(dom5,
		       GuardLayers<5>(2),
		       DistributedTag());
  }
  // with loc divisor
  {
    UniformGridLayout<5> GL5R(dom5,Loc<5>(2),ReplicatedTag());
    UniformGridLayout<5> GL5D(dom5,Loc<5>(2),DistributedTag());
  }
  
  {
    UniformGridLayout<5> GL5R(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),ReplicatedTag());
    UniformGridLayout<5> GL5D(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),DistributedTag());
  }
    
  {
    UniformGridLayout<5> GL5R(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),
		       GuardLayers<5>(2),
		       ReplicatedTag());
    UniformGridLayout<5> GL5D(dom5,
		       Loc<5>(2),
		       GuardLayers<5>(2),
		       GuardLayers<5>(2),
		       DistributedTag());

    UniformGridLayout<5> foo;
    foo.initialize(dom5,
		   Loc<5>(2),
		   GuardLayers<5>(2),
		   GuardLayers<5>(2),
		   ReplicatedTag());

    tester.out() << " UGL<5> initialzed "<<std::endl;
    tester.out() << foo << std::endl;

  }
 

 
  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("UniformGridLayout operations");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: uniformgridlayout_test.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
