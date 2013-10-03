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
// DynamicArray test 1: Create/destroy operations.
//-----------------------------------------------------------------------------

// include files 

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/DynamicArrays.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": DynamicArray create/destroy ops." << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  // Create simple dynamic and multipatch dynamic arrays.

  tester.out() << "Creating DynamicArray ..." << std::endl;
  Interval<1> D(6);
  DynamicArray<int> b2(D);

  tester.out() << "Creating MultiPatch DynamicArray ..." << std::endl;
  Loc<1> blocks(3);
  GridPartition<1> gpar(blocks);
  LocalMapper<1> cmap(gpar);  
  DynamicLayout dynlayout(D,gpar,cmap);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > c2(dynlayout);

  // Initialize the arrays with scalars.
  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  int sum2 = 0;

  tester.out() << "Initializing DynamicArrays ..." << std::endl;
  for (int i=0; i < D.size(); ++i)
    {
      b2(i) = 10 + i;
      c2(i) = 10 + i;
      sum2 += (10 + i);
    }
  tester.out() << "Initialization complete, sum2 = " << sum2 << std::endl;
  tester.out() << "DynamicArray b2 = " << b2 << std::endl;
  tester.out() << "MP DynamicArray c2 = " << c2 << std::endl;
  tester.check("DynamicArray initial sum",
	       sum(b2) == sum2);
  tester.check("MP DynamicArray initial sum",
	       sum(c2) == sum2);

  // Now create some elements at the end of the arrays, and make sure they
  // have the proper length.

  tester.out() << "Creating elements ..." << std::endl;
  b2.create(2);
  c2.create(2);
  b2.sync();
  c2.sync();
  tester.out() << "Domain of b2 is now = " << b2.domain() << std::endl;
  tester.out() << "Domain of c2 is now = " << c2.domain() << std::endl;
  tester.check("DynamicArray size after create",
	       b2.domain().size() == D.size() + 2);
  tester.check("MP DynamicArray size after create",
	       c2.domain().size() == D.size() + 2);

  // Initialize the new elements
  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  for (int i=D.size(); i < b2.domain().size(); ++i)
    {
      b2(i) = 10 + i;
      c2(i) = 10 + i;
      sum2 += (10 + i);
    }
  tester.out() << "New initialization complete, sum2 = " << sum2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.out() << "c2 = " << c2 << std::endl;
  tester.check("DynamicArray sum after create",
	       sum(b2) == sum2);
  tester.check("MP DynamicArray sum after create",
	       sum(c2) == sum2);

  // Delete the third element in each, using ShiftUp

  tester.out() << "Deleting third element of each w/ ShiftUp ..." << std::endl;
  int elem = 2;
  b2.destroy(Interval<1>(elem,elem), ShiftUp());
  c2.destroy(Interval<1>(elem,elem), ShiftUp());
  b2.sync();
  c2.sync();
  sum2 -= (10 + elem);
  tester.out() << "ShiftUp delete complete, sum2 = " << sum2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.out() << "c2 = " << c2 << std::endl;
  tester.check("DynamicArray sum after ShiftUp delete",
	       sum(b2) == sum2);
  tester.check("MP DynamicArray sum after ShiftUp delete",
	       sum(c2) == sum2);

  // Delete the first element in each, using BackFill

  tester.out() << "Deleting 1st element of each w/ BackFill ..." << std::endl;
  elem = 0;
  b2.destroy(Interval<1>(elem,elem), BackFill());
  c2.destroy(Interval<1>(elem,elem), BackFill());
  b2.sync();
  c2.sync();
  sum2 -= 10;
  tester.out() << "BackFill delete complete, sum2 = " << sum2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.out() << "c2 = " << c2 << std::endl;
  tester.check("DynamicArray sum after BackFill delete",
	       sum(b2) == sum2);
  tester.check("MP DynamicArray sum after BackFill delete",
	       sum(c2) == sum2);

  // Perform create/destroy ops on a DynamicArray with non-zero origin

  tester.out() << "Creating non-zero-offset DynamicArrays ..." << std::endl;
  Interval<1> D4(5,10);
  DynamicArray<int> b4(D4);
  Loc<1> blocks2(2);
  GridPartition<1> gpar2(blocks2);
  LocalMapper<1> cmap2(gpar2);  
  DynamicLayout dynlayout2(D4,gpar2,cmap2);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > c4(dynlayout2);

  // Initialize the arrays with scalars.
  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  sum2 = 0;

  for (int i=0; i < D4.size(); ++i)
    {
      b4(i + D4.first()) = 100 + i;
      c4(i + D4.first()) = 100 + i;
      sum2 += (100 + i);
    }
  tester.out() << "DynamicArray b4 = " << b4 << std::endl;
  tester.out() << "MP DynamicArray c4 = " << c4 << std::endl;
  tester.out() << "initial sum2 = " << sum2 << std::endl;
  tester.check("DynamicArray 2 initial sum",
	       sum(b4) == sum2);
  tester.check("MP DynamicArray 2 initial sum",
	       sum(c4) == sum2);

  // Destroy the middle two elements using a global domain.

  Interval<1> D5(7,8);
  tester.out() << "Deleting elements " << D5 << " of domain ";
  tester.out() << b4.domain() << " w/ ShiftUp ..." << std::endl;
  sum2 -= sum(b4(D5));
  b4.destroy(D5, ShiftUp());
  c4.destroy(D5, ShiftUp());
  b4.sync();
  c4.sync();
  tester.out() << "ShiftUp delete complete, sum2 = " << sum2 << std::endl;
  tester.out() << "b4 = " << b4 << std::endl;
  tester.out() << "c4 = " << c4 << std::endl;
  tester.check("DynamicArray 2 sum after ShiftUp delete",
	       sum(b4) == sum2);
  tester.check("MP DynamicArray 2 sum after ShiftUp delete",
	       sum(c4) == sum2);

  // Destroy second remaining element using a global domain and backfill

  Interval<1> D6(6,6);
  tester.out() << "Deleting elements " << D6 << " of domain ";
  tester.out() << b4.domain() << " w/ BackFill ..." << std::endl;
  sum2 -= sum(b4(D6));
  b4.destroy(D6, BackFill());
  c4.destroy(D6, BackFill());
  b4.sync();
  c4.sync();
  tester.out() << "BackFill delete complete, sum2 = " << sum2 << std::endl;
  tester.out() << "b4 = " << b4 << std::endl;
  tester.out() << "c4 = " << c4 << std::endl;
  tester.check("DynamicArray 2 sum after BackFill delete",
	       sum(b4) == sum2);
  tester.check("MP DynamicArray 2 sum after BackFill delete",
	       sum(c4) == sum2);

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("DynamicArray create/destroy");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_array_test1.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
