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
// DynamicArray test 3: Create/destroy for MP Arrays with shared layouts
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
  tester.out() << argv[0] << ": MP DynamicArray w/ shared layouts."
               << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  // Create some Interval objects to create and index into Array's with

  tester.out() << "Creating Interval<1> objects ..." << std::endl;
  Interval<1> D1(3);
  tester.out() << "D1 = " << D1 << std::endl;

  // Create MultiPatch dynamic arrays that share a layout.

  tester.out() << "Creating MP DynamicArray using domain D1 ... " << std::endl;
  Loc<1> blocks(3);
  GridPartition<1> gpar(blocks);
  LocalMapper<1> cmap(gpar);  
  DynamicLayout dynlayout(D1,gpar,cmap);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > a2(dynlayout);
  tester.check("a2 size", a2.domain().size() == D1.size());
  tester.check("a2 patches", a2.layout().sizeLocal() == 3);

  tester.out() << "Creating MP DynamicArray w/ same layout ..." << std::endl;
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > b2(a2.layout());
  tester.check("b2 size", b2.domain().size() == D1.size());
  tester.check("b2 patches", b2.layout().sizeLocal() == 3);

  // Test looping over layout nodes

  tester.out() << "DynamicArray< MultiPatch<DynamicTag,Dynamic> > layout:\n";
  tester.out() << a2.layout() << std::endl;

  // Initialize dynamic arrays with scalars.

  a2 = 30;
  b2 = 40;
  Pooma::blockAndEvaluate();
  tester.out() << "Initialized MP DynamicArray's to 30, 40:" << std::endl;
  tester.out() << "a2 = " << a2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.check("a2 initial sum", sum(a2) == (a2.domain().size() * 30));
  tester.check("b2 initial sum", sum(b2) == (b2.domain().size() * 40));

  // Create elements in the shared-layout MPE arrays

  tester.out() << "Creating 2 elements at end of a2 and b2 ..." << std::endl;
  a2.create(2);
  a2.sync();
  a2(3)=a2(4)=-50;
  b2(3)=b2(4)=-60;

  a2(a2.engine().domain().last()-1)=0;
  a2(a2.engine().domain().last())=0;
  
  tester.out() << "a2 = " << a2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.check("a2 size after create", a2.domain().size() == (D1.size() + 2));
  tester.check("b2 size after create", b2.domain().size() == (D1.size() + 2));

  // Delete an element in the shared-layout MPE arrays

  tester.out() << "Deleting 2nd element of a2 & b2 w/backfill ..."<<std::endl;
  b2.destroy(Interval<1>(1,1), BackFill());
  b2.sync();
  tester.out() << "a2 = " << a2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.check("a2 size after BackFill",a2.domain().size() == (D1.size() + 1));
  tester.check("b2 size after BackFill",b2.domain().size() == (D1.size() + 1));

  // Copy values from the beginning of a2 and b2 to their end

  tester.out() << "Copying first three elements of a2 and b2 ..." << std::endl;
  a2.copy(Interval<1>(3));
  a2.sync();
  tester.out() << "a2 = " << a2 << std::endl;
  tester.out() << "b2 = " << b2 << std::endl;
  tester.check("a2 size after copy", a2.domain().size() == (D1.size() + 4));
  tester.check("b2 size after copy", b2.domain().size() == (D1.size() + 4));

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("MP DynamicArray w/ shared layouts");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_array_test3.cpp,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:16:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
