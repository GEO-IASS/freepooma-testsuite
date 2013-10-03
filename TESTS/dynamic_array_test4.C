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
// DynamicArray test 4: Create/Destroy ops on arrays that are views
// on the same data.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/DynamicArrays.h"
#include "Pooma/GMPArrays.h"


int main(int argc, char *argv[])
{
  int i;

  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": DynamicArray view dynamic ops." << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  // Create simple single-patch dynamic arrays.

  tester.out() << "Creating DynamicArray objects ..." << std::endl;
  Interval<1> D(6);
  DynamicArray<int> b2(D);
  b2 = 4;
  DynamicArray<int> b3(b2);
  Pooma::blockAndEvaluate();
  tester.out() << "Created DynamicArray b2 = " << b2 << std::endl;
  tester.out() << "Created DynamicArray b3 = " << b3 << std::endl;
  for (i=D.first(); i <= D.last(); ++i)
    {
      tester.check(b2.read(i) == b3.read(i));
      tester.check(b3.read(i) == 4);
    }
  tester.check("b2 sum matches b3", sum(b2) == sum(b3));

  // Create multi-patch dynamic arrays.

  tester.out() << "Creating MP DynamicArray objects ...";
  tester.out() << std::endl;
  Loc<1> blocks(3);
  GridPartition<1> gpar(blocks);
  LocalMapper<1> cmap(gpar);  
  DynamicLayout dynlayout(D,gpar,cmap);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > c2(dynlayout);
  c2 = 7;
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > c3(c2);
  Pooma::blockAndEvaluate();
  tester.out() << "Created MP DynamicArray c2 = " << c2 << std::endl;
  tester.out() << "Created MP DynamicArray c3 = " << c3 << std::endl;
  for (i=D.first(); i <= D.last(); ++i)
    {
      tester.check(c2.read(i) == c3.read(i));
      tester.check(c3.read(i) == 7);
    }
  tester.check("c2 sum matches c3", sum(c2) == sum(c3));

  // Change the array value and see that they are both affected.

  tester.out() << "Changing Array's to be equal to 2 ..." << std::endl;
  b3 = 2;
  c3 = 2;
  Pooma::blockAndEvaluate();
  tester.out() << "New b2 = " << b2 << std::endl;
  tester.out() << "New b3 = " << b3 << std::endl;
  tester.out() << "New c2 = " << c2 << std::endl;
  tester.out() << "New c3 = " << c3 << std::endl;
  for (i=D.first(); i <= D.last(); ++i)
    {
      tester.check(b2.read(i) == b3.read(i));
      tester.check(b2.read(i) == 2);
      tester.check(c2.read(i) == c3.read(i));
      tester.check(c2.read(i) == 2);
    }
  tester.check("b2 sum matches b3", sum(b2) == sum(b3));
  tester.check("c2 sum matches c3", sum(c2) == sum(c3));

  // Create new elements in c2; this should change c3 as well.

  tester.out() << "Creating 3 new elements in c2, set to 3 ..." << std::endl;
  c2.create(3);
  c2.sync();
  c2(i++) = 3;
  c2(i++) = 3;
  c2(i++) = 3;
  tester.out() << "New c2 = " << c2 << std::endl;
  tester.out() << "New c3 = " << c3 << std::endl;
  for (i=D.first(); i <= (D.last() + 3); ++i)
    {
      tester.check(c2.read(i) == c3.read(i));
      tester.check(c2.read(i) == (i <= D.last() ? 2 : 3));
      c2(i) = i - D.first();
    }

  // Create a view of the first three odd elements, and then do
  // a dynamic op.

  tester.out() << "Creating view of the 0, 2, 4th elements ..." << std::endl;
  Array< 1, int, DynamicView> v2(b2(Range<1>(0,4,2)));
  Array< 1, int, MultiPatchView<DynamicTag,Dynamic,1> > v3(c2(Range<1>(0,4,2)));
  tester.out() << "Current b2 = " << b2 << std::endl;
  tester.out() << "Current b3 = " << b3 << std::endl;
  tester.out() << "Current v2 = " << v2 << std::endl;
  tester.out() << "Current c2 = " << c2 << std::endl;
  tester.out() << "Current c3 = " << c3 << std::endl;
  tester.out() << "Current v3 = " << v3 << std::endl;
  tester.out() << "Deleting first four elements of original array c2 ...";
  tester.out() << std::endl;
  c2.destroy(Interval<1>(4), ShiftUp());
  c2.sync();
  tester.out() << "New c2 = " << c2 << std::endl;
  tester.out() << "New c3 = " << c3 << std::endl;
  tester.out() << "New v3 = " << v3 << std::endl;
  for (i=D.first(); i <= (D.last() - 1); ++i)
    {
      int offi = i - D.first();
      tester.check(c2.read(i) == c3.read(i));
      tester.check(c2.read(i) == offi + 4);
      if (offi < 3)
        tester.check(v3.read(offi) == c2.read(2*offi + D.first()));
    }

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("DynamicArray view dynamic ops");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_array_test4.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
