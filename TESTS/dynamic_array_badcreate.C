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
// DynamicArray bad create: Try create on an DynamicArray with existing view
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/DynamicArrays.h"

// a test function that takes as arguments:
//   1. A DynamicArray
//   2. An array which is a view on the DynamicArray

template<class T, class CA>
bool
testview(Pooma::Tester &tester, DynamicArray<T,Dynamic> &da,
	 const CA &daview)
{
  tester.out() << "In testview:" << std::endl;
  tester.out() << "    da = " << da << std::endl;
  tester.out() << "daview = " << daview << std::endl;

  // the following should crash
  tester.out() << "Trying to create values within da ..." << std::endl;
  bool result = false;
#if POOMA_EXCEPTIONS
  try {
    da.create(3);
    tester.out() << "Ack! create call didn't throw!!!" << std::endl;
    result = false;
  }
  catch(const Pooma::Assertion &a)
  {
    tester.out() << "Caught assertion - it worked!" << std::endl;
    result = true;
  }
#else
  da.create(3);
  tester.out() << "Ack! Program should have aborted and never gotten here!"
               << std::endl;
#endif

  return result;
}


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0];
  tester.out() << ": DynamicArray dynamic ops w/views." << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  // Create an Interval object to create and index into an Array with

  tester.out() << "Creating an Interval<1> object ..." << std::endl;
  Interval<1> D1(3);
  tester.out() << "D1 = " << D1 << std::endl;

  // Create simple dynamic array.

  tester.out() << "Creating DynamicArray using domain ..." << std::endl;
  DynamicArray<int,Dynamic> a(D1);
  tester.check(a.domain().size() == D1.size());

  // Initialize dynamic array with scalar.

  a = 3;
  tester.out() << "Initialized DynamicArray to 3:" << std::endl;
  tester.out() << "a = " << a << std::endl;
  tester.check(sum(a) == (a.domain().size() * 3));

  // Create elements in the array.

  tester.out() << "Creating 2 elements at end of a ..." << std::endl;
  a.create(2);
  a.sync();
  tester.out() << "a = " << a << std::endl;
  tester.check(a.domain().size() == (D1.size() + 2));

#if 0
  // This test is bogous, as the comment in Engine/DynamicEngine.cpp::create()
  // tells (the check for shared data is commented out):
  //   "It would be nice to assert that no-one else is looking at the engine
  //    when we perform dynamic operations, but all the particle swap operations
  //    take place inside iterates which means the engine is a copy of another
  //    engine, so the data is shared."

  // Call a function which takes a view and the original DynamicArray

  std::cout << "The program should abort in the next operation when it\n";
  std::cout << "tries to create elements in an array with an existing view.";
  std::cout << std::endl;
  tester.out() << "Calling testview with a and a(1,3) ..." << std::endl;
  tester.check(testview(tester, a, a(Interval<1>(1,3))));
#endif

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("DynamicArray dynamic ops w/views");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_array_badcreate.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
