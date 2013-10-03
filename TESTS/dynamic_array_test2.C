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
// DynamicArray test 2: Assignment with Array objects, using Brick-like engines
//-----------------------------------------------------------------------------

// include files 

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/DynamicArrays.h"
#include "Pooma/BrickArrays.h"


int main(int argc, char *argv[])
{
  int i;

  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0];
  tester.out() << ": DynamicArray <--> Array assignment." << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  // Create some Interval objects to create and index into Array's with

  tester.out() << "Creating Interval<1> objects ..." << std::endl;
  Interval<1> D1(3);
  Interval<1> D2(4);
  tester.out() << "D1 = " << D1 << std::endl;
  tester.out() << "D2 = " << D2 << std::endl;

  // Create simple single-patch dynamic arrays.

  tester.out() << "Creating DynamicArray objects ..." << std::endl;
  DynamicArray<int> a(D1);

  // Create simple Brick-based regular arrays.

  tester.out() << "Creating regular Array objects ..." << std::endl;
  Array<1,long>      b(D2);

  // Initialize dynamic array with scalar.

  a = 3;
  tester.out() << "Initialized DynamicArray a to the value 3." << std::endl;
  tester.out() << "a = " << a << std::endl;
  tester.check("Initially DynamicArray", sum(a) == (a.domain().size() * 3));

  // Initialize the regular array with scalars.
  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();

  tester.out() << "Initializing regular Array objects ..." << std::endl;
  for (i=0; i < b.domain().size(); ++i)
    {
      b(i) = i + 11;
    }
  tester.out() << "b = " << b << std::endl;

  // Resize a to the same size as b, and do operations.

  int oldsum = sum(a);
  tester.out() << "Resizing a to domain " << b.domain() << std::endl;
  a.create(1);
  a.sync();
  a(a.domain().size() - 1) = 1000;
  tester.out() << "a = " << a << std::endl;
  tester.check("Resize a sum", sum(a) == (oldsum + 1000));

  int suma = sum(a);
  long sumb = sum(b);
  tester.out() << "Trying a += b:" << std::endl;
  a += b;
  tester.out() << "a = " << a << std::endl;
  tester.out() << "b = " << b << std::endl;
  tester.check("a += b", sum(a) == (suma + sumb));

  tester.out() << "Trying b = a:" << std::endl;
  b = a;
  tester.out() << "a = " << a << std::endl;
  tester.out() << "b = " << b << std::endl;
  tester.check("b = a", sum(a) == sum(b));

  tester.out() << "Trying a = (b + b):" << std::endl;
  a = (b + b);
  tester.out() << "a = " << a << std::endl;
  tester.out() << "b = " << b << std::endl;
  tester.check("a = (b + b)", sum(a) == (sum(b) + sum(b)));

  tester.out() << "Trying a = (a + a) - b" << std::endl;
  suma = sum(a);
  a = (a + a) - b;
  tester.out() << "a = " << a << std::endl;
  tester.out() << "b = " << b << std::endl;
  tester.check("a = (a + a) - b", sum(a) == (2*suma - sum(b)));

  tester.out() << "Trying b = (a * b) + (b - a)" << std::endl;
  sumb = sum((a * b) + (b - a));
  b = (a * b) + (b - a);
  tester.out() << "a = " << a << std::endl;
  tester.out() << "b = " << b << std::endl;
  tester.check("b = (a * b) + (b - a)", sum(b) == sumb);

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("DynamicArray <--> Array expressions");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_array_test2.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
