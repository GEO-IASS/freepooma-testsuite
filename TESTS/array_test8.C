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
// Array test 8: conformance tests.
//
// Note: exceptions are not thread safe so this program may not work in
//       in parallel.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Array<1> x(7), y(7), z(6);
  
  y = 0;
  z = 0;

  // This test depends on bounds checking throwing an
  // exception.  For it to work, exceptions must be turned
  // on, bounds checking must be on and we can't be running
  // in parallel.

#if ! POOMA_THREADS
#if POOMA_EXCEPTIONS
#if POOMA_BOUNDS_CHECK

  bool worked = false;
  try {
    x = y + z;
  }
  catch (const Pooma::Assertion &) {
    worked = true;
  }
  tester.check(worked);

  Array<3> a(4,5,6), b(4,5,6), c(4,4,6), d(4,5,6);
  
  b = 0;
  c = 0;
  d = 0;
  
  worked = false;
  try {
    a = b + 3 * b + c - sin(d);
  }
  catch (const Pooma::Assertion &) {
    worked = true;
  }

  tester.check(worked);

#endif
#endif
#endif

  int ret = tester.results( "array_test8" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test8.cpp,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
