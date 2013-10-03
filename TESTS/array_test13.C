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
// Array test 13: bounds checking
//
// Note: exceptions are not thread safe so this program may not work in
//       in parallel.
//-----------------------------------------------------------------------------

// Make sure bounds checking is on.

#ifdef POOMA_BOUNDS_CHECK
#undef POOMA_BOUNDS_CHECK
#endif

#ifdef POOMA_BOUNDS_CHECK_DEFAULT
#undef POOMA_BOUNDS_CHECK_DEFAULT
#endif

#define POOMA_BOUNDS_CHECK POOMA_YES
#define POOMA_BOUNDS_CHECK_DEFAULT POOMA_TRUE

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"

#include <cmath>


int main(int argc, char* argv[])
{
  // Initialize Pooma.
  
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  // This test depends on exceptions being present.

#if ! POOMA_THREADS
#if POOMA_EXCEPTIONS

  int cnt = 0;

  int n = 10;
  Array<3> a(n, n, n);
  Array<2, Vector<3> > b(n, n);
  
  try
    {
      a(-1, 0, 0) = 3.0;
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }
    
  try
    {
      double d = a(-1, 0, 0);
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }
    
  try
    {
      a(0, n, 0) = 3.0;
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }

  try
    {
      double d = a(0, n, 0);
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }

  try
    {
      Interval<1> I(n+1);
      a(I, 0, 0) = 3.0;
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }

  try
    {
      Interval<1> I(n+1);
      Array<1> v(I);
      v = a(I, 0, 0);
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }

  try
    {
      Interval<1> I(0,0);
      b(I,I).comp(4) = 3.0;
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }

  try
    {
      Interval<1> I(0,0);
      Array<2> v(I,I);
      v = b(I,I).comp(4);
    }
  catch(Pooma::Assertion &as)
    {
      tester.exceptionHandler( as );
      cnt++;
    }

  if (cnt != 8)
    tester.set(false);

#endif
#endif

  int ret = tester.results("array_test13");
  Pooma::finalize();
  return ret;

}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test13.cpp,v $   $Author: richi $
// $Revision: 1.12 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
