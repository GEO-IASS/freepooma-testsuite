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
// Array test 6: global reductions.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/AllDomain.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Array<2, int> a(2, 2);
  Array<2, int> b(2, 2);
  bool OK = true;

  b = 0;
  
  Pooma::blockAndEvaluate();
  a(0, 0) = 1; a(0,1) = 2; a(1, 0) = 3; a(1, 1) = 4;
  b(1,0) = 1;
  
  OK = (OK && sum(a) == 10);
  OK = (OK && min(a(AllDomain<1>(), 1)) == 2);
  OK = (OK && max(3.0 * a) == 12);
  OK = (OK && prod(a + b) == 32);
  OK = (OK && all(b) == false);
  OK = (OK && any(b) == true);
  OK = (OK && bitOr(a) == 7);
  OK = (OK && bitAnd(a) == 0);

  tester.check(OK);

  int ret = tester.results( "array_test6" );
  Pooma::finalize();
  return ret; 
}
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test6.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
