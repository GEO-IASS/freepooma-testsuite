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
// ----------------------------------------------------------------------
// Assertion test #1.
//-----------------------------------------------------------------------------

#include "Pooma/Arrays.h"
#include "Utilities/Tester.h"

#include <iostream>

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  Pooma::Tester tester(argc, argv);


  Array<1> x(7), y(7);    
  Array<1> z(6);

#if POOMA_EXCEPTIONS
  try {
    x = y + z;
  } 
  catch (const Pooma::Assertion & err) {
    err.print(tester.out()); tester.out() << std::endl;
  }
#endif
  
  int res = tester.results("assert_test1 " );
  Pooma::finalize();

  return res;
}
