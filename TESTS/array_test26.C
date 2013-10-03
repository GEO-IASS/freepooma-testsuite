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
// Array test 26: 2-argument min/max functions.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  int i, j;
  Array<2, int> a(4, 4), b(4, 4), c(4, 4);

  Pooma::blockAndEvaluate();

  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
      {
	a(i, j) = -2 + i + j;
	b(i, j) = 4 - i - j;
      }

  tester.out() << a << '\n' << std::endl;
  tester.out() << b << '\n' << std::endl;

  tester.out() << min(a, b) << '\n' << std::endl;
  tester.out() << max(a, b) << std::endl;

  c = min(a, b);

  Pooma::blockAndEvaluate();

  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
      tester.check(c(i, j) == (a(i, j) < b(i, j) ? a(i, j) : b(i, j)));

  c = max(a, b);

  Pooma::blockAndEvaluate();

  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
      tester.check(c(i, j) == (a(i, j) > b(i, j) ? a(i, j) : b(i, j)));

  int retval = tester.results("array_test26");

  Pooma::finalize();

  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test26.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
