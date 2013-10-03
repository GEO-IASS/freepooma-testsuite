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
// test of RNGSimple
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Functions/RNGSimple.h"

#include <iostream>

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  RNGSimple r0, r1, r2;

  tester.out() << "some random numbers" << std::endl;

  int i;
  for (i = 0; i < 10; ++i)
  {
    tester.out() << i << ": "
		 << r0.value() << ","
		 << r1.value() << ","
		 << r2.value() << std::endl; 
    r0.advance();
    r1.advance();
    r2.advance();
  }

  tester.out() << r0.value() - r1.value() << std::endl;

  tester.check("same values", abs(r0.value() - r1.value()) < 0.00001);

  tester.out() << "different seeds!" << std::endl;
  r0.advance(1);
  r1.advance(2);
  r2.advance(3);

  for (i = 0; i < 10; ++i)
  {
    tester.out() << i << ": "
		 << r0.value() << ","
		 << r1.value() << ","
		 << r2.value() << std::endl;
    r0.advance();
    r1.advance();
    r2.advance();
  }

  tester.out() << r0.value() - r1.value() << std::endl;

  tester.check("different values", r0.value() != r1.value());

  int ret = tester.results("rngSimple");

  Pooma::finalize();

  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: rngSimple.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:51 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
