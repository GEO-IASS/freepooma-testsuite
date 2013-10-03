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
// Array test 24: elementwise tensor and vector operations.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Tiny/Tensor.h"

#include <iostream>
#include <cmath>

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Elementwise tensor and vector tests.."
	       << std::endl;
  tester.out() << "------------------------------------------------"
	       << std::endl;

  Array<2, Tensor<2> > a(4, 4), b(4, 4);
  Tensor<2> t(1.0, 2.0, 3.0, 4.0);
  
  b = t;
  a = sqrt(b - 1) * sqrt(b + 1);
  
  tester.out() << a << std::endl;

  tester.out() << "------------------------------------------------"
	       << std::endl;

  int retval = tester.results("array_test24");

  Pooma::finalize();

  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test24.cpp,v $   $Author: richi $
// $Revision: 1.6 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
