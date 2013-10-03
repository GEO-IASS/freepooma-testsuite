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
// Array test 25: initialize function.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Array<2> a(4, 4), b;
  
  a = 3.0;
  b.initialize(a);

  Pooma::blockAndEvaluate();

  b(2, 2) = -1.0;
  
  tester.out() << a << std::endl;
  tester.out() << b << std::endl;
  
  tester.check("simple", a(2,2), -1.0);

  Interval<1> x(0,5),y(0,5);
  Interval<2> dom(x,y);
  Array<2> xy(dom);
  xy = 0;

  Pooma::blockAndEvaluate();
  
  xy(3,3) = 303;

  b.initialize(xy.engine());

  tester.check("engine-init", b.domain(), xy.domain());

  Interval<1> x1(47,57);
  Interval<1> y1(5);
  Interval<1> z1(10);
  Interval<3> xyz(x1,y1,z1);
  Array<3,double> foo3;
  foo3.initialize(x1,y1,z1);

  tester.check("domain-init", foo3.domain(), xyz);

  int retval = tester.results("array_test25");

  Pooma::finalize();

  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test25.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
