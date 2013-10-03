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
// Array test 22: miscellaneous bugs that were reported.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Layout/DomainLayout.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"
#include "Tiny/Tensor.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> dom(1,20);
  DomainLayout<1> layout(dom);

  Array<1,double,Brick> a(dom), c(layout);
  typedef Array<1, Vector<2, double>, Brick> Array_t;
  typedef ComponentView<Loc<1>, Array_t>::Type_t CView_t; 
  Array_t b(dom);
  CView_t d(b.comp(1));

  a = 2.0;
  //  b = Vector<2,double>(2.0,3.0);
  b = Vector<2,double>(2.0, 1.0);
  d = 3.0;

  c = (a+2*b).comp(1);

  Pooma::blockAndEvaluate();

  tester.out() << "Created arrays:" << std::endl;
  tester.out() << "  a = " << a << std::endl;
  tester.out() << "  b = " << b << std::endl;
  tester.out() << "  c = " << c << std::endl;

  // Make sure that a particular element from c is OK

  tester.check("c(2) == 8", c(2) == 8.0);

  // check that assignment of a scalar to a vector field
  // compiles.

  b = 1.0;
  Pooma::blockAndEvaluate();

  tester.check("assigning scalar", b(2) == Vector<2,double>(1.0,1.0));

  typedef Array<1, Tensor<2, double, Antisymmetric>,
    Brick> Array2_t;
  typedef ComponentView<Loc<2>, Array2_t>::Type_t CView2_t; 
  Array2_t aa(dom);
  CView2_t bb(aa.comp(0, 1));

  //  aa.comp(0, 1) = 2.0;
  bb = 2.0;

  tester.out() << aa << std::endl;

  tester.check("antisymmetry", aa(3)(1, 0) == -2.0);
  
  int ret = tester.results("array_test22");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test22.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
