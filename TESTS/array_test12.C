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
// Array test 12: where()
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Array/ArrayOperators.h"

#include <cmath>

template <int Dim>
void check(Pooma::Tester& tester)
{
  tester.out() << Dim << "-dimensional tests:\n";
  Interval<Dim> I;
  for (int i=0; i<Dim; ++i)
    I[i] = Interval<1>(10);
  Array<Dim> a(I), b(I);
  a = 1.0;
  b = 0.0;
  b = where(a == 1.0, a);
  tester.check("2-arg where with array rhs", all(b == 1.0));
  b = 0.0;
  b = where(a == 1.0, 5.0);
  tester.check("2-arg where with scalar rhs", all(b == 5.0));
  b = 0.0;
  b = where(a == 1.0, a, a);
  tester.check("3-arg where with array/array rhs", all(b == 1.0));
  b = 0.0;
  b = where(a == 1.0, a, 3.0);
  tester.check("3-arg where with array/scalar rhs", all(b == 1.0));
  b = 0.0;
  b = where(a == 1.0, 3.0, a);
  tester.check("3-arg where with scalar/array rhs", all(b == 3.0));
  b = 0.0;
  b = where(a == 1.0, 1.0, 3.0);
  tester.check("3-arg where with scalar/scalar rhs", all(b == 1.0));
} 

int main(int argc, char* argv[])
{
  // Initialize Pooma.
  
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int n = 20;
  double pi = 3.1415926535897932;
  Array<1> a(n), b(n), c(n), d(n);
  int i;
  for (i = 0; i < n; ++i) 
  {
    a(i) = sin(0.1 * pi * i);
    b(i) = cos(0.1 * pi * i);
    c(i) = i * 1.0;
  }

  d = c;
  d = where( GT(a, b), where( LT(c, 6.5), 1.0, 0.0), d + 2.0);

  double compare[20] =
  {
    2.0,3.0,4.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,
    0.0,0.0,0.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0
  };
  double d2 = 0;
  double diff;
  
  Pooma::blockAndEvaluate();

  tester.out() << "Created arrays:" << std::endl;
  tester.out() << "  a = " << a << std::endl;
  tester.out() << "  b = " << b << std::endl;
  tester.out() << "  c = " << c << std::endl;
  tester.out() << "  d = " << d << std::endl;

  for (i = 0; i < n; ++i)
  {
    diff = d(i) - compare[i];
    d2 += diff * diff;
  }

  tester.out() << "Computed difference^2 from expected result = " << d2;
  tester.out() << std::endl;
  tester.check("d2 < 0.000001", d2 < 0.000001);

  d /= where( GT(d, 2.5), d );

  Pooma::blockAndEvaluate();

  for (i = 0; i < n; ++i)
  {
    if (compare[i] > 2.5)
    {
      diff = d(i) - 1.0;
    }
    else
    {
      diff = d(i) - compare[i];
    }

    d2 += diff * diff;
  }

  tester.out() << "Computed difference^2 from expected result = " << d2;
  tester.out() << std::endl;
  tester.check("d2 < 0.000001", d2 < 0.000001);

  int cnt = sum(where(d == 0.0, 1));
  tester.check("counting zeros with where reduction", cnt == 6);

  tester.check("where reduction", prod(where(d == 0.0, d)) == 0.0);

  // generic 2/3-arg where with array/scalar rhs

  check<1>(tester);
  check<2>(tester);
  check<3>(tester);
  
  int ret = tester.results("array_test12");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test12.cpp,v $   $Author: richi $
// $Revision: 1.17 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
