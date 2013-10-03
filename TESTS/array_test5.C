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
// Array test 5: complex array elements.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Arrays.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"

#include <complex>
#include <cmath>

template<class T>
inline void check(const T &ans, const T &correct, Pooma::Tester &tester)
{
  tester.check(ans == correct);
}

template<class T>
inline void floatCheck(const T &ans, const T &correct, Pooma::Tester &tester)
{
  tester.check(std::abs(ans - correct) < 1.0e-08);
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  complex<double> x(1.0, 2.0);
  Array<2> a(2, 2, modelElement(7.0));
  Array<2, complex<double> > b(2, 2, modelElement(x)), c(2, 2);
  Array<2> d(2,2);
  int i, j;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        floatCheck(a(i,j), 7.0, tester);
        floatCheck(b(i,j), x, tester);
      }

  c = a + 2.0 * b;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        floatCheck(c(i,j), complex<double>(9.0,4.0), tester);
      }

  complex<double> y(-3, -4);
  c += a + y * conj(b);
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        floatCheck(c(i,j), complex<double>(5.0,6.0), tester);
      }

  d = norm(a + y * conj(b));

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        floatCheck(d(i,j), 20.0, tester);
      }

  d = real(y * pow(b, 2));

  bool OK = true;
  Pooma::blockAndEvaluate();  
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        OK = (OK && (fabs(d(i,j) - 25.0) < 1e-6));
      }
  check(OK, true, tester);
  
  Array<1, complex<double> > e(2);
  Array<1, Vector<2, complex<double> > > f(2), g(2);
  Vector<2, complex<double> > v(complex<double>(1.0, 2.0),
    complex<double>(3.0, 4.0));
  Vector<2, complex<double> > v1(complex<double>(-3.0, -1.0),
    complex<double>(-7.0, -1.0));
  e = complex<double>(-1.0, 1.0);
  f = v;
  g = f * e;

  tester.check(all(g == v1));
  
  int ret = tester.results( "array_test5" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test5.cpp,v $   $Author: richi $
// $Revision: 1.16 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
