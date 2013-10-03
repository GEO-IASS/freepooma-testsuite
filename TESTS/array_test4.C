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
// Array test 4: TinyMatrix array elements.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Tiny/TinyMatrix.h"
#include "Tiny/Vector.h"
#include "Tiny/VectorTinyMatrix.h"

#include <cmath>

template<class T>
inline void check(const T &ans, const T &correct, Pooma::Tester &tester)
{
  tester.check(fabs(ans - correct) <= 1e-6);
}

template<int D1, int D2, class T, class E>
inline void check(const TinyMatrix<D1,D2,T,E> &ans, 
  const TinyMatrix<D1,D2,T,E> &correct, Pooma::Tester &tester)
{
  for (int i = 0; i < D1; i++)
    for (int j = 0; j < D2; j++)
      tester.check(fabs(ans(i,j)-correct(i,j)) <= 1e-6);
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  TinyMatrix<3,3> x(0.0, 1.0, 2.0, 0.1, 1.1, 2.1, 0.2, 1.2, 2.2);
  Array<2, int> a(2, 2, modelElement(7));
  Array<2, TinyMatrix<3,3> > b(2, 2, modelElement(x)), c(2, 2);
  Array<2> d(2,2);
  int i, j;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
	check((double)a(i,j), 7.0, tester);
        check(b(i,j), x, tester);
      }

  b.comp(1,2) = 6.0;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(b(i,j), TinyMatrix<3,3>(0.0,1.0,2.0,0.1,1.1,2.1,0.2,6,2.2), tester);
      }

  b.comp(0,1) = a + b.comp(1,0) + b.comp(2,1);

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(b(i,j), TinyMatrix<3,3>(0.0,1.0,2.0,10.1,1.1,2.1,0.2,6,2.2), tester);
      }

  c = a + 2 * b;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(c(i,j), TinyMatrix<3,3>(7.0,9.0,11.0,27.2,9.2,11.2,7.4,19.0,11.4),
	      tester);
      }

  TinyMatrix<3,3> y(-1,-2,-3,1,2,3,-1,-2,-3);
  c = a + y * b;
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(c(i,j), TinyMatrix<3,3>(7.0,5.0,1.0,17.1,9.2,13.3,6.8,-5,0.4), tester);
      }

  Vector<3> z(3,4,5);
  d = a + dot(z,dot(b,z));
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(d(i,j), 407.8, tester);
      }

  int ret = tester.results( "array_test4" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test4.cpp,v $   $Author: richi $
// $Revision: 1.18 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
