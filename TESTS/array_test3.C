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
// Array test 3: Vector array elements.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"

#include <iostream>

static int checkNum = 1;

template<class T>
inline void check(const T &ans, const T &correct, Pooma::Tester &tester)
{
  bool ok = (ans == correct);
  tester.out() << "Check #" << checkNum << std::endl;
  tester.out() << "Answer:  " << ans << std::endl;
  if (!ok)
    {
      tester.out() << "Correct: " << correct << std::endl;
    }

  checkNum++;      
  tester.check(ok);
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Array<2, int> a(2, 2, modelElement(7));
  Array<2, Vector<3> > b(2, 2, modelElement(Vector<3>(1.0, 2.0, 3.0))), 
    c(2, 2);
  Array<2> d(2,2);
  int i, j;

  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(a(i,j), 7, tester);
        check(b(i,j), Vector<3>(1.0, 2.0, 3.0), tester);
      }
  
  b.comp(1) = 6.0;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(b(i,j), Vector<3>(1.0, 6.0, 3.0), tester);
      }

  b.comp(0) = a + b.comp(1) + b.comp(2);

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(b(i,j), Vector<3>(16.0, 6.0, 3.0), tester);
      }

  c = a + b;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(c(i,j), Vector<3>(23.0, 13.0, 10.0), tester);
      }

  c = a + 2 * b;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(c(i,j), Vector<3>(39.0, 19.0, 13.0), tester);
      }

  Vector<3> x(-1,-2,-3);
  c = a + x * b;
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(c(i,j), Vector<3>(-9.0, -5.0, -2.0), tester);
      }

  c = a + b * x;
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(c(i,j), Vector<3>(-9.0, -5.0, -2.0), tester);
      }

  d = a + dot(x,b);
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(d(i,j), -30.0, tester);
      }

  d = a - dot(c,b);
  
  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(d(i,j), 187.0, tester);
      }

  b.comp(0) = a + b.comp(1) + b.comp(2) - 1;

  Pooma::blockAndEvaluate();
  for (j = 0; j < 2; j++)
    for (i = 0; i < 2; i++)
      {
        check(b(i,j), Vector<3>(15.0, 6.0, 3.0), tester);
      }
         
  int ret = tester.results( "array_test3" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test3.cpp,v $   $Author: richard $
// $Revision: 1.19 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
