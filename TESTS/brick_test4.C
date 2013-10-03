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
// BrickEngine test code
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Engine/BrickEngine.h"

#include <iomanip>

typedef Engine<1,double,Brick> Array_t;
typedef Engine<2,double,Brick> Array2_t;
typedef Engine<1,double,BrickView>  View1_t;
typedef Engine<1,double,BrickView> View2_t;

void print(const Array2_t &, Pooma::Tester &tester);

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTesting BrickEngine with Fortran-like offsets."
              << std::endl;

    Interval<1> I(1,4);
    Interval<2> II(I,I);

    Array2_t a(II), b(II), c(II);

    for (int i = 1; i <=4; i++)
      for (int j = 1; j <= 4; j++)
	{
	  a(i,j) = i;
	  b(i,j) = j;
	}

    for (int i = 1; i <=4; i++)
      for (int j = 1; j <= 4; j++)
	c(i,j) = 0.0;

    for (int i = 1; i <=4; i++)
      for (int j = 1; j <= 4; j++)
	for (int k = 1; k <= 4; k++)
	  {
	    c(i,j) += a(i,k)*b(k,j);
	  }

    tester.out() << "\na = " << std::endl;

    print(a,tester);

    tester.out() << "\nb = " << std::endl;

    print(b,tester);

    tester.out() << "\nc = matmul(a,b) = " << std::endl;

    print(c,tester);
#if POOMA_EXCEPTIONS
  }
  catch(const char *err) 
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
  catch(const Pooma::Assertion &err)
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
#endif    
  int ret = tester.results("brick_test4");
  Pooma::finalize();
  return ret;
}

void print(const Array2_t &a, Pooma::Tester &tester)
{
  for (int i = 1; i <= 4; i++)
    {
      for (int j = 1; j <=4; j++)
	{
	  tester.out() << std::setw(4) << a(i,j) << " ";
	}
      tester.out() << std::endl;
    }
}

	     
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brick_test4.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
