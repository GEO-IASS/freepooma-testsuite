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
typedef Engine<1,double,BrickView>  View1_t;
typedef Engine<1,double,BrickView> View2_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTesting BrickEngine." << std::endl;

    Interval<1> I(10);

    Array_t A(I+5);

    for (int i = 5; i < 15; i++)
      {
	int ii = i-5;
	A(i) = 2.0 + ii - ii*ii;
      }
    
    for (int i = 5; i < 15; i++)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;

    Interval<1> J(2,5);

    View1_t B(A,J+5);

    for (int i = 0; i < 4; i++)
      tester.out() << B(i) << " ";
    tester.out() << std::endl;

    Range<1> K(1,9,2);
    View2_t C(A,K+5);
    
    for (int i = 0; i < 5; i++)
      tester.out() << C(i) << " ";
    tester.out() << std::endl;

    Array_t AC = A;

    AC(7) = -999;
    tester.out() << "AC(2) = " << AC(7) << std::endl;
    tester.out() << "A(2) = " << A(7) << std::endl;

    AC.makeOwnCopy();

    AC(12) = -111;
    tester.out() << "AC(2) = " << AC(12) << std::endl;
    tester.out() << "A(2) = " << A(12) << std::endl;

    tester.out() << "\nTesting BrickEngine<double,3>." << std::endl;

    Interval<3> III(I,I,I);
    Engine<3,double,Brick> AAA(III);

    int Imax = I.length();

    for (int i = 0; i < Imax; i++)
      for (int j = 0; j < Imax; j++)
	for (int k = 0; k < Imax; k++)
	  {
	    AAA(i,j,k) = i + j + k;
	  }

    for (int i = 0; i < Imax; i++)
      {
	tester.out() << "Slice i = " << i << std::endl;
	for (int j = 0; j < Imax; j++)
	  {
	    for (int k = 0; k < Imax; k++)
	      {
		tester.out() << std::setw(3) << AAA(i,j,k) << " ";
	      }
	    tester.out() << std::endl;
	  }
	tester.out() << std::endl;
      }

    Range<1> J2(2,8,2);

    Range<3> JJJ(J2,J2,J2);

    Engine<3,double,BrickView> AV(AAA,JJJ);

    int Jmax = J2.length();

    for (int i = 0; i < Jmax; i++)
      for (int j = 0; j < Jmax; j++)
	for (int k = 0; k < Jmax; k++)
	  {
	    AV(i,j,k) = -1;
	  }

    for (int i = 0; i < Imax; i++)
      {
	tester.out() << "Slice i = " << i << std::endl;
	for (int j = 0; j < Imax; j++)
	  {
	    for (int k = 0; k < Imax; k++)
	      {
		tester.out() << std::setw(3) << AAA(i,j,k) << " ";
	      }
	    tester.out() << std::endl;
	  }
	tester.out() << std::endl;
      }

    Range<1> J3(0,2,2);
    Range<1> J0(3);
    Range<3> JJJJ(J0,J3,J0);

    Engine<3,double,BrickView> AVV(AV,JJJJ);
     
    tester.out() << "Domain of AAA = " << std::endl << std::endl;
    tester.out() << AAA.domain() << std::endl << std::endl;

    Interval<3> avdom = AV.domain();

    tester.out() << "Domain of AV  = " << std::endl << std::endl;
    tester.out() << avdom << std::endl << std::endl;
    tester.out() << avdom[0].length() << std::endl << std::endl;

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
  int ret = tester.results("brick_test3");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brick_test3.cpp,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
