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
// The POOMA Framework
//
// This program was prepared by the Regents of the University of California
// at Los Alamos National Laboratory (the University) under Contract No.
// W-7405-ENG-36 with the U.S. Department of Energy (DOE). The University has
// certain rights in the program pursuant to the contract and the program
// should not be copied or distributed outside your organization. All rights
// in the program are reserved by the DOE and the University. Neither the U.S.
// Government nor the University makes any warranty, express or implied, or
// assumes any liability or responsibility for the use of this software
//
// Visit http://www.acl.lanl.gov/POOMA for more details
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
#include "Utilities/Tester.h"

typedef Engine<1,double,CompressibleBrick> Array_t;
typedef Engine<1,double,CompressibleBrickView> VArray_t;
typedef Engine<1,double,BrickView>  View1_t;
typedef Engine<1,double,BrickView> View2_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  Interval<1> II(10);
  Array_t AA(II);
  AA.compressedReadWrite() = 4.14;

    for (int i=AA.layout().domain().first();
       i<AA.layout().domain().last();++i)
    tester.out() << AA.read(i)<< " ";
  tester.out()<<std::endl;

#if 0
  tester.out()<<" before create " << std::endl;

  AA.create(5);

  tester.out() << AA.layout().domain()<<std::endl;
#endif
  
  for (int i=AA.layout().domain().first();
       i<AA.layout().domain().last();++i)
    tester.out() << AA.read(i)<< " ";
  tester.out()<<std::endl;

  AA(6) = 9.9;
   
  tester.out() << " after modifying one element "
	       << "of the CompressibleBrick"
	       << std::endl;

  for (int i=AA.layout().domain().first();
       i<AA.layout().domain().last();++i)
    tester.out() << AA.read(i)<< " ";
  tester.out()<<std::endl;

#if 0
  AA.create(5);

  tester.out() << " AA layout domain " 
	       << AA.layout().domain()<<std::endl;
#endif

  for (int i=AA.layout().domain().first();
       i<AA.layout().domain().last();++i)
    tester.out() << AA.read(i)<< " ";
  tester.out()<<std::endl;
#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTesting CompressibleBrickEngine." 
		 << std::endl;

    Interval<1> I(10);
    Array_t A(I);
    Interval<1> J(2,5);
    VArray_t B(A,J);
    Range<1> K(1,9,2);
    VArray_t C(A,K);
    Range<1> L(1,4,3);
    VArray_t D(C,L);

    A.compressedReadWrite() = 3.14;
    tester.out() << A.read(3) << std::endl;
    tester.out() << B.read(2) << std::endl;
    tester.out() << C.read(3) << std::endl;

    int i;
    for (i = 0; i < 10; i++)
      A(Loc<1>(i)) = 2.0 + i - i*i;
    
    tester.out() << "A: ";
    for (i = 0; i < 10; i++)
      tester.out() << A.read(Loc<1>(i)) << " ";
    tester.out() << std::endl;
    tester.out() << "B: ";
    for (i = 0; i < 3; i++)
      tester.out() << B.read(i) << " ";
    tester.out() << std::endl;
    tester.out() << "C: ";
    for (i = 0; i < 5; i++)
      tester.out() << C.read(i) << " ";
    tester.out() << std::endl;
    tester.out() << "D: ";
    for (i = 0; i < 2; i++)
      tester.out() << D.read(Loc<1>(i)) << " ";
    tester.out() << std::endl;

    //    for (int i = 0; i < 4; i++)
    //      tester.out() << B(Loc<1>(i)) << " ";
    //    tester.out() << std::endl;
    
    Array_t AC = A;

    AC(3) = -999;
    tester.out() << "AC(3) = " << AC(3) << std::endl;
    tester.out() << "A(3) = " << A(3) << std::endl;
    tester.out() << "B(1) = " << B(1) << std::endl;
    tester.out() << "C(1) = " << C(1) << std::endl;
    tester.out() << "D(0) = " << D(0) << std::endl;

    AC.makeOwnCopy();

    AC(7) = -111;
    tester.out() << "AC(7) = " << AC(7) << std::endl;
    tester.out() << "A(7) = " << A(7) << std::endl;
    tester.out() << "C(3) = " << C(3) << std::endl;

    Array_t E(I);
    for (i = 0; i < 10; i++)
      E(i) = i;

    tester.out() << "E: ";
    for (i = 0; i < 10; i++)
      tester.out() << E(i) << " ";
    tester.out() << std::endl;

    Array_t F = E;
    
    tester.out() << "F == E" << std::endl;
    tester.out() << "F: ";
    for (i = 0; i < 10; i++)
      tester.out() << F(i) << " ";
    tester.out() << std::endl;

    {
      View1_t G(A);
      for (i = 0; i < 10; i++)
	G(i) = 3.4;
      tester.out() << "A.compressed(): " 
		   << (A.compressed() ? "true" : "false")
		   << std::endl;
    }

    tester.out() << "A.compressed(): " 
		 << (A.compressed() ? "true" : "false")
		 << std::endl;

    tester.out() << "C: ";
    for (i = 0; i < 5; i++)
      tester.out() << C.read(i) << " ";
    tester.out() << std::endl;
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
  int ret = tester.results("compbrick_test1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: compbrick_test1.cpp,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
