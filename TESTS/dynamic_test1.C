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
// DynamicEngine test code
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Engine/DynamicEngine.h"

typedef Engine<1,double,Dynamic> Array_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTesting DynamicEngine." << std::endl;

    Interval<1> I(10);
    Array_t A(I);

    for (int i = 0; i < 10; i++)
      A(Loc<1>(i)) = 2.0 + i - i*i;
    
    for (int i = 0; i < 10; i++)
      tester.out() <<  A(Loc<1>(i)) << " ";
    tester.out() << std::endl;

    Interval<1> J(2,5);
    NewEngine<Array_t,Interval<1> >::Type_t B(A,J); 
    
    for (int i = 0; i < 4; i++)
      {
        tester.check(B.read(i) == A.read(i + 2));
        tester.out() <<  B.read(Loc<1>(i)) << " ";
      }
    tester.out() << std::endl;

    Range<1> K(1,9,2);
    NewEngine<Array_t, Range<1> >::Type_t C(A,K); 
    
    tester.check(C.stride() == 2);
    
    for (int i = 0; i < 5; i++)
      {
        tester.check(
          C(i) == A(i * 2 + 1) );
        tester.out() << C.read(i) << " ";
      }
    tester.out() << std::endl;

    A(3) = -444;
    tester.check(A(3) == -444);
    tester.check(B(1) == -444);
    tester.check(C(1) == -444);
    
    Range<1> KV(0,4,2);
    typedef NewEngine<Array_t, Range<1> >::Type_t View1_t;
    
    NewEngine<View1_t,Range<1> >::Type_t CV(C,KV); // A(1), A(5), A(9)
                                                   // C(0), C(2), C(4)
                                                   //       B(3)
    tester.check(CV(0) == A(1));
    tester.check(CV(1) == B(3));

    Interval<1> IV(0,4);
    NewEngine<View1_t,Range<1> >::Type_t CV2(C,IV); 

    tester.check(CV2(0) == A(1));
    tester.check(CV2(1) == A(3));

    Array_t AC = A;

    AC(Loc<1>(2)) = -999;
    
    tester.check(AC(2) == -999);
    tester.check(A(2)  == -999);
    tester.check(AC.read(2) == A.read(2));
    tester.check(B(0) == -999);
    tester.out() << "AC(2) = " << AC(Loc<1>(2)) << std::endl;
    tester.out() << "A(2) = " << A(Loc<1>(2)) << std::endl;

    tester.check(A.isShared());
    tester.check(AC.isShared());
    
    AC.makeOwnCopy();
    tester.check(A.isShared());
    tester.check(!AC.isShared());
    
    double save = AC.read(7);
    A(7) = -111;
    
    tester.check(A.read(7) == -111);
    tester.check(AC.read(7)  == save);
    tester.check(C(3) == -111);
    tester.out() << "AC(2) = " << AC(7) << std::endl;
    tester.out() << "A(2) = " << A(7) << std::endl;

    Array_t E(I);
    for (int i = 0; i < 10; i++)
      E(i) = i;

    tester.out() << "E: ";
    for (int i = 0; i < 10; i++)
      tester.out() << E.read(i) << " ";
    tester.out() << std::endl;

    Array_t F;
    tester.check(F.domain().size() == 0);
    
    F = E;
    
    tester.check(F.isShared());
    tester.check(E.isShared());
    
    tester.out() << "F == E" << std::endl;
    tester.out() << "F: ";
    for (int i = 0; i < 10; i++)
      {
        tester.out() << F(i) << " ";
        tester.check(F(i) == E(i));
      }
    tester.out() << std::endl;

    Array_t G(I);
    for (int i = 0; i < 10; i++)
      G(i) = i*i;

    tester.out() << "G: ";
    for (int i = 0; i < 10; i++)
      tester.out() << G(i) << " ";
    tester.out() << std::endl;
    
    tester.check(!G.isShared());
    
    E = G;

    tester.check(E.isShared());
    tester.check(G.isShared());
    tester.check(!F.isShared());
    
    tester.out() << "E = G;" << std::endl;
    tester.out() << "E: ";
    for (int i = 0; i < 10; i++)
      {
        tester.check(E(i) == G(i));
        tester.out() << E(i) << " ";
      }
    tester.out() << std::endl;

    AC(Loc<1>(2)) = -222;
    tester.check(AC(2) == -222);
    tester.check(A(2)  == -999);
    tester.check(B(0)  == -999);
    tester.out() << "AC(2) = " << AC(Loc<1>(2)) << std::endl;
    tester.out() << "A(2) = " << A(Loc<1>(2)) << std::endl;

    A.makeOwnCopy();
    
    tester.check(!A.isShared());
    tester.check(B.dataBlock().isShared());
    tester.check(C.dataBlock().isShared());
    
    B(1) = -888;
    tester.check(A(3) == -444);
    tester.check(B(1) == -888);
    tester.check(C(1) == -888);
    
    B(3) = -555;
    tester.check(C(2) == -555);
    tester.check(B(3) == -555);
    tester.check(A(5) != -555);
    tester.check(CV(1) == -555);

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
  int ret = tester.results("dynamic_test1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_test1.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
