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
// Dynamic operations on Dynamic engines
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Engine/DynamicEngine.h"
#include "Utilities/Tester.h"

#include <algorithm>

int main(int argc, char *argv[])
{
  int i;

  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  tester.out() << "\nTesting Dynamic Engine." << std::endl;

  Interval<1> I(20);
  tester.out() << " Testing the dynamic aspects of DynamicEngine " 
               << std::endl;
  tester.out() << " Arrays are defined on the interval " << I << std::endl;

  typedef Engine<1,double,Dynamic> Array_t;

  tester.out() << "\n Testing the destroy function " << std::endl;
  {
    Array_t A(I);
    Array_t B(I);
    Array_t C(I);
    Array_t D(I);

    for (i = 0; i < I.length(); ++i)
      {
	A(i) = i;
	B(i) = i;
	C(i) = i;
	D(i) = i;
      }
  
    tester.out()<< " Array A is: " << std::endl;
    for (i = 0; i < I.length(); ++i)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;

    Range<1> kill_list(3,9,2); // 3,5,7,9
    A.destroy(kill_list,BackFill());

    tester.out() << "A's new length = " << A.domain().length() << std::endl;
    tester.check(A.domain().length() == 16);

    tester.out() << " Array A after destroying " << kill_list
		 << " and BackFill() " << std::endl;

    for (i = 0; i < A.domain().last()+1; ++i)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;

    double Atst[] = {0, 1, 2, 16, 4, 17, 6, 18, 8, 19, \
                     10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
                        
    bool test = std::equal(Atst,Atst+A.domain().length(),&A(0));
    tester.check(test);
    
    B.destroy(kill_list,ShiftUp());

    tester.out() << "B's new length = " << B.domain().length() << std::endl;
    tester.check(B.domain().length() == 16);
    
    tester.out() << " Array B after destroying " << kill_list
		 << " and ShiftUp() " << std::endl;
    for (i = 0; i < B.domain().last()+1; ++i)
      tester.out() << B(i) << " ";
    tester.out() << std::endl;

    double Btst[] = {0, 1, 2, 4, 6, 8, 10, 11, 12, 13, 
                     14, 15, 16, 17, 18, 19};

    test = std::equal(Btst,Btst+B.domain().length(),&B(0));
    tester.check(test);
    
    int kill_array[8] = {0, 1, 5, 6, 7, 14, 18, 19};
    
    C.destroy(kill_array, kill_array+8, BackFill());
    
    tester.out() << "C's new length = " << C.domain().length() << std::endl;
    tester.check(C.domain().length() == 12);
    
    tester.out() << " Array C after destroying [";
    
    InformIterator<int> Out(tester.out(), ", ");
    std::copy(kill_array, kill_array+7, Out);
    tester.out() << kill_array[7] << "], and BackFill() " << std::endl;
    
    for (i = 0; i < C.domain().last()+1; ++i)
      tester.out() << C(i) << " ";
    tester.out() << std::endl;
    
    double Ctst[] = {12, 13, 2, 3, 4, 17, 15, 16, 8, 9, 
                     10, 11};

    test = std::equal(Ctst,Ctst+C.domain().length(),&C(0));
    
    D.destroy(kill_array, kill_array+8, ShiftUp());
    
    tester.out() << "D's new length = " << D.domain().length() << std::endl;
    tester.check(D.domain().length() == 12);
    
    tester.out() << " Array D after destroying [";
    std::copy(kill_array, kill_array+7, Out);
    tester.out() << kill_array[7] << "], and ShiftUp() " << std::endl;
    
    for (i = 0; i < D.domain().last()+1; ++i)
      tester.out() << D(i) << " ";
    tester.out() << std::endl;
    
    double Dtst[] = {2, 3, 4, 8, 9, 10, 11, 12, 13, 15, 16, 17};

    test = std::equal(Dtst,Dtst+D.domain().length(),&D(0));
    

#if 0	
  // DynamicEngine requires kill list to be in increasing order!
    Range<1> backwardskilllist(9,3,-2);
    tester.out() << " destroy the elements specified by: " << backwardskilllist
		 << " and BackFill()" << std::endl;
    C.destroy(backwardskilllist,BackFill());
	 
    for (i = 0; i < C.domain().last()+1; ++i)
      tester.out() << C(i) << " ";
    tester.out() << std::endl;

    tester.out() << " destroy the elements specified by: " << backwardskilllist
		 << " and ShiftUp()" << std::endl;
    D.destroy(backwardskilllist,ShiftUp());
  
    for (i = 0; i < D.domain().last()+1; ++i)
      tester.out() << D(i) << " ";
    tester.out() << std::endl;
#endif
  }
	
  tester.out() << "\n Testing the create(int num)" << std::endl;
  {
    Array_t A(I);
    Array_t B(I);
    Array_t C(I);
    Array_t D(I);

    for (i = 0; i < I.length(); ++i)
    {
      A(i) = i;
      B(i) = i;
      C(i) = i;
      D(i) = i;
    }

    tester.out() << " A.create(3) " << std::endl;

    A.create(3);
    
    tester.check(A.domain().length() == 23);
    
    for (i = 0; i < 3; ++i)
      A(Loc<1>(I.length() + i)) = -(I.length() + i);

    for (i = 0; i < A.domain().last()+1; ++i)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;

    double Atst[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                     15, 16, 17, 18, 19, -20, -21, -22};
                        
    bool test = std::equal(Atst,Atst+A.domain().length(),&A(0));
    tester.check(test);

    Range<1> rend(4,16,4); // 4, 8, 12, 16

    tester.out() << " in A, destroy " << rend << " BackFill()  " << std::endl;
    A.destroy(rend,BackFill());

    tester.check(A.domain().length() == 19);
    
    for (i = 0; i < A.domain().last()+1; ++i)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;
    
    double Btst[] = {0, 1, 2, 3, 19, 5, 6, 7, -20, 9, 10, 11, 
                     -21, 13, 14, 15, -22, 17, 18};

    test = std::equal(Btst,Btst+A.domain().length(),&A(0));
    tester.check(test);
    
    Interval<1> middle(5,7);
    
    A.destroy(middle,BackFill());
    
    tester.check(A.domain().length() == 16);
    
    for (i = 0; i < A.domain().last()+1; ++i)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;
    
    double Ctst[] = {0, 1, 2, 3, 19, -22, 17, 18, -20, 9, 10, 11, 
                     -21, 13, 14, 15};

    test = std::equal(Ctst,Ctst+A.domain().length(),&A(0));
    tester.check(test);
    
    A.destroy(middle,ShiftUp());
    
    tester.check(A.domain().length() == 13);
    
    for (i = 0; i < A.domain().last()+1; ++i)
      tester.out() << A(i) << " ";
    tester.out() << std::endl;
    
    double Dtst[] = {0, 1, 2, 3, 19, -20, 9, 10, 11, 
                     -21, 13, 14, 15};

    test = std::equal(Dtst,Dtst+A.domain().length(),&A(0));
    tester.check(test);
    
  }

  int ret = tester.results("dynamic_test2");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_test2.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
