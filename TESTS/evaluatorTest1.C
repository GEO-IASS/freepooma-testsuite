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
// Compressed Evaluation.
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Engine/BrickEngine.h"
#include "Pooma/Arrays.h"
#include "Utilities/Clock.h"
#include "Utilities/Tester.h"

#include <iostream>

static bool printStuff = false;

template<class A1,class A2>
bool verify(int line1, const A1 &a1, int line2, const A2 &a2,Pooma::Tester & tester)
{
  Pooma::blockAndEvaluate();
  bool passed    =  isSmall(a1-a2);
  if (printStuff)
  {
    if (!passed)
    {
      tester.out() << "Failure: line #" << line1 << " != line #" << line2
                << std::endl;
    }
  }
  return (passed);
}


template<class Array>
bool isSmall(const Array &a)
{
  double epsilon = 0.000001;
  int i;
  int first = a.domain()[0].first();
  int last  = a.domain()[0].last();
  double sum = 0.0;
  for (i = first; i <= last; ++i)
  {
    sum += a.read(i) * a.read(i);
  }
  return (sum < epsilon);
}

typedef Array<1,double,Brick> AB_t;
typedef Array<1,double,CompressibleBrick> AC_t;
typedef Array<1,Vector<2,double>,CompressibleBrick> ACV_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  // if the user runs this code with options then print stuff.
  printStuff = (argc > 1);

  int return_status = 0;

  bool worked = true;
  bool check;

  int from = 1;
  int to = 20000;
  int fromInterior = from+1;
  int toInterior = to-1;

  Interval<1> dom(from,to);
  Interval<1> I(fromInterior,toInterior);
  int i;
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
  double t11,t12,t13,t14,t15,t16,t17;

  AB_t  ba(dom),bb(dom),b1(dom),b2(dom),b3(dom);
  AC_t  ca(dom),cb(dom),c1(dom),c2(dom),c3(dom);
  AC_t  c4(dom),c5(dom),c6(dom);
  ACV_t da(dom),db(dom),d1(dom),d2(dom),d3(dom);
  ACV_t d4(dom),d5(dom);
  AC_t d6(dom);

  for (i = from; i <= to; ++i)
  {
    ba(i) = i+2;
    bb(i) = 3;
    ca(i) = i+2;
    da(i) = Vector<2,double>(i+2,i+2);
  }
  cb = 0.0;
  db = Vector<2,double>(0.0,0.0);

  Pooma::blockAndEvaluate();
  cb.engine().compressedReadWrite() = 3.0;
  db.engine().compressedReadWrite() = Vector<2,double>(3.0,3.0);
  Pooma::blockAndEvaluate();

  // line #1
  t1 = Pooma::Clock::value();
  b1 = ba+bb;
  Pooma::blockAndEvaluate();

  // line #2
  t2 = Pooma::Clock::value();
  b2 = ba+cb;
  Pooma::blockAndEvaluate();

  // line #3
  t3 = Pooma::Clock::value();
  b3 = ca+cb;
  Pooma::blockAndEvaluate();

  // line #4
  t4 = Pooma::Clock::value();
  c1 = ba+bb;
  Pooma::blockAndEvaluate();

  // line #5
  t5 = Pooma::Clock::value();
  c2 = ba+cb;
  Pooma::blockAndEvaluate();

  // line #6
  t6 = Pooma::Clock::value();
  c3 = ca+cb;
  Pooma::blockAndEvaluate();

  // line #7
  t7 = Pooma::Clock::value();
  c4 = bb+bb;
  Pooma::blockAndEvaluate();

  // line #8
  t8 = Pooma::Clock::value();
  c5 = bb+cb;
  Pooma::blockAndEvaluate();

  // line #9
  t9 = Pooma::Clock::value();
  c6 = cb+cb;
  Pooma::blockAndEvaluate();

  t10 = Pooma::Clock::value();
  Pooma::blockAndEvaluate();

  // line #10
  t11 = Pooma::Clock::value();
  d1.comp(0) = ba+bb;
  Pooma::blockAndEvaluate();

  // line #11
  t12 = Pooma::Clock::value();
  d2.comp(0) = ba+db.comp(0);
  Pooma::blockAndEvaluate();

  // line #12
  t13 = Pooma::Clock::value();
  d3 = da+db;
  Pooma::blockAndEvaluate();

  // line #13
  t14 = Pooma::Clock::value();
  d4.comp(0) = bb+bb;
  Pooma::blockAndEvaluate();

  // line #14
  t15 = Pooma::Clock::value();
  d5.comp(0) = bb+db.comp(0);
  Pooma::blockAndEvaluate();

  // line #15
  t16 = Pooma::Clock::value();
  d6 = db.comp(0)+db.comp(0);
  Pooma::blockAndEvaluate();

  t17 = Pooma::Clock::value();
  Pooma::blockAndEvaluate();

  bool check1 = c6.engine().compressed();
  worked = worked && check1;
  if (printStuff && !check1)
  {
    tester.out() << "c6 is not compressed!" << std::endl;
  }

  check1 = d6.engine().compressed();
  worked = worked && check1;
  if (printStuff && !check1)
  {
    tester.out() << "d6 is not compressed!" << std::endl;
  }

  check1 = engineFunctor((db.comp(0)).engine(),Compressed());
  worked = worked && check1;
  if (printStuff && !check1)
  {
    tester.out() << "db.comp(0) is not compressed!" << std::endl;
  }

  check = verify(1,b1,4,c1,tester);
  worked = worked && check;
  check = verify(2,b2,5,c2,tester);
  worked = worked && check;
  check = verify(3,b3,6,c3,tester);
  worked = worked && check;
  check = verify(7,c4,8,c5,tester);
  worked = worked && check;
  check = verify(9,c6,8,c5,tester);
  worked = worked && check;
 
  check = verify(1,b1,11,d1.comp(0),tester);
  worked = worked && check;
  /*
  check = verify(2,b2,5,c2);
  worked = worked && check;
  check = verify(3,b3,6,c3);
  worked = worked && check;
  check = verify(7,c4,8,c5);
  worked = worked && check;
  check = verify(9,c6,8,c5);
  worked = worked && check;
  */
 
  double cval = (t10-t9)/(t7-t6);
  if (cval > 0.1)
  {
    tester.out() << "warning! compressed eval took " << cval
	      << " times ordinary eval" << std::endl;
  }

  cval = (t17-t16)/(t14-t13);
  if (cval > 0.1)
  {
    tester.out() << "warning! compressed eval took " << cval
	      << " times ordinary eval" << std::endl;
  }

  if (worked)
    tester.out() << "PASSED" << std::endl;
  else
    tester.out() << "FAILED" << std::endl;

  tester.results("evaluatorTest1 " );


  Pooma::finalize();


  return return_status;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: evaluatorTest1.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:41 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
