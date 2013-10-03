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
// Big Expression Test Code
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/AllDomain.h"
#include "Engine/BrickEngine.h"
#include "Engine/CompressibleBrick.h"
#include "Engine/IndirectionEngine.h"
#include "Tiny/Vector.h"
#include "Array/Array.h"
#include "Array/tests/ExpressionTest.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int from = 1;
  int to = 20;
  int fromInterior = 2;
  int toInterior = 19;

  Interval<1> dom(from,to);
  Interval<1> I(fromInterior,toInterior);
  int loc1 = 4;
  int loc2 = 12;

  Array<1,double,Brick>
    a1(dom),a2(dom),a3(dom),a4(dom), initial(dom);

  initial = 0.0;

  Pooma::blockAndEvaluate();

  initial(loc1) = 2.0;
  initial(loc2) = 3.0;

  test1(tester, 1, a1, a2, a3, a4, initial, I);
  test2(tester, 2, a1, a2, a3, a4, initial, I);
  test3(tester, 3, a1, a2, a3, a4, initial, I);
  test4(tester, 4, a1, a2, a3, a4, initial, I);

  Array<1,Vector<2,double>,Brick>
    av1(dom),av2(dom),av3(dom),av4(dom), initialv(dom);

  initialv = Vector<2,double>(0.0,0.0);

  Pooma::blockAndEvaluate();

  initialv(4) = Vector<2,double>(2.0,3.0);
  initialv(12) = Vector<2,double>(3.0,-1.0);

  test5(tester,5,av1,av2,av3,av4, initialv,I);

  Array<1,double,CompressibleBrick>
    ac1(dom),ac2(dom),ac3(dom),ac4(dom);

  test1(tester, 6, ac1, ac2, ac3, ac4, initial, I);
  test2(tester, 7, ac1, ac2, ac3, ac4, initial, I);
  test4(tester, 9, ac1, ac2, ac3, ac4, initial, I);

  Array<1,Vector<2,double>,CompressibleBrick>
    avc1(dom),avc2(dom),avc3(dom),avc4(dom);

  test5(tester, 10, avc1, avc2, avc3, avc4, initialv, I);
  test6(tester, 11, avc1, avc2, avc3, avc4, initialv, I);

  //---------------------------------------------------------------------
  // simple indirection test - rotate some values

  test7(tester, 12, a1, a2, a3, a4, initial, dom);

  //---------------------------------------------------------------------
  // slices  

  Interval<2> dom2(dom,dom);

  Array<2,double,Brick>
    a21(dom2),a22(dom2),a23(dom2),a24(dom2);

  AllDomain<1> all;

  test1(tester,13,a21(all,3),a22(all,3),a23(all,3),a24(all,3), initial,I);
  test2(tester,14,a21(all,3),a22(all,3),a23(all,3),a24(all,3), initial,I);
  test3(tester,15,a21(all,3),a22(all,3),a23(all,3),a24(all,3), initial,I);
  test4(tester,16,a21(all,3),a22(all,3),a23(all,3),a24(all,3), initial,I);

  int ret = tester.results("array_test15");
  Pooma::finalize();
  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test15.cpp,v $   $Author: richard $
// $Revision: 1.26 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
