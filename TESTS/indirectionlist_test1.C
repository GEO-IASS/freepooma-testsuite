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

#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include "Domain/Region.h"
#include "Domain/Touches.h"
#include "Domain/Contains.h"
#include "Domain/Split.h"
#include "Domain/Intersect.h"
#include "Engine/BrickEngine.h"
#include "Engine/DynamicEngine.h"
#include "Array/Array.h"
#include "DynamicArray/DynamicArray.h"
#include "Domain/IndirectionList.h"
#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"

#include <iostream>

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  tester.out() << "Starting IndirectionList test." << std::endl << std::endl;

  Interval<1> foo(0,6);

  Array<1,int,Brick> klist(foo);

  klist = 1;
  Pooma::blockAndEvaluate();

  for(int i=1;i<7;i++)
    klist(i) = klist(i-1)+i;
  klist(2)=3;
  klist(5)=12;
  klist(6)=20;
  
  tester.out() << klist << std::endl;

  Interval<1> fff(0,20);

  DynamicArray<double,Dynamic> goo(fff),roo(fff);
  for(int i=0;i<=goo.domain().last();i++)
    goo(i)=roo(i)=i;

  IndirectionList<int> iklist(klist);

  tester.out()  << " iklist.first() = " << iklist.first() << std::endl;
  tester.out()  << " iklist.last() = " << iklist.last() << std::endl;
  tester.out()  << " iklist.size() = " << iklist.size() << std::endl;
  
  tester.out() << "DynamicArray to be altered" << goo << std::endl;

  Range<1> sss(0,3);
	
  goo.destroy(iklist,ShiftUp());

  tester.out() << "after destroy with ShiftUp" << std::endl;
  tester.out() << goo << std::endl;


  roo.destroy(iklist,BackFill());
  tester.out() << "after destroy with BackFill" << std::endl;
  tester.out() << roo << std::endl;

  tester.out() << "Finished IndirectionList test." << std::endl << std::endl;

  int res = tester.results("indirectionlist_test1");
  Pooma::finalize();
  return res;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: indirectionlist_test1.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
