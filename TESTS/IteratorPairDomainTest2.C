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

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/IteratorPairDomain.h"

#include "Array/PrintArray.h"
#include "DynamicArray/DynamicArray.h"
#include "Engine/DynamicEngine.h"

#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  using std::vector;

  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  PrintArray printer(2,2);

  tester.out() << "Starting IteratorPairDomain test.\n" << endl;

  vector<int> klist(7); // klist(0) ... klist(6)

  klist.assign(klist.size(),1);

  for (int i = 1; i < 7; ++i)
    klist[i] = klist[i-1] + i;
  klist[2] = 3;
  klist[5] = 12;
  klist[6] = 20;
  
  typedef vector<int>::iterator Iter_t;
  typedef Pooma::IteratorPairDomain<Iter_t> IPDomain_t;
  
  IPDomain_t dom(klist.begin(),klist.end());

  Interval<1> fff(0,20);

  DynamicArray<double,Dynamic> goo(fff), roo(fff);

  for(int i = 0; i < goo.domain().size(); ++i)
    goo(i) = roo(i) = i;

  tester.out() << "DynamicArray to be altered  : ";
  printer.print(tester.out(),goo);

  tester.out() << "Elements to be deleted      : " << dom << endl;
	
  goo.destroy(dom,ShiftUp());

  tester.out() << "After destroy with ShiftUp  : ";
  printer.print(tester.out(),goo);

  roo.destroy(dom,BackFill());

  tester.out() << "After destroy with BackFill : ";
  printer.print(tester.out(),roo);

  tester.out() << "Finished IteratorPairDomain test 2.\n" << endl;

  int res = tester.results("IteratorPairDomainTest2");
  Pooma::finalize();
  return res;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IteratorPairDomainTest2.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
