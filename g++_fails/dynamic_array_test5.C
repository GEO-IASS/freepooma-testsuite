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
// DynamicArray test 5: Create/destroy for MP Arrays with shared
// layouts and iterators instead of domains to describe the kill list.
//-----------------------------------------------------------------------------

// include files 

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/DynamicArrays.h"

#include <vector>
using std::vector;

#include <iostream>
using std::endl;

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": MP DynamicArray w/ shared layouts."
               << endl;
  tester.out() << "-------------------------------------------" << endl;

  // Create some Interval objects to create and index into Array's with

  tester.out() << "Creating Interval<1> objects ..." << endl;
  Interval<1> D1(3);
  tester.out() << "D1 = " << D1 << endl;

  // Create MultiPatch dynamic arrays that share a layout.

  tester.out() << "Creating MP DynamicArray using domain D1 ... " << endl;
  GridPartition<1> gpar(Loc<1>(3));
  LocalMapper<1> cmap(gpar);  
  DynamicLayout dynlayout(D1,gpar,cmap);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > a1(dynlayout);
  tester.check("a1 size", a1.domain().size() == D1.size());
  tester.check("a1 patches", a1.layout().sizeLocal() == 3);

  tester.out() << "Creating MP DynamicArray w/ same layout ..." << endl;
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > b1(a1.layout());
  tester.check("b1 size", b1.domain().size() == D1.size());
  tester.check("b1 patches", b1.layout().sizeLocal() == 3);

  // Test looping over layout nodes

  tester.out() << "DynamicArray< MultiPatch<DynamicTag,Dynamic> > layout:\n";
  tester.out() << a1.layout() << endl;

  // Initialize dynamic arrays with scalars.

  a1 = 30;
  b1 = 40;
  Pooma::blockAndEvaluate();
  tester.out() << "Initialized MP DynamicArray's to 30, 40:" << endl;
  tester.out() << "a1 = " << a1 << endl;
  tester.out() << "b1 = " << b1 << endl;
  tester.check("a1 initial sum", sum(a1) == (a1.domain().size() * 30));
  tester.check("b1 initial sum", sum(b1) == (b1.domain().size() * 40));

  // Create elements in the shared-layout MPE arrays

  tester.out() << "Creating 2 elements at end of a1 and b1 ..." << endl;
  a1.create(2);
  a1.sync();
  a1(3)=a1(4)=-50;
  b1(3)=b1(4)=-60;

  a1(a1.engine().domain().last()-1)=0;
  a1(a1.engine().domain().last())=0;
  
  tester.out() << "a1 = " << a1 << endl;
  tester.out() << "b1 = " << b1 << endl;
  tester.check("a1 size after create", a1.domain().size() == (D1.size() + 2));
  tester.check("b1 size after create", b1.domain().size() == (D1.size() + 2));

  // Delete an element in the shared-layout MPE arrays

  tester.out() << "Deleting 2nd element of a1 & b1 w/backfill ..."<<endl;
  b1.destroy(Interval<1>(1,1), BackFill());
  b1.sync();
  tester.out() << "a1 = " << a1 << endl;
  tester.out() << "b1 = " << b1 << endl;
  tester.check("a1 size after BackFill",a1.domain().size() == (D1.size() + 1));
  tester.check("b1 size after BackFill",b1.domain().size() == (D1.size() + 1));

  // Copy values from the beginning of a1 and b1 to their end

  tester.out() << "Copying first three elements of a1 and b1 ..." << endl;
  a1.copy(Interval<1>(3));
  a1.sync();
  tester.out() << "a1 = " << a1 << endl;
  tester.out() << "b1 = " << b1 << endl;
  tester.check("a1 size after copy", a1.domain().size() == (D1.size() + 4));
  tester.check("b1 size after copy", b1.domain().size() == (D1.size() + 4));

  // Delete elements using int* iterator:

  int killList[3] = {0, 3, 4};

  // No destroy method specified --- should use BackFill() by default.
  b1.destroy(killList,killList+3);
  b1.sync();
  tester.out() << "a1 = " << a1 << endl;
  tester.out() << "b1 = " << b1 << endl;

  // Create some larger multi-patch arrays.

  tester.out() 
    << "Creating dynamic arrays with initial domain of 50 and 10 patches."
    << endl;

  Interval<1> D2(50);
  GridPartition<1> gpar2(Loc<1>(10));
  LocalMapper<1> cmap2(gpar2);  
  DynamicLayout dynlayout2(D2,gpar2,cmap2);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > a2(dynlayout2);
  DynamicArray< int, MultiPatch<DynamicTag,Dynamic> > b2(dynlayout2);

  tester.out() << "Domain = " << D1 << endl;
  tester.out() << "Layout = " << dynlayout2 << endl;
  tester.check("a2 size", a2.domain().size() == D2.size());
  tester.check("a2 patches", a2.layout().sizeLocal() == 10);
  tester.check("b2 size", b2.domain().size() == D2.size());
  tester.check("b2 patches", b2.layout().sizeLocal() == 10);

  // Assign some values:

  for (int i = 0; i < a2.domain().size(); ++i)
    {
      a2(i) = i;
      b2(i) = -i;
    }

  PrintArray printer(2,3);

  tester.out() << "a2 = ";
  printer.print(tester.out(),a2);
  tester.out() << endl;

  tester.out() << "b2 = ";
  printer.print(tester.out(),b2);
  tester.out() << endl;

  vector<int> klist2(10);
  for (int i = 0; i < 10; ++i)
    {
      klist2[i] = i*5;
    }

  Pooma::IteratorPairDomain<int*> kdom(&klist2[0],&klist2[9]+1);
  tester.out() << "Kill domain = " << kdom << endl;

  tester.out() << "Destroying elements..." << endl;

  a2.destroy((const int*)&klist2[0],(const int*)&klist2[9]+1,ShiftUp());
  a2.sync();

  tester.out() << "a2.domain() = " << a2.domain() << endl;
  tester.out() << "a2.layout() = " << a2.layout() << endl;

  tester.out() << "a2 = ";
  printer.print(tester.out(),a2);
  tester.out() << endl;

  tester.out() << "b2 = ";
  printer.print(tester.out(),b2);
  tester.out() << endl;

  // Next we want to delete elements from a couple of individual
  // patches.

  int kplist[] = {0, 2, 3};
  Pooma::IteratorPairDomain<const int*> kpdom(&kplist[0], &kplist[0]+3);
  tester.out() << "Deleting elements from patch 3, kill domain = " 
               << kpdom << endl;

  a2.destroy(&kplist[0],&kplist[0]+3,3,BackFill());

  tester.out() << "Deleting same domain from patch 6..." << endl;

  vector<int> kfoo(&kplist[0], &kplist[0]+3);

  // No destroy method specified --- should use BackFill() by default.

  // NOTE: This may not compile!!! It depends on vector<int>::iterator
  // being assignable to "const int *", which is not guaranteed.  The
  // problem is that the patch-based destroy operation results in a
  // virtual function call to all observers of the underlying layout,
  // with the type of the event passed via an event base class. The
  // handler code must then determine the event type, which is
  // different for every possible domain type. The problem is that
  // vector<int>::const_iterator is usually a typedef for "const int
  // *" (and vector<int>::iterator is a typedef for "int *"), and thus
  // this code will work (and furthermore it is impossible to
  // enumerate them as different types, since they're not). But on
  // some systems this isn't the case. It would be nice to support
  // vector<int>::iterator even when it is a separate type, but I
  // don't know how to do this within the language. If someone else
  // wants to conditionally compile a seperate execution path for
  // this, and add the necessary configuration options, that would
  // probably be the best solution. For now, I've put #if 1...#endif
  // around the call and if it fails to compile, just change it to 
  // #if 0. [JAC]

#if 1
  b2.destroy(kfoo.begin(),kfoo.end(),6);
#endif
  
  a2.sync();

  tester.out() << "a2.domain() = " << a2.domain() << endl;
  tester.out() << "a2.layout() = " << a2.layout() << endl;

  tester.out() << "a2 = ";
  printer.print(tester.out(),a2);
  tester.out() << endl;

  tester.out() << "b2 = ";
  printer.print(tester.out(),b2);
  tester.out() << endl;

  

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "-------------------------------------------" << endl;
  int retval = tester.results("MP DynamicArray w/ shared layouts");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dynamic_array_test5.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
