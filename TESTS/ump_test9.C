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
// Dirty flag test.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/UMPArrays.h"
#include "Layout/GuardLayers.h"
#include "Layout/Node.h"
#include "Utilities/Tester.h"

#include <iterator>
#include <vector>
using std::vector;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);
  
  // Useful typedefs, constants, declarations, &c.
  
  typedef Brick                      PTag_t;
  typedef BrickView                  PVTag_t;
  typedef UniformTag                 LTag_t;
  typedef MultiPatch<LTag_t,PTag_t>  MPTag_t;
  typedef Engine<2,int,MPTag_t>      UMPEngine_t;
  typedef Array<2,int,MPTag_t>       UMPArray_t;
  
  typedef Engine<2,int,PTag_t>       PatchEngine_t;
  typedef Engine<2,int,PVTag_t>      PatchViewEngine_t;
  typedef Array<2,int,PTag_t>        PatchArray_t;
  typedef Array<2,int,PVTag_t>       PatchViewArray_t;
  typedef Array<2,int,Brick>         BrickArray_t;
  
  typedef UniformGridLayout<2>       Layout_t; // Change if LTag_t is changed.
  typedef UniformGridLayoutView<2,2> ViewLayout_t;
  typedef Node<Interval<2> >         Node_t;
  
  typedef MultiPatchView<LTag_t,PTag_t,2> VTag_t;
  typedef Engine<2,int,VTag_t>            ViewEngine_t;
  typedef Array<2,int,VTag_t>             ViewArray_t;
  
  using std::endl;
  
  // Run parameters...
  
  const int size            = 9;       // Array will be size x size
  const int nblocks         = 3;       // Will have nblocks x nblocks patches
  const int internal_guards = 2;
  const int external_guards = 1;
  
  // Create the total domain.
  
  Interval<1> D(size);
  Interval<2> domain(D, D);
  
  // Create the block sizes.
  
  Loc<2> blocks(nblocks,nblocks);
  
  // OK, let's try some guard cells.
  
  GuardLayers<2> igcs(internal_guards), egcs(external_guards);

  // Create the partitioners.
  
  UniformGridPartition<2> partition(blocks,igcs,egcs);
  
  // Create the layout.
  
  UniformGridLayout<2> layout(domain, partition, ReplicatedTag());
  
  tester.out() << "\nCreating array a and assigning to it." << endl;

  UMPArray_t a(layout);

  a = 1;  // The ultimate test of whether POOMA is working 8-).

  tester.out() << "a's dirty flag is " << a.engine().isDirty() << endl;
  tester.check(a.engine().isDirty() == true);
  
  UMPArray_t b = a;

  tester.out() << "b's dirty flag is " << b.engine().isDirty() << endl;
  tester.check(b.engine().isDirty() == true);

  a.engine().fillGuards();

  tester.out() << "\nFilled a's guards." << endl;

  tester.out() << "a's dirty flag is " << a.engine().isDirty() << endl;
  tester.check(a.engine().isDirty() == false);
  
  tester.out() << "b's dirty flag is " << b.engine().isDirty() << endl;
  tester.check(b.engine().isDirty() == false);  
  
  // Create the view domain.
  
  CTAssert(size-2 > 2);
  Interval<1> DV(2,size-2);
  Interval<2> viewDomain(DV, DV);

  tester.out() << "\nCreating a view of a and assigning to it..." << endl;

  ViewArray_t av(a,viewDomain);
  
  av = 2;

  tester.out() << "a's dirty flag is " << a.engine().isDirty() << endl;
  tester.out() << "b's dirty flag is " << b.engine().isDirty() << endl;
  tester.check(a.engine().isDirty() == true);
  tester.check(b.engine().isDirty() == true);

  tester.out() << "av's dirty flag is " << av.engine().isDirty() << endl;
  tester.check(av.engine().isDirty() == true);

  tester.out() << "\nFilling av's guards..." << endl;

  av.engine().fillGuards();

  tester.out() << "a's dirty flag is " << a.engine().isDirty() << endl;
  tester.out() << "b's dirty flag is " << b.engine().isDirty() << endl;
  tester.check(a.engine().isDirty() == false);
  tester.check(b.engine().isDirty() == false);

  tester.out() << "av's dirty flag is " << av.engine().isDirty() << endl;
  tester.check(av.engine().isDirty() == false);

  int retval = tester.results("ump_test9: dirty flag test.");

  Pooma::finalize();    
  return retval;
}
    
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ump_test9.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
