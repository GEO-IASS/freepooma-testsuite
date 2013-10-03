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
// ump_test6
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
  typedef UniformTag                 LTag_t;
  typedef MultiPatch<LTag_t,PTag_t>  MPTag_t;
  typedef Engine<2,int,MPTag_t>      UMPEngine_t;
  typedef Array<2,int,MPTag_t>       UMPArray_t;
  
  typedef Engine<2,int,PTag_t>       PatchEngine_t;
  typedef Array<2,int,PTag_t>        PatchArray_t;
  typedef Array<2,int,Brick>         BrickArray_t;
  
  typedef UniformGridLayout<2>       Layout_t;  // Change if LTag_t is changed.
  typedef Node<Interval<2> >         Node_t;
  
  using std::endl;
  
  // Run parameters...
  
  const int size            = 9;       // Array will be size x size
  const int nblocks         = 3;       // Will have nblocks x nblocks patches
  const int internal_guards = 2;
  const int external_guards = 1;
  const int badval          = -77777;
  
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
  
  // Make a UMP engine
    
  UMPEngine_t a(layout);
  UMPArray_t aa(a);

  aa = badval;
    
  tester.out() << aa << endl << endl;

  Layout_t::iterator niter;
  int i, j;
  
  // Print out the patches.

  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter; 
    tester.out() << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));

    // The guards are uninitialized, so this will cause UMRs.

    tester.out() << pa << endl << endl;
    ++niter;
  }
  
  // Zero the guards and print the patches again.
  
  a.setGuards(0);
  
  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter; 
    tester.out() << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    ++niter;
  }

  // Check that engine indexing is working.

  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      aa(i,j) = i + j;
      
  tester.out() << aa << endl << endl;
    
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      {
        tester.check(aa(i,j) == i + j);
        tester.check(aa.read(i,j) == i + j);
      }
  
  // Set the guards to badval and check the patches again.
  
  a.setGuards(badval);
  
  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter; 
    tester.out() << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    ++niter;
  }

  // Finally, check guard cell accumulation.
  
  aa = 0;
  
  a.setGuards(1);
  
  tester.out() << aa << endl << endl;
  
  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter; 
    tester.out() << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    ++niter;
  }
  
  // Now accumulate from the guards and see what we get. 
  
  a.accumulateFromGuards();
  
  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter; 
    tester.out() << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    ++niter;
  }

  // The total result should simply be the number of guard cells
  // overlapping any particular position.

  tester.out() << aa << endl << endl;

  int retval = tester.results("ump_test6: guard cell test.");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ump_test6.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
