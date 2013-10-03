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
// Grid-based Multi-Patch Array's test 4
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Layout/GridLayout.h"
#include "Engine/MultiPatchEngine.h"
//#include "Evaluator/MultiPatchEval.h"
#include "Domain/Grid.h"
#include "Layout/GuardLayers.h"
#include "Layout/Node.h"
#include "Utilities/Tester.h"

#include <iostream>
#include <iterator>
#include <vector>


using std::vector;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);
  
  // Useful typedefs, constants, declarations, &c.
  
  typedef Brick                      PTag_t;
  typedef GridTag                    LTag_t;
  typedef MultiPatch<LTag_t,PTag_t>  MPTag_t;
  typedef Engine<2,int,MPTag_t>      GMPEngine_t;
  typedef Array<2,int,MPTag_t>       GMPArray_t;
  
  typedef Engine<2,int,PTag_t>       PatchEngine_t;
  typedef Array<2,int,PTag_t>        PatchArray_t;
  typedef Array<2,int,Brick>         BrickArray_t;
  
  typedef GridLayout<2>              Layout_t;  // Change if LTag_t is changed.
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
  
  GridPartition<2> partition(blocks,igcs,egcs);

  Range<1> range(D[0].first(), D[0].last()+1, D[0].size()/nblocks);
  Grid<2> grid(range,range);
  GridPartition<2> gpartition(grid,igcs,egcs);


  // Create the layout.
  
  GridLayout<2> layout(domain, partition, ReplicatedTag());

  GridLayout<2> glayout(domain, gpartition, ReplicatedTag());

  // Make a GMP engine
    
  GMPEngine_t a(layout);
  GMPEngine_t ga(glayout);
  
  int i, j;
  
  tester.out()<< " GP with num block partitioner, layout "<<endl;

  Layout_t::iterator niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter; 
    tester.out() << std::endl << std::endl;
    PatchArray_t pa(a.globalPatch(*niter));
    pa = badval;
    int res = sum((pa - badval)*(pa - badval));
    tester.check(res == 0);
    tester.out() << pa << endl << endl;
    ++niter;
  }
  tester.out()<< " GP with Grid based partitioner, layout "<<endl;

  Layout_t::iterator gniter = glayout.beginGlobal();
  while (gniter != glayout.endGlobal())
  {
    tester.out() << *gniter; 
    tester.out() << std::endl << std::endl;
    PatchArray_t gpa(ga.globalPatch(*gniter));
    gpa = badval;
    int res = sum((gpa - badval)*(gpa - badval));
    tester.check(res == 0);
    tester.out() << gpa << endl << endl;
    ++gniter;
  }

  
  // Check that engine indexing is working.

  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      a(i,j) = i + j;
      
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      {
        tester.check(a(i,j) == i + j);
        tester.check(a.read(i,j) == i + j);
      }

  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      ga(i,j) = i + j;
      
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      {
        tester.check(ga(i,j) == i + j);
        tester.check(ga.read(i,j) == i + j);
      }
  // Now do it for an array with the same engine.

  GMPArray_t aa(a);
  GMPArray_t gaa(ga);
  
  tester.out() << " Blocked GP, array aa " <<endl;
  tester.out() << aa << endl << endl;

  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      {
        tester.check(aa(i,j) == i + j);
        tester.check(aa.read(i,j) == i + j);
      }
  tester.out() << " Grid'ed GP, array gaa " <<endl;
  tester.out() << gaa << endl << endl;

  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      {
        tester.check(gaa(i,j) == i + j);
        tester.check(gaa.read(i,j) == i + j);
      }

  // Now look at the patches:
  
  tester.out()<< " Blocked GP patches " <<endl;

  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    BrickArray_t ans(niter->allocated());
    ans = badval;
    ans(niter->domain()) = aa(niter->domain());
    int res = sum((ans-pa)*(ans-pa));
    tester.check(res == 0);
    ++niter;
  }

  tester.out()<< " Grid'ed GP patches " <<endl;

  gniter = glayout.beginGlobal();
  while (gniter != glayout.endGlobal())
  {
    tester.out() << *gniter << endl << endl;
    PatchArray_t gpa(ga.globalPatch(*gniter));
    tester.out() << gpa << endl << endl;
    BrickArray_t gans(gniter->allocated());
    gans = badval;
    gans(gniter->domain()) = gaa(gniter->domain());
    int res = sum((gans-gpa)*(gans-gpa));
    tester.check(res == 0);
    ++gniter;
  }
  // Not easy to set up checks for this. Just look at them 
  // when -v is used.
  
  tester.out()<< " Block GP patch fill list " << endl<<endl;

  Layout_t::FillIterator_t fiter = layout.beginFillList();
  while (fiter != layout.endFillList())
  {
    tester.out() << "From: "  << fiter->ownedID_m 
                 << ", To: "  << fiter->guardID_m
                 << ", Dom: " << fiter->domain_m << endl << endl;
    ++fiter;
  }
  
  tester.out()<< " Grid'ed GP patch fill list " << endl<<endl;
  fiter = glayout.beginFillList();
  while (fiter != glayout.endFillList())
  {
    tester.out() << "From: "  << fiter->ownedID_m 
                 << ", To: "  << fiter->guardID_m
                 << ", Dom: " << fiter->domain_m << endl << endl;
    ++fiter;
  }
  


  // Fill the guard cells. 
  
  a.fillGuards();
  
  ga.fillGuards();

  // Look at the patches again. 
  
  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    BrickArray_t ans(niter->allocated());

    tester.out()<< " niter->allocated() " << niter->allocated()<<endl;

    ans(niter->allocated()) = aa(niter->allocated());
    int res = sum((ans-pa)*(ans-pa));
    tester.check(res == 0);
    ++niter;
  }


  gniter = glayout.beginGlobal();
  while (gniter != glayout.endGlobal())
  {
    tester.out() << *gniter << endl << endl;
    PatchArray_t gpa(ga.globalPatch(*gniter));
    tester.out() << gpa << endl << endl;
    BrickArray_t gans(gniter->allocated());

    tester.out()<< " gniter->allocated() " << gniter->allocated()<<endl;

    gans(gniter->allocated()) = gaa(gniter->allocated());
    int res = sum((gans-gpa)*(gans-gpa));
    tester.check(res == 0);
    ++gniter;
  }

  // Test the touches calculations.
  
  vector<Node_t> domains;
  Interval<2> look;
  look[0] = Interval<1>(5);
  look[1] = Interval<1>(2);
  
  tester.out() << "Owned domains intersecting " << look << endl;
  tester.out() << "========================================================" << endl;

  layout.touches(look,std::back_inserter(domains));
  
  vector<Node_t>::const_iterator ni;
  for (ni = domains.begin(); ni != domains.end(); ++ni)
    tester.out() << *ni << endl;
  
  tester.out() << "========================================================" << endl;
  tester.out() << endl;
  
  domains.clear();
  
  tester.out() << "Allocated domains intersecting " << look << endl;
  tester.out() << "========================================================" << endl;

  layout.touchesAlloc(look,std::back_inserter(domains));

  for (ni = domains.begin(); ni != domains.end(); ++ni)
    tester.out() << *ni << endl;
  
  tester.out() << "========================================================" << endl;
  tester.out() << endl;
  
  domains.clear();
  
  look[0] = Interval<1>(3,3);
  look[1] = Interval<1>(3,3);
  
  tester.out() << "Owned domains intersecting " << look << endl;
  tester.out() << "========================================================" << endl;

  layout.touches(look,std::back_inserter(domains));
  
  for (ni = domains.begin(); ni != domains.end(); ++ni)
    tester.out() << *ni << endl;
  
  tester.out() << "========================================================" << endl;
  tester.out() << endl;
  
  domains.clear();
  
  tester.out() << "Allocated domains intersecting " << look << endl;
  tester.out() << "========================================================" << endl;

  layout.touchesAlloc(look,std::back_inserter(domains));

  for (ni = domains.begin(); ni != domains.end(); ++ni)
    tester.out() << *ni << endl;
  
  tester.out() << "========================================================" << endl;
  tester.out() << endl;
  
  domains.clear();
  
  look[0] = Interval<1>(3,5);
  look[1] = Interval<1>(3,5);
  
  tester.out() << "Owned domains intersecting " << look << endl;
  tester.out() << "========================================================" << endl;

  layout.touches(look,std::back_inserter(domains));
  
  for (ni = domains.begin(); ni != domains.end(); ++ni)
    tester.out() << *ni << endl;
  
  tester.out() << "========================================================" << endl;
  tester.out() << endl;
  
  domains.clear();
  
  tester.out() << "Allocated domains intersecting " << look << endl;
  tester.out() << "========================================================" << endl;

  layout.touchesAlloc(look,std::back_inserter(domains));

  for (ni = domains.begin(); ni != domains.end(); ++ni)
    tester.out() << *ni << endl;
  
  tester.out() << "========================================================" << endl;
  tester.out() << endl;
  
  domains.clear();
  
  // One more check that things didn't get messed up.
  
  for (i = 0; i < size; ++i)
    for (j = 0; j < size; ++j)
      {
        tester.check(aa(i,j) == i + j);
        tester.check(aa.read(i,j) == i + j);
      }
   
  tester.out() << aa << endl << endl;
  
  // Finally, check if we can write into the guards directly through the array
  
  if (external_guards > 0)
    {
      for (i = -external_guards; i < size + external_guards; ++i)
        for (j = -external_guards; j < size + external_guards; ++j)
          {
            aa(i,j) = i + j;
          }
      
      tester.out() << aa << endl << endl;
  
      for (i = -external_guards; i < size + external_guards; ++i)
        for (j = -external_guards; j < size + external_guards; ++j)
          {
            tester.check(aa(i,j) == i + j);
            tester.check(aa.read(i,j) == i + j);
          }
    }  

  int ret = tester.results("gmp_test4 (Guard cell fill test)");
  Pooma::finalize();
  return ret; 
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: gmp_test4.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
