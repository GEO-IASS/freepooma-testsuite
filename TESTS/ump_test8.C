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
// ump_test8
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
  const int v[3]            = { 3, 8, 1 };
  const int nblocks         = 3;       // Will have nblocks x nblocks patches
  const int internal_guards = 2;
  const int external_guards = 1;
  const int badval          = -77777;
  
  // Create the total domain.
  
  Interval<1> D(size);
  Interval<2> domain(D, D);
  
  Interval<1> VD(v[0], v[1], v[2]);
  Interval<2> vdom(VD, VD);
  
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
  
  int i, j;
  
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
  
  // Now do it for an view of an array with the same engine.

  UMPArray_t aa(a);
//  ViewArray_t av(aa,vdom);
//  ViewArray_t av = aa(vdom);
  View1<UMPArray_t,Interval<2> >::Type_t av = aa(vdom);
  
  tester.out() << av << endl << endl;
  
  int iv, jv;
  for (i = vdom[0].first(), iv = 0; i <= vdom[0].last(); ++i, ++iv)
    for (j = vdom[1].first(), jv = 0; j <= vdom[1].last(); ++j, ++jv)
      {
        tester.check(av(iv,jv) == i + j);
        tester.check(av.read(iv,jv) == i + j);
      }
   
  // Now look at the patches:
  
  ViewLayout_t vlayout = av.engine().layout();
  ViewEngine_t vengine = av.engine();
  
  niter = vlayout.beginGlobal();
  while (niter != vlayout.endGlobal())
  {
    tester.out() << *niter << endl << endl;
    PatchViewArray_t pa(vengine.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    BrickArray_t ans(niter->domain());
    ans = badval;
    ans(niter->domain()) = av(niter->domain());
    int res = sum((ans() - pa) * (ans() - pa));
    tester.check(res == 0);
    ++niter;
  }
  
  // Fill the guard cells. 
  
  vengine.fillGuards();
  
  // Look at the patches again. 
  
  niter = vlayout.beginGlobal();
  while (niter != vlayout.endGlobal())
  {
    tester.out() << *niter << endl << endl;
    PatchViewArray_t pa(vengine.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    BrickArray_t ans(niter->domain());
    Range<2> bd = Pooma::NoInit();
    ans(niter->domain()) = aa(vlayout.localToBase(niter->domain(),bd));
    int res = sum((ans() - pa) * (ans() - pa));
    tester.check(res == 0);
    ++niter;
  }
  
  // Look at the patches from the viewed engine.
  
  niter = layout.beginGlobal();
  while (niter != layout.endGlobal())
  {
    tester.out() << *niter << endl << endl;
    PatchArray_t pa(a.globalPatch(*niter));
    tester.out() << pa << endl << endl;
    BrickArray_t ans(niter->allocated());
    ans(niter->allocated()) = aa(niter->allocated());
    int res = sum((ans() - pa()) * (ans() - pa()));
    tester.check(res == 0);
    ++niter;
  }
  
  
  int retval = tester.results("ump_test8");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ump_test8.cpp,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
