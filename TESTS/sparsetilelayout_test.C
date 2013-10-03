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
// SparseTileLayout test: Create and use SparseTileLayout objects
//-----------------------------------------------------------------------------
#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Pooma/GMPArrays.h"
#include "Layout/SparseTileLayout.h"
#include "Partition/TilePartition.h"
#include "Utilities/Tester.h"
#include "Pooma/Arrays.h"
#include <iostream>
#include <iterator>

typedef Brick                                       PatchEngineTag_t;
typedef  MultiPatch<SparseTileTag,PatchEngineTag_t> mp_t;
typedef  Array<2,double,mp_t>                       Array_t;


int main(int argc, char *argv[]) 
{
  Pooma::initialize(argc, argv);


  // Initialize POOMA and output stream, using Tester class


  Pooma::Tester tester(argc, argv); 
  tester.out() << argv[0] << ": SparseTileLayout operations." << std::endl;
  tester.out() << "----------------------------------------" << std::endl;

  typedef SparseTileLayout<2>::Domain_t Domain_t;
  Domain_t f(Interval<1>(0,9),Interval<1>(0,9));

  typedef SparseTileLayout<2>::PatchList_t  PatchList_t;
  PatchList_t plist;

  plist.push_back(Interval<2>(Interval<1>(0,4),Interval<1>(0,4)));
  plist.push_back(Interval<2>(Interval<1>(5,9),Interval<1>(0,4)));
  plist.push_back(Interval<2>(Interval<1>(0,4),Interval<1>(5,9)));
  plist.push_back(Interval<2>(Interval<1>(5,9),Interval<1>(5,9)));

  GuardLayers<2> igl(2),egl(2);

  TilePartition<2> tp(plist,igl,egl);

  SparseTileLayout<2> stl_pl(f,plist,ReplicatedTag());

  stl_pl.print(tester.out());

  stl_pl.syncPatch();

  SparseTileLayout<2> pp(f,tp,ReplicatedTag());

  pp.syncPatch();

  tester.out()<< std::endl;
  tester.out()<< " printing out the sparse tile layout " <<std::endl;
  tester.out()<< "   4 equal size patches tile the domain " <<std::endl;
  tester.out()<< "   this is equivalent to a "<<std::endl;
  tester.out()<< "    UGL(domain,Loc<2>(2),GuardLayers<2>(2))"<<std::endl;

  pp.print(tester.out());

  //*************
  
  PatchList_t pplist;
  
  pplist.push_back(Interval<2>(Interval<1>(0,4),Interval<1>(3,9)));
  pplist.push_back(Interval<2>(Interval<1>(5,9),Interval<1>(0,7)));

  TilePartition<2> tp2(pplist,igl,egl);
  SparseTileLayout<2> twopatch(f,tp2,ReplicatedTag());

  Array_t STa(twopatch);

  STa = 1.1;


  // this bit assigns border regions. 

  SparseTileLayout<2>::BorderFillIterator_t
    tp_b_start = twopatch.beginBorderFillList();

  SparseTileLayout<2>::BorderFillIterator_t
    tp_b_end = twopatch.endBorderFillList();
  
  tester.out() << " testing assigning into border regions " << std::endl;

  tester.out() << " layout is: " << std::endl;
  tester.out() << twopatch << std::endl;


  for ( ; tp_b_start != tp_b_end ; ++tp_b_start )
    {
      tester.out() << " domain "<< tp_b_start->domain() 
		   << " Patch id "<< tp_b_start->patchID()<<std::endl;

      tester.out() << " domain of the patch is "<<
	STa.patchLocal( tp_b_start->patchID()).domain()<<std::endl;

      STa.patchLocal( tp_b_start->patchID() ) ( tp_b_start->domain() ) = 2.2;

      tester.out() << " patch "<<tp_b_start->patchID() <<std::endl;
      tester.out() <<  STa.patchLocal( tp_b_start->patchID()) <<std::endl;
      
      tester.out() << "view of the same patch "<<tp_b_start->patchID() <<std::endl;
      tester.out() <<  STa.patchLocal( tp_b_start->patchID())() <<std::endl;
    }


  // SparseTL based arrays work in expressions...
 
  Array_t STb(twopatch);
  Array_t STc(twopatch);
 

  STa = 1.0;        
  STb = 9.9;

  STc = STa + (STb - STa);

  tester.out()<< " print out the sparse tile layout based array "<<std::endl;

  tester.out()<<STc <<std::endl;

  tester.out()<< " print an expression using a STlayout "<<std::endl;

  tester.out()<<  17 + (STc+4)*3 <<std::endl;

 
 
  // Now assign to views that include undefined areas.


  Interval<2> view(Interval<1>(2,6),Interval<1>(0,9));
  STc(view) = .5*STa(view);

  tester.out()<<STc <<std::endl;

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("SparseTileLayout operations");
  Pooma::finalize();
  return retval;
} 

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: sparsetilelayout_test.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
