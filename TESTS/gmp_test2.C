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
// Grid-based Multi-Patch Array's test 2
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Layout/GridLayout.h"
#include "Engine/MultiPatchEngine.h"
//#include "Evaluator/MultiPatchEval.h"
#include "Domain/Grid.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  int i;
  
  // Create the total domain.
  
  Interval<1> domain(12);
  
  // Create the block sizes.
  
  Loc<1> blocks(3), blocks2(4);

  // Create the partitioners.
  
  GridPartition<1> partition(blocks), partition2(blocks2);
  
  Range<1> r(domain[0].first(), domain[0].last()+1,
	     domain[0].size()/blocks[0].first());
  Grid<1> g(r);
  Range<1> r2(domain[0].first(), domain[0].last()+1,
	      domain[0].size()/blocks2[0].first());
  Grid<1> g2(r2);
  GridPartition<1> gp(g),gp2(g2);

  // Create the layouts.
  
  GridLayout<1> layout(domain, partition, ReplicatedTag());
  GridLayout<1> layout2(domain, partition2, ReplicatedTag());
  GridLayout<1> layout3(g, ReplicatedTag());
  GridLayout<1> layout4(g2, ReplicatedTag());
  
  tester.out() << layout << std::endl;

  tester.out() << layout2 << std::endl;

  tester.out() << layout3 << std::endl;

  tester.out() << layout4 << std::endl;

  // Make some GMP arrays and fill them.
  
  Array<1, double, Brick > a(12), ans(12);
  Array<1, double, Brick > ga(12),gans(12);
  Array<1, double, MultiPatch<GridTag,Brick> > bb(layout), cc(layout2);
  Array<1, double, MultiPatch<GridTag,Brick> > gbb(layout3), gcc(layout4);

  for (i = 0; i < 12; i++)
    {
      bb(i) = 1.0 + i;
      cc(i) = -2.3 * i;
      ans(i) = bb(i) + 3.0 * cc(i);

      gbb(i) = 1.0 + i;
      gcc(i) = -2.3 * i;
      gans(i) = gbb(i) + 3.0 * gcc(i);
    }
  
  a = bb + 3.0 * cc;

  ga = gbb + 3.0 * gcc;

  Pooma::blockAndEvaluate();

  for (i = 0; i < 12; i++)
    {
      tester.check(a(i) == ans(i));
    }
      
  int ret = tester.results("gmp_test2");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: gmp_test2.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
