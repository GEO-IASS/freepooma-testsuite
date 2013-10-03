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
// Grid-based Multi-Patch Array's test 3
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Arrays.h"
#include "Domain/Grid.h"
#include "Layout/GridLayout.h"
#include "Engine/MultiPatchEngine.h"
//#include "Evaluator/MultiPatchEval.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  // Create the total domain.
  
  const int N = 40;
  Interval<1> D(1, N);
  Interval<2> domain(D, D);
  Interval<1> I(2, N-1), J(2, N-1), X(N/2 - 4, N/2 + 4);
  
  // Create the block sizes.
  
  Loc<2> blocks(2,2);

  // Create the partitioners.
  
  GridPartition<2> partition(blocks);
  
  Range<1> range(D.first(), D.last()+1, D.size()/2);
  Grid<2> grid(range,range);
  GridPartition<2> gpartition(grid);

  // Create the layout.
  
  GridLayout<2> layout(domain, partition, ReplicatedTag());
  GridLayout<2> glayout(domain,gpartition, ReplicatedTag());
  
  tester.out() << layout << std::endl;
  tester.out() << glayout << std::endl;

  // Make some GMP arrays and fill them.
  
  Array<2, double, MultiPatch<GridTag,Brick> > a(layout), b(layout);
  Array<2, double, MultiPatch<GridTag,Brick> > ga(layout), gb(layout);

  a = 0;
  b = 0;
  ga = 0;
  gb = 0;

  b(N/2, N/2) = 1000;
  gb(N/2, N/2) = 1000;

  tester.out() << a(X,X) << std::endl;
  tester.out() << b(X,X) << std::endl;

  tester.out() << ga(X,X) << std::endl;
  tester.out() << gb(X,X) << std::endl;

  a(I,J) = (1.0 / 9.0) *
    (b(I+1,J+1) + b(I+1,J  ) + b(I+1,J-1) +
     b(I  ,J+1) + b(I  ,J  ) + b(I  ,J-1) +
     b(I-1,J+1) + b(I-1,J  ) + b(I-1,J-1));
     

  ga(I,J) = (1.0 / 9.0) *
    (gb(I+1,J+1) + gb(I+1,J  ) + gb(I+1,J-1) +
     gb(I  ,J+1) + gb(I  ,J  ) + gb(I  ,J-1) +
     gb(I-1,J+1) + gb(I-1,J  ) + gb(I-1,J-1));

  tester.out() << a(X,X) << std::endl;
  tester.out() << ga(X,X) << std::endl;

  int ret = tester.results("gmp_test3");
  Pooma::finalize();
  return ret; 
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: gmp_test3.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
