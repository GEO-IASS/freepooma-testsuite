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
// ump_test3
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/UMPArrays.h"
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
  
  UniformGridPartition<2> partition(blocks);
  
  // Create the layout.
  
  UniformGridLayout<2> layout(domain, partition, ReplicatedTag());
  
  // Make some UMP arrays and fill them.
  
  Array<2, double, MultiPatch<UniformTag,Brick> > a(layout), b(layout);
  a = 0;
  b = 0;
  b(N/2, N/2) = 1000;
  
  tester.out() << a(X,X) << std::endl;
  tester.out() << b(X,X) << std::endl;

  a(I,J) = (1.0 / 9.0) *
    (b(I+1,J+1) + b(I+1,J  ) + b(I+1,J-1) +
     b(I  ,J+1) + b(I  ,J  ) + b(I  ,J-1) +
     b(I-1,J+1) + b(I-1,J  ) + b(I-1,J-1));
     
  tester.out() << a(X,X) << std::endl;

  int ret = tester.results("ump_test3");
  Pooma::finalize();
  return ret;
}
	     
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ump_test3.cpp,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
