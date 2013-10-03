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
// ump_test2
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/UMPArrays.h"
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
  
  UniformGridPartition<1> partition(blocks), partition2(blocks2);
  
  // Create the layouts.
  
  UniformGridLayout<1> layout(domain, partition, ReplicatedTag());
  UniformGridLayout<1> layout2(domain, partition2, ReplicatedTag());
  
  // Make some UMP arrays and fill them.
  
  Array<1, double, Brick > a(12), ans(12);
  Array<1, double, MultiPatch<UniformTag,Brick> > bb(layout), cc(layout2);
  for (i = 0; i < 12; i++)
    {
      bb(i) = 1.0 + i;
      cc(i) = -2.3 * i;
      ans(i) = bb(i) + 3.0 * cc(i);
    }
  
  a = bb + 3.0 * cc;

  Pooma::blockAndEvaluate();

  for (i = 0; i < 12; i++)
    {
      tester.check(a(i) == ans(i));
    }
  
  int ret = tester.results("ump_test2");
  Pooma::finalize();
  return ret;
}
	     
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ump_test2.cpp,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
