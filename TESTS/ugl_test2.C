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
// UniformGridLayout test: Create and use UniformGridLayout objects
//-----------------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Pooma/UMPArrays.h"
#include "Partition/ContextMapper.h"
#include "Partition/SpatialPartition.h"
#include "Partition/DistributedMapper.h"
#include "Utilities/Tester.h"
#include <iostream>
#include <iterator>


int main(int argc, char *argv[]) 
{

  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Testing default constructors & assignment." 
               << std::endl;
  tester.out() << "----------------------------------------" << std::endl;


  Interval<1> I1(0,999);
  Interval<2> I2(I1,I1);
  
  GuardLayers<2> gl(2);
  Loc<2> blocks(5,5);
  
  UniformGridPartition<2> partition(blocks, gl, gl);
  
  UniformGridLayout<2> tgl(I2, partition, DistributedTag());
  
  tester.out() << "Here's the original layout:\n" 
               << tgl << std::endl;

  UniformGridLayout<2> tgl2;
  
  tester.out() << "Here's an empty layout:\n" 
               << tgl2 << std::endl;

  tgl2 = tgl;
  
  tester.out() << "Here's the second layout after assignment:\n" 
               << tgl2 << std::endl;

  Interval<1> IV1a(500,900);
  Interval<1> IV1b(600,700);
  Interval<2> IV2(IV1a,IV1b);
  
  UniformGridLayoutView<2,2> tglv(tgl,IV2);
  
  tester.out() << "Here's a non-slice view of the original layout:\n"
               << tglv << std::endl;
  
  UniformGridLayoutView<2,2> tglv2(tgl2,IV2);

  tester.out() << "Here's the same view of the second layout:\n"
               << tglv2 << std::endl;

  UniformGridLayoutView<2,2> tglv0;
  tglv0 = tglv2;
  
  tester.out() << "Here's the last layout after assignment from the second view:\n"
               << tglv0 << std::endl;
  
  
  SliceRange<2,1> S1(I2,IV1a,400);
  
  UniformGridLayoutView<1,2> tgls1(tgl, S1);
  
  tester.out() << "Here's a slice view:\n"
               << tgls1 << std::endl;
  
  UniformGridLayoutView<1,2> tgls0;
  
  tgls0 = tgls1;
  
  tester.out() << "Here's the default constructed slice after assignment:\n"
               << tgls0 << std::endl;
  
  
  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("UniformGridLayout Test 2");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ugl_test2.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
