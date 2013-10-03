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
// BrickEngine test code
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Engine/BrickEngine.h"
#include "Engine/IndirectionEngine.h"
#include "Array/Array.h"

#include <iostream>


typedef Array<1,double,Brick> Array_t;
typedef Array<1,int,Brick> ArrayIn_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> in(1,20),sh(5);
  Array_t a(in);
  ArrayIn_t h(sh);

  int i;
  for (i=1;i<=20;++i)
  {
    a(i) = i;
  }
  for (i=0;i<5;++i)
  {
    h(i) = 2*i+3;
  }

  a(h) += 4;

  Pooma::blockAndEvaluate();

  bool worked = true;

  for (i=0;i<5;++i)
  {
    worked = worked && ( a(h(i)) == h(i) + 4 );
  }

  tester.check(worked);

  int ret = tester.results("array_test14");
  Pooma::finalize();
  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test14.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
