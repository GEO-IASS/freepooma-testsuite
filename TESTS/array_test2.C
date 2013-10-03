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
// Array test 2: indexing.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Array/Array.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  
  Array<3, int> a(6, 6, 6);
  
  // Block since we're starting scalar code.
    
  Pooma::blockAndEvaluate();
  
  for (int i2 = 0; i2 < 6; i2++)
    for (int i1 = 0; i1 < 6; i1++)
      for (int i0 = 0; i0 < 6; i0++)
	a(i0,i1,i2) = i2+10*(i1+10*i0);

  tester.check("a(int,int,int)", a(1, 2, 3), 123);
  tester.check("a.read(int,int,int)", a.read(1, 2, 3), 123);
  tester.check("a(int,long,int)", a(1, 2L, 3), 123);
  tester.check("a.read(int,long,int)", a.read(1, 2L, 3), 123);
  tester.check("a(int,int,unsigned)", a(1, 2, 3U), 123);
  tester.check("a.read(int,int,unsigned)", a.read(1, 2, 3U), 123);
  tester.check("a(unsigned long,int,unsigned)", a(1UL, 2, 3U), 123);
  tester.check("a.read(unsigned long,int,unsigned)", a.read(1UL, 2, 3U), 123);

  int ret = tester.results( "array_test2" );
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test2.cpp,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
