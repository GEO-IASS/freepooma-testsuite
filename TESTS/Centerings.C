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
// Test of the new Centerings class.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Field/FieldCentering.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Explicitly obtain the centerings.
  CanonicalCentering<2> centering2;
  CanonicalCentering<3> centering3;

  Centering<2> cell = centering2(CellType, Continuous);

  Centering<3> allFace = centering3(FaceType, Continuous, XDim | YDim);

  tester.out() << cell << std::endl;
  tester.out() << allFace << std::endl;

  // Use the functional interface to obtain the centerings.
  tester.out() << canonicalCentering<2>(CellType, Continuous) << std::endl;
  tester.out() << canonicalCentering<3>(FaceType, Discontinuous, XDim | YDim) << std::endl;

  // Briefly test the comparison operators.
  tester.check(cell == canonicalCentering<2>(CellType, Continuous));
  tester.check(cell != centering2(FaceType, Continuous, XDim | YDim));
  tester.check(allFace == centering3(FaceType, Continuous, XDim | YDim));
  tester.check(centering3(FaceType, Continuous, XDim | YDim) == allFace);

  int ret = tester.results("Centerings");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Centerings.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
