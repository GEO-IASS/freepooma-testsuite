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
// DiskMetaTest1: Test the DiskMeta class
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

// Class to test

#include "IO/DiskMeta.h"

// Basic infrastructure

#include "Utilities/Tester.h"
#include "Domain/Interval.h"
using namespace Pooma;

// Standard includes

#include <string>
#include <iostream>
#include <fstream>

using std::string;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
// Test data
//-----------------------------------------------------------------------------

// This array was taken from a dump of VolFrac.meta, supplied by
// Jean Marshall.

char testdata[] = "\n\
                   # This is some test data for creating a .meta file \n\
                   Type =           unknown # unknown OK\n\
                   Dim =            3\n\
                   Domain =         0 3 1\n\
                   Domain =         0 4 1 \n\
                   Domain =         0 5 1\n\
                   Fields =         2\n\
                   Records =        1\n\
                   SMPs =           1\n\
                   VnodesInRecord =  4\n\
                   VnodeTally=     0\n";


//-----------------------------------------------------------------------------
// Main program
//  -- create disk files from the testdata array.
//  -- run the test.
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::Tester tester(argc, argv);

  // Write the test file.

  std::ofstream testfile;
  testfile.open("TestData.meta");
  testfile.write(&testdata[0], sizeof(testdata));
  testfile.close();

  // Create a DiskMeta object and read the data file.

  DiskMeta metareader("TestData");

  bool success = metareader.open();
  tester.check(success);

  success = metareader.read();
  tester.check(success);

  // Check that the file was read correctly.

  tester.check(metareader.filename() == "TestData.meta");
  tester.check(metareader.type() == "unknown");
  tester.check(metareader.dimension() == 3);

  typedef Interval<1> Dom_t;
  tester.check(metareader.domain(0) == Dom_t(4));
  tester.check(metareader.domain(1) == Dom_t(5));
  tester.check(metareader.domain(2) == Dom_t(6));

  tester.check(metareader.fieldsPerRecord() == 2);
  tester.check(metareader.numRecords() == 1);
  tester.check(metareader.numFileSets() == 1);

  tester.out() << "PatchesPerRecord size = " 
               << metareader.patchesPerRecord().size() << std::endl;
  tester.out() << "PatchesPerRecord : ";
  unsigned int i;
  for (i = 0; i < metareader.patchesPerRecord().size(); ++i)
    {
      tester.out() << metareader.patchesPerRecord()[i] << " ";
    }
  tester.out() << std::endl;

  tester.check(metareader.patchesPerRecord().size() == 1);
  tester.check(metareader.patchesPerRecord()[0] == 4);

  tester.out() << "PatchTally size = " 
               << metareader.patchTally().size() << std::endl;
  tester.out() << "PatchTally : ";
  for (i = 0; i < metareader.patchTally().size(); ++i)
    {
      tester.out() << metareader.patchTally()[i] << " ";
    }
  tester.out() << std::endl;

  tester.check(metareader.patchTally().size() == 1);
  tester.check(metareader.patchTally()[0] == 0);

  int ret = tester.results("DiskMetaTest1");
  
  return ret;
}

