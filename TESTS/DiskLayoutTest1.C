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
// DiskLayoutTest1: Test the DiskLayout functionality for a single
//   fileset, including the ability to dynamically detect the need to
//   fix byte ordering. 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

// Class to test

#include "IO/DiskLayout.h"

// Basic infrastructure

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"

// Standard includes

#include <fstream>

#if !POOMA_NO_IOS_HEADER
# include <ios>
#endif

//-----------------------------------------------------------------------------
// Test data
//-----------------------------------------------------------------------------

#include "VolFracLayoutData.h"

//-----------------------------------------------------------------------------
// Test routine
//-----------------------------------------------------------------------------

// Test that the above data is properly interpreted.

void testDiskLayout(const char *basename, Pooma::Tester &tester)
{
  typedef DiskLayout<3> DiskLayout_t;

  DiskLayout_t dl(basename);

  // Read the first layout in the file.

  bool success = dl.open();
  tester.check(success);
  
  if (success)
    {
      success = dl.read();
      tester.check(success);

      // Iterate through the "Nodes" and check that all of the values are
      // correct. 

      const DiskLayout_t::NodeList_t &nlist = dl.allNodes();

      tester.check(nlist.size() == 4);

      unsigned int i;

      for (i = 0; i < nlist.size(); ++i) 
        {
          tester.check(nlist[i].context_m == 0);
        }
  
      tester.check(nlist[0].domain_m[0] == Interval<1>(0,1));
      tester.check(nlist[0].domain_m[1] == Interval<1>(0,1));
      tester.check(nlist[0].domain_m[2] == Interval<1>(0,5));
    
      tester.check(nlist[1].domain_m[0] == Interval<1>(0,1));
      tester.check(nlist[1].domain_m[1] == Interval<1>(2,4));
      tester.check(nlist[1].domain_m[2] == Interval<1>(0,5));
    
      tester.check(nlist[2].domain_m[0] == Interval<1>(2,3));
      tester.check(nlist[2].domain_m[1] == Interval<1>(0,1));
      tester.check(nlist[2].domain_m[2] == Interval<1>(0,5));
    
      tester.check(nlist[3].domain_m[0] == Interval<1>(2,3));
      tester.check(nlist[3].domain_m[1] == Interval<1>(2,4));
      tester.check(nlist[3].domain_m[2] == Interval<1>(0,5));

      Interval<3> testdomain(Interval<1>(4),Interval<1>(5),Interval<1>(6));

      tester.out().setOutputContext(-1);
      tester.out() << "Global domain = " << dl.domain() << std::endl;
      tester.check(dl.domain() == testdomain);

      // Print out the nodelist if we're in verbose mode.

      for (i = 0; i < nlist.size(); ++i)
        {
          tester.out() << "Node " << i 
                       << ": context = " << nlist[i].context_m
                       << ", domain = " << nlist[i].domain_m << std::endl;
        }

      // Try to read another layout - this should fail gracefully

      success = dl.read();
    }

  tester.check(!success);
}

//-----------------------------------------------------------------------------
// Main program
//  -- create disk files from the binary testdata and rtestdata
//     arrays.
//  -- run the test on the files.
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Write the test file.

  if (Pooma::context() == 0)
    {
      // Only write the test file on context 0

      std::ofstream testfile;
      testfile.open("TestData.layout",std::ios::binary);
      testfile.write(&VolFrac_layout_dump[0], sizeof(VolFrac_layout_dump));
      testfile.close();

      testfile.open("RTestData.layout",std::ios::binary);
      testfile.write(&VolFrac_layout_dump_reversed[0], 
		     sizeof(VolFrac_layout_dump_reversed));
      testfile.close();
    }

  // Run the test on each file.

  tester.out() << "Testing with big-endian data..." << std::endl;

  testDiskLayout("TestData", tester);

  tester.out() << "\nTesting with little-endian data..." << std::endl;

  testDiskLayout("RTestData", tester);

  int ret = tester.results("DiskLayoutTest1");
  
  Pooma::finalize();

  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DiskLayoutTest1.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
