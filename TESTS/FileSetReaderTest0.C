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
// FileSetReaderTest0: This doesn't actually test anything - rather it
// creates a file set to be read by FileSetReaderTest1 and 2. This is
// a binary dump of a file set created with r1 on nirvana. 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

// Pooma includes

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

// Dumps of Jean Marshall's VolFrac.{layout,offset,data} files and the
// text from VolFrac.meta.

#include "VolFracLayoutData.h"
#include "VolFracOffsetData.h"
#include "VolFracDataData.h"

char VolFrac_meta_dump[] = "\n\
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
// Utility routine
//-----------------------------------------------------------------------------

// Write the data above to files for testing

void setup()
{
  // Write the test files.

  std::ofstream testfile;

  testfile.open("TestData.layout",std::ios::binary);
  testfile.write(&VolFrac_layout_dump[0], sizeof(VolFrac_layout_dump));
  testfile.close();

  testfile.open("TestData.offset",std::ios::binary);
  testfile.write(&VolFrac_offset_dump[0], sizeof(VolFrac_offset_dump));
  testfile.close();

  testfile.open("TestData.data",std::ios::binary);
  testfile.write(&VolFrac_data_dump[0], sizeof(VolFrac_data_dump));
  testfile.close();

  testfile.open("TestData.meta"); // not binary
  testfile.write(&VolFrac_meta_dump[0], sizeof(VolFrac_meta_dump));
  testfile.close();


}

//-----------------------------------------------------------------------------
// Main program
//  -- Set up files for FileSetReaderTest1,2
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Set up the input files

  setup();

  int ret = tester.results("FileSetReaderTest0");
  
  Pooma::finalize();

  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetReaderTest0.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
