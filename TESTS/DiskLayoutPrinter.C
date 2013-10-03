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
// DiskLayoutPrinter: Test/Utility that prints out a textual
// representation of the data in a "DiscField" .layout file.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

// DiskLayout class

#include "IO/DiskLayout.h"

// Basic infrastructure

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"

// Standard includes

#include <iostream> 

//-----------------------------------------------------------------------------
// Utility routine
//-----------------------------------------------------------------------------

// List the contents of the file.

void listDiskLayout(const char *basename)
{
  typedef DiskLayout<3> DiskLayout_t;

  DiskLayout_t dl(basename);

  bool success = dl.open();

  if (!success)
    {
      std::cout << "Could not open file!" << std::endl;
    }
  else
    {
      std::cout << "Reading layout from " << basename << ".layout\n";
      if (dl.bytesReversed()) std::cout << "Layout has bytes reversed\n";

      // Read the first layout in the file.

      int rec = 1;
      while (dl.read())
        {
          if (rec == 1) 
            std::cout << "Global domain = " << dl.domain() << std::endl;

          std::cout << 
            "---------------------------------------------------------"
            "--------\n";
          std::cout << "Record " << rec++ << std::endl;

          // Iterate through the "Nodes" and print out the data.

          const DiskLayout_t::NodeList_t &nlist = dl.allNodes();

          std::cout << "Number of nodes in this layout: " << nlist.size() 
                    << std::endl;

          unsigned int i;

          if (Pooma::context() == 0)
            {
              // Print out the nodelist if we're in verbose mode.

              for (i = 0; i < nlist.size(); ++i)
                {
                  std::cout << "Node " << i 
                            << ": context = " << nlist[i].context_m
                            << ", domain = " << nlist[i].domain_m << std::endl;
                }
            }
        }
      std::cout << 
        "-----------------------------------------------------------------"
                << std::endl;
    }
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

  if (argc == 2 && *argv[1] != '-') 
    {
      listDiskLayout(argv[1]);
    }
  else // If no argument was specified, just print "PASSED" message
    {  // (If -v was passed, print a usage message as well.)

      Pooma::Tester tester(argc,argv);
      tester.out() << "Usage: DiskLayoutPrinter basename" << std::endl;
      int ret = tester.results("DiskLayoutPrinter");
      return ret;
    }

  Pooma::finalize();

  return 0;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DiskLayoutPrinter.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
