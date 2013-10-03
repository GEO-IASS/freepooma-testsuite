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
//
// The POOMA Framework
// 
// This program was prepared by the Regents of the University of
// California at Los Alamos National Laboratory (the University) under
// Contract No.  W-7405-ENG-36 with the U.S. Department of Energy (DOE).
// The University has certain rights in the program pursuant to the
// contract and the program should not be copied or distributed outside
// your organization.  All rights in the program are reserved by the DOE
// and the University.  Neither the U.S.  Government nor the University
// makes any warranty, express or implied, or assumes any liability or
// responsibility for the use of this software
//
// Visit http://www.acl.lanl.gov/POOMA for more details
//
//-----------------------------------------------------------------------------
// Pooma general interface test program.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"

int main(int argc, char *argv[])
{
  int i;

  // Print out argc, argv before initialization.

  std::cerr << "Before initialize: argc = " << argc << ", " << std::endl;
  for (i=0; i < argc; ++i)
    std::cerr << "  argv[" << i << "] = '" << argv[i] << "'" << std::endl;

  // Initialize POOMA

  std::cerr << "Initializing POOMA ..." << std::endl;
  Pooma::initialize(argc, argv);

  std::cerr << "After initialize: argc = " << argc << ", " << std::endl;
  for (i=0; i < argc; ++i)
    std::cerr << "  argv[" << i << "] = '" << argv[i] << "'" << std::endl;

  // Print out some results of POOMA calls to the different POOMA streams

  POOMA_PRINT(Pooma::pinfo, "POOMA version = " << Pooma::version() << std::endl);

  POOMA_PRINT(Pooma::pwarn, "POOMA build date = " << Pooma::buildDate()<<std::endl);

  POOMA_PRINT(Pooma::perr, "POOMA major ver = "<< Pooma::majorVersion()<<std::endl);
  POOMA_PRINT(Pooma::perr, "POOMA minor ver = "<< Pooma::minorVersion()<<std::endl);

  // Start logging output to a file

  Pooma::logMessages("pooma.out");

  POOMA_PRINT(Pooma::pinfo, "Now logging messages to file 'pooma.out'."<<std::endl);

  POOMA_PRINT(Pooma::pwarn, "My context = " << Pooma::context() << std::endl);
  POOMA_PRINT(Pooma::perr, "Total contexts = " << Pooma::contexts() << std::endl);

  // Do some debugging statements

  POOMA_PRINT(Pooma::pinfo, "About to start printing debug messages." << std::endl);

  POOMA_DEBUG(0, "This is a level-0 debug message." << std::endl);

  POOMA_DEBUG(1, "This is a level-1 debug message." << std::endl);

  POOMA_DEBUG(3, "This is a level-3 debug message." << std::endl);

  POOMA_DEBUG(5, "This is a level-5 debug message." << std::endl);

  // Ready with basic testing.
  Pooma::finalize();

  return 0;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: pooma.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:17:06 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
