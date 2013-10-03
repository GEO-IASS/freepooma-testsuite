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

#include <signal.h>
#include <stdlib.h>

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"

Pooma::Tester *tester;
bool handler_ok = false;

void newAbortHandler()
{
  tester->out() << "Running newly installed abort handler." << std::endl;
  handler_ok = true;
}

// This function is registered as the signal handler for SIGABRT.

void abortSignalHandler(int)
{
  // Exit with a `0' exit status to indicate that all went well.
  // This test is *expected* to abort.
  tester->check(handler_ok);
  int res = tester->results("pAbort");
  Pooma::finalize();
  exit(res);
}


int main(int argc, char *argv[])
{
  // Initialize POOMA

  Pooma::initialize(argc, argv);

  tester = new Pooma::Tester(argc, argv);

  // Shut down POOMA

  tester->out() << "Shutting down POOMA with abort()..." << std::endl;

  // Register a signal handler so that when Pooma::pAbort calls the
  // abort standard library function, this program does not exit with
  // a nonzero exit code.
  signal(SIGABRT, abortSignalHandler);

  Pooma::abortHandler(newAbortHandler);
  Pooma::pAbort("This is the abort message.", 2);

  // If we get here, the call to Pooma::pAbort did not work.
  int res = tester->results("pAbort");
  Pooma::finalize();
  return res;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: pabort.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:17:06 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
