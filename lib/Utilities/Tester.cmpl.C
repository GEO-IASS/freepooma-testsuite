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
// Classes:
// Tester
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/Tester.h"
#include "Utilities/Inform.h"
#include "Utilities/PAssert.h"
#include <string>
#include <iostream>


//-----------------------------------------------------------------------------
// Create a default Tester object
//-----------------------------------------------------------------------------

Pooma::Tester::Tester()
  : ok_m(true), quiet_m(false), inform_m("Pooma"), 
    verbose_m(false), abort_m(false)
{
}


//-----------------------------------------------------------------------------
// Create a Tester by specifying argc,argv values.
//-----------------------------------------------------------------------------

Pooma::Tester::Tester(int argc, char **argv)
  : ok_m(true), quiet_m(false), inform_m("Pooma"), 
    verbose_m(false), abort_m(false)
{

  // Parse the command-line paramters

  parse(argc, argv);
}


//-----------------------------------------------------------------------------
// Destructor.
//-----------------------------------------------------------------------------

Pooma::Tester::~Tester()
{
}


//-----------------------------------------------------------------------------
// Print out a message to cout about the current status
//-----------------------------------------------------------------------------

int Pooma::Tester::results(const char *msg) const
{
  // Only print out a message if we're not supposed to be quiet

  if (!quiet_m)
    {
      Inform coutmsg;
      coutmsg << (ok_m ? "PASSED" : "FAILED");
      if (msg != 0)
	coutmsg << " ... " << msg;
      coutmsg << std::endl;
    }

  // Return the current value for the error code

  return returnValue();
}


//-----------------------------------------------------------------------------
// Print notification of exception.
//-----------------------------------------------------------------------------

void Pooma::Tester::exceptionHandler(const char *msg)
{
  if (verbose_m)
  {
    Inform exout("EXCEPTION");

    exout << "### Exception handled by Tester. ###\n"
	  << "### Exception message:\n";
  
    if (msg)
      exout << msg;
    else
      exout << "[none]";
  
    exout << std::endl;
  }
}

void Pooma::Tester::exceptionHandler(const Pooma::Assertion &asrt)
{
  if (verbose_m)
  {
    Inform exout("EXCEPTION");

    exout << "### POOMA Assertion Failure ###\n"
	  << "### " << asrt.what() << '\n'
	  << "### File " << asrt.file() << "; Line " << asrt.line()
	  << std::endl;
  }
}


//-----------------------------------------------------------------------------
// Parse the given command-line arguments in argc,argv, and use them to
// change the tester's behavior.  This will not modify argc or argv.
//
// Tester arguments:
//   -v       : turn on verbose output from test program
//   -p <str> : change prefix of test program messages to <str>
//   -q       : do not print out anything at all, just have program
//              return 0 or 1
//-----------------------------------------------------------------------------

void Pooma::Tester::parse(int argc, char **argv)
{
  int i = 1;
  while (i < argc)
    {
      std::string word(argv[i]);

      // Check the word for one we recognize, and set the proper variable.

      if (word == "-v")
	{
	  verbose_m = true;
	  quiet_m = false;
	}
      else if (word == "-p")
	{
	  if (i < (argc + 1))
	    out().setPrefix(argv[++i]);
	}
      else if (word == "-q")
	{
	  verbose_m = false;
	  quiet_m = true;
	}
      else if (word == "-abort")
        {
          abort_m = true;
        }

      // Move on to the next argument.

      ++i;
    }

  // Use the values from the flags to set up our state.

  // If we're not verbose, turn off the inform stream.

  if (!verbose_m || quiet_m)
    out().setOutputLevel(Inform::off);
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Tester.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
