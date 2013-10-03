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
// Options
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/Options.h"
#include "Utilities/Inform.h"
#include "Pooma/Configuration.h"
#include <string>
#include <iostream>
#include <stdlib.h>

namespace Pooma {

//-----------------------------------------------------------------------------
// Create an Options object with default settings.
//-----------------------------------------------------------------------------

Options::Options()
{
  // Just reset the options to their default values.

  reset();
}


//-----------------------------------------------------------------------------
// Create an Options object with a set of argc,argv values.  This
// will parse the arguments and use them to override the current
// settings.  POOMA-specific arguments are stripped out.
//-----------------------------------------------------------------------------

Options::Options(int &argc, char ** argv)
{
  // First, reset all the options to their default values.

  reset();

  // Then parse the command-line arguments and use them to set values.
  // Strip out pooma-specific args.

  parse(argc, argv);
}


//-----------------------------------------------------------------------------
// Copy constructor.
//-----------------------------------------------------------------------------

Options::Options(const Options &opts)
{
  // Just use operator= for this.

  *this = opts;
}

Options &Options::operator=(const Options &opts)
{
  concurrency_m = opts.concurrency();

  info_m = opts.printInfo();
  warn_m = opts.printWarnings();
  err_m  = opts.printErrors();
  
  logfile_m = opts.logfile();

  stats_m = opts.printStats();
  debug_m = opts.debug();
  
  neverCompress_m = opts.neverCompress();
  deferredFills_m = opts.deferredGuardFills();
  
  hardinit_m            = opts.hardInit();
  hardrun_m             = opts.hardRun();
  lockthreads_m         = opts.lockThreads();
  blockingExpressions_m = opts.blockingExpressions();

  return *this;
}


//-----------------------------------------------------------------------------
// Destructor.
//-----------------------------------------------------------------------------

Options::~Options()
{
}


//-----------------------------------------------------------------------------
// Print out a usage summary for POOMA arguments to cerr.
//-----------------------------------------------------------------------------

void Options::usage()
{
  Inform msg("POOMA Usage", std::cerr);
  msg << ">>>-----------------------------------<<<\n";
  msg << ">>> POOMA command-line option summary <<<\n";
  msg << ">>>-----------------------------------<<<\n";
  msg << "Standard options:\n";
  msg << "--pooma-threads <N> ......... set concurrency level (N >= 1)\n";
  msg << "--pooma-info ................ turn on output of info messages\n";
  msg << "--pooma-warn ................ turn on output of warning messages\n";
  msg << "--pooma-err ................. turn on output of error messages\n";
  msg << "--pooma-log <file> .......... turn on logging of output to <file>\n";
  msg << "--pooma-stats ............... turn on output of stats at end\n";
  msg << "--pooma-nocompress .......... disable compression of\n";
  msg << "                              compressible brick-engines\n";
// Not used yet, so don't mention in the usage message:
//  msg << "--pooma-nodeferred-guardfills disable deferred guard fills\n";
  msg << "--pooma-help ................ print out this summary\n";
  msg << "Developer options:\n";
  msg << "--pooma-debug <N> ........... set debug output level to <N>\n";
  msg << "--pooma-smarts-hardinit\n";
  msg << "--pooma-smarts-hardrun\n";
  msg << "--pooma-smarts-lockthreads\n";
  msg << "--pooma-blocking-expressions\n";
  msg << "All options exist as \"yes\" and \"no\" pairs.\n";
  msg << "For example --pooma-info and --pooma-noinfo.\n";
  msg << "The \"no\" versions listed above imply that \"yes\" is the default.";
  msg << std::endl;
}


//-----------------------------------------------------------------------------
// Reset all the options to their default values.
//-----------------------------------------------------------------------------

void Options::reset()
{
  concurrency_m = 1;

  info_m = true;
  warn_m = true;
  err_m  = true;
  
  logfile_m = "";
  stats_m   = false;
  debug_m   = Inform::off;
  
  neverCompress_m = false;
  deferredFills_m = true;
  
  // These defaults can be changed with configure.

  hardinit_m            = POOMA_DEFAULT_SMARTS_HARDINIT;
  hardrun_m             = POOMA_DEFAULT_SMARTS_HARDRUN;
  lockthreads_m         = POOMA_DEFAULT_SMARTS_LOCKTHREADS;
  blockingExpressions_m = POOMA_DEFAULT_BLOCKING_EXPRESSIONS;
}


//-----------------------------------------------------------------------------
// Parse the given command-line arguments in argc,argv, and use them to
// change the options values.  If an error occurs, print a message and
// abort.  POOMA-specific arguments are stripped out.
//-----------------------------------------------------------------------------

void Options::parse(int &argc, char ** &argv)
{
  // If there are no arguments beyond the first, just return.

  if (argc < 2)
    return;

  // Create items to store the return number of arguments and their pointers.

  int retargc = 1;
  char **retargv = new char *[argc];
  retargv[0] = argv[0];

  // Now, scan through the arguments, changing values as necessary.
  // Strip out pooma-specific arguments.

  int i = 1;
  while (i < argc)
    {
      bool argok = true;
      bool argvalerr = false;
      std::string word(argv[i]);

      // Check the word for one we recognize, and set the proper variable.

      if (word == "--pooma-threads")
	{
	  argok = intArgument(argc, argv, i+1, concurrency_m);
	  argvalerr = (concurrency_m < 1);
	  ++i;
	}
      else if (word == "--pooma-nothreads")
	{
	  concurrency_m = 1;
	}
      else if (word == "--pooma-info" || word == "--pooma-noinfo")
	{
	  info_m = (word == "--pooma-info");
	}
      else if (word == "--pooma-nocompress" || word == "--pooma-compress")
	{
	  neverCompress_m = (word == "--pooma-nocompress");
	}
      else if (word == "--pooma-nodeferred-guardfills" || 
               word == "--pooma-deferred-guardfills")
	{
	  deferredFills_m = (word == "--pooma-deferred-guardfills");
	}
      else if (word == "--pooma-warn" || word == "--pooma-nowarn")
	{
	  warn_m = (word == "--pooma-warn");
	}
      else if (word == "--pooma-err" || word == "--pooma-noerr")
	{
	  err_m = (word == "--pooma-err");
	}
      else if (word == "--pooma-stats" || word == "--pooma-nostats")
	{
	  stats_m = (word == "--pooma-stats");
	}
      else if (word == "--pooma-log")
	{
	  argok = stringArgument(argc, argv, i+1, logfile_m);
	  ++i;
	}
      else if (word == "--pooma-debug")
	{
	  argok = intArgument(argc, argv, i+1, debug_m);
	  ++i;
	}
      else if (word == "--pooma-nodebug")
	{
	  debug_m = Inform::off;
	}
      else if (word == "--pooma-smarts-hardinit")
	{
	  hardinit_m = true;
	}
      else if (word == "--pooma-smarts-nohardinit")
	{
	  hardinit_m = false;
	}
      else if (word == "--pooma-smarts-hardrun")
	{
	  hardrun_m = true;
	}
      else if (word == "--pooma-smarts-nohardrun")
	{
	  hardrun_m = false;
	}
      else if (word == "--pooma-smarts-lockthreads")
	{
	  lockthreads_m = true;
	}
      else if (word == "--pooma-smarts-nolockthreads")
	{
	  lockthreads_m = false;
	}
      else if (word == "--pooma-blocking-expressions")
	{
	  blockingExpressions_m = true;
	}
      else if (word == "--pooma-noblocking-expressions")
	{
	  blockingExpressions_m = false;
	}
      else if (word == "--pooma-help")
	{
	  usage();
	  exit(0);
	}
      else
	{
	  // This is an unrecognized option, so we should return this
	  // back to the caller.

	  retargv[retargc++] = argv[i];
	}

      // Check if there was an error ...
      if (!argok)
	{
	  // There was, print a message, a usage summary, and exit.

	  std::cerr << "\nERROR: Bad format for POOMA command-line option '";
	  std::cerr << word.c_str() << "'.\n" << std::endl;
	  usage();
	  exit(0);
	}

      // Check if the argument had an illegal value ...
      if (argvalerr)
	{
	  // It did, print a message, a usage summary, and exit.

	  std::cerr << "\n";
	  std::cerr << "ERROR: Illegal value for POOMA command-line option '";
	  std::cerr << word.c_str() << "'.\n" << std::endl;
	  usage();
	  exit(0);
	}

      // Otherwise, move on to the next argument.

      ++i;
    }

  // If we're here, the parsing was successful, so return back the changed
  // argc,argv values.

  argc = retargc;

  for (i = 0; i < retargc; ++i) argv[i] = retargv[i];

  delete [] retargv;
}


//-----------------------------------------------------------------------------
// A function to check for an integer argument in the pos position,
// and set the last argument equal to it.  If it does not exist, or the
// next argument is not an integer, return false and do not change the
// last argument's value
//-----------------------------------------------------------------------------

bool intArgument(int argc, char **argv, int pos, int &val)
{
  // Make sure there is an argument available

  if (pos >= argc)
    return false;

  // Make sure the 'pos' argument is a number.  If it starts with a number
  // or -number, it is OK.

  char firstchar = argv[pos][0];
  if (firstchar < '0' || firstchar > '9')
    {
      // first char is not a number.  Is the second, with the first a '-/+'?

      if ((firstchar != '-' && firstchar != '+') || argv[pos][1] == 0 ||
	  (argv[pos][1] < '0' || argv[pos][1] > '9'))
	return false;
    }

  // Get the value and return it in the last argument

  val = atoi(argv[pos]);
  return true;
}


//-----------------------------------------------------------------------------
// A function to check for a string argument in the pos position, and
// set the last argument equal to it.  If it does not exist, return false
// and leave the last argument unchanged.
//-----------------------------------------------------------------------------

bool stringArgument(int argc, char **argv, int pos, std::string &val)
{
  // Make sure there is an argument available

  if (pos >= argc)
    return false;

  // Get the value and return it in the last argument

  val = argv[pos];
  return true;
}


//-----------------------------------------------------------------------------
// A function to check for a floating-point argument in the pos position.
//-----------------------------------------------------------------------------

bool doubleArgument(int argc, char **argv, int pos, double &val)
{
  // Make sure there is an argument available

  if (pos >= argc)
    return false;

  // Get the value and return it in the last argument

  val = atof(argv[pos]);
  return true;
}

} // namespace Pooma

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Options.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
