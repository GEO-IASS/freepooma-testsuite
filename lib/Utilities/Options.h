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

#ifndef POOMA_UTILITIES_OPTIONS_H
#define POOMA_UTILITIES_OPTIONS_H

//-----------------------------------------------------------------------------
// Classes:
// Options
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * A simple container class that holds information about how
 * POOMA should be used when it starts up.
 *
 * The user can set up an Options
 * object with the settings they want, and give that to Pooma::initialize,
 * instead of passing argc, argv values.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <string>


namespace Pooma {

/**
 * Options stores the list of run-time configurable options for POOMA.
 * Internally, Pooma stores an Options instance that holds the values for
 * these run-time configurable values.  When the user calls Pooma::initialize,
 * they can either provide an argc, argv pair with command-line options,
 * or they can provide an Options object directly.  In the former case,
 * Pooma::initialize will create a new Options object that will parse the
 * command-line options to get the settings.  In the latter, Pooma::initialize
 * will just use the settings in the given Options object.
 *
 * An Options object just stores values for what should be used as the
 * settings, it does not call any other POOMA routines to actually affect
 * the changes.  It is meant to be used as:
 *   -# An alternative to using argc, argv if the user does not want to;
 *   -# A way for POOMA to store the run-time configuration settings.
 *
 * Options has the following types of methods:
 *   - A default constructor that just sets the options to their default values
 *   - A copy constructor
 *   - A constructor taking argc,argv, that parses the args for --pooma-*
 *     arguments.  Pooma-specific arguments are stripped out.
 *   - Query or set the value of the options.
 *   - Reset the options to their default values.
 *   - Print out a usage summary.
 */

class Options
{
public:
  //============================================================
  // Constructors
  //============================================================

  // Create an Options object with default settings.

  Options();

  // Create an Options object with a set of argc,argv values.  This
  // will parse the arguments and use them to override the current
  // settings.  POOMA-specific arguments are stripped out.

  Options(int &argc, char ** argv);

  // Copy constructor.

  Options(const Options &opts);

  // Assignment operator.

  Options &operator=(const Options &opts);


  //============================================================
  // Destructors
  //============================================================

  // Clean up any allocated storage for options.

  ~Options();


  //============================================================
  // Option accessors: routines to get or set option values.
  //============================================================

  // Return or set the current value for the concurrency

  int concurrency() const { return concurrency_m; }

  void concurrency(int c)
    {
      PAssert(c >= 1);
      concurrency_m = c;
    }

  // Return or set whether info messages should be printed or not.

  bool printInfo() const     { return info_m; }

  void printInfo(bool p)     { info_m = p; }

  // Return or set whether warning messages should be printed or not.

  bool printWarnings() const { return warn_m; }

  void printWarnings(bool p) { warn_m = p; }

  // Return or set whether error messages should be printed or not.

  bool printErrors() const   { return err_m; }

  void printErrors(bool p)   { err_m = p; }

  // Return or set the name of a log file for messages.
  // If this is an empty string, no logging is requested.

  const std::string &logfile() const { return logfile_m; }

  void logfile(const std::string &s) { logfile_m = s; }

  // Return or set whether statistics should be printed at the end.

  bool printStats() const { return stats_m; }

  void printStats(bool p) { stats_m = p; }

  // Return or set the debug output level.

  int debug() const { return debug_m; }

  void debug(int p) { debug_m = p; }

  // Return or set the "compressible" status flag.
  
  bool neverCompress() const { return neverCompress_m; }
  
  void neverCompress(bool p) { neverCompress_m = p; }
  
  // Return or set the deferred guard fill flag.
  
  bool deferredGuardFills() const { return deferredFills_m; }
  
  void deferredGuardFills(bool p) { deferredFills_m = p; }
  
  // Return or set whether hard-initialization should be used with smarts.

  bool hardInit() const { return hardinit_m; }

  void hardInit(bool p) { hardinit_m = p; }

  // Return or set whether hard run affinity should be used with smarts.

  bool hardRun() const { return hardrun_m; }

  void hardRun(bool p) { hardrun_m = p; }

  // Return or set whether threads should be locked to a processor with smarts.

  bool lockThreads() const { return lockthreads_m; }

  void lockThreads(bool p) { lockthreads_m = p; }

  // Should a block 'n evaluate be done after each expresson?

  bool blockingExpressions() const { return blockingExpressions_m; }

  void blockingExpressions(bool p) { blockingExpressions_m = p; }


  //============================================================
  // Option operations.
  //============================================================

  // Print out a usage summary for POOMA arguments to cerr.

  void usage();

  // Reset all the options to their default values.

  void reset();

  // Parse the given command-line arguments in argc,argv, and use them to
  // change the options values.  If an error occurs, print a message and
  // abort.  POOMA-specific arguments are stripped out.

  void parse(int &argc, char ** &argv);

private:
  //============================================================
  // Private data members.
  //============================================================

  // The initial level of concurrency requested.

  int concurrency_m;

  // Whether to turn on or off display of messages.

  bool info_m;
  bool warn_m;
  bool err_m;

  // A filename for logging output; if empty, do not log.

  std::string logfile_m;

  // Should we print out statistics at the end?

  bool stats_m;

  // What level of debug output should we use (if debug output is
  // enabled)?  Default is Inform::off

  int debug_m;
  
  // Flag to allow user to disable compression of Compressible-Brick engines.
  
  bool neverCompress_m;
  
  // By default, filling guards is deferred until someone tries to read
  // from the guards. Setting this flag to false will disable this behavior,
  // causing the guards to always be filled when an array is modified.
  
  bool deferredFills_m;

  // Should hard-initialization be used with smarts data?

  bool hardinit_m;

  // Should hard run affinity be used with smarts?

  bool hardrun_m;

  // Should threads be locked to a processor with smarts?

  bool lockthreads_m;

  // Should a block 'n evaluate be done after each expresson?
  
  bool blockingExpressions_m;
};

/// @name Utility functions.
//@{
/// These used to be private methods in the Options class, but they 
/// are generally useful for parsing options, so they're now in the Pooma
/// namespace.

/// A function to check for an integer argument in the pos position,
/// and set the last argument equal to it.  If it does not exist, or the
/// next argument is not an integer, return false and do not change the
/// last argument's value

bool intArgument(int argc, char **argv, int pos, int &val);

/// A function to check for a string argument in the pos position, and
/// set the last argument equal to it.  If it does not exist, return false
/// and leave the last argument unchanged.

bool stringArgument(int argc, char **argv, int pos, std::string &val);

/// A function to check for a floating-point argument in the pos position.

bool doubleArgument(int argc, char **argv, int pos, double &val);

//@}

} // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_OPTIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Options.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
