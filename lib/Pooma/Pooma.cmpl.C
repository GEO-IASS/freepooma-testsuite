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
// Pooma general interface definitions.
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "Pooma/Pooma.h"
#include "Threads/PoomaMutex.h"
#include "Threads/PoomaSmarts.h"
#include "Tulip/Messaging.h"
#include "Tulip/ReduceOverContexts.h"
#include "Utilities/Inform.h"
#include "Utilities/Options.h"
#include "Utilities/PAssert.h"
#include "Utilities/Statistics.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

#if POOMA_MESSAGING
# include "Tulip/Messaging.h"
#endif

//-----------------------------------------------------------------------------
// The POOMA library global state objects, accessible by users
//-----------------------------------------------------------------------------

namespace Pooma {

// An Inform stream for informative messages.

Inform pinfo("Pooma");

// An Inform stream for warning messages.

Inform pwarn("Warning", std::cerr, Inform::allContexts);

// An Inform stream for error messages.

Inform perr("Error", std::cerr, Inform::allContexts);

// An Inform stream for debug messages.

Inform pdebug("** Debug **", std::cerr, Inform::allContexts);

#if POOMA_CHEETAH

Cheetah::Controller *controller_g = 0;

#endif

//-----------------------------------------------------------------------------
// The POOMA library global state internal objects, only used here
//-----------------------------------------------------------------------------

// The context objects are accessed through Pooma::context() and
// Pooma::contexts().  Since these functions may be called dozens of times
// per iterate generated, they need to be inlined, which is why this data
// isn't included in the following localized scope.

// Our context number (from 0 ... # contexts - 1)

Context_t myContext_g = 0;

// The total number of contexts

int numContexts_g = 1;

// The current expression number

int expression_g = 0;

  namespace {                   
    // anonymous - these things have file scope! 
    // In order to make these stand out, I've added an "_s" suffix
    // to everything defined here (JAC).

    // Has the initialize() function been called yet?

    bool initialized_s = false;

    // Did we initialize the RTS and/or arch, or did the user?

    bool weInitializedRTS_s = false, weInitializedArch_s = false;

    // The Options object used to store information about the current options.

    Options options_s;

    // The POOMA threads scheduler.  The type for Scheduler_t is determined
    // in Threads/PoomaSmarts.h.

    Scheduler_t mainScheduler_s;
    
    // The POOMA statistics object.
    
    Statistics statistics_s;

    // An output stream that might be used to log output messages.

    std::ofstream *logstream_s = 0;
    
    Inform::ID_t pinfoLogID_s;
    Inform::ID_t pwarnLogID_s;
    Inform::ID_t perrLogID_s;
    Inform::ID_t pdebugLogID_s;

    // A default abort handler function, that just does nothing.

    void defAbortHandler_s()
    {
      std::cerr << "In default abort handler." << std::endl;
    }

    // The abort handler function to call if pAbort() is called.

    AbortHandler_t currentAbortHandler_s = defAbortHandler_s;

    // A filter function to perform a reduction on context zero and return
    // the result.
    
    long reductionFilter_s(long val)
    {
      ReduceOverContexts<long, OpAddAssign> reduce(val, 0);
      return reduce;
    }
    
    // A simple routine used here to clean up the Pooma state before the
    // program exits.

    void cleanup_s()
    {
      // If we're supposed to print out statistics, do so now.
      
      if (printStats())
        statistics_s.print(pinfo, reductionFilter_s);

      // Close the log file, if necessary.

      logMessages(0);

      // Turn off all POOMA streams.

      infoMessages(false);
      warnMessages(false);
      errorMessages(false);
      debugLevel(Inform::off);
    }

  } // anonymous namespace

//-----------------------------------------------------------------------------
// Initialize statistics here.
//
// Don't forget to add a corresponding 
// POOMA_DECLARE_STATISTIC call to Pooma.h.
// Include detailed desciptions below of what a statistic measures
// along with which source file(s) it resides.
//-----------------------------------------------------------------------------

// Evaluator/Evaluator.h
// The number of times Evaluator<MainEvaluatorTag>::evaluate() is called.

POOMA_INIT_STATISTIC(NumExpressions, 
  "Number of expressions evaluated")

// Evaluator/Evaluator.h
// The number of times Evaluator<MainEvaluatorTag>::evaluateZeroBased() 
// is called.

POOMA_INIT_STATISTIC(NumZBExpressions, 
  "Number of zero-based expressions evaluated")

// Evaluator/Evaluator.h
// The number of times Evaluator<MultiPatchEvaluatorTag>::evaluate() is called.

POOMA_INIT_STATISTIC(NumMultiPatchExpressions, 
  "Number of multi-patch expressions evaluated")

// Evaluator/CompressibleEval.h
// The number of times a fully compressed assignment (single number) is
// performed. Appears two places inside of 
// KernelEvaluator<CompressibleKernelTag>.

POOMA_INIT_STATISTIC(NumCompressedAssigns, 
  "Number of fully compressed assignments")

// Evaluator/CompressibleEval.h
// The number of times a either the RHS or LHS must be uncompressed to
// do an assignment. Appears in KernelEvaluator<CompressibleViewKernelTag>::
// evaluate().

POOMA_INIT_STATISTIC(NumAssignsRequiringUnCompression, 
  "Number of assignments requiring uncompression")

// Evaluator/InlineEvaluator.h
// The number of times a KernelEvaluator<InlineKernelTag>::evaluate(). 
// is called (both versions).

POOMA_INIT_STATISTIC(NumInlineEvaluations, 
  "Number of assignments using the inline evaluator")

// Evaluator/Evaluator.h
// The number of patches sent to single patch evaluators from the multi-
// patch version. Appears in Evaluator<MultiPatchEvaluatorTag>::evaluate().

POOMA_INIT_STATISTIC(NumLocalPatchesEvaluated, 
  "Number of local patches evaluated")

// Evaluator/Reduction.h
// The times Reduction<MainEvaluatorTag>::evaluate() is called.

POOMA_INIT_STATISTIC(NumReductions, 
  "Number of reductions performed")

// Engine/CompressibleBlock.h
// The number of times a compressible block uncompresses.

POOMA_INIT_STATISTIC(NumUnCompresses, 
  "Number of times a compressible block uncompresses")

// Engine/CompressibleBlock.h
// The number of times a tryCompress() call fails.

POOMA_INIT_STATISTIC(NumUnsuccessfulTryCompresses, 
  "Number of times a compression attempt fails")

// Engine/CompressibleBlock.h
// The number of times a tryCompress() call succeeds.

POOMA_INIT_STATISTIC(NumSuccessfulTryCompresses, 
  "Number of times a compression attempt succeeds")

// Pooma/Pooma.h
// The number of times Pooma::poll() is called.

POOMA_INIT_STATISTIC(NumPolls, 
  "Number of calls to Pooma::poll()")
    

//-----------------------------------------------------------------------------
// Starts up POOMA, does all required initializations, and initializes the
// run-time system.  Special command-line arguments can be used to indicate
// how the run-time system should operate.  These POOMA arguments are stripped
// out of the list.  Return success.
//-----------------------------------------------------------------------------

bool initialize(int &argc, char ** &argv, bool initRTS, bool getCLArgsArch,
  bool initArch)
{
  // Do special command line argument initialization, if necessary.

  if (getCLArgsArch)
    Pooma::Arch::getCommandLineArguments(argc, argv);

  // Initialize Cheetah, if necessary.  Command-line arguments are used
  // to initialize Cheetah.  Eventually, we should have a way to
  // initialize Cheetah with pooma-related command-line options.  Then
  // we can do this in the other initialize routine by querying for
  // the Cheetah options from the Options object.

#if POOMA_MPI
# ifdef _OPENMP
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  PInsist(provided >= MPI_THREAD_FUNNELED, "No MPI support for OpenMP");
# else
  MPI_Init(&argc, &argv);
# endif
#elif POOMA_CHEETAH
  controller_g = new Cheetah::Controller(argc, argv);
#endif

  // Just create an Options object for this argc, argv set, and give that
  // to the other version of initialize.  Options will parse the arguments.

  Options opts(argc, argv);

  return initialize(opts, initRTS, initArch);
}


//-----------------------------------------------------------------------------
// Starts up POOMA, does all required initializations, and initializes the
// run-time system.  Special command-line arguments can be used to indicate
// how the run-time system should operate.  These POOMA arguments are stripped
// out of the list.  Return success.
//-----------------------------------------------------------------------------

bool initialize(Options &opts, bool initRTS, bool initArch)
{  
  // If we've already been initialized, it is an error.

  PInsist(!initialized_s, "You can only call Pooma::initialize once.");

  // Indicate for future reference that we've done initialization and whether
  // we starts the RTS or not.

  initialized_s = true;
  weInitializedRTS_s = initRTS;
  weInitializedArch_s = initArch;

  // Initialize special for the architecture.
  
  if (initArch)
    Pooma::Arch::initialize();

  // Set the debug level.

  debugLevel(opts.debug());

  // Save the options_s for future reference.

  options_s = opts;

  // Now, initialize the Run-Time System, if requested and we're compiled
  // with parallelism.

  if (initRTS)
  {
    // set the Run Time System concurrency

    Smarts::concurrency(opts.concurrency());

    // PUT OTHER INIT RTS CODE HERE
  }

  // Set myContext_s and numContexts_s to the context numbers.

#if POOMA_MESSAGING

#if POOMA_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myContext_g);
  MPI_Comm_size(MPI_COMM_WORLD, &numContexts_g);
  for (int i=0; i<Smarts::SystemContext::max_requests; ++i)
    Smarts::SystemContext::free_requests_m.insert(i);
#elif POOMA_CHEETAH
  PAssert(controller_g != 0);

  myContext_g   = controller_g->mycontext();
  numContexts_g = controller_g->ncontexts();
#endif

  initializeCheetahHelpers(numContexts_g);

#else

  myContext_g   = 0;
  numContexts_g = 1;

#endif

  // Enable logging to a file, if requested.

  logstream_s = 0;
  logMessages(opts.logfile().c_str());

  // Enable or disable the output streams.

  infoMessages(opts.printInfo());
  warnMessages(opts.printWarnings());
  errorMessages(opts.printErrors());

  // This barrier is here so that Pooma is initialized on all contexts
  // before we continue.  (Another context could invoke a remote member
  // function on us before we're initialized... which would be bad.)

#if POOMA_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#elif POOMA_CHEETAH
  controller_g->barrier();
#endif

  // Initialize the Inform streams with info on how many contexts we
  // have and what context we're running on.  We cannot do this until
  // now since the run-time system has not been fully initialized.

  Inform::setContext(myContext_g);
  Inform::setNumContexts(numContexts_g);

  return true;
}


//-----------------------------------------------------------------------------
// Shut down POOMA parallism, and perform any at-exit actions that POOMA
// is requested to do (such as print statistics or finish I/O).  This will
// shut down the run-time system, but only if POOMA started it.  Return
// success.
//-----------------------------------------------------------------------------

bool finalize()
{
  return finalize(weInitializedRTS_s, weInitializedArch_s);
}


//-----------------------------------------------------------------------------
// Shut down POOMA. Return success.
//-----------------------------------------------------------------------------

bool finalize(bool quitRTS, bool quitArch)
{
  Pooma::blockAndEvaluate();

  if (initialized_s)
  {
    // Wait for threads to finish.

    Smarts::wait();

    // Do other POOMA cleanup tasks.

    cleanup_s();

#if POOMA_MESSAGING
    // Clean up the Cheetah helpers.

    finalizeCheetahHelpers();
#endif

    // Shut down the RTS, if necessary.

    if (quitRTS)
    {
#if POOMA_MESSAGING

      // Deleting the controller shuts down the cross-context communication
      // if this is the last thing using this controller.  If something
      // else is using this, Cheetah will not shut down until that item
      // is destroyed or stops using the controller.

#if POOMA_MPI
      MPI_Finalize();
#elif POOMA_CHEETAH
      if (controller_g != 0)
	delete controller_g;
#endif

#endif
    }
  }

  // Perform special shutdown tasks.

  if (quitArch)
    Pooma::Arch::finalize();
  
  return true;
}


//-----------------------------------------------------------------------------
// Quit POOMA, with the given error code.  The abort handler for POOMA
// will be invoked right before exit.
//-----------------------------------------------------------------------------

void pAbort(int errorcode)
{
  pAbort("Pooma::pAbort called.", errorcode);
}


//-----------------------------------------------------------------------------
// Quit POOMA, with the given error code, and print out the given message.
// The abort handler for POOMA will be invoked right before exit.
//-----------------------------------------------------------------------------

void pAbort(const char *msg, int)
{
  // Print out the message, to cerr.

  if (msg != 0)
    std::cerr << msg << std::endl;

  // Invoke the abort handler.

  currentAbortHandler_s();

  // Do other POOMA cleanup tasks.
  // (Only if POOMA is initialized.)

  if (initialized_s)
  {
    cleanup_s();
  }

  // And that's all she wrote ...
  ::abort();
}


//-----------------------------------------------------------------------------
// Return the current abort hander function pointer.  The abort hander
// function should simply be a function with no return value or arguments.
//-----------------------------------------------------------------------------

AbortHandler_t abortHandler()
{
  return currentAbortHandler_s;
}


//-----------------------------------------------------------------------------
// Set the POOMA abort handler function pointer to the given value, and
// return the previous abort handler function pointer that was changed.
//-----------------------------------------------------------------------------

AbortHandler_t abortHandler(AbortHandler_t ah)
{
  AbortHandler_t oldah = currentAbortHandler_s;
  currentAbortHandler_s = ah;
  return oldah;
}


//-----------------------------------------------------------------------------
// Reset the POOMA abort handler to the default function, and return the
// previous abort handler function pointer that was changed.
//-----------------------------------------------------------------------------

AbortHandler_t resetAbortHandler()
{
  return abortHandler(defAbortHandler_s);
}


//-----------------------------------------------------------------------------
// Return a string with the POOMA version.
//-----------------------------------------------------------------------------

const char *version()
{
  PAssert(initialized_s);
  return POOMA_VERSION_STRING;
}


//-----------------------------------------------------------------------------
// Return the major version number, as an integer.
//-----------------------------------------------------------------------------

int majorVersion()
{
  PAssert(initialized_s);
  return POOMA_MAJOR_VERSION;
}


//-----------------------------------------------------------------------------
// Return the minor version number, as an integer.
//-----------------------------------------------------------------------------

int minorVersion()
{
  PAssert(initialized_s);
  return POOMA_MINOR_VERSION;
}


//-----------------------------------------------------------------------------
// Return the build date, as a string
//-----------------------------------------------------------------------------

const char *buildDate()
{
  PAssert(initialized_s);
  return POOMA_BUILD_DATE;
}


//-----------------------------------------------------------------------------
// Return or set whether statistics should be printed when POOMA exits.
//-----------------------------------------------------------------------------

bool printStats()
{
  PAssert(initialized_s);
  return options_s.printStats();
}

void printStats(bool on)
{
  PAssert(initialized_s);
  options_s.printStats(on);
}


//-----------------------------------------------------------------------------
// Return or set whether informative messages from POOMA will be printed.
//-----------------------------------------------------------------------------

bool infoMessages()
{
  PAssert(initialized_s);
  return (pinfo.outputLevel() >= 0);
}

void infoMessages(bool on)
{
  PAssert(initialized_s);
  pinfo.setOutputLevel(on ? Inform::on : Inform::off);
}


//-----------------------------------------------------------------------------
// Return or set whether warning messages from POOMA will be printed.
//-----------------------------------------------------------------------------

bool warnMessages()
{
  PAssert(initialized_s);
  return (pwarn.outputLevel() >= 0);
}

void warnMessages(bool on)
{
  PAssert(initialized_s);
  pwarn.setOutputLevel(on ? Inform::on : Inform::off);
}


//-----------------------------------------------------------------------------
// Return or set whether error messages from POOMA will be printed.
//-----------------------------------------------------------------------------

bool errorMessages()
{
  PAssert(initialized_s);
  return (perr.outputLevel() >= 0);
}

void errorMessages(bool on)
{
  PAssert(initialized_s);
  perr.setOutputLevel(on ? Inform::on : Inform::off);
}


//-----------------------------------------------------------------------------
// Tell POOMA to echo all messages to the given log file.  Only messages
// from context 0 will be logged.  If a null string is provided, this
// will turn off logging if it has previously been enabled.
//-----------------------------------------------------------------------------

void logMessages(const char *filename)
{
  PAssert(initialized_s);

  // If we're currently logging to a file, stop doing so.

  if (logstream_s != 0)
    {
      // Remove this stream from the Inform objects.

      pinfo.close(pinfoLogID_s);
      pwarn.close(pwarnLogID_s);
      perr.close(perrLogID_s);
      pdebug.close(pdebugLogID_s);
      
      // Close the file.

      delete logstream_s;
      logstream_s = 0;
    }

  // If a filename is given, start logging to that file.

  if (filename != 0 && *filename != 0)
    {
      // Create a new file

      logstream_s = new std::ofstream(filename, std::ios::out);

      // Tell the main streams to start logging to this file.

      pinfoLogID_s   = pinfo.open(*logstream_s);
      pwarnLogID_s   = pwarn.open(*logstream_s);
      perrLogID_s    = perr.open(*logstream_s);
      pdebugLogID_s  = pdebug.open(*logstream_s);

      // Set the output level for all the streams properly

      pinfo.setOutputLevel(pinfo.outputLevel());
      pwarn.setOutputLevel(pwarn.outputLevel());
      perr.setOutputLevel(perr.outputLevel());
      pdebug.setOutputLevel(pdebug.outputLevel());
    }
}


//-----------------------------------------------------------------------------
// Return or set the debug message output level.  This is only useful
// if POOMA is compiled with the printdebug option.
//-----------------------------------------------------------------------------

int debugLevel()
{
  PAssert(initialized_s);
  return pdebug.outputLevel();
}

void debugLevel(int level)
{
  PAssert(initialized_s);
  pdebug.setOutputLevel(level);
}


//-----------------------------------------------------------------------------
// Return or set the "non-compressible" status flag.
//-----------------------------------------------------------------------------
 
bool neverCompress() 
{ 
  PAssert(initialized_s);
  return options_s.neverCompress();
}
  
void neverCompress(bool p) 
{ 
  PAssert(initialized_s);
  options_s.neverCompress(p); 
}
  

//-----------------------------------------------------------------------------
// Return or set the deferred guard fill flag.
//-----------------------------------------------------------------------------
  
bool deferredGuardFills() 
{ 
  PAssert(initialized_s);
  return options_s.deferredGuardFills(); 
}
  
void deferredGuardFills(bool p) 
{ 
  PAssert(initialized_s);
  options_s.deferredGuardFills(p); 
}
  

//-----------------------------------------------------------------------------
// Return a reference to the RTS scheduler used by r2.
// (The actual scheduler is local to this file.)
//-----------------------------------------------------------------------------

Scheduler_t &scheduler()
{
  PAssert(initialized_s);
  return mainScheduler_s;
}


//-----------------------------------------------------------------------------
// Wait for all the expressions to be evaluated by the RTS.  It's
// necessary to call this function before attempting serial access to
// arrays that were modified by the expressions.
//-----------------------------------------------------------------------------

void blockAndEvaluate()
{
  PAssert(initialized_s);

#if POOMA_CHEETAH

# if POOMA_SMARTS_SCHEDULER_SERIALASYNC

  typedef Smarts::SystemContext SystemContext_t;

  while (Pooma::incomingMessages() || SystemContext_t::workReady())
  {
    controller_g->poll();
    SystemContext_t::runSomething();
  }

# else // we're using the serial scheduler, so we only need to get messages

  while (Pooma::incomingMessages())
  {
    controller_g->poll();
  }

# endif // schedulers

#else // !POOMA_CHEETAH

  mainScheduler_s.blockingEvaluate();

#endif // !POOMA_CHEETAH
}


//-----------------------------------------------------------------------------
// Return or set whether SMARTS hard initialization of data should be used.
//-----------------------------------------------------------------------------

bool hardInit()
{
  PAssert(initialized_s);
  return options_s.hardInit();
}

void hardInit(bool on)
{
  PAssert(initialized_s);
  options_s.hardInit(on);
}


//-----------------------------------------------------------------------------
// Return whether SMARTS hard affinity for running iterates should be used.
//-----------------------------------------------------------------------------

bool hardRun()
{
  PAssert(initialized_s);
  return options_s.hardRun();
}

void hardRun(bool on)
{
  PAssert(initialized_s);
  options_s.hardRun(on);
}

//-----------------------------------------------------------------------------
// Return whether threads should be locked to processors.
//-----------------------------------------------------------------------------

bool lockThreads()
{
  PAssert(initialized_s);
  return options_s.lockThreads();
}

void lockThreads(bool on)
{
  PAssert(initialized_s);
  options_s.lockThreads(on);
}

//-----------------------------------------------------------------------------
// Return whether threads should be locked to processors.
//-----------------------------------------------------------------------------

bool blockingExpressions()
{
  PAssert(initialized_s);
  return options_s.blockingExpressions();
}

void blockingExpressions(bool on)
{
  PAssert(initialized_s);
  options_s.blockingExpressions(on);
}

} // namespace Pooma


//-----------------------------------------------------------------------------
// A handy place for setting break points when the debugger is buggy.
//-----------------------------------------------------------------------------

extern "C" void pooma_stop_here();

void Pooma::stopHere() 
{ 
  ::pooma_stop_here();
}

extern "C" void pooma_stop_here() 
{ 
  // This space intentionally left blank.
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Pooma.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.44 $   $Date: 2004/11/01 18:17:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
