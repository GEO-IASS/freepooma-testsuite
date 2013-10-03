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

#ifndef POOMA_POOMA_POOMA_H
#define POOMA_POOMA_POOMA_H

//-----------------------------------------------------------------------------
// Classes:
//   none
//
// Functions:
//   Pooma::initialize
//   Pooma::finalize
//   Pooma::pAbort
//   Pooma::abortHandler
//   Pooma::resetAbortHandler
//
//   Pooma::version
//   Pooma::majorVersion
//   Pooma::minorVersion
//   Pooma::buildDate
//
//   Pooma::printStats
//   Pooma::debugLevel
//   Pooma::infoMessages
//   Pooma::warnMessages
//   Pooma::errorMessages
//   Pooma::logMessages
//
//   Pooma::context
//   Pooma::contexts
//   Pooma::scheduler
//   Pooma::beginGeneration
//   Pooma::endGeneration
//   Pooma::blockAndEvaluate
//   Pooma::hardInit
//   Pooma::hardRun
//   Pooma::lockThreads
//   Pooma::blockingExpressions
//   Pooma::controller
//   Pooma::poll
//
// Global objects:
//   Inform Pooma::pinfo
//   Inform Pooma::pwarn
//   Inform Pooma::perr
//   Inform Pooma::pdebug
//
// Macros:
//   POOMA_PRINT(stream, text)
//   POOMA_DEBUG(level, text)
//   POOMA_INFO(text)
//   POOMA_WARN(text)
//   POOMA_ERROR(text)
//   POOMA_DECLARE_STATISTIC(var)
//   POOMA_INIT_STATISTIC(var, name)
//   POOMA_INIT_STATISTIC_WITH(var, name, initval)
//   POOMA_INCREMENT_STATISTIC(var)
//   POOMA_INCREMENT_STATISTIC_BY(var, val)
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Pooma
 * @brief
 * Pooma.h includes all the declarations of the basic POOMA library
 * interface functions.  These general routine are used to initialize,
 * query, and shut down the POOMA library environment, including the
 * underlying run-time system.  This is generally included at the top of
 * an application source file by just using
 *
 *   #include "Pooma/Pooma.h"
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include "Threads/PoomaSmarts.h"
#include "Utilities/Inform.h"
#include "Utilities/Options.h"

#if POOMA_MESSAGING
#include "Tulip/Messaging.h"
#endif


//-----------------------------------------------------------------------------
// Macro definitions
//-----------------------------------------------------------------------------

// Print out a message to the given stream, making sure to set locks.

#define POOMA_PRINT(stream, text)     stream << text

// Print out a message to the POOMA debug output stream, but only if the
// POOMA_PRINTDEBUG option has been selected during configuration.

#if POOMA_PRINTDEBUG
# define POOMA_DEBUG(level, text) \
   Pooma::pdebug.setMessageLevel(level) << text;
#else
# define POOMA_DEBUG(level, text)
#endif

// Print out a message to the POOMA info, warn, or err output streams.

#define POOMA_INFO(text)  POOMA_PRINT(Pooma::pinfo, text)
#define POOMA_WARN(text)  POOMA_PRINT(Pooma::pwarn, text)
#define POOMA_ERROR(text) POOMA_PRINT(Pooma::perr, text)

// Set up a statistic. Declare statistics using POOMA_DECLARE_STATISTIC
// in this file and initialize them using POOMA_INIT_STATISTIC or
// POOMA_INIT_STATISTIC_WITH. A typical usage is:
//
//   POOMA_DECLARE_STATISTIC(NumExpressions)
//   POOMA_INIT_STATISTIC(NumExpressions, 
//                        "Number of expressions evaluated")
//
// Then, increment the statistic at the appropriate place using
//
//   POOMA_INCREMENT_STATISTIC(NumExpressions)
//
// or, for a non-unity change
//
//   POOMA_INCREMENT_STATISTIC_BY(NumExpressions, n)
//
// Note there are no semi-colons after the macro invocations.

#define POOMA_DECLARE_STATISTIC(var)                                       \
  void increment##var(long val = 1);

#define POOMA_INIT_STATISTIC(var,name)                                     \
  namespace {                                                              \
    Pooma::StatisticsData *stat##var##_s = statistics_s.add(name);         \
  }                                                                        \
  void increment##var(long val)                                            \
  {                                                                        \
    static Pooma::Mutex_t mutex;                                           \
    mutex.lock();                                                          \
    stat##var##_s->increment(val);                                         \
    mutex.unlock();                                                        \
  }

#define POOMA_INIT_STATISTIC_WITH(var,name,ival)                           \
  namespace {                                                              \
    Pooma::StatisticsData *stat##var##_s = statistics_s.add(name, ival);   \
  }                                                                        \
  void increment##var(long val)                                            \
  {                                                                        \
    static Pooma::Mutex_t mutex;                                           \
    mutex.lock();                                                          \
    stat##var##_s->increment(val);                                         \
    mutex.unlock();                                                        \
  }

#define POOMA_INCREMENT_STATISTIC(var)                                     \
  Pooma::increment##var();

#define POOMA_INCREMENT_STATISTIC_BY(var,val)                              \
  Pooma::increment##var(val);


//-----------------------------------------------------------------------------
// Pooma:: global functions and objects
//-----------------------------------------------------------------------------

namespace Pooma {

  //------------------------------------------------------
  // General typedefs and enumerations for POOMA.
  //------------------------------------------------------

  // The type for the POOMA abort handler routine.

  typedef void (*AbortHandler_t)();

  // The type used to refer to contexts

  typedef int Context_t;

  // The type used to refer to patches

  typedef int PatchID_t;


  //------------------------------------------------------
  // architecture-specific functions
  //------------------------------------------------------

  // Some architectures may need do special stuff to get command line
  // arguments, initialize/finalize the architecture, or pass time. This
  // can be accomplished by defining POOMA_ARCH_SPECIFIC_FUNCTIONS to 
  // POOMA_YES in a configuration file along with custom versions of
  // the functions below.
  
#if !POOMA_ARCH_SPECIFIC_FUNCTIONS

  namespace Arch {
    inline void dawdle() { }
    inline void getCommandLineArguments(int &, char** &) { }
    inline void initialize() { }
    inline void finalize() { }
  }

#endif // !POOMA_ARCH_SPECIFIC_FUNCTIONS


  //------------------------------------------------------
  // statistics
  //------------------------------------------------------

  // Here is where we declare statistics.
  //
  // Don't forget to add a corresponding 
  // POOMA_INIT_STATISTIC call to Pooma.cmpl.cpp.
  // Detailed desciptions of what a statistic measures
  // along with which source file(s) it resides in are
  // in Pooma.cmpl.cpp.
  
  POOMA_DECLARE_STATISTIC(NumExpressions)
  POOMA_DECLARE_STATISTIC(NumMultiPatchExpressions)
  POOMA_DECLARE_STATISTIC(NumZBExpressions)
  POOMA_DECLARE_STATISTIC(NumCompressedAssigns)
  POOMA_DECLARE_STATISTIC(NumAssignsRequiringUnCompression)
  POOMA_DECLARE_STATISTIC(NumInlineEvaluations)
  POOMA_DECLARE_STATISTIC(NumLocalPatchesEvaluated)
  POOMA_DECLARE_STATISTIC(NumReductions)
  POOMA_DECLARE_STATISTIC(NumUnCompresses)
  POOMA_DECLARE_STATISTIC(NumUnsuccessfulTryCompresses)
  POOMA_DECLARE_STATISTIC(NumSuccessfulTryCompresses)
  POOMA_DECLARE_STATISTIC(NumPolls)


  //------------------------------------------------------
  // Global objects
  //------------------------------------------------------

  // An Inform stream for informative messages.

  extern Inform pinfo;

  // An Inform stream for warning messages.

  extern Inform pwarn;

  // An Inform stream for error messages.

  extern Inform perr;

  // An Inform stream for debug messages.

  extern Inform pdebug;


  //------------------------------------------------------
  // Initialization, cleanup, and abort functions.
  //------------------------------------------------------

  // Initialize POOMA, using the given argc, argv values.  POOMA_specific
  // arguments will be removed.  If the 3rd argument is true, also initialize
  // the run-time system.  If the 4th argument is true, use arch-specific routine
  // to get command line arguments. If the 5th argument is true, call 
  // arch-specific initialize().
  // Return success.

  bool initialize(int &argc, char ** &argv, 
    bool initRTS = true, bool getCLArgsArch = true, bool initArch = true);

  // Initialize POOMA, using the given Options container instead of argc,argv.
  // If the 2nd argument is true, also initialize the run-time system. 
  // If the 3rd argument is true, call arch-specific initialize().
  // Return success. 

  bool initialize(Pooma::Options &opts, bool initRTS = true, 
    bool initArch = true);

  // Shut down POOMA parallism, and perform any at-exit actions that POOMA
  // is requested to do (such as print statistics or finish I/O).  This will
  // shut down the run-time system and/or call the arch-specific finalize(), 
  // but only if POOMA started the RTS and/or called the arch-specific
  // initialize().  
  // Return success.

  bool finalize();

  // Shut down POOMA. If the 1st argument is true, also shut-down the RTS.
  // If the 2nd argument is true, call arch-specific finalize().
  // Return success.
  
  bool finalize(bool quitRTS, bool quitArch);

  // Quit POOMA, with the given error code.  The abort handler for POOMA
  // will be invoked right before exit.

  void pAbort(int errorcode = 0) POOMA_ATTRIBUTE_NORETURN;

  // Quit POOMA, with the given error code, and print out the given message.
  // The abort handler for POOMA will be invoked right before exit.

  void pAbort(const char *msg, int errorcode = 0) POOMA_ATTRIBUTE_NORETURN;

  // A handy place to set break points when the debugger is not being
  // cooperative.

  void stopHere();

  // Return the current abort hander function pointer.  The abort hander
  // function should simply be a function with no return value or arguments.

  AbortHandler_t abortHandler();

  // Set the POOMA abort handler function pointer to the given value, and
  // return the previous abort handler function pointer that was changed.

  AbortHandler_t abortHandler(AbortHandler_t);

  // Reset the POOMA abort handler to the default function, and return the
  // previous abort handler function pointer that was changed.

  AbortHandler_t resetAbortHandler();


  //------------------------------------------------------
  // access general information about POOMA
  //------------------------------------------------------

  // Return a string with the POOMA version.

  const char *version();

  // Return the major version number, as an integer.

  int majorVersion();

  // Return the minor version number, as an integer.

  int minorVersion();

  // Return the build date, as a string

  const char *buildDate();


  //------------------------------------------------------
  // functions to query or modify the state of the POOMA flags
  // and I/O streams pinfo, pwarn, perr.  Streams can be turned
  // on or off, and can be logged to a file.
  //------------------------------------------------------

  // Return or set whether statistics should be printed when POOMA exits.
  // This is only useful if POOMA is compiled with the profile option.

  bool printStats();

  void printStats(bool on);

  // Return or set whether informative messages from POOMA will be printed.

  bool infoMessages();

  void infoMessages(bool on);

  // Return or set whether warning messages from POOMA will be printed.

  bool warnMessages();

  void warnMessages(bool on);

  // Return or set whether error messages from POOMA will be printed.

  bool errorMessages();

  void errorMessages(bool on);

  // Tell POOMA to echo all messages to the given log file.  Only messages
  // from context 0 will be logged.  If a null string is provided, this
  // will turn off logging if it has previously been enabled.

  void logMessages(const char *filename);

  // Return or set the debug message output level.  This is only useful
  // if POOMA is compiled with the printdebug option.

  int debugLevel();

  void debugLevel(int val);
  
  // Return or set the "compressible" status flag.
  
  bool neverCompress();
  
  void neverCompress(bool p);
  
  // Return or set the deferred guard fill flag.
  
  bool deferredGuardFills();
  
  void deferredGuardFills(bool p);
  

  //------------------------------------------------------
  // functions used to query and modify the parallelism
  //------------------------------------------------------

  extern Context_t myContext_g;
  extern int numContexts_g;
  extern int expression_g;

  // Return the context number for this context (0 ... # contexts - 1).
  
  inline Context_t context() { return myContext_g; }

  // Return the total number of contexts in use.

  inline int contexts() { return numContexts_g; }

  // Return a reference to the main RTS scheduler used by POOMA.

  Scheduler_t &scheduler();

  // Wait for all the expressions to be evaluated by the RTS.  It's
  // necessary to call this function before attempting serial access to
  // arrays that were modified by the expressions.

  void blockAndEvaluate();

  // Return or set whether SMARTS-style hard initialization of 
  // data should be used.

  bool hardInit();

  void hardInit(bool on);

  // Return whether SMARTS-style hard affinity for running iterates 
  // should be used.

  bool hardRun();

  void hardRun(bool on);

  // Return whether threads should be locked to processors.

  bool lockThreads();

  void lockThreads(bool on);

  // Return whether a blockAndEvaluate should be done after each expression

  bool blockingExpressions();

  void blockingExpressions(bool on);
  
  // begin a new expression

  inline void beginExpression()
  {
    scheduler().beginGeneration();
  }

  // end an expression

  inline void endExpression()
  {
    scheduler().endGeneration();
    expression_g++;

    if (blockingExpressions()) 
    {
      blockAndEvaluate();
    }
  }

  // give the current expression number

  inline int expression()
  {
    return expression_g;
  }

#if POOMA_CHEETAH

  extern Cheetah::Controller *controller_g;

  inline Cheetah::Controller *controller()
  {
    PAssert(controller_g != 0);
    return controller_g;
  }

#endif

// poll() is used to push the messaging system along when
// we're waiting on some condition (in the counting semaphore,
// or waiting for a remote value in RemoteProxy).  Because message
// can depend on iterates,  we also need to push iterates.
// (For example, you might be waiting for a message from another
// context, but that context is waiting for some data that an iterate
// on this context is going to produce.)

inline void poll()
{
#if POOMA_CHEETAH
  controller()->poll();
#endif
#if POOMA_SMARTS_SCHEDULER_SERIALASYNC
  Smarts::SystemContext::runSomething();
#endif

    POOMA_INCREMENT_STATISTIC(NumPolls)
}
  
} // namespace Pooma


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_POOMA_POOMA_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Pooma.h,v $   $Author: richi $
// $Revision: 1.36 $   $Date: 2004/11/10 21:51:50 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
