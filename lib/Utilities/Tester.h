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

#ifndef POOMA_UTILITIES_TESTER_H
#define POOMA_UTILITIES_TESTER_H

//-----------------------------------------------------------------------------
// Classes:
// Tester
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * A simple class with some utility methods used to write test
 * programs.
 *
 * It includes an Inform class used to print messages and
 * results, a method to update an 'OK' boolean variable about whether the
 * test is working correctly, and the ability to search for a few common
 * command-line arguments to control how the test program should function.
 * In POOMA test codes, Tester object should be created at the beginning,
 * after initializing POOMA.  It can be destroyed after finalizing POOMA.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/Inform.h"
#include "Utilities/PAssert.h"
#include <cstdlib>
#include <cmath>


namespace Pooma {


/**
 * Tester is used to help make it easier to write simple test programs.
 * It provides:
 *   - A built-in Inform stream
 *   - A boolean flag with an "OK"/"Not OK" current setting
 *   - Methods to set or update this OK flag
 *   - Methods to print out messages to the stream
 *   - Exception handlers for use in catch() blocks
 *   - The ability to parse some simple command-line options to control
 *     a test code's behavior.
 *
 * You use Tester in this way:
 * -# Initialize Pooma as normal.
 * -# Create a Tester instance, giving it the argc, argv from main:
 *        Pooma::Tester tester(argc, argv);
 * -# For each test you want to do, you can print out messages to tester,
 *    then call the "check" method with a boolean argument:
 *        tester << "This is the first test." << endl;
 *        tester.check(testval == true);
 *    If the argument to "check" is false, it will clear the flag.  Once
 *    cleared, it will stay cleared even if further "check" calls are done
 *    with an argument of true (so once one test fails, the whole test code
 *    will fail).
 * -# When you're done, you can print out a result message, like this:
 *        int retval = tester.results("Test code description string");
 *    This will print out the message:
 *        PASSED ..... message
 *    or
 *        FAILED ..... message
 *    depending on the current value of the OK flag.  It returns an exit
 *    code value, 1 if the test FAILED, 0 if the test PASSED, that you can
 *    use to return from main.  For example:
 *
 *        int retval = tester.results("My current test.");
 *        Pooma::finalize();
 *        return retval;
 *
 * The optional command-line arguments that Tester understands are: (these
 * will NOT be stripped from the argc, argv given to Tester in its constructor)
 *   - -v       : turn on verbose output from test program
 *   - -p <str> : change prefix of test program messages to <str>
 *   - -q       : do not print out anything at all, just have program
 *                return 0 or 1
 */

class Tester
{
public:
  //============================================================
  // Constructors
  //============================================================

  // Create a default Tester object.

  Tester();

  // Create a Tester object which parses argc and argv.

  Tester(int argc, char **argv);


  //============================================================
  // Destructors
  //============================================================

  ~Tester();


  //============================================================
  // Testing accessors
  //============================================================

  // Return the Inform stream used to print out test messages

  inline Inform &out()
    {
      return inform_m;
    }

  // Return the current state of the status flag

  inline bool ok() const
    {
      return ok_m;
    }

  // Return the proper main() return value.  This should be used
  // at the end of a test program to return the proper error code.

  inline int returnValue() const
    {
      return (ok() ? 0 : 1);
    }


  //============================================================
  // Testing operations
  //============================================================

  // Take the provided boolean value, and if it is false, set our
  // ok status to false.  If it is true, do not change the ok status.
  // Return the value of the boolean we are checking.

  inline bool check(bool val)
    {
      ok_m = (ok_m && val);
      if (!ok_m && abort_m) 
        {
          PInsist(0,"Check failed!");
        }
      return val;
    }

  // Take the provided boolean value (2nd arg), and if it is false, set our
  // ok status to false.  If it is true, do not change the ok status.
  // Return the value of the boolean we are checking.  This version 
  // of check will also print out a message of the form:
  //   "Checking <string>: check = <val>, updated status = <status>"

  inline bool check(const char *str, bool val)
    {
      check(val);
      if (str != 0)
	out() << "Checking " << str;
      else
	out() << "Checking";
      out() << ": check = " << val << ", updated status = " << ok_m;
      out() << std::endl;
      return val;
    }

  template<class T>    
  bool check(const char *str, const T &val, const T &correct)
    {
      bool res = check(val == correct);
      if (str != 0)
	out() << "Checking " << str;
      else
	out() << "Checking";
      out() << ": val = " << val << ", correct = " << correct
            << ", updated status = " << ok_m;
      out() << std::endl;
      return res;
    }

  template<class T>    
  bool check(const char *str, const T &val, const T &correct, 
    const T &tol)
    {
      bool res = check(std::abs(val - correct) < tol);
      if (str != 0)
	out() << "Checking " << str;
      else
	out() << "Checking";
      out() << ": val = " << val << ", correct = " << correct
            << ", updated status = " << ok_m;
      out() << std::endl;
      return res;
    }

  // Just set the OK flag to the given value

  inline void set(bool val)
    {
      ok_m = val;
    }

  // Setters for the quiet_m member and output prefix (needed at least for now
  // on Windows with CodeWarrior, where command-line arguments don't exist.

  inline void setQuiet(bool quiet) { 
    quiet_m = quiet;
    // If we're not verbose, turn off the inform stream.
    if (!verbose_m || quiet_m) {
      out().setOutputLevel(Inform::off);
    } else {
      out().setOutputLevel(Inform::on);
    }
  }
  inline void setVerbose(bool verbose) { 
    verbose_m = verbose;
    // If we're not verbose, turn off the inform stream.
    if (!verbose_m || quiet_m) {
      out().setOutputLevel(Inform::off);
    } else {
      out().setOutputLevel(Inform::on);
    }
  }
  
  inline bool verbose() { return verbose_m; }
  
  inline void setPrefix(char *prefix) { out().setPrefix(prefix); }
  

  // Print out a message to cout about the current status.  If a
  // string is given, print that message out on the same line.
  // This will either print PASSED or FAILED as the first word.
  // This will also return the current main() return code setting.

  int results(const char *msg = 0) const;

  // Exception handlers for use in catch() blocks.  The handlers print
  // notifications of exceptions.

  void exceptionHandler(const char *msg = 0);
  void exceptionHandler(const Assertion &asrt);

private:
  //============================================================
  // Private data members.
  //============================================================

  // The status of the test

  bool ok_m;

  // Should we be quiet about printing everything out?

  bool quiet_m;

  // An Inform stream used to print out messages.

  Inform inform_m;

  // Turn on/off the Inform output stream inform_m:

  bool verbose_m;

  // If this is set, failure will cause a PInsist to be called.
  
  bool abort_m;
  
  //============================================================
  // Private methods.
  //============================================================

  // Parse the arguments and set up our testing state.

  void parse(int argc, char **argv);
};

} // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_TESTER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Tester.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
