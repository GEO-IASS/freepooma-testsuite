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

/** @file
 * @ingroup Utilities
 * @brief
 * Benchmark framework.
 */

#ifndef POOMA_UTILITIES_BENCHMARK_H
#define POOMA_UTILITIES_BENCHMARK_H

//-----------------------------------------------------------------------------
// Classes: 
//   Implementation
//   Benchmark
//-----------------------------------------------------------------------------

#include "Utilities/Inform.h"

#include <cmath>
#include <string>
#include <vector>


/**
 * Implementation provides a framework for implementing a benchmark in a
 * specific way. It is a virtual base class. Users must override almost
 * all of the member functions.
 */

class Implementation {
public:

  //---------------------------------------------------------------------------
  // Make the destructor virtual since we will usually be deleting 
  // Implementations through a base-class pointer.

  virtual ~Implementation() { };
  
  //---------------------------------------------------------------------------
  // Returns the type of this implementation.
  // For example, C, C++Tran, etc.

  virtual const char *type() const = 0;

  //---------------------------------------------------------------------------
  // Returns the a qualification for the type of this implementation.
  // For example, UMP, Opt, etc.

  virtual const char *qualification() const { return ""; }

  //---------------------------------------------------------------------------
  // Performs initialization for the specified problem size.

  virtual void initialize(int n) = 0;

  //---------------------------------------------------------------------------
  // Runs the benchmark.

  virtual void run() = 0;

  //---------------------------------------------------------------------------
  // Runs a virtual function that computes overhead of timing and any setup.

  virtual void runSetup() { }

  //---------------------------------------------------------------------------
  // Returns a value to give the user a sliver of belief that the benchmark
  // ran correctly. 

  virtual double resultCheck() const = 0;

  // --------------------------------------------------------------------------
  // Return the op count of the kernel being benchmarked.

  virtual double opCount() const = 0;  
  
  // --------------------------------------------------------------------------
  // Return whether the implementation has internal Pooma::Clock::value()
  // calls assigning stop and start double values or (default) not:

  virtual bool internalClockCalls() const { return false; }

  // --------------------------------------------------------------------------
  // If our run() method has internal clock calls, this method can
  // be used to retrieve the timing result.

  virtual double internalTimingResult() const { return 0.0; }

  // --------------------------------------------------------------------------
  // Return whether the implementation specifies forcing the running of
  // multiple iterations even if highSpeed timers are in use:

  virtual bool forceMultipleIterations() const { return false; }

  //---------------------------------------------------------------------------
  // Some canned implementation types.
  
  static const char *CType() { return "C"; }
  static const char *CppType() { return "C++"; }
  static const char *P2Type() { return "PoomaII"; }
  static const char *CppTranType() { return "CppTran"; }
  static const char *F77Type() { return "Fortran77"; }
  static const char *F90Type() { return "Fortran90"; }
  static const char *BlitzType() { return "Blitz++"; }

};
 

/**
 * Benchmark provides a user interface for creating, running, and tabulating
 * results for benchmarks.
 */

class Benchmark {
public:

  // --------------------------------------------------------------------------
  // Constructor. Parses arguments to configure benchmark. Also can specify a
  // benchmark variation name here. Finally, the outputContext, context, and
  // nContexts arguments specify the context information used by the Inform
  // class (Benchmark uses an Inform object for printing output).

  Benchmark(int argc, char *argv[], const char *varName = "",
            Inform::Context_t outputContext = 0);

  //---------------------------------------------------------------------------
  // Destructor. Deletes the implementations. 

  virtual ~Benchmark();

  //---------------------------------------------------------------------------
  // Adds an implementation to this benchmark.

  void addImplementation(Implementation *impl);
  
  //---------------------------------------------------------------------------
  // Runs the benchmark and optionally prints out results.

  void run();

  //---------------------------------------------------------------------------
  // Sets the number of run kernel iterations,
  // assuming the value has not been previously set on the command line.

  void setIterations(long iters);

  //---------------------------------------------------------------------------
  // Sets the sampling parameters assuming these values have not been 
  // previously set by the command line.

  void setSamplingParameters(int startSize, int numDecades, int numPoints);

  //---------------------------------------------------------------------------
  // Sets the default number of patches for multi-patch arrays to use
  // assuming this value hasn't been set by the command line.
  
  void setNumPatches(int numPatches);
  
  //---------------------------------------------------------------------------
  // Returns the default number of patches for multi-patch arrays to use.
  
  int numPatches() const { return numPatches_m; }
  
  // --------------------------------------------------------------------------
  // Point the Inform object pointer to a user-specified Inform object:
  
  void setInform(Inform *inform) { inform_m = inform; }
  
private:

  //---------------------------------------------------------------------------
  // State variables that tell whether or not we've specified various things
  // from the command line.
  
  bool setIterations_m;
  bool setParams_m;
  bool setNumPatches_m;
  bool setSamples_m;
  
  //---------------------------------------------------------------------------
  // If true, we are supposed to display results or show diagnostic
  // output or print running time, not Mflops.

  bool print_m;
  bool diags_m;
  bool report_time_m;

  // --------------------------------------------------------------------------
  // Points to the Inform object pointer used for printing output:
  
  Inform *inform_m;

  //---------------------------------------------------------------------------
  // If true, we are supposed to test for validity.

  bool test_m;
  
  //---------------------------------------------------------------------------
  // If true, only list the implementations and then exit.
  
  bool listimpls_m;
  
  //---------------------------------------------------------------------------
  // The number of times that we are supposed to run an implementation.

  long iters_m;
    
  //---------------------------------------------------------------------------
  // The number of decades of problem size to sample, the number of points
  // per decade, and the starting point for the sample.
  
  int decades_m, points_m, start_m;

  //---------------------------------------------------------------------------
  // The default number of patches for multi-patch arrays to use.
  
  int numPatches_m;

  //---------------------------------------------------------------------------
  // The default number of samples to use.
  
  int samples_m;
  
  //---------------------------------------------------------------------------
  // The name of this benchmark and any variation.

  std::string name_m, variation_m;
    
  //---------------------------------------------------------------------------
  // A container giving the index of the Implementations to actually run in
  // this benchmark.
  
  std::vector<int> implsToRun_m;
    
  //---------------------------------------------------------------------------
  // A container giving the names of the Variations to actually run in
  // this executable.
  
  std::vector<std::string> varsToRun_m;
    
  //---------------------------------------------------------------------------
  // A container holding the Implementations owned by this benchmark.
  
  std::vector<Implementation *> impls_m;
    
  //---------------------------------------------------------------------------
  // A container holding the Implementations that actually ran during the
  // immediately preceding benchmark.
  
  std::vector<Implementation *> implsRan_m;
    
  //---------------------------------------------------------------------------
  // A container holding the times for each sample for each implementation
  // that actually ran during the immediately preceding benchmark.
  
  std::vector< std::vector<double> > times_m;

  //---------------------------------------------------------------------------
  // Reinitializes result data in preparation for another run.
  
  void getReadyToRun();
  
  //---------------------------------------------------------------------------
  // Decides which implementations to run and runs them.

  void runIt();
  
  //---------------------------------------------------------------------------
  // Prints results for all of the implementations that ran.
  
  void printResults();
  
  //---------------------------------------------------------------------------
  // Runs a specific implementation, tests the results, and stores timing and
  // validity data.

  void runImplementation(Implementation *impl, int sample);

  //---------------------------------------------------------------------------
  // Computes the number of points for the ith trial.
  
  int trialPoints(int i)
    {
      return int(start_m * pow(10.0, double(i) / double(points_m)));
    }
    
  //---------------------------------------------------------------------------
  // Computes the total number of sample points.

  int numPoints() const
    {
      return decades_m * points_m + 1;
    }
    
  //---------------------------------------------------------------------------
  // Returns the total number of samples.

  int numSamples() const
    {
      return samples_m;
    }
    
  //---------------------------------------------------------------------------
  // Prints an ordered list of implementations

  void printImplementations();
  
  //---------------------------------------------------------------------------
  // Prints out a usage message.

  static void usage(const char *name);
  
  
  //---------------------------------------------------------------------------
  // Prints an input error and the usage message.
  
  static void inputError(const char *msg, const char *progname, int err = 1);
  
};

#endif // POOMA_UTILITIES_BENCHMARK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Benchmark.h,v $   $Author: richi $
// $Revision: 1.32 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
