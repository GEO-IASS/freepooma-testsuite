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
// Benchmark non-template definitions.
//-----------------------------------------------------------------------------

#include "Utilities/Benchmark.h"
#include "Utilities/Clock.h"
#include "Utilities/Options.h"
#include "Utilities/PAssert.h"
#include "Pooma/Configuration.h"
#include "Pooma/Pooma.h"

using Pooma::intArgument;
using Pooma::stringArgument;

#include <algorithm>

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>
#include <iomanip>
#include <ios>

///////////////////////////////////////////////////////////////////////////////
//
// Member functions for class Benchmark.
//
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// Benchmark constructor.  Parses the command line arguments to configure the
// benchmark. These arguments will override values subseqently set in member
// functions.  Member functions can be used to override values not set in the
// command line and the defaults set here.
// ----------------------------------------------------------------------------

Benchmark::Benchmark(int argc, char *argv[], const char *varName, 
                     Inform::Context_t outputContext)
{
  // As of now, we have not set anything.

  setIterations_m = false;
  setParams_m = false;
  setNumPatches_m = false;
  setSamples_m = false;
  listimpls_m = false;

  // By default, we will print results and show diagnostics.

  print_m = true;
  diags_m = true;
  report_time_m = false;

  // Default Inform object has null prefix and only prints from context 0:

  inform_m = new Inform("", outputContext);

  // By default, we test.

  test_m = true;

  // By default, we will run the benchmark kernel for 10 iterations, unless
  // we have a slow clock.

  iters_m = 10;

  // By default, we'll run 0 decades with 1 point per decade starting at 100.
  
  decades_m = 0;
  points_m = 1;
  start_m = 100;
  
  // By default, we'll use one sample.
  
  samples_m = 1;

  // By default, we'll use 10 patches.
  
  numPatches_m = 10;
  
  // Set the name of the benchmark from the name of the executable.

  name_m = argv[0];

  // Get the variation name, if any.

  variation_m = varName;

  // Go through all of the arguments. Setting parameters apprpriately.

  bool convertOK = true;
  
  int i = 1;
  while (i < argc)
    {
      if (strcmp("--run-impls", argv[i]) == 0)
	{
	  // Specify implementations to run. This should be a list of
	  // names with none of them starting with a "-".
	  
	  int j = i + 1;
	  if (j == argc) 
	    inputError("--run-impls needs at least one argument", argv[0]);
	  
	  while (j < argc)
	    {
	      if (argv[j][0] == '-') 
	        {
	          if (j == i + 1)
	            inputError("--run-impls needs at least one argument",
	                       argv[0]);
	          break;
	        }
              int num;
	      convertOK = intArgument(argc, argv, j, num);
	      if (!convertOK) 
	        inputError("Bad argument to --run-impls", argv[0]);
              implsToRun_m.push_back(num);
              ++j;
	    }
	  i = j;
	}
      else if (strcmp("--run-vars", argv[i]) == 0)
	{
	  // Specify variations to run. This should be a list of
	  // names with none of them starting with a "-".
	  
	  int j = i + 1;
	  if (j == argc) inputError("--run-vars needs at least one argument", 
	                            argv[0]);

	  while (j < argc)
	    {
	      if (argv[j][0] == '-') 
	        {
	          if (j == i + 1)
	            inputError("--run-vars needs at least one argument",
	                       argv[0]);
	          break;
	        }
	      std::string s;
	      convertOK = stringArgument(argc, argv, j, s);
	      if (!convertOK) inputError("Bad argument to --run-vars", argv[0]);
              varsToRun_m.push_back(s);
	    }
	  i = j;
	}
      else if (strcmp("--samples", argv[i]) == 0)
	{
	  // Specify the number of samples we'll run.

	  setSamples_m = true;
	  convertOK = intArgument(argc, argv, i+1, samples_m);
	  if (!convertOK) inputError("Bad argument to --samples", argv[0]);

	  i += 2;
	}
      else if (strcmp("--iters", argv[i]) == 0)
	{
	  // Specify the number of iterations we'll run.

	  setIterations_m = true;
	  int tmp;
	  convertOK = intArgument(argc, argv, i+1, tmp);
	  if (!convertOK) inputError("Bad argument to --iters", argv[0]);
          iters_m = tmp;
	  i += 2;
	}
      else if (strcmp("--sim-params", argv[i]) == 0)
	{
	  // Specify the starting value, number of decades,
	  // and number of points.

	  setParams_m = true;
	  if (i+3 >= argc) inputError("--sim-params requires three arguments",
	                              argv[0]);
  	  
  	  convertOK = intArgument(argc, argv, i+1, start_m);
	  if (!convertOK) inputError("Bad first argument to --sim-params", 
	                             argv[0]);
  	  
  	  convertOK = intArgument(argc, argv, i+2, decades_m);
	  if (!convertOK) inputError("Bad second argument to --sim-params", 
	                             argv[0]);

  	  convertOK = intArgument(argc, argv, i+3, points_m);
	  if (!convertOK) inputError("Bad third argument to --sim-params", 
	                             argv[0]);
  	  
	  i += 4;
	}
      else if (strcmp("--no-diags", argv[i]) == 0)
	{
	  diags_m = false;
	  ++i;
	}
      else if (strcmp("--no-print", argv[i]) == 0)
	{
	  print_m = false;
	  ++i;
	}
      else if (strcmp("--report-time", argv[i]) == 0)
	{
	  report_time_m = true;
	  ++i;
	}
      else if (strcmp("--num-patches", argv[i]) == 0)
	{
	  setNumPatches_m = true;
  	  convertOK = intArgument(argc, argv, i+1, numPatches_m);
	  if (!convertOK) inputError("Bad argument to --num-patches", 
	                             argv[0]);
	  i += 2;
	}
      else if (strcmp("--benchmark-help", argv[i]) == 0)
	{
	  usage(argv[0]);
	  exit(0);
	}
      else if (strcmp("--list-impls", argv[i]) == 0)
	{
	  listimpls_m = true;
	  ++i;
	}
      else
	{
	  Pooma::perr << "Unknown option: " 
	              << argv[i] << "." << std::endl;
	  usage(argv[0]);
	  exit(1);
	}
      PAssert(convertOK);
    }
}


//-----------------------------------------------------------------------------
// Prints usage information.
//-----------------------------------------------------------------------------

void Benchmark::usage(const char *name)
{
  Pooma::perr.setPrefix();
  Pooma::perr << name << " options:\n"
    << "--benchmark-help...................print this message.\n"
    << "--sim-params N D P.................run a series of cases\n"
    << "                                   starting with problem size N\n"
    << "                                   through size = N * 10^D\n"
    << "                                   with P points per decade.\n"
    << "--list-impls.......................prints enumerated list of\n"
    << "                                   available implementations.\n"
    << "--run-impls N1, N2, etc............run the series of implementations\n"
    << "                                   N1, N2, etc., where N1, N2, ...\n"
    << "                                   are the numbers listed by\n"
    << "                                   --run-impls\n"
    << "--run-vars V1, V2, etc.............run the series of variations\n"
    << "                                   V1, V2, etc.\n"
    << "--no-print.........................don't print anything (useful if\n"
    << "                                   profiling using an external tool).\n"
    << "--no-diags.........................suppress diagnostic output.\n"
    << "--report-time......................print time, not Mflops.\n"
    << "--iters N..........................run benchmark for N iterations\n"
    << "                                   (no effect if using SGI timers).\n"
    << "--samples N........................repeat runs N time.\n"
    << "--num-patches N....................run UMP cases with N patches in\n"
    << "                                   each dimension."
    << std::endl;
}

//-----------------------------------------------------------------------------
// Prints an error message and usage information and then exits.
//-----------------------------------------------------------------------------

void Benchmark::inputError(const char *msg, const char *progname, int errcode)
{
  Pooma::perr << msg << std::endl;
  usage(progname);
  exit(errcode);
}

//-----------------------------------------------------------------------------
// Prints a list of available implementations and exits
//-----------------------------------------------------------------------------

void Benchmark::printImplementations()
{
  Pooma::perr.setPrefix();
  if (variation_m != std::string(""))
    {
      Pooma::perr << "Variation " << variation_m << std::endl;
      Pooma::perr.setPrefix("  ");
    }
  int size = impls_m.size();
  for (int i = 0; i < size; ++i)
    {
      Pooma::perr << i << "\t  " << impls_m[i]->type();
      if (strlen(impls_m[i]->qualification()) != 0)
        Pooma::perr << " " << impls_m[i]->qualification();
      Pooma::perr << std::endl;
    }
}

//-----------------------------------------------------------------------------
// Delete the implementations.
//-----------------------------------------------------------------------------

Benchmark::~Benchmark()
{
  std::vector<Implementation *>::iterator i = impls_m.begin();
  while (i != impls_m.end()) delete *i++;
}


//-----------------------------------------------------------------------------
// Adds an Implementation to the ones we are supposed to run. The Benchmark
// class takes responsibility for deleting the implementation.
//-----------------------------------------------------------------------------

void Benchmark::addImplementation(Implementation *impl)
{
  impls_m.push_back(impl);
}


//-----------------------------------------------------------------------------
// Sets the number of iterations assuming the value has not been previously set
// by the command line.
//-----------------------------------------------------------------------------

void Benchmark::setIterations(long iters)
{
  if (!setIterations_m)
    iters_m = iters;
}


//-----------------------------------------------------------------------------
// Sets the default number of patches assuming the value has not been
//  previously set by the command line.
//-----------------------------------------------------------------------------

void Benchmark::setNumPatches(int numPatches)
{
  if (!setNumPatches_m)
    numPatches_m = numPatches;
}


//-----------------------------------------------------------------------------
// Sets the sampling parameters assuming these values have not been previously 
// set by the command line.
//-----------------------------------------------------------------------------

void Benchmark::setSamplingParameters(int startVal, int numDecades, 
  int numPoints)
{
  if (!setParams_m)
    {
      start_m = startVal;
      decades_m = numDecades;
      points_m = numPoints;
    }
}


//-----------------------------------------------------------------------------
// Initializes the benchmarking machinary, runs the benchmark for some or all
// of the implementations, and optionally prints results.
//-----------------------------------------------------------------------------

void Benchmark::run()
{
  // If the user just wants a list of available implementations, we
  // print it and exit:
  
  if (listimpls_m) 
    {
      printImplementations();
      return;
    }
    
  // If the user has specified a variation list and we have a name but
  // are not on it, simply return.

  if (variation_m != std::string("") && varsToRun_m.size() != 0 &&
      std::find(varsToRun_m.begin(), varsToRun_m.end(), variation_m) ==
      varsToRun_m.end())
    return;

  if (print_m && diags_m)
    {
      int len = name_m.length() + 10;
      *inform_m << "\n" << name_m.c_str() << " Benchmark";
      if (variation_m != std::string(""))
	{
	  len += 14 + variation_m.length();
          *inform_m << ", variation \"" << variation_m.c_str() << "\"";
        }
      
      std::string dash(len, '-' );
      *inform_m << "\n" << dash.c_str() << std::endl;
    }

#if POOMA_EXCEPTIONS
  try 
    {
      getReadyToRun();
      
      runIt();
    
      printResults();
    }
  catch (Pooma::Assertion &a)
    {
      a.print(std::cerr);
      std::cerr << std::endl;
    }
  catch (const char *msg)
    {
      std::cerr << "Caught exception: " << msg << std::endl;
    }
  catch(...)
    {
      std::cerr << "Unknown exception." << std::endl;
    }
#else
  getReadyToRun();
  
  runIt();
  
  printResults();
#endif
}


//-----------------------------------------------------------------------------
// Reinitializes result data in preparation for another run.
//-----------------------------------------------------------------------------

void Benchmark::getReadyToRun()
{
  implsRan_m.clear();
  times_m.clear();
}


//-----------------------------------------------------------------------------
// Decides which implementations to run and runs them.
//-----------------------------------------------------------------------------

void Benchmark::runIt()
{
  // Loop through all of the implementations we have.

  for (int i = 0; i < impls_m.size(); ++i)
    {
      // If we've specified the names of the implementations that we 
      // are going to run, only run those.

      bool go1 = implsToRun_m.size() == 0;      
      bool go2 = std::find(implsToRun_m.begin(), implsToRun_m.end(), i) 
                           != implsToRun_m.end();
                           
      if (go1 || go2)
	{
	  // We're going to run this implementation.
          Implementation *impl = impls_m[i];
          for (int j = 0; j < samples_m; j++)
            {
	      runImplementation(impl, j);
            }
          implsRan_m.push_back(impl);
	}
    }
}


//-----------------------------------------------------------------------------
// Runs a specific implementation, tests the results, and stores timing and
// validity data.
//-----------------------------------------------------------------------------

void Benchmark::runImplementation(Implementation *impl, int sample)
{
  // Create space to store the times and test results.

  int numRuns = numPoints();
  std::vector<double> times(numRuns);
  std::vector<bool> valid(numRuns, true);

  // If we're printing, let the user know what's happening.

  if (print_m && diags_m)
    {
      *inform_m << "Running sample #" << sample+1 << " for " << impl->type();
      if (strlen(impl->qualification()) != 0)
        *inform_m << " " << impl->qualification();
      *inform_m << " Implementation:" << std::endl;
    }

  for (int i = 0; i < numRuns; i++)
    {
      // Compute the size of the trial.

      int npts = trialPoints(i);

      // If we're printing, print out the sample size.
      
      if (print_m && diags_m)
        *inform_m << "  N = " << npts << "..." << std::endl;

      // Initialize. 

      impl->initialize(npts);
        
      // Run the benchmark.  Compute overhead as you go.
      
      double start, stop, total = 0.0;
      double subTotal = 0.0;
      long iters = 0;

      // If we have the high-speed timers, just run one iteration.

      if (Pooma::Clock::highSpeed && !(impl->forceMultipleIterations())) 
        {
          iters = 1;
          if (impl->internalClockCalls()) 
            {
              // Run benchmark once and get timing result.
              
              impl->run();
	      total = impl->internalTimingResult();
            }
          else 
            {
              // Run and time benchmark once.
          
              start = Pooma::Clock::value();
              impl->run();
              stop = Pooma::Clock::value();
              total = stop - start;
          
              // Subtract out the looping overhead.
          
              start = Pooma::Clock::value();
              impl->runSetup();
              stop = Pooma::Clock::value();
              total -= (stop - start);
            }
        }
      else 
        { 
          if (impl->internalClockCalls()) 
            {
              for (iters = 0; iters < iters_m; ++iters) 
                {
                  // Run benchmark and accumulate timing result
                  
                  impl->run();
                  total += impl->internalTimingResult();
                }
            }
          else 
            {
              // Compute elapsed time here rather than accumulated time
              // to avoid getting a result of zero for low-resolution timers.
              
              long iterMax = iters_m;
              while (true)
                {
                  start = Pooma::Clock::value();
                  for (iters = 0; iters < iterMax; ++iters) 
                    {
                      // Run benchmark
            
                      impl->run();
	            }
              
                  // Compute elapsed total time.
	      
	          stop = Pooma::Clock::value();
	          total = stop - start;
	          
	          // Did we see any go by? (Can happen with short kernels.)
	          
	          if (total != 0.0)
	            break;
	          else
	            {
	              // If we're running with multiple contexts, we're hosed.
	              // Otherwise, increase the number of iterations.
	              
	              PAssert(Pooma::contexts() == 1);
	              iterMax *= 10;
	            }
	        }
          
              // Make an attempt to subtract out the looping overhead.
              
              start = Pooma::Clock::value();
              for (iters = 0; iters < iterMax; ++iters) 
                {
                  // Run setup
                  
                  impl->runSetup();
	        }
          
              // Compute elapsed overhead time and subtract from total.
	      
	      stop = Pooma::Clock::value();
	      subTotal = stop - start;
              total -= subTotal;
            }
        }

      // Compute run time per iteration in seconds.

      double timeper = total / double(iters);

      // Either store the running time or the MOps.

      if (report_time_m)
	times[i] = total;
      else
	times[i] = impl->opCount() / timeper / 1.0e6;

      // If we're testing results and we're printing, do this now.

      if (test_m && print_m && diags_m)
        *inform_m << "    Correctness test value for N = "
	  << npts << " is " << impl->resultCheck() << "." << std::endl;
    }

  // Store all of the test info.

  times_m.push_back(times);
}


//-----------------------------------------------------------------------------
// Prints results for all of the implementations that ran.
//-----------------------------------------------------------------------------

void Benchmark::printResults()
{
  // If we're not supposed to print, don't.
  
  if (!print_m)
    return;

  // Print out header line.

  std::vector<Implementation *>::iterator i = implsRan_m.begin();
  *inform_m << "        ";

  while (i != implsRan_m.end())
    {
      const char *qualCStr = (*i)->qualification();
      PInsist(strlen(qualCStr) < 13, 
	      "Benchmark::printResults: qualification string too long!\n"
              "                         Must be 12 or fewer characters.\n");
      std::string qu(qualCStr);
      std::string in;
      if (qu.size() > 0)
        in = (*i)->type();
      int n1 = (14 - in.size()) / 2, n2 = 14 - in.size() - n1;
      std::string sp1(n1, ' ' ), sp2(n2, ' ' );
      *inform_m << sp1.c_str() << in.c_str() << sp2.c_str();
      ++i;
    }
  *inform_m << '\n';

  i = implsRan_m.begin();
  *inform_m << "N       ";
  while (i != implsRan_m.end())
    {
      const char *qualCStr = (*i)->qualification();
      PInsist(strlen(qualCStr) < 13, 
	      "Benchmark::printResults: qualification string too long!\n"
              "                         Must be 12 of fewer characters.\n");
      std::string in(qualCStr);
      if (in.size() == 0)
        in = (*i)->type();
      int n1 = (14 - in.size()) / 2, n2 = 14 - in.size() - n1;
      std::string sp1(n1, ' ' ), sp2(n2, ' ' );
      *inform_m << sp1.c_str() << in.c_str() << sp2.c_str();
      ++i;
    }
  *inform_m << '\n';

  // Print out the values, one for each implementation.
      
  for (int j = 0; j < numPoints(); j++)
    {
      std::vector<double> pt(implsRan_m.size());
      for (int ki = 0; ki < implsRan_m.size(); ki++)
	pt[ki] = 0.0;
      for (int js = 0; js < numSamples(); js++)
        {
          for (int k = 0; k < implsRan_m.size(); k++)
            {
	      if (times_m[k * numSamples() + js][j] > pt[k])
		pt[k] = times_m[k * numSamples() + js][j];
	    }
	}
#if POOMA_MISSING_IOMANIPS
      *inform_m << std::setw(7) << trialPoints(j);
      for (int k = 0; k < implsRan_m.size(); k++)
	{
	  *inform_m << "    " << std::setw(6) << std::setprecision(4)
	            << pt[k] << "    ";
	}
#else
      *inform_m << std::setw(7) << std::left << trialPoints(j);
      for (int k = 0; k < implsRan_m.size(); k++)
	{
	  *inform_m << "    " << std::setw(6) << std::right << std::fixed 
	            << std::setprecision(2) << pt[k] << "    ";
	}
#endif
      *inform_m << '\n';
    }
  *inform_m << std::flush;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Benchmark.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.47 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
