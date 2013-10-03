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
 * Statistics keeps statistics about a given POOMA job, and can report on
 * a summary of these statistics when asked.
 */

#ifndef POOMA_UTILITIES_STATISTICS_H
#define POOMA_UTILITIES_STATISTICS_H

//-----------------------------------------------------------------------------
// Classes: 
//   Statistics
//   StatisticsData: helper class
//-----------------------------------------------------------------------------

#include <string>
#include <vector>
#include <utility>

class Inform;


namespace Pooma {

class StatisticsData
{
  friend class Statistics;
  
public:

  const std::string &description() const { return data_m.first; }
  long value() const { return data_m.second; }
  
  void increment(long val = 1) { data_m.second += val; }
  
private:

  // Make sure riff-raff can't instantiate these puppies.

  StatisticsData(const char *description, long initialValue = 0)
  : data_m(description, initialValue)
  { }
  
  ~StatisticsData() { }

  std::pair<std::string, long> data_m;
};

/**
 * Statistics keeps statistics about a given POOMA job, and can report on
 * a summary of these statistics when asked.
 *
 * This interface is extensible ... you can add new types of statistics
 * by calling 'add(description, initval)' with a string description
 * of the stat, and the initial value it should have.
 *
 * The StatisticsData class is a helper that contains the description and
 * data.
 */

class Statistics {
private:

  //---------------------------------------------------------------------------
  // A NO-OP filter function. Up here because CW requires it.
  
  static long defaultFilter(long val);

public:
  
  //---------------------------------------------------------------------------
  // Constructor: initialize statistics.

  Statistics();

  //---------------------------------------------------------------------------
  // Destructor

  ~Statistics();

  //---------------------------------------------------------------------------
  // Print out the statistics to the given Inform object

  void print(Inform &, long (*filter)(long) = defaultFilter);

  //---------------------------------------------------------------------------
  // Add a statistics object to our list of stats ... return an integer
  // which is the index of the stat in our list, which can be used with
  // the 'incrementStat' and 'decrementStat' methods to change that statistic.

  StatisticsData *add(const char *description, long initval = 0) 
  {
    StatisticsData *sd = new StatisticsData(description, initval);
    statList_m.push_back(sd);
    return sd;
  }

private:

  //---------------------------------------------------------------------------
  // A vector of statistics data objects, which will be used to print
  // out the results at the end.

  std::vector<StatisticsData *> statList_m;
};

} // namespace Pooma

#endif // POOMA_UTILITIES_STATISTICS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Statistics.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
