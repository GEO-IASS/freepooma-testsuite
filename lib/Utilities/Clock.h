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
 * Running timer.
 */

#ifndef POOMA_UTILITIES_CLOCK_H
#define POOMA_UTILITIES_CLOCK_H

#if POOMA_CLOCK_USES_GETTIMEOFDAY
#include <sys/time.h> // gettimeofday for faster timer
#endif

#include <time.h>

namespace Pooma {

/**
 * Clock provides a running timer, utilizing high-speed SGI timers if 
 * available.
 */

class Clock 
{
public:

  //---------------------------------------------------------------------
  // Set a static const that tells whether or not this class is utilizing
  // high-speed timers.
  
#if defined(CLOCK_SGI_CYCLE)
  enum { highSpeed = true };
#else
  enum { highSpeed = false };
#endif  

  //---------------------------------------------------------------------
  // Return the current value of the timer [sec].
  //
  // NOTE: some of these timers return CPU time and some return
  // "real" time. You need to know which you're using to understand
  // your timing results, particularly in parallel. 
  
  inline static double value()
  {
#if POOMA_CLOCK_USES_CLOCK_GETTIME

    // If clock_gettime is available, it should be used
    // as it has up-to nanosecond resolution.

    timespec ts;

# if defined(CLOCK_SGI_CYCLE)

    // This provides direct access to SGI's hardware performance
    // registers.

    clock_gettime(CLOCK_SGI_CYCLE, &ts);     // CPU-time???

# else

    clock_gettime(CLOCK_REALTIME, &ts);      // clock time

# endif

    return ts.tv_sec + 1e-9 * ts.tv_nsec;

#else
# if POOMA_CLOCK_USES_GETTIMEOFDAY

    // gettimeofday has up to microsecond resolution.

    timeval tv;
    gettimeofday(&tv, 0);                    // clock time
    return tv.tv_sec + 1e-6 * tv.tv_usec;

# else

    // Don't believe CLOCKS_PER_SEC - on many systems (i.e. linux) it
    // is set to 1000000, but actually only has a resolution
    // determined by the timer interrupt, which is about ever 10 ms
    // under Linux.

    return double(clock()) / CLOCKS_PER_SEC; // CPU-time

# endif
#endif
  }
};

} // namespace Pooma


//////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_CLOCK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Clock.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
