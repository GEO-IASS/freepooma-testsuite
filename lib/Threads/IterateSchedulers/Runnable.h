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

#ifndef _Runnable_h_
#define _Runnable_h_

/** @file
 * @ingroup IterateSchedulers
 * @brief
 * Base class for a schedulable object or function to be executed
 * by the scheduler asynchronously.
 */ 

#include <string.h>

namespace Smarts {

/**
 * Runnable is the base class for system classes "Thread" and
 * "Iterate".  However, the user may define his/her own
 * sub-class. Any class derived from Runnable, is an object that
 * the scheduler understands and therefore is the mechanism to
 * have something executed in parallel by the scheduler on behalf
 * of the user.
 */

class Runnable
{
public:

  Runnable()
  {
    priority_m = 0;
  }

  /// The parameter to this constructor is the CPU id for
  /// hard affinity.

  Runnable(int)
  {
    priority_m = 0;
  }

  virtual ~Runnable() {}

  /// Accessor function to priority.

  inline int
  priority() { return priority_m; }

  /// Set priority of this runnable relative to other runnables 
  /// being scheduled.

  inline void
  priority(int _priority) { priority_m = _priority; }

  virtual void execute()
  {
    run();
  }

protected:
  virtual void run() {}

private:
  int priority_m;
};

typedef Runnable *RunnablePtr_t;

/// Schedulers need to implement this function to add
/// a runnable to the execution queue.

inline void add(RunnablePtr_t);

} // namespace Smarts

#endif
