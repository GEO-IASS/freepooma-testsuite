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

#ifndef POOMA_THREADS_POOMACSEM_H
#define POOMA_THREADS_POOMACSEM_H

//-----------------------------------------------------------------------------
// Class:
//   Pooma::CountingSemaphore
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview: 
// Pooma::CountingSemaphore is a counting semaphore object that can be used for
// blocking the parse thread until a condition has been met.  (Typically that
// a given set of iterates has run.)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include "Utilities/PAssert.h"
#include "Utilities/PurifyConstructors.h"

//-----------------------------------------------------------------------------
//
// CountingSemaphore:
//
// Pooma uses CountingSemaphore objects to implement a counting semaphore
// capability.
// A counting semaphore is initialized with a limit, and then various
// entities can increment this count until a counter hits the limit.  A
// thread can then "wait" on the semaphore to reach the limit, typically to
// implement a "join" operation.
//
// The Pooma CountingSemaphore wraps the Smarts counting semaphore if we're
// using Smarts, because the Smarts counting semaphore is a thread-specific
// concept that is not aware of iterates.  We need to tell the scheduler to
// release iterates before waiting, otherwise we get a deadlock because the
// parse thread is waiting for iterates to increment the semaphore and the
// iterates may be waiting on in holding queues waiting for the parse thread
// to hit the next generation.
//
// The SerialAsync version also need iterate-aware code.  When the semaphore
// wait() is called in the parse thread, we make calls to run iterates off
// the Smarts ready queues until the semaphore height has been reached.
// Also, the SerialAsync version makes call to the Cheetah poll() method
// if we're using messaging, because in the cross-context code, iterates
// can be waiting on messages.
//
// The POOMA wrapper assumes the following interface to the counting semaphore:
//
//   constructor:
//     CountingSemaphore() ... initialize the limit to zero.
//
//   methods:
//     void wait() ... wait for the semaphore to reach the limit
//     int count() ... return the current count
//     int height() ... return the limit "height"
//     void height(int d) ... change the limit "height" value
//     void raise_height(int d) ... increase the limit "height" value
//     void incr() ... increment the count by one
//     int operator+=(int d) ... add to the count, and return new count
//
// In the serial case, CountingSemaphore is a dummy class defined below.
//
//-----------------------------------------------------------------------------

#if POOMA_THREADS 

// NOTE: This is probably not correct if real PThreads are being used,
// unless we've adjusted the path or something to pick up a version of
// CSem.h that wraps a PThread semaphore.

#include "CSem.h"

namespace Pooma {

class CountingSemaphore
{
public:

  CountingSemaphore() : csem_m() { }

  // Causes the calling thread to wait until the internal counter has
  // reached the threshold.
  // For this serial version, it is an error if the limit has not
  // already been reached.
    
  void wait()
  {
    Pooma::scheduler().releaseIterates();
    csem_m.wait();
  }
  
  // Returns the value of the current internal counter.

  inline
  int count()
  {
    return csem_m.count();
  }

  // Returns the threshold(limit) of the counting semaphore

  int height()
  {
    return csem_m.height();
  }

  // Changes the threshold(limit) of the counting semaphore

  void height(int d)
  {
    csem_m.height(d);
  }

  // Adds specified to the current threshold(limit)

  void raise_height(int d)
  {
    csem_m.raise_height(d);
  }

  // Increments the threshold(limit) by 1

  void incr()
  {
    csem_m.incr();
  }

  // Overload pre-increment. We don't do post-increment
  // as it doesn't really make sense to make a copy of the
  // the semaphore. 
    
  CountingSemaphore &operator++()
  {
    incr();
    return *this;
  }
    
  // Adds specified to the current threshold(limit), same as raise_height()

  int operator+=(int d)
  {
    return (csem_m += d);
  }


private:

  Smarts::CSem csem_m;
};

}

#elif POOMA_SMARTS_SCHEDULER_SERIALASYNC

#include "Threads/IterateSchedulers/SerialAsync.h"

namespace Pooma {

/*------------------------------------------------------------------------
CLASS
	CSem_Serial_Async_version

	Counting semaphore that works with the Serial Async scheduler.

KEYWORDS
	Data-parallelism, Native-interface, Counting-semaphore.

DESCRIPTION

        SerialAsync needs a special counting semaphore since the context
        behaves differently than the multithreaded case.  When it waits,
        we need to run iterates off the queue until the count has been
        increased sufficiently.
-----------------------------------------------------------------------------*/

class CountingSemaphore 
{
public:
  inline
  CountingSemaphore(int limit=0)
  {
    height_m = limit;
    count_m = 0;
  }

  inline
  ~CountingSemaphore()
  {
    PAssert(count_m == height_m);
  }

  // The SerialAsync version of wait() just calls Pooma::poll() which
  // pushes both messages and iterates through queues.

  void wait()
  {
    PAssert(count_m <= height_m);
    while (count_m < height_m)
    {
      Pooma::poll();
    }
  }

  ///////////////////////////
  // Returns the value of the current internal counter.
  //
  inline
  int count()
  {
    return count_m;
  }

  ///////////////////////////
  // Returns the threshold(limit) of the counting semaphore
  //
  inline
  int height()
  {
    return height_m;
  }

  ///////////////////////////
  // Changes the threshold(limit) of the counting semaphore
  //
  inline
  void height(int d)
  {
    height_m = d;
  }

  ///////////////////////////
  // Adds specified to the current threshold(limit)
  //
  inline
  void raise_height(int d)
  {
    height_m += d;
  }

  ///////////////////////////
  // Increments the threshold(limit) by 1
  //
  inline
  void incr()
  {
    count_m += 1;
  }

  inline
  void operator++(void)
  {
    incr();
  }

private:
  // current count
  int count_m;

  // height of the semaphore
  int height_m;
};

}

#else

// Dummy CSem class

namespace Pooma {

class CountingSemaphore
{
public:

#if defined(NOPAssert)

  // Constructors and copy assignment.
  // Hopefully these will nix warnings from Purify.

  CountingSemaphore() { }
  CountingSemaphore(const CountingSemaphore &) { }
  CountingSemaphore &operator=(const CountingSemaphore &) { return *this; }

  // CSem interface:
  // Optimized version is a no-op.

  void wait() const { }
  int count() const { return 0; }
  int height() const { return 0; }
  void height(int) { }
  void raise_height(int) { }
  void incr() { }
  CountingSemaphore &operator++() { incr(); return *this; }
  int operator+=(int) { return 0; }

#else

  // Testing version.  Stores a limit and count, and checks when
  // wait() is called that the limit has been reached.
  // that are not in sync.

  // Constructors and copy assignment.
  // Hopefully these will nix warnings from Purify.

  CountingSemaphore() 
    : count_m(0), height_m(0)
  { }

  CountingSemaphore(const CountingSemaphore &model) 
    : count_m(model.count_m), height_m(model.height_m)
  { }

  CountingSemaphore &operator=(const CountingSemaphore &model) 
  { 
    count_m = model.count_m;
    height_m = model.height_m;
    return *this; 
  }

  // Causes the calling thread to wait until the internal counter has
  // reached the threshold.
  // For this serial version, it is an error if the limit has not
  // already been reached.
    
  void wait() const { PAssert(count_m == height_m); }
  
  // Returns the value of the current internal counter.

  int count() const { return count_m; }

  // Returns the threshold(limit) of the counting semaphore

  int height() const { return height_m; }

  // Changes the threshold(limit) of the counting semaphore

  void height(int d) { height_m = d; }

  // Adds specified to the current threshold(limit)

  void raise_height(int d) { height_m += d; }

  // Increments the threshold(limit) by 1

  void incr() { count_m += 1; }

  // Overload pre-increment. We don't do post-increment
  // as it doesn't really make sense to make a copy of the
  // the semaphore. 
    
  CountingSemaphore &operator++() { incr(); return *this; }
    
  // Adds specified to the current threshold(limit), same as raise_height()

  int operator+=(int d) { count_m += d; return count_m; }

  private:
  // Current count

  int count_m;

  // Current limit height

  int height_m;

#endif

};


}

#endif

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_THREADS_POOMACSEM_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PoomaCSem.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:07 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
