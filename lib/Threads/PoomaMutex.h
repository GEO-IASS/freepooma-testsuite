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

#ifndef POOMA_THREADS_POOMAMUTEX_H
#define POOMA_THREADS_POOMAMUTEX_H

//-----------------------------------------------------------------------------
// Class:
//   Pooma::Mutex_t
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview: 
// Pooma::Mutex_t is a typedef for the appropriate type of mutex
// object to be used when running in parallel. 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Pooma uses Mutex_t objects to protect various pieces of data from
// simultaneous access by multiple threads. The actual type of the
// mutex is determined here. Usually this will just be a smarts
// mutex. However, if we're using a different threading system, we may 
// need to include a wrapper for that system's mutex. Pooma assumes
// that the mutex has lock() and unlock() member functions, so any
// type with this interface can be typedef'd to Pooma::Mutex_t below.
//
// In the serial case, Mutex_t is a typedef for the dummy Mutex class
// defined below.
//
//-----------------------------------------------------------------------------

#if POOMA_THREADS && !POOMA_SMARTS_SCHEDULER_SERIALASYNC

// NOTE: This is probably not correct if real PThreads are being used,
// unless we've adjusted the path or something to pick up a version of
// Mutex.h that wraps a PThread mutex.

#include "Mutex.h"

// NOTE: This is probably not correct if real PThreads are being used,
// unless we've adjusted the path or something to pick up a version of
// Mutex.h that wraps a PThread mutex.

namespace Pooma {
  typedef Smarts::Mutex Mutex_t;
}

#else

// Dummy Mutex class

namespace Pooma {

  class DummyMutex
  {
  public:

#if defined(NOPAssert)

    // Constructors and copy assignment.
    // Hopefully these will nix warnings from Purify.

    inline DummyMutex() { }
    inline DummyMutex(const DummyMutex &) { }
    inline DummyMutex &operator=(const DummyMutex &) { return *this; }

    // Mutex interface:
    // Optimized version is a noop.

    inline void lock() { }
    inline void unlock() { }

#else

    // Testing version. Checks for lock/unlock calls
    // that are not in sync.

    // Constructors and copy assignment.
    // Hopefully these will nix warnings from Purify.

    inline DummyMutex() 
      : locked_m(false) 
    { }

    inline DummyMutex(const DummyMutex &model) 
      : locked_m(model.locked_m) 
    { }

    inline DummyMutex &operator=(const DummyMutex &model) 
    { 
      locked_m = model.locked_m;
      return *this; 
    }

    inline void lock() 
    {
      PAssert(!locked_m);
      locked_m = true;
    }

    inline void unlock() 
    {
      PAssert(locked_m);
      locked_m = false;
    }

  private:

    bool locked_m;

#endif

  };

  typedef DummyMutex Mutex_t;

}

#endif

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_THREADS_POOMAMUTEX_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PoomaMutex.h,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:17:07 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
