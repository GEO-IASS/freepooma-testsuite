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
// Class:
// SmartsTag_t
//-----------------------------------------------------------------------------

#ifndef POOMA_THREADS_SCHEDULER_H
#define POOMA_THREADS_SCHEDULER_H

/** @file
 * @ingroup Threads
 * @brief
 * Scheduler multiplexing based on configuration.
 *
 * This file exist to wrap the correct includes from Smarts based on the
 * scheduler that we've selected.  If we're running in serial then we include
 * the a stub file.  This file defines a single typedef: SmartsTag_t, a policy
 * tag which is used to select the appropriate smarts data object etc.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"

//-----------------------------------------------------------------------------
// First include the appropriate Smarts scheduler file and define 
// the scheduler policy tag and identification string appropriate
// to the particular choice of scheduler.
//-----------------------------------------------------------------------------

#if POOMA_THREADS

// Set up to use the proper scheduler.  One of the following blocks will
// end up being included.

# if POOMA_SMARTS_SCHEDULER_ASYNC

#  include "Smarts.h"
#  include "IterateSchedulers/FastAsyncScheduler.h"

namespace Pooma {
  typedef Smarts::FastAsync SmartsTag_t;
}

# elif POOMA_SMARTS_SCHEDULER_MCVE_MULTIQ

#   include "Smarts.h"
#   include "IterateSchedulers/MCVE_MultiQ.h"

namespace Pooma {
  typedef Smarts::MCVE_MultiQ SmartsTag_t;
}

# else

#  error "You have not selected a scheduler"

# endif

#else

# if POOMA_SMARTS_SCHEDULER_SERIALASYNC

#  include "Threads/IterateSchedulers/SerialAsync.h"

namespace Pooma {
  typedef Smarts::SerialAsync SmartsTag_t;
}

# else

// Set up a stub version of the Smarts interface for use in just serial runs

#  include "SmartsStubs.h"

namespace Pooma {
  typedef Smarts::Stub SmartsTag_t;
}

# endif

#endif

#endif     // POOMA_THREADS_SCHEDULER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Scheduler.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:07 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
