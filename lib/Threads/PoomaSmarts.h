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

#ifndef POOMA_THREADS_POOMA_SMARTS_H
#define POOMA_THREADS_POOMA_SMARTS_H

/** @file
 * @ingroup Threads
 * @brief
 * The POOMA wrapper around defines, includes, and typedefs for the Smarts
 * run-time evaluation system.
 *
 * Based on the settings of POOMA_THREADS and
 * the selected scheduler, define several typedefs and include the necessary
 * files.  If we're compiling only for serial runs, use a stub version of
 * the Smarts interface instead.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/Scheduler.h"

// Given the scheduler policy tag, we can now define the types of the
// scheduler, dataobject, and iterate:

namespace Pooma {

static const char schedulerVersion[] = POOMA_SCHEDULER_NAME;

typedef Smarts::IterateScheduler<SmartsTag_t> Scheduler_t;
typedef Smarts::DataObject<SmartsTag_t>       DataObject_t;
typedef Smarts::Iterate<SmartsTag_t>          Iterate_t;

typedef Smarts::Runnable                      Runnable_t;
inline void addRunnable(Runnable_t *runnable)
{
  Smarts::add(runnable);
}

} // namespace Pooma


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_THREADS_POOMA_SMARTS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PoomaSmarts.h,v $   $Author: richard $
// $Revision: 1.18 $   $Date: 2004/11/01 18:17:07 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
