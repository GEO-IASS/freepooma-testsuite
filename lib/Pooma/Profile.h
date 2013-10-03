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

#ifndef POOMA_POOMA_PROFILE_H
#define POOMA_POOMA_PROFILE_H

/** @file
 * @ingroup Pooma
 * @brief
 * TAU Profiling. In the absence of profiling or tracing, don't include TAU 
 * headers.  Otherwise, include the standard TAU profiling header.
 *
 *   Wrapper around the decision on whether to include the full Tau profiling
 * headers, or just define empty macros.  This will include the full Tau
 * files if POOMA_PROFILE or POOMA_TRACE are set to non-zero values.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"


#if POOMA_PROFILE || POOMA_TRACE

// if we're using profiling, include the standard Tau header
# include "Profile/Profiler.h"

#else 

// In the absence of profiling, define the macros as null
# define TYPE_STRING(profileString, str)
# define PROFILED_BLOCK(name, type)

# define TAU_TYPE_STRING(profileString, str)
# define TAU_PROFILE(name, type, group)
# define TAU_PROFILE_TIMER(var, name, type, group)
# define TAU_PROFILE_START(var)
# define TAU_PROFILE_STOP(var)
# define TAU_PROFILE_STMT(stmt)
# define TAU_PROFILE_EXIT(msg)
# define TAU_PROFILE_INIT(argc, argv)
# define TAU_PROFILE_SET_NODE(node)
# define TAU_PROFILE_SET_CONTEXT(context)
# define TAU_PROFILE_CALLSTACK()

# define TAU_REGISTER_EVENT(event, name)
# define TAU_EVENT(event, data)
# define TAU_EVENT_DISABLE_MIN(event)
# define TAU_EVENT_DISABLE_MAX(event)
# define TAU_EVENT_DISABLE_MEAN(event)
# define TAU_EVENT_DISABLE_STDDEV(event)
# define TAU_STORE_ALL_EVENTS
# define TAU_REPORT_STATISTICS()
# define TAU_REPORT_THREAD_STATISTICS()
# define TAU_REGISTER_THREAD()

# define CT(obj)

# define TAU_TRACE_SENDMSG(type, destination, length)
# define TAU_TRACE_RECVMSG(type, source, length)

#endif // POOMA_PROFILE || POOMA_TRACE


#endif // POOMA_UTILITIES_POOMAPROFILE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Profile.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:17:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
