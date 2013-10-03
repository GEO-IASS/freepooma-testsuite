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
// NotifyEngineWrite
//
// function
// notifyEngineWrite()
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_NOTIFYENGINEWRITE_H
#define POOMA_ENGINE_NOTIFYENGINEWRITE_H

//////////////////////////////////////////////////////////////////////

/** @file 
 * @ingroup Engine
 * @brief
 * NotifyEngineWrite is a general wrapper class that is used to tell
 * an engine that we're going to write to it.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/WrappedInt.h"


/**
 * NotifyEngineWrite is a general wrapper class that is used to tell
 * an engine that we're going to write to it.
 *
 * Multipatch engines will want to fill their guard cells after being written
 * to.  This class allows us to notify the engines that need to be notified,
 * and do nothing to other engines.
 *
 * You must specialize NotifyEngineWrite for multiPatch engines.
 */

template<class Engine>
struct NotifyEngineWrite
{
  NotifyEngineWrite(){}
 ~NotifyEngineWrite(){}

 inline static void
  notify(const Engine &)
  {
    // Engines that are multipatch must specialize this functor.
    CTAssert(!(Engine::multiPatch));
  }
};

/// This helper function simplifies use of the NotifyEngineWrite functor.

template<class Engine>
inline
void notifyEngineWrite(const Engine &e)
{
  NotifyEngineWrite<Engine>::notify(e);
}

/// This function lets us skip the notification at compile time.  (If we're
/// actually reading from the engine, for example.)

template<class Engine>
inline
void notifyEngineWrite(const Engine &, const WrappedInt<false> &)
{
}

/// This function lets us skip the notification at compile time.  (If we're
/// actually reading from the engine, for example.)

template<class Engine>
inline
void notifyEngineWrite(const Engine &e, const WrappedInt<true> &)
{
  NotifyEngineWrite<Engine>::notify(e);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_NOTIFYENGINEWRITE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NotifyEngineWrite.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
