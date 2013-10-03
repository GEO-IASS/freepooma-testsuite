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
// Global TagGenerator
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Tulip/Messaging.h"

#include "Tulip/CollectFromContexts.h"
#include "Tulip/ReduceOverContexts.h"
#include "Tulip/RemoteProxy.h"
#include "Tulip/PatchSizeSyncer.h"
#include "Tulip/SendReceive.h"

#if POOMA_CHEETAH
int  ReduceOverContextsBase::tagBase_m = 0;
int  CollectFromContextsBase::tagBase_m = 0;
bool RemoteProxyBase::ready_m;
int  RemoteProxyBase::tag_m = 0;
#endif


//-----------------------------------------------------------------------------
// Tag generator creates a set of tags for global use in r2.  There is a
// separate tag number for each context we communicate with.
//-----------------------------------------------------------------------------

TagGenerator tagGenerator_g;

namespace Pooma {

//-----------------------------------------------------------------------------
// expectedMessages
//-----------------------------------------------------------------------------

int expectedMessages_g = 0;

#if POOMA_CHEETAH
Cheetah::MatchingHandler *collectionHandler_g   = 0;
Cheetah::MatchingHandler *indexHandler_g        = 0;
Cheetah::MatchingHandler *reductionHandler_g    = 0;
Cheetah::MatchingHandler *remoteEngineHandler_g = 0;
Cheetah::MatchingHandler *particleSwapHandler_g = 0;
#endif

void initializeCheetahHelpers(int contexts)
{
  tagGenerator_g     = TagGenerator(contexts);
  expectedMessages_g = 0;
#if POOMA_CHEETAH
  collectionHandler_g        = new Cheetah::MatchingHandler(*Pooma::controller());
  indexHandler_g             = new Cheetah::MatchingHandler(*Pooma::controller());
  reductionHandler_g         = new Cheetah::MatchingHandler(*Pooma::controller());
  remoteEngineHandler_g      = new Cheetah::MatchingHandler(*Pooma::controller());
  particleSwapHandler_g      = new Cheetah::MatchingHandler(*Pooma::controller());
#endif
}

void finalizeCheetahHelpers()
{
  PAssert(expectedMessages_g == 0);
#if POOMA_CHEETAH
  if (collectionHandler_g != 0)
    delete collectionHandler_g;
  if (indexHandler_g != 0)
    delete indexHandler_g;
  if (reductionHandler_g != 0)
    delete reductionHandler_g;
  if (remoteEngineHandler_g != 0)
    delete remoteEngineHandler_g;
  if (particleSwapHandler_g != 0)
    delete particleSwapHandler_g;
#endif
}

int sendTag(int context)
{
  return tagGenerator_g.send(context);
}

int receiveTag(int context)
{
  return tagGenerator_g.receive(context);
}

} // namespace Pooma

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Messaging.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
