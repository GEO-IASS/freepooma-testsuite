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
// SendKernel
// ReceiveKernel
// SendReceive
//-----------------------------------------------------------------------------

#ifndef POOMA_CHEETAH_SEND_RECEIVE_H
#define POOMA_CHEETAH_SEND_RECEIVE_H

//////////////////////////////////////////////////////////////////////

/** @file 
 * @ingroup Tulip
 * @brief
 * SendKernel and ReceiveKernel are special iterates that interact with cheetah
 * to send and receive data that gets used in expressions.
 *
 * SendReceive is a wrapper class that contains a send() and receive()
 * function that encapsulate generating the necessary tag and launching the
 * SendKernel and ReceiveKernel iterates.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tulip/Messaging.h"
#include "Pooma/Pooma.h"
#include "Evaluator/InlineEvaluator.h"
#include "Evaluator/RequestLocks.h"
#include "Engine/DataObject.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

#if POOMA_CHEETAH

/** 
 * A SendIterate requests a read lock on a piece of data.  When that read lock
 * is granted, we call a cheetah matching handler to send the data to the
 * appropriate context.  We construct the SendIterate with a tag that is used
 * to match the appropriate ReceiveIterate on the remote context.
 */

template<class View>
class SendIterate
  : public Pooma::Iterate_t
{
public:
  SendIterate(const View &view, int toContext, int tag)
    : Pooma::Iterate_t(Pooma::scheduler()),
      toContext_m(toContext),
      tag_m(tag),
      view_m(view)
  {
    PAssert(toContext >= 0);
   
    hintAffinity(engineFunctor(view_m,
			       DataObjectRequest<BlockAffinity>()));

    // Priority interface was added to r2 version of serial async so that
    // message iterates would run before any other iterates.
    priority(-1);

    DataObjectRequest<WriteRequest> writeReq(*this);
    DataObjectRequest<ReadRequest> readReq(writeReq);
    engineFunctor(view_m, readReq);
  }

  virtual ~SendIterate()
  {
    DataObjectRequest<WriteRelease> writeReq;
    DataObjectRequest<ReadRelease> readReq(writeReq);
    engineFunctor(view_m, readReq);
  }

  virtual void run()
  {
    Pooma::remoteEngineHandler()->send(toContext_m, tag_m, view_m);
  }

private:

  // Context we're sending the data to.

  int toContext_m;

  // A tag used to match the sent data with the right receive.

  int tag_m;

  // The data we're sending (typically a view of an array).

  View view_m;
};


/**
 * ReceiveIterate requests a write lock on a piece of data.  When that lock
 * is granted, we register the data with the cheetah matching handler which
 * will fill the block when a message arrives.  The write lock is released
 * by the matching handler.
 */

template<class View, class IncomingView>
class ReceiveIterate
  : public Pooma::Iterate_t
{
public:

  typedef ReceiveIterate<View, IncomingView> This_t;

  ReceiveIterate(const View &view, int fromContext, int tag)
    : Pooma::Iterate_t(Pooma::scheduler()),
      fromContext_m(fromContext),
      tag_m(tag),
      view_m(view)
  {
    PAssert(fromContext >= 0);
    
    hintAffinity(engineFunctor(view,
			       DataObjectRequest<BlockAffinity>()));

    // Priority interface was added to r2 version of serial async so that
    // message iterates would run before any other iterates.
    priority(-1);

    DataObjectRequest<WriteRequest> writeReq(*this);
    engineFunctor(view, writeReq);

    Pooma::addIncomingMessage();
  }

  virtual ~ReceiveIterate()
  {
  }

  // If we're using cheetah, but don't support out-of-order execution, then
  // the run method of this iterate must block until the message has been
  // received.  Unlike typical iterates, the work implied by this iterate
  // isn't actually performed in the run method.  The run method merely
  // registers a method that gets handled by cheetah when the appropriate
  // message arrives.

  virtual void run()
  {
    Pooma::remoteEngineHandler()->request(fromContext_m, tag_m,
					  This_t::apply, view_m);
  }

private:

  static void apply(const View &viewLocal, IncomingView &viewMessage)
  {
    // For now, we just copy the message into the brick accepting the data.

    KernelEvaluator<InlineKernelTag>::evaluate(viewLocal, OpAssign(),
					       viewMessage);
    
    // Release the received block:
    DataObjectRequest<WriteRelease> writeReq;
    engineFunctor(viewLocal, writeReq);
    
    Pooma::gotIncomingMessage();
  }

  // Context we're sending the data to.

  int fromContext_m;

  // A tag used to match the sent data with the right send.

  int tag_m;

  // The place to put the data we're receiving (typically a view of the
  // engine).;

  View view_m;
};

/**
 * SendReceive contains two static functions, send(view, context) and
 * receive(view, context).  These functions encapsulate generating matching
 * tags for the send and receive and launching the iterates to perform the
 * send and receive.
 */

struct SendReceive
{
  template<class View>
  static
  void send(const View &view, int toContext)
  {
    int tag = Pooma::sendTag(toContext);
    Pooma::scheduler().handOff(new SendIterate<View>(view, toContext, tag));
  }
};

template<class IncomingView>
struct Receive
{
  template<class View>
  static
  void receive(const View &view, int fromContext)
  {
    PAssert(fromContext >= 0);
    int tag = Pooma::receiveTag(fromContext);
    Pooma::scheduler().handOff(new ReceiveIterate<View, IncomingView>
					(view, fromContext, tag));
  }
};


#elif POOMA_MPI


/** 
 * A SendIterate requests a read lock on a piece of data.  When that read lock
 * is granted, we call a cheetah matching handler to send the data to the
 * appropriate context.  We construct the SendIterate with a tag that is used
 * to match the appropriate ReceiveIterate on the remote context.
 */

template<class View>
class SendIterate
  : public Pooma::Iterate_t
{
public:
  SendIterate(const View &view, int toContext, int tag)
    : Pooma::Iterate_t(Pooma::scheduler()),
      toContext_m(toContext),
      tag_m(tag),
      view_m(view)
  {
    PAssert(toContext >= 0);
   
    hintAffinity(engineFunctor(view_m,
			       DataObjectRequest<BlockAffinity>()));

    // Priority interface was added to r2 version of serial async so that
    // message send iterates would run before any other iterates.
    priority(-1);

    DataObjectRequest<WriteRequest> writeReq(*this);
    DataObjectRequest<ReadRequest> readReq(writeReq);
    engineFunctor(view_m, readReq);
  }

  virtual void run()
  {
    typedef Cheetah::Serialize<Cheetah::CHEETAH, View> Serialize_t;

    // serialize and send buffer
    int length = Serialize_t::size(view_m);
    buffer_m = new char[length];
    Serialize_t::pack(view_m, buffer_m);
    MPI_Request *request = Smarts::SystemContext::getMPIRequest(this);
    int res = MPI_Isend(buffer_m, length, MPI_CHAR, toContext_m, tag_m,
			MPI_COMM_WORLD, request);
    PAssert(res == MPI_SUCCESS);

    // release locks
    DataObjectRequest<WriteRelease> writeReq;
    DataObjectRequest<ReadRelease> readReq(writeReq);
    engineFunctor(view_m, readReq);
  }

  virtual ~SendIterate()
  {
    // cleanup temporary objects.
    delete[] buffer_m;
  }

private:

  // Context we're sending the data to.

  int toContext_m;

  // A tag used to match the sent data with the right receive.

  int tag_m;

  // Communication buffer.

  char *buffer_m;

  // The data we're sending (typically a view of an array).

  View view_m;
};


/**
 * ReceiveIterate requests a write lock on a piece of data.  When that lock
 * is granted, we register the data with the cheetah matching handler which
 * will fill the block when a message arrives.  The write lock is released
 * by the matching handler.
 */

template<class View, class IncomingView>
class ReceiveIterate
  : public Pooma::Iterate_t
{
public:

  typedef ReceiveIterate<View, IncomingView> This_t;

  ReceiveIterate(const View &view, int fromContext, int tag)
    : Pooma::Iterate_t(Pooma::scheduler()),
      fromContext_m(fromContext),
      tag_m(tag), buffer_m(NULL),
      view_m(view)
  {
    PAssert(fromContext >= 0);
    
    hintAffinity(engineFunctor(view,
			       DataObjectRequest<BlockAffinity>()));

    // Priority interface was added to r2 version of serial async so that
    // message receive iterates would run after any other iterates.
    priority(-1);

    DataObjectRequest<WriteRequest> writeReq(*this);
    engineFunctor(view, writeReq);

    Pooma::addIncomingMessage();

    // pre-allocate incoming buffer and issue async receive
    // we may hog on requests here - so maybe we need to conditionalize
    // this a bit on request availability?
    if (Smarts::SystemContext::haveLotsOfMPIRequests()) {
      int length = Cheetah::Serialize<Cheetah::CHEETAH, View>::size(view_m);
      buffer_m = new char[length];
      MPI_Request *request = Smarts::SystemContext::getMPIRequest(this);
      int res = MPI_Irecv(buffer_m, length, MPI_CHAR, fromContext_m, tag_m,
			  MPI_COMM_WORLD, request);
      PAssert(res == MPI_SUCCESS);
    }
  }

  virtual void run()
  {
    // nothing - work is done in destructor, if we had enough requests free
    if (!buffer_m) {
      int length = Cheetah::Serialize<Cheetah::CHEETAH, View>::size(view_m);
      buffer_m = new char[length];
      MPI_Request *request = Smarts::SystemContext::getMPIRequest(this);
      int res = MPI_Irecv(buffer_m, length, MPI_CHAR, fromContext_m, tag_m,
			  MPI_COMM_WORLD, request);
      PAssert(res == MPI_SUCCESS);
    }
  }

  virtual ~ReceiveIterate()
  {
    typedef Cheetah::Serialize<Cheetah::CHEETAH, View> Serialize_t;

    // de-serialize into target view directly
    Serialize_t::unpack(view_m, buffer_m);

    // cleanup temporary objects
    delete[] buffer_m;

    // release locks
    DataObjectRequest<WriteRelease> writeReq;
    engineFunctor(view_m, writeReq);

    Pooma::gotIncomingMessage();
  }

private:

  // Context we're sending the data to.

  int fromContext_m;

  // A tag used to match the sent data with the right send.

  int tag_m;

  // Communication buffer.

  char *buffer_m;

  // The place to put the data we're receiving (typically a view of the
  // engine).;

  View view_m;
};

/**
 * SendReceive contains two static functions, send(view, context) and
 * receive(view, context).  These functions encapsulate generating matching
 * tags for the send and receive and launching the iterates to perform the
 * send and receive.
 */

struct SendReceive
{
  template<class View>
  static
  void send(const View &view, int toContext)
  {
    int tag = Pooma::sendTag(toContext);
    Pooma::scheduler().handOff(new SendIterate<View>(view, toContext, tag));
  }
};

template<class IncomingView>
struct Receive
{
  template<class View>
  static
  void receive(const View &view, int fromContext)
  {
    PAssert(fromContext >= 0);
    int tag = Pooma::receiveTag(fromContext);
    Pooma::scheduler().handOff(new ReceiveIterate<View, IncomingView>
					(view, fromContext, tag));
  }
};


#else // not POOMA_MESSAGING


/**
 * The non-cheetah versions of send and receive are empty and should never
 * actually be used, since a remote view should only happen when the data
 * lives on another context.
 */

struct SendReceive
{
  template<class View>
  static
  void send(const View &view, int toContext)
  {
    PAssert(false);
  }
};

template<class IncomingView>
struct Receive
{
  template<class View>
  static
  void receive(const View &view, int fromContext)
  {
    PAssert(false);
  }
};


#endif // not POOMA_MESSAGING

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CHEETAH_SEND_RECEIVE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SendReceive.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
