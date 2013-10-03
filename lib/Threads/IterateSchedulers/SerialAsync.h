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
// This program was prepared by the Regents of the University of California
// at Los Alamos National Laboratory (the University) under Contract No. 
// W-7405-ENG-36 with the U.S. Department of Energy (DOE). The University has 
// certain rights in the program pursuant to the contract and the program 
// should not be copied or distributed outside your organization. All rights
// in the program are reserved by the DOE and the University. Neither the U.S.
// Government nor the University makes any warranty, express or implied, or
// assumes any liability or responsibility for the use of this software
//-----------------------------------------------------------------------------
// Class:
// IterateScheduler<SerialAsync>
// Iterate<SerialAsync>
// DataObject<SerialAsync>
//-----------------------------------------------------------------------------

#ifndef _SerialAsync_h_
#define _SerialAsync_h_

/** @file
 * @ingroup IterateSchedulers
 * @brief
 * Smarts classes for times when you want no threads but you do want
 * dataflow evaluation.
 *
 * SerialAsync IterateScheduler is a policy template to create a
 * dependence graphs and executes the graph respecting the
 * dependencies without using threads.
 * There is no (thread level) parallelism, but Iterates may be executed
 * out-of-order with respect to the program text. Also this scheduler is
 * used for message based parallelism in which case asyncronous execution
 * leads to reduced communication latencies.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include <list>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <stack>
#include "Pooma/Configuration.h"
#if POOMA_MPI
# include <mpi.h>
#endif
#include "Threads/IterateSchedulers/IterateScheduler.h"
#include "Threads/IterateSchedulers/Runnable.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//inline ostream& lock(ostream& s) { return s; }
//inline ostream& unlock(ostream& s) { return s; }

namespace Smarts {

/**
 * Tag class for specializing IterateScheduler, Iterate and DataObject.
 */

struct SerialAsync 
{
  enum Action { Read, Write };
};


/**
 * Iterate<SerialAsync> is used to implement the SerialAsync
 * scheduling policy.
 *
 * An Iterate is a non-blocking unit of concurrency that is used
 * to describe a chunk of work. It inherits from the Runnable
 * class and as all subclasses of Runnable, the user specializes
 * the run() method to specify the operation.
 * Iterate<SerialAsync> is a further specialization of the
 * Iterate class to use the SerialAsync Scheduling algorithm to
 * generate the data dependency graph for a data-driven
 * execution.
 */

template<>
class Iterate<SerialAsync> : public Runnable
{
  friend class IterateScheduler<SerialAsync>;
  friend class DataObject<SerialAsync>;

public:

  typedef DataObject<SerialAsync> DataObject_t;
  typedef IterateScheduler<SerialAsync> IterateScheduler_t;


  /// The Constructor for this class takes the IterateScheduler and a
  /// CPU affinity.  CPU affinity has a default value of -1 which means
  /// it may run on any CPU available.

  inline Iterate(IterateScheduler<SerialAsync> & s, int affinity=-1)
    : scheduler_m(s), notifications_m(1), generation_m(-1), togo_m(1)
  {}

  /// The dtor is virtual because the subclasses will need to add to it.

  virtual ~Iterate() {}

  /// The run method does the core work of the Iterate.
  /// It is supplied by the subclass.

  virtual void run() = 0;

  //@name Stubs for the affinities
  /// There is no such thing in serial.
  //@{

  inline int affinity() const {return 0;}

  inline int hintAffinity() const {return 0;}

  inline void affinity(int) {}

  inline void hintAffinity(int) {}

  //@}

  /// Notify is used to indicate to the Iterate that one of the data
  /// objects it had requested has been granted. To do this, we dec a
  /// dependence counter which, if equal to 0, the Iterate is ready for
  /// execution.

  void notify()
  {
    if (--notifications_m == 0)
      add(this);
  }

  /// How many notifications remain?

  int notifications() const { return notifications_m; }

  void addNotification() { notifications_m++; }

  int& generation() { return generation_m; }

  int& togo() { return togo_m; }

protected:

  /// What scheduler are we working with?
  IterateScheduler<SerialAsync> &scheduler_m;

  /// How many notifications should we receive before we can run?
  int notifications_m;

  /// Which generation we were issued in.
  int generation_m;

  /// How many times we need to go past a "did something" to be ready
  /// for destruction?
  int togo_m;

};


struct SystemContext
{
  void addNCpus(int) {}
  void wait() {}
  void concurrency(int){}
  int concurrency() { return 1; }
  void mustRunOn() {}

  // We have a separate message queue because they are
  // higher priority.
  typedef Iterate<SerialAsync> *IteratePtr_t;
  static std::list<RunnablePtr_t> workQueueMessages_m;
  static std::list<RunnablePtr_t> workQueue_m;
#if POOMA_MPI
  static const int max_requests = 1024;
  static MPI_Request requests_m[max_requests];
  static std::map<int, IteratePtr_t> allocated_requests_m;
  static std::set<int> free_requests_m;
#endif


#if POOMA_MPI

  /// Query, if we have lots of MPI_Request slots available

  static bool haveLotsOfMPIRequests()
  {
    return free_requests_m.size() > max_requests/2;
  }

  /// Get a MPI_Request slot, associated with an iterate

  static MPI_Request* getMPIRequest(IteratePtr_t p)
  {
    PInsist(!free_requests_m.empty(), "No free MPIRequest slots.");
    int i = *free_requests_m.begin();
    free_requests_m.erase(free_requests_m.begin());
    allocated_requests_m[i] = p;
    p->togo()++;
    return &requests_m[i];
  }

  static void releaseMPIRequest(int i)
  {
    IteratePtr_t p = allocated_requests_m[i];
    allocated_requests_m.erase(i);
    free_requests_m.insert(i);
    if (--(p->togo()) == 0)
      delete p;
  }

  static bool waitForSomeRequests(bool mayBlock)
  {
    if (allocated_requests_m.empty())
      return false;

    int last_used_request = allocated_requests_m.rbegin()->first;
    int finished[last_used_request+1];
    MPI_Status statuses[last_used_request+1];
    int nr_finished;
    int res;
    if (mayBlock)
      res = MPI_Waitsome(last_used_request+1, requests_m,
			 &nr_finished, finished, statuses);
    else
      res = MPI_Testsome(last_used_request+1, requests_m,
			 &nr_finished, finished, statuses);
    PAssert(res == MPI_SUCCESS || res == MPI_ERR_IN_STATUS);
    if (nr_finished == MPI_UNDEFINED || nr_finished == 0)
      return false;

    // release finised requests
    while (nr_finished--) {
      if (res == MPI_ERR_IN_STATUS) {
	if (statuses[nr_finished].MPI_ERROR != MPI_SUCCESS) {
	  char msg[MPI_MAX_ERROR_STRING+1];
	  int len;
	  MPI_Error_string(statuses[nr_finished].MPI_ERROR, msg, &len);
	  msg[len] = '\0';
	  PInsist(0, msg);
	}
      }
      releaseMPIRequest(finished[nr_finished]);
    }
    return true;
  }

#else

  static bool waitForSomeRequests(bool mayBlock)
  {
    return false;
  }

#endif


  /// This function lets you check if there are iterates that are
  /// ready to run.

  static bool workReady()
  {
    return !(workQueue_m.empty()
	     && workQueueMessages_m.empty()
#if POOMA_MPI
	     && allocated_requests_m.empty()
#endif
	     );
  }

  /// Run an iterate if one is ready.  Returns if progress
  /// was made.

  static bool runSomething(bool mayBlock = true)
  {
    // do work in this order to minimize communication latency:
    // - process finished messages
    // - issue all messages
    // - do some regular work
    // - wait for messages to complete

    if (waitForSomeRequests(false))
      return true;

    RunnablePtr_t p = NULL;
    if (!workQueueMessages_m.empty()) {
      p = workQueueMessages_m.front();
      workQueueMessages_m.pop_front();
    } else if (!workQueue_m.empty()) {
      p = workQueue_m.front();
      workQueue_m.pop_front();
    }

    if (p) {
      p->execute();
      Iterate<SerialAsync> *it = dynamic_cast<IteratePtr_t>(p);
      if (it) {
	if (--(it->togo()) == 0)
	  delete it;
      } else
	delete p;
      return true;

    } else
      return waitForSomeRequests(mayBlock);
  }

};

/// Adds a runnable to the appropriate work-queue.

inline void add(RunnablePtr_t rn)
{
  if (rn->priority() == -1)
  {
    SystemContext::workQueueMessages_m.push_front(rn);
  }
  else
  {
    SystemContext::workQueue_m.push_front(rn);
  }
}

inline  void concurrency(int){}
inline  int concurrency() {return 1;}
inline  void wait() {}
inline  void mustRunOn(){}


/**
 * Implements a asynchronous scheduler for a data driven execution.
 * Specializes a IterateScheduler.
 *
 * The SerialAsync IterateScheduler, Iterate and DataObject
 * implement a SMARTS scheduler that does dataflow without threads.
 * What that means is that when you hand iterates to the
 * IterateScheduler it stores them up until you call
 * IterateScheduler::blockingEvaluate(), at which point it evaluates
 * iterates until the queue is empty.
 */

template<>
class IterateScheduler<SerialAsync>
{
  friend class DataObject<SerialAsync>;
  friend class Iterate<SerialAsync>;
public:

  typedef DataObject<SerialAsync> DataObject_t;
  typedef Iterate<SerialAsync> Iterate_t;

  IterateScheduler()
    : generation_m(0)
  {}

  ~IterateScheduler() {}

  void setConcurrency(int) {}

  /// Tells the scheduler that the parser thread is starting a new
  /// data-parallel statement.  Any Iterate that is handed off to the
  /// scheduler between beginGeneration() and endGeneration() belongs
  /// to the same data-paralllel statement and therefore has the same
  /// generation number.
  /// Nested invocations are handled as being part of the outermost
  /// generation.

  void beginGeneration()
  {
    // Ensure proper overflow behavior.
    if (++generation_m < 0)
      generation_m = 0;
    generationStack_m.push(generation_m);
  }

  /// Tells the scheduler that no more Iterates will be handed off for
  /// the data parallel statement that was begun with a
  /// beginGeneration().

  void endGeneration()
  {
    PAssert(inGeneration());
    generationStack_m.pop();

#if POOMA_MPI
    // this is a safe point to block until we have "lots" of MPI Requests
    if (!inGeneration())
      while (!SystemContext::haveLotsOfMPIRequests())
	SystemContext::runSomething(true);
#endif
  }

  /// Wether we are inside a generation and may not safely block.

  bool inGeneration() const
  {
    return !generationStack_m.empty();
  }

  /// What the current generation is.

  int generation() const
  {
    if (!inGeneration())
      return -1;
    return generationStack_m.top();
  }

  /// The parser thread calls this method to evaluate the generated
  /// graph until all the nodes in the dependence graph has been
  /// executed by the scheduler.  That is to say, the scheduler
  /// executes all the Iterates that has been handed off to it by the
  /// parser thread.

  void blockingEvaluate()
  {
    if (inGeneration()) {
      // It's not safe to block inside a generation, so
      // just do as much as we can without blocking.
      while (SystemContext::runSomething(false))
	;

    } else {
      // Loop as long as there is anything in the queue.
      while (SystemContext::workReady())
        SystemContext::runSomething(true);
    }
  }

  /// The parser thread calls this method to ask the scheduler to run
  /// the given Iterate when the dependence on that Iterate has been
  /// satisfied.

  void handOff(Iterate<SerialAsync>* it)
  {
    // No action needs to be taken here.  Iterates will make their
    // own way into the execution queue.
    it->generation() = generation();
    it->notify();
  }

  void releaseIterates() { }

private:

  typedef std::list<Iterate_t*> Container_t;
  typedef Container_t::iterator Iterator_t;

  static std::stack<int> generationStack_m;
  int generation_m;

};


/**
 * Implements a asynchronous scheduler for a data driven execution.
 *
 * The DataObject Class is used introduce a type to represent
 * a resources (normally) blocks of data) that Iterates contend
 * for atomic access. Iterates make request for either a read or
 * write to the DataObjects. DataObjects may grant the request if
 * the object is currently available. Otherwise, the request is
 * enqueue in a queue private to the data object until the
 * DataObject is release by another Iterate. A set of read
 * requests may be granted all at once if there are no
 * intervening write request to that DataObject.
 * DataObject<SerialAsync> is a specialization of DataObject for
 * the policy template SerialAsync. 
 *
 * There are two ways data can be used: to read or to write.
 * Don't change this to give more than two states;
 * things inside depend on that.
 */

template<>
class DataObject<SerialAsync>
{
  friend class IterateScheduler<SerialAsync>;
  friend class Iterate<SerialAsync>;
public:

  typedef IterateScheduler<SerialAsync> IterateScheduler_t;
  typedef Iterate<SerialAsync> Iterate_t;


  /// Construct the data object with an empty set of requests
  /// and the given affinity.

  DataObject(int affinity=-1)
    : released_m(queue_m.end()), notifications_m(0)
  {
    // released_m to the end of the queue (which should) also be the
    // beginning.  notifications_m to zero, since nothing has been
    // released yet.
  }
    
  /// for compatibility with other SMARTS schedulers, accept
  /// Scheduler arguments (unused)

  inline DataObject(int affinity, IterateScheduler<SerialAsync>&)
    : released_m(queue_m.end()), notifications_m(0)
  {}

  /// Stub out affinity because there is no affinity in serial.

  int affinity() const { return 0; }

  /// Stub out affinity because there is no affinity in serial.

  void affinity(int) {}

  /// An iterate makes a request for a certain action in a certain
  /// generation.

  inline void request(Iterate<SerialAsync>&, SerialAsync::Action);

  /// An iterate finishes and tells the DataObject it no longer needs
  /// it.  If this is the last release for the current set of
  /// requests, have the IterateScheduler release some more.

  void release(SerialAsync::Action)
  {
    if (--notifications_m == 0)
      releaseIterates();
  }

private:

  /// If release needs to let more iterates go, it calls this.
  inline void releaseIterates();

  /**
   * The type for a request.
   */
  class Request
  {
  public:
    Request() {}
    Request(Iterate<SerialAsync>& it, SerialAsync::Action act) 
      : iterate_m(&it), act_m(act) {}

    Iterate<SerialAsync>& iterate() const { return *iterate_m; }
    SerialAsync::Action act() const { return act_m; }
  private:
    Iterate<SerialAsync>* iterate_m;
    SerialAsync::Action act_m;
  };

  /// The type of the queue and iterator.
  typedef std::list<Request> Container_t;
  typedef Container_t::iterator Iterator_t;

  /// The list of requests from various iterates.
  /// They're granted in FIFO order.
  Container_t queue_m;

  /// Pointer to the last request that has been granted.
  Iterator_t released_m;

  /// The number of outstanding notifications.
  int notifications_m;
};

/// void DataObject::releaseIterates(SerialAsync::Action)
/// When the last released iterate dies, we need to
/// look at the beginning of the queue and tell more iterates
/// that they can access this data.

inline void
DataObject<SerialAsync>::releaseIterates()
{
  // Get rid of the reservations that have finished.
  queue_m.erase( queue_m.begin() , released_m );

  // Next we'll see if we can release some new ones.
  // Get the begin and end iterators.
  released_m = queue_m.begin();
  DataObject<SerialAsync>::Iterator_t end;
  end = queue_m.end();

  // If there are any in the queue, we'll release something.
  if ( released_m != end )
    {
      // Release the first one whatever it is.
      released_m->iterate().notify();
      ++notifications_m;

      // Record what action that one will take
      // and record its generation number
      SerialAsync::Action act = released_m->act();

      // Look at the next iterate.
      ++released_m;

      // If the first one was a read, release more.
      if ( act == SerialAsync::Read )
	{

        // As long as we aren't at the end and we have more reads... 
        while ((released_m != end) && 
               (released_m->act()==SerialAsync::Read))
          {
            // Release it...
            released_m->iterate().notify();
            ++notifications_m;

            // And go on to the next.
            ++released_m;
          }

	}

    }
}

/// void DataObject::request(Iterate&, action)
/// An iterate makes a reservation with this DataObject for a given
/// action in a given generation.  The request may be granted
/// immediately.

inline void
DataObject<SerialAsync>::request(Iterate<SerialAsync>& it, 
                                 SerialAsync::Action act) 

{
  // The request can be granted immediately if:
  // The queue is currently empty, or
  // the request is a read and everything in the queue is a read,
  // or (with relaxed conditions), everything is the same generation.

  // Set notifications dynamically and automatically 
  //     every time a request is made by the iterate 
  it.notifications_m++;

  bool allReleased = (queue_m.end() == released_m);
  bool releasable =  queue_m.empty() ||
      ((act == SerialAsync::Read) &&
       (queue_m.begin()->act() == SerialAsync::Read) && 
       allReleased);

  // Push the request on the stack.
  queue_m.push_back( Request(it, act) );

  // If it's releasable, release it and record the release.
  if (releasable)
    {
      it.notify();
      ++notifications_m;
    }

  // Otherwise, if this is the first nonreleasable iterate,
  // make released_m point to the last element of the queue.
  else if (allReleased)
    {
      --released_m;
    }
}


} // namespace Smarts

//////////////////////////////////////////////////////////////////////

#endif     // _SerialAsync_h_

/***************************************************************************
 * $RCSfile: SerialAsync.h,v $   $Author: richard $
 * $Revision: 1.13 $   $Date: 2004/11/01 18:17:08 $
 ***************************************************************************/
