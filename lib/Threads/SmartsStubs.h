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

#ifndef POOMA_THREADS_SMARTSSTUBS_H
#define POOMA_THREADS_SMARTSSTUBS_H

/** @file
 * @ingroup IterateSchedulers
 * @brief
 * Stub scheduler for serial in-order evaluation.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/IterateSchedulers/IterateScheduler.h"
#include "Threads/IterateSchedulers/Runnable.h"

namespace Smarts {

//----------------------------------------------------------------------
// The tag class we'll use for the template parameter.
//----------------------------------------------------------------------

class Stub
{
public:
  enum Action { Read, Write };
};

template<> class Iterate<Stub>;
template<> class IterateScheduler<Stub>;
template<> class DataObject<Stub>;


////////////////////////////////////////////////////////////////////////
//
// The specialization of Iterate for Stub
//
////////////////////////////////////////////////////////////////////////

template<>
class Iterate<Stub> : public Runnable
{
public:
    // Construct the Iterate with:
    //   The scheduler it will inform when it is ready.
    //   Its affinity.
  inline Iterate(IterateScheduler<Stub> & scheduler, int affinity=-1);

  // The dtor is virtual because the subclasses will need to add to it.
  virtual ~Iterate() {}
  virtual void run() = 0;
  int affinity() { return 0; }
  int hintAffinity() { return 0; }

  void affinity(int) {}
  void hintAffinity(int) {}
  void generation(int gen) { generation_m=gen; }
  int generation() { return generation_m; }
  
protected:
  IterateScheduler<Stub> &scheduler_m;

private:
  int generation_m;

};


////////////////////////////////////////////////////////////////////////////
// The specialization of IterateScheduler for Stub                        //
////////////////////////////////////////////////////////////////////////////

template<>
class IterateScheduler<Stub>
{
public:
  //---------------------------------------------------------------------------
  // Some type definition used in all Schedulers
  //---------------------------------------------------------------------------
  // Tell the scheduler that we are beginning a new generation.
  inline
  void beginGeneration() { generation_m++; }

  // Tell the scheduler that we are finishing the current generation.
  inline
  void endGeneration() {}

  inline
  void blockingEvaluate() {}

  inline
  void releaseIterates() { }

  //---------------------------------------------------------------------------
  // Constructors & Destructors;
  //---------------------------------------------------------------------------
  
  inline
  IterateScheduler() 
    : generation_m(0) 
  { }

  inline
  ~IterateScheduler() {}

  // Return the current generation.
  inline int generation() const {  return generation_m; }

  // This class simply runs an iterate as soon as it is handed off.

  inline
  void handOff(Iterate<Stub>* it)
  {
    it->run();
    delete it;
  }
  
protected:
private:
  // Record the current generation.
  int generation_m;

}; // class IterateScheduler<Stub>

inline Iterate<Stub>::Iterate(IterateScheduler<Stub> & scheduler, int)
  : scheduler_m(scheduler) 
{
  generation(scheduler.generation());
}


////////////////////////////////////////////////////////////////////////
//
// The specialization of DataObject for Stub                        //// 
// Concrete class for holding the access requests to a user object. ////
//
////////////////////////////////////////////////////////////////////////

template<>
class DataObject<Stub>
{
public:

    // There are two ways data can be used: to read or to write.
    // Don't change this to give more than two states:
    // things inside depend on that.

    // Construct the data object with an empty set of requests
    // and the given affinity.
    inline DataObject(int=-1) {}

    // Get the affinity.
    inline int affinity() const { return 0; };

    // Set the affinity.
    inline void affinity(int) {}

    // An iterate makes a request for a certain action in a certain
    // generation.
    inline void request(Iterate<Stub>&, Stub::Action) {}

    // An iterate finishes and tells the DataObject it no longer needs
    // it.  If this is the last release for the current set of requests,
    // have the IterateScheduler release some more.
    inline void release(Stub::Action) {}

protected:
private:

};

inline void concurrency(int)
{
}

inline int concurrency()
{
  return 1;
}

inline void wait()
{
}

inline void add(Runnable *runnable)
{
  runnable->execute();
  delete runnable;
}

} // namespace Smarts


#endif // POOMA_THREADS_SMARTSSTUBS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SmartsStubs.h,v $   $Author: richard $
// $Revision: 1.18 $   $Date: 2004/11/01 18:17:07 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
