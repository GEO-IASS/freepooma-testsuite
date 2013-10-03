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

#ifndef POOMA_UTILITIES_REFCOUNTED_H
#define POOMA_UTILITIES_REFCOUNTED_H

//-----------------------------------------------------------------------------
// Classes:
//   RefCounted
//   Shared<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * RefCounted and Shared classes.
 *   - RefCounted: Mix-in class that encapsulates the reference-counting 
 *               of an object.
 *   - Shared<T>:  A template that simply inherits from RefCounted,
 *               and has a member that returns the contained data.
 *               The data member is protected, so this can be used
 *               via inheritance.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Threads/PoomaMutex.h"

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

/**
 * RefCounted is a mix-in class that supports reference counting
 * of an object. It encapsulates the count and provides an interface
 * for manipulating and checking the count.
 *
 * When running in a threaded environment, RefCounted protects the
 * reference count with a mutex, so this class is thread safe.
 */

class RefCounted 
{
public:
  
  //============================================================
  // Constructors, Destructor, Assignment...
  //============================================================

  // Default constructor.
  // Starts count at zero. The client that creates the RefCounted
  // object is responsible for calling addReference().

  RefCounted() 
    : count_m(0) 
    { }

  // Copy constructor.
  // This may appear a bit odd. If a RefCounted object is
  // copied, then this creates a NEW RefCounted object
  // that must be reference counted separately from the
  // old one. Thus its count must be initialized to zero.
  // Ordinarily RefCounted objects aren't copied. However,
  // clients may wish to implement a clone() operation
  // that does explicitly make a (deep) copy.

  RefCounted(const RefCounted &) 
    : count_m(0) 
    { }

  // Trivial destructor:

  ~RefCounted() { }

  //============================================================
  // Accessors
  //============================================================

  bool isShared() const;

  //============================================================
  // Mutators
  //============================================================

  // These simply increment and decrement the reference count.

  void addReference();
  void removeReference();

  // Remove reference and check if it leaves garbage.

  bool removeRefAndCheckGarbage();
  
  // Expose lock and unlock:

  void lock() const
  {
    mutex_m.lock();
  }

  void unlock() const
  {
    mutex_m.unlock();
  }

  // Return the current value of the reference count.

  int count() const;
  
  // Ditto, but without locking the mutex while copying it.
  
  int countUnlocked() const;

private:

  // Assignment is private and unimplemented.

  RefCounted & operator=(const RefCounted &); // UNIMPLEMENTED

  int count_m;

  // Mutex is declared mutable since we must be able to lock and
  // unlock the mutex in accessors that don't change the logical state
  // of the RefCounted object.

  mutable Pooma::Mutex_t mutex_m;

};


//-----------------------------------------------------------------------------
// Inline functions
//-----------------------------------------------------------------------------

inline bool 
RefCounted::isShared() const 
{ 
  mutex_m.lock();
  bool test = count_m > 1;
  mutex_m.unlock();
  return test; 
}


inline void 
RefCounted::addReference()
{
  mutex_m.lock();
  ++count_m;
  mutex_m.unlock();
}

inline void 
RefCounted::removeReference()
{
  mutex_m.lock();
  --count_m;
  PAssert(count_m >= 0);
  mutex_m.unlock();
}

inline bool
RefCounted::removeRefAndCheckGarbage() 
{ 
  mutex_m.lock();
  PAssert(count_m > 0); 
  bool test = --count_m == 0;
  mutex_m.unlock();
  return test;
}

inline int
RefCounted::count() const
{
  mutex_m.lock();
  int count = count_m;
  mutex_m.unlock();
  return count;
}

inline int
RefCounted::countUnlocked() const
{
  return count_m;
}


/**
 *  Simple template class encapsulating a single data item and
 *  inheriting from RefCounted. 
 */

template <class T>
class Shared : public RefCounted
{
public:

  //============================================================
  // Constructors, Destructor, Assignment...
  //============================================================

  Shared(const T &d) : data_m(d) {};

  Shared(const Shared<T> & model) 
    : data_m(model.data_m) 
  { }

  Shared<T> & operator=(const Shared<T> &model)
  { 
    if (&model == this) return *this;
    data_m = model.data_m;
    return *this;
  }

  Shared<T> & operator=(const T & d) 
  { 
    data_m = d; 
    return *this; 
  }

  //============================================================
  // Accessors
  //============================================================

  inline
  T &data() { return data_m; }

  inline
  const T &data() const { return data_m; }

  bool operator==(const Shared<T> &rhs) const
    { return data_m == rhs.data_m; }

  bool operator!=(const Shared<T> &rhs) const
    { return data_m != rhs.data_m; }

protected:

  T data_m;
};



// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_REFCOUNTED_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RefCounted.h,v $   $Author: richard $
// $Revision: 1.19 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
