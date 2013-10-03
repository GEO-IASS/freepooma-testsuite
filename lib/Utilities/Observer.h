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

#ifndef POOMA_UTILITIES_OBSERVER_H
#define POOMA_UTILITIES_OBSERVER_H

//-----------------------------------------------------------------------------
// Classes:
// Observer<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * An abstract base class for objects which need to use the
 * observer pattern.
 *
 * Observer objects register themselves with one or
 * more Observable objects, and are informed of events by the Observable
 * through their virtual notify() method.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/ObserverEvent.h"


//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * The Observer class, along with the Observable class, are used to implement
 * the observer pattern.  In this pattern, there are two sets of objects:
 *   -# Observable<T> objects, which contain a list of Observer<T>  pointers.
 *   -# Observer<T> objects, which check in as observers of any number of
 *      Observable objects.
 *
 * When the Observer<T> is initialized, it should call the 'attach' method of
 * all Observable<T> objects it needs to watch.  When the Observable changes
 * in some way, for example when it changes state is or is deleted, the
 * Observable will call the 'notify' method of all the Obserers registered
 * with it.  An Observer<T> can stop watching an object, by calling the
 * 'detach' method of that Observable.
 *
 * Observable<T> is actually y a base or wrapper class for the actual class
 * which is to be observed.  The template parameter T is the type of object
 * which is being observed.  This is the same template parameter for Observer.
 * Thus, Observer is templated on the type of object it will observe.
 * Observer<T> can attach as an observer of Observable<T> objects.  When
 * Observable<T>  notifies its observers of events, it calls the method:
 *
 *    virtual void notify(T &observable, int event);
 *
 * in each Observer<T>.  An Observer<T> can attach to more than one
 * Observable<T>, and can distinguish which one is notifying it by the
 * first argument to notify.  'event' is an argument which can contain
 * an integer code to tell the Observer what is happening.  It is up to
 * the class implementing this Observer interface to know how to interpret
 * the event value, based on the type of object it is observing.
 *
 * Since Observer is an abstract base class, objects must inherit from
 * Observer and supply an implementation of the virtual 'notify' method.
 * Derived classes are also responsible for calling 'attach' and 'detach'
 * for the necessary Observable's.
 *
 * When an Observable is deleted, it notifies each registered Observer that
 * it is being deleted by using the reserved event code '0'.  When an
 * Observer gets a notification of this, it should NOT try to call 'detach' 
 * for that Observable; it should just remember that that Observable is no
 * longer available and assume that it has been 'detached' already.
 */

template<class T>
class Observer
{
public:
  //============================================================
  // Constructors
  //============================================================

  // The constructor and destructor for Observer do nothing.  It is
  // up to the class deriving from Observer to call 'attach' and 'detach'
  // when necessary.

  Observer()
    {
    }


  //============================================================
  // Destructors
  //============================================================

  virtual ~Observer()
    {
    }


  //============================================================
  // Observer operations
  //============================================================

  // The one virtual public interface method for Observer.
  // notify is called by an Observable when it needs to tell attached
  // Observers that some event has occurred.  It is up to the derived
  // class, for the type T, to be able to interpret the meaning of the
  // integer event code in the provided ObserverEvent object
  // (or to ignore it, if it needs to).
  // notify is called with a reference to the object being observed
  // and the event which occurred.
  // Note that event code '0' is special; it means that the given Observable
  // is being destroyed, so this Observer should just note that it is no
  // longer attached to that Observable.

  virtual void notify(T &observed, const ObserverEvent &event) = 0;

  // A non-virtual notify that just wrapps the given integer in
  // an event object and passes that on.

  inline void notify(T &observed, int event)
    {
      notify(observed, ObserverEvent(event));
    }
};


template<class T>
class SingleObserver
{
public:
  //============================================================
  // Constructors
  //============================================================

  // The constructor and destructor for Observer do nothing.  It is
  // up to the class deriving from Observer to call 'attach' and 'detach'
  // when necessary.
  SingleObserver() { }


  //============================================================
  // Destructors
  //============================================================

  virtual ~SingleObserver() { }


  //============================================================
  // Observer operations
  //============================================================

  // The one virtual public interface method for Observer.
  // notify is called by an Observable when it needs to tell attached
  // Observers that some event has occurred.  It is up to the derived
  // class, for the type T, to be able to interpret the meaning of the
  // integer event code in the provided ObserverEvent object
  // (or to ignore it, if it needs to).
  // notify is called with a reference to the object being observed
  // and the event which occurred.
  // Note that event code '0' is special; it means that the given Observable
  // is being destroyed, so this Observer should just note that it is no
  // longer attached to that Observable.

  virtual void notify(const T &observed, const ObserverEvent &event) = 0;

  // A non-virtual notify that just wrapps the given integer in
  // an event object and passes that on.

  inline void notify(const T &observed, int event)
    {
      notify(observed, ObserverEvent(event));
    }
};

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_OBSERVER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Observer.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
