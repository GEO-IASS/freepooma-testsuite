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

#ifndef POOMA_UTILITIES_OBSERVABLE_H
#define POOMA_UTILITIES_OBSERVABLE_H

//-----------------------------------------------------------------------------
// Classes:
// Observable<T>
// SingleObservable<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * A base or wrapper class for an object of type T that
 * needs to allow other objects to 'observe' it.
 *
 * Observable, with Observer,
 * is used to implement the observer pattern.  Observer<T> objects will
 * register themselves with the Observable, and the Observable will notify
 * them of changes to the observed object.
 *
 * SingleObservable<T>: An optimized observable that can only be viewed 
 * by a single observer.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Threads/PoomaMutex.h"
#include "Utilities/Observer.h"
#include "Utilities/PAssert.h"
#include <vector>


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {


/**
 * The Observable class, along with the Observer class, are used to implement
 * the observer pattern.  In this pattern, there are two sets of objects:
 *   -# Observable<T> objects, which contain a list of Observer<T> pointers.
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
 * Observable<T> is actually a base or wrapper class for the actual class
 * that is to be observed.  The template parameter T is the type of object
 * that is being observed.  This is the same template parameter for Observer.
 * Observer<T> can attach as an observer of Observable<T> objects.  When
 * Observable<T> notifies its observers of events, it calls the method:
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
 * Observable contains a list of Observers, methods 'attach(Observer<T>)' and
 * 'detach(Observer<T>)' to add and remove Observes, and a 'notify(int event)'
 * method which will pass on the event to to the notify method of all
 * attached Observer's.  Observer can be used as a pass class for some
 * object of type T, as in:
 *
 *   class A : public Observable<A> { ...}
 *
 * or as a wrapper class for some class A.
 *
 * When an Observable is deleted, it notifies each registered Observer that
 * it is being deleted by using the reserved event code '0'.  When an
 * Observer gets a notification of this, it should NOT try to call 'detach'
 * for that Observable; it should just remember that that Observable is no
 * longer available and assume that it has been 'detached' already.
 */

template<class T>
class Observable
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  enum { deleteEvent = 0 };


  //============================================================
  // Constructors
  //============================================================

  // The constructor for Observable initializing the reference to the
  // object being observed, and sets up an empty list of Observers.

  Observable(T &o) : observed_m(o), count_m(0)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  // When destructed, an Observable informs all registered objects
  // that it is going away.  It does this by calling notify with the
  // special reserved event code '0'.

  ~Observable()
    {
      notify(deleteEvent);
    }


  //============================================================
  // Observable accessors
  //============================================================

  // Return the number of observers we have registered

  int observers() const
    {
      return count_m;
    }

  //============================================================
  // Observable operations
  //============================================================

  // Allow an Observer to register with this Observable.  This does not
  // check for duplicates, so if the same object attaches twice, it will
  // be notified twice.

  void attach(Observer<T> *o)
    {
      mutex_m.lock();
      observers_m.push_back(o);
      count_m += 1;
      mutex_m.unlock();
    }

  void attach(Observer<T> &o)
    {
      attach(&o);
    }

  // Allow an Observer to indicate it no longer wants to be informed of
  // events from this Observable.  If the object is not currently registered,
  // this results in an assertion failure.

  void detach(Observer<T> *o)
    {
      mutex_m.lock();
      for (int i=0; i < count_m; ++i) {
	if (observers_m[i] == o) {
	  count_m -= 1;
	  observers_m.erase(observers_m.begin() + i);
	  break;
	}
      }
      mutex_m.unlock();
    }
 
  void detach(Observer<T> &o)
    {
      detach(&o);
    }

  // When called, notify calls the notify method in each attached Observer,
  // passing on which observed object this is referring to and what the
  // event code is.

  inline void notify(int event)
    {
      //      mutex_m.lock();
      for (int i=0; i < count_m; ++i)
	observers_m[i]->notify(observed_m, event);
      //      mutex_m.unlock();
    }

  // When called, notify calls the notify method in each attached Observer,
  // passing on which observed object this is referring to and what the
  // event code is.

  inline void notify(const ObserverEvent &event)
    {
      //      mutex_m.lock();
      for (int i=0; i < count_m; ++i)
	observers_m[i]->notify(observed_m, event);
      //      mutex_m.unlock();
    }

private:
  // A reference to the object being observed.  This is passed on to the
  // Observers in the 'notify' method.

  T &observed_m;

  // The list of currently attached observers.  We store pointers since
  // we will be calling the virtual function 'notify' in each.

  std::vector<Observer<T> *> observers_m;

  int count_m;

  Pooma::Mutex_t mutex_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Observable();
  Observable(const Observable<T> &);
  Observable<T> &operator=(const Observable<T> &);
};


/**
 * SingleObservable is an optimized observable that can only be observed
 * by one observer.
 */

template<class T>
class SingleObservable
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  enum { deleteEvent = 0 };


  //============================================================
  // Constructors
  //============================================================

  // The constructor for Observable initializing the reference to the
  // object being observed, and sets up an empty Observer.

  SingleObservable() : observer_m(0)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  // When destructed, an Observable informs all registered objects
  // that it is going away.  It does this by calling notify with the
  // special reserved event code '0'.

  ~SingleObservable()
    {
      notify(T(),0);
    }


  //============================================================
  // Observable operations
  //============================================================

  // Allow an Observer to register with this Observable.

  void attach(SingleObserver<T> *o)
    {
      PAssert(observer_m == 0); // only one observer is legal
      observer_m = o;
    }

  void attach(SingleObserver<T> &o)
    {
      attach(&o);
    }

  // Allow an Observer to indicate it no longer wants to be informed of
  // events from this Observable.

  void detach()
    {
      observer_m = 0;
    }

  // When called, notify calls the notify method in each attached Observer,
  // passing on which observed object this is referring to and what the
  // event code is.

  inline void notify(const T& value, int event)
    {
      if (observer_m != 0)
	observer_m->notify(value, event);
    }

  inline void notify(const T& value, const ObserverEvent &event)
    {
      if (observer_m != 0)
	observer_m->notify(value, event);
    }


private:
  // The currently attached observer.

  SingleObserver<T> *observer_m;

  // The copy constructors are made private and undefined
  // since they should not be used

  SingleObservable(const SingleObservable<T> &);
  SingleObservable<T> &operator=(const SingleObservable<T> &);
};

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_OBSERVABLE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Observable.h,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
