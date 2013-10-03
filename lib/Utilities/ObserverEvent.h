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

#ifndef POOMA_UTILITIES_OBSERVER_EVENT_H
#define POOMA_UTILITIES_OBSERVER_EVENT_H

//-----------------------------------------------------------------------------
// Classes:
// ObserverEvent
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * ObserverEvent class - a base class for all events that will be passed
 * on to observers from observables.
 *
 * It includes one integer data member
 * used to indicate to observer subclasses what kind of event it is.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/Unique.h"


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
 * ObserverEvent is the type of object passed to the notify method.
 * It contains an integer indicating the event "code", that classes can
 * examine and use to downcast the event if necessary.  There is also a
 * version of notify that just takes an integer; this is wrapped in
 * an ObserverEvent, and passed on.
 *
 * If you have an event that requires some more information beyond just
 * an event code, make a subclass of ObserverEvent and have the Observer's
 * that get that event cast the event object to the proper type.
 *
 * ObserverEvent's also have a unique ID value, obtained via the
 * ID() method, with type ObserverEvent::ID_t
 */

class ObserverEvent
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  typedef Unique::Value_t ID_t;


  //============================================================
  // Constructors
  //============================================================

  // Constructor: Must specify the integer code

  ObserverEvent(int event)
    : event_m(event), ID_m(Unique::lockedGet())
    {
    }

  // Copy constructor

  ObserverEvent(const ObserverEvent &oe)
    : event_m(oe.event_m), ID_m(oe.ID_m)
    {
    }

  // Assignment operator

  ObserverEvent &operator=(const ObserverEvent &oe)
    {
      event_m = oe.event();
      ID_m = oe.ID();
      return *this;
    }


  //============================================================
  // Destructors
  //============================================================

  // Nothing to do here.  We make this virtual so that we can do
  // dynamic casts.

  virtual ~ObserverEvent()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  // Return our current event code

  inline int event() const
    {
      return event_m;
    }

  // Return our ID value

  inline ID_t ID() const
    {
      return ID_m;
    }

  // Return a value which indicates a "null ID", meaning one that does
  // not refer to any particular event.  This is useful for initializing
  // event values in constructors, etc.  It is static so that you do
  // not need to create a particular ObserverEvent instance to get this
  // value.

  static inline ID_t nullID()
    {
      return (-1);
    }

private:
  // The integer event code

  int event_m;

  // The unique ID vlaue

  ID_t ID_m;
};


/// checkDynamicID(obj, ID) is a specializable function that is used
/// by some classes to check the dynamic ID value stored in the first
/// argument by some means.  If it is the same as the given ID, this
/// returns false.  If it is not the same, it should return true and
/// change the state of obj to indicate that it has "seen" the given ID.
///
/// The default version of this just returns true, generally meaning,
/// "this ID has not been seen, proceed".

template<class Obj>
inline bool checkDynamicID(Obj &, ObserverEvent::ID_t)
{
  return true;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_OBSERVER_EVENT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ObserverEvent.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo

