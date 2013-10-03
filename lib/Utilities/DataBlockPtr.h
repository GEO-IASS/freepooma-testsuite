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

#ifndef POOMA_UTILITIES_DATABLOCKPTR_H
#define POOMA_UTILITIES_DATABLOCKPTR_H

//-----------------------------------------------------------------------------

// Classes: 
//   DataBlockPtr<T>
//   DataBlockController
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * DataBlockPtr acts like a RefCountedBlockPtr that has two
 * additional pieces of functionality:
 *   -# it contains a pointer to a Smarts DataObject used for 
 *      constructing and running the Smarts data-flow graph. 
 *   -# it can notify an observer when the destructor is called
 *      (i.e. when views of an engine go away).
 *
 * DataBlockController is an extension of the RefBlockController
 * used by RefCountedBlockPtr. This is where the data object and
 * observable actually reside.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/PoomaSmarts.h"
#include "Utilities/Observable.h"
#include "Utilities/ObserverEvent.h"
#include "Utilities/RefCountedBlockPtr.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template <class T> class SingleObserver;

//-----------------------------------------------------------------------------
//
// Full Description:
//  
// DataBlockController
//
//-----------------------------------------------------------------------------

template <class T> 
class DataBlockController 
  : public RefBlockController<T>
{
public:

  //============================================================
  // Public typedefs
  //============================================================

  typedef Pooma::DataObject_t DataObject_t;

  typedef SingleObservable<int> Observable_t;

  typedef RefBlockController<T> Base_t;
  
  // Required typedef
  
  typedef typename RefBlockController<T>::NoInitTag NoInitTag;

  // A typedef used to store a dynamic-operation ID value.  This will
  // be used to identify when a dynamic operation has been performed on
  // a block of data.  We only want this to happen once, so we will
  // need to indicate what has been done here recently so we can avoid
  // repeating it.  See Utilities/ObserverEvent.h for this typedef.

  typedef ObserverEvent::ID_t DynamicID_t;

  //============================================================
  // Constructors
  //============================================================

  // These simply call the base class constructor.
  // By default, the DataBlockController owns the allocated DataObject.
  // Observable's default constructor builds an unobserved observable.
  // We initialize the dynamicID value to -1, which indicates that it
  // does not refer to any particular dynamic event (they will all be
  // values >= 0).

  explicit 
  DataBlockController(size_t size)
    : Base_t(size), dataObjectPtr_m(new DataObject_t(-1)), owned_m(true),
      dynamicID_m(ObserverEvent::nullID())
  { }

  DataBlockController(size_t size, const T & model)
    : Base_t(size,model), dataObjectPtr_m(new DataObject_t(-1)), owned_m(true),
      dynamicID_m(ObserverEvent::nullID())
  { }

  DataBlockController(T *p, size_t size)
    : Base_t(p,size), dataObjectPtr_m(new DataObject_t(-1)), owned_m(true),
      dynamicID_m(ObserverEvent::nullID())
  { }

  DataBlockController(size_t size, const NoInitTag &tag)
    : Base_t(size,tag), dataObjectPtr_m(new DataObject_t(-1)), owned_m(true),
      dynamicID_m(ObserverEvent::nullID())
  { }

  // This one is new as it sets the affinity for the DataObject.
  // It would be nice to do away with the tag here, but
  // this would be ambiguous for DataBlockController<int>.

  struct WithAffinity 
  { 
    WithAffinity() { }
    WithAffinity(const WithAffinity&) { }
    WithAffinity &operator=(const WithAffinity &) { return *this; }
  };


  DataBlockController(size_t size, int affinity, const WithAffinity &)
    : Base_t(size), dataObjectPtr_m(new DataObject_t(affinity)), owned_m(true),
      dynamicID_m(ObserverEvent::nullID())
  { }
  
  DataBlockController(size_t size, int affinity, const WithAffinity &,
                      const NoInitTag &tag)
    : Base_t(size,tag), dataObjectPtr_m(new DataObject_t(affinity)),
      owned_m(true),
      dynamicID_m(ObserverEvent::nullID())
  { }

  // This one takes a specified DataObject. This is for use by clients
  // that need to maintain ownership of the DataObject, like 
  // CompressibleBlockController.

  DataBlockController(size_t size, DataObject_t &dobj)
    : Base_t(size), dataObjectPtr_m(&dobj), owned_m(false),
      dynamicID_m(ObserverEvent::nullID())
  { }

  DataBlockController(size_t size, const T& model, DataObject_t &dobj)
    : Base_t(size,model), dataObjectPtr_m(&dobj), owned_m(false),
      dynamicID_m(ObserverEvent::nullID())
  { }

  DataBlockController(size_t size, DataObject_t &dobj, const NoInitTag &tag)
    : Base_t(size,tag), dataObjectPtr_m(&dobj), owned_m(false),
      dynamicID_m(ObserverEvent::nullID())
  { }

  // Copy Constructor
  //
  // From RefBlockController:
  //
  //   Ordinarily, this will NOT be used. However, if one
  //   wants to have a RefCountedBlockPtr<T1> where T1 itself
  //   is or contains a RefCountedBlockPtr<T2>, then this
  //   may occasionally be used. When it IS used, a DEEP
  //   copy is required. The RefCounted base class's copy
  //   constructor properly sets the count of the new
  //   class to zero.
  //
  // In addition, we have to deal with the DataBlock and
  // Observable. Since the copy represents a NEW object, we do NOT
  // want to copy the DataObject and Observable.  Rather we create a
  // new DataObject with the same affinity as the old, and a new
  // Observable with the default constructor.
   
  DataBlockController(const DataBlockController &model)
    : Base_t(model),
      dataObjectPtr_m(model.dataObjectPtr_m ? 
		      new DataObject_t(model.affinity()) : 0),
      owned_m(model.dataObjectPtr_m ? true : false),
      observable_m(),
      dynamicID_m(ObserverEvent::nullID())
  { }
  
  // Similar to above, but let the client specify the new DataObject.
  // (They specified it so they own it.)

  DataBlockController(const DataBlockController &model, DataObject_t &dobj)
    : Base_t(model), observable_m(), owned_m(false),
      dataObjectPtr_m(&dobj),
      dynamicID_m(ObserverEvent::nullID())
  { }


  //============================================================
  // Destructor
  //============================================================

  // The base class destructor will take care of the actual block
  // memory, so all we have to do is delete the data-object, if we
  // own it.
  
  ~DataBlockController()
  { 
    if (owned_m) delete dataObjectPtr_m;
  }

  
  //============================================================
  // Accessor and Mutator functions
  //============================================================

  // Attach observers to our observable:

  void attach(SingleObserver<int> *o)
  {
    observable_m.attach(o);
  }

  // Detach the observer from our observable:

  void detach()
  {
    observable_m.detach();
  }

  // Access the smarts data object:

  inline DataObject_t* dataObject() const
  {
    return dataObjectPtr_m;
  }

  inline void dataObject(DataObject_t *obj)
  {
    // If you let people set the dataObject, then you don't
    // own it.  Delete the old one.

    if (owned_m) delete dataObjectPtr_m;
    owned_m = false;
    dataObjectPtr_m = obj;
  }

  // return the affinity for smarts:

  inline int affinity() const
  {
    return dataObjectPtr_m->affinity();
  }

  // set the affinity for smarts:

  inline void affinity(int affin)
  {
    dataObjectPtr_m->affinity(affin);
  }


  //============================================================
  // Notify observer that an action has occurred
  //
  // Event codes sent to observer:
  //  addViewEvent:    Inform the CBC that there is another view.
  //  removeViewEvent: Inform the CBC that view is going away.
  //
  // CBC keeps its own reference count. This seems redundant since
  // this is a RefCounted object already, but there are problems
  // dealing with the underlying reference count in a thread safe
  // manner.
  //============================================================
  
  enum Notifier { addViewEvent, removeViewEvent };

  inline void notifyOnDestruct()
  {
    observable_m.notify(0,removeViewEvent);
  }

  inline void notifyOnConstruct()
  {
    observable_m.notify(0,addViewEvent);
  }


  //============================================================
  // DynamicID operations
  //
  // When there are several different objects using a single
  // DataBlockPtr, and those objects can perform dynamic
  // operations on the data, you must be careful to avoid
  // doing the same operation to a single DataBlockPtr more than
  // once.  We avoid this by having a "dynamic ID" value in
  // the single thing shared by all the other objects, namely
  // this DataBlockPtr.  When those objects try to do a dynamic
  // op involving this object, they first check the dynamic ID.
  // If it matches the ID of the dynamic op they are trying to
  // do, then the operation is skipped for that object.  If it
  // does NOT match, then the operation must be a new one.  In
  // that case, the other objects using this DataBlockPtr can
  // set the dynamic ID here (using setDynamicID), then perform
  // the op.  This could potentially be used for things other
  // than dynamic operations, basically any kind of operation
  // you want to do that might be done to this more than once and
  // you want to do it only once.
  //============================================================

  // Return the ID value for the most recent dynamic operation

  DynamicID_t dynamicID() const
  {
    return dynamicID_m;
  }

  // Change the ID value for the most recent dynamic operation

  void setDynamicID(DynamicID_t id)
  {
    dynamicID_m = id;
  }

private:

  //============================================================
  // Data
  //============================================================
  
  mutable DataObject_t *dataObjectPtr_m;

  // owned is true if we own the dataObjectPtr.
  
  bool owned_m;

  Observable_t observable_m;

  // An identifier for the most recent dynamic operation (see comments
  // above for dynamicID() and setDynamicID() routines for why we
  // have this here).

  DynamicID_t dynamicID_m;
};


/**
 * DataBlockPtr is a customized RefCountedBlockPtr for use in various
 * Pooma Engines. In particular, it adds these things:
 *   -# it contains a pointer to a Smarts DataObject used for 
 *      constructing and running the Smarts data-flow graph. 
 *   -# It is also an observable, for the purpose of notifying
 *      an observer whenever a destructor is called (i.e. when 
 *      views of an engine go away).
 *
 * To accomodate the new data, DataBlockPtr adds the following to the
 * RefCountedBlockPtr interface:
 *
 * Constructors:
 *   - DataBlockPtr(size,affinity) 
 *                  construct a block with affinity and size
 *
 * Members:
 *  - dataObject()   return the Smarts data object
 *  - affinity()     returns the affinity of the data object
 *  - attach()       attach a SingleObserver to our observable
 *  - detach()       detach a SingleObserver to our observable
 */

template <class T, 
          bool BoundsChecked=POOMA_BOUNDS_CHECK_DEFAULT>
class DataBlockPtr 
  : public RefCountedBlockPtr<T,BoundsChecked,DataBlockController<T> >
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef T Pointee_t;
  typedef T Element_t;

  // Convenience typedefs

  typedef DataBlockPtr<T,BoundsChecked> This_t;

  typedef Pooma::DataObject_t DataObject_t;

  typedef SingleObservable<int> Observable_t;

  // The "sister" type, with the other bounds-checking polarity.

  typedef DataBlockPtr<T,!BoundsChecked> That_t;

  // The RefCountedBlockPtr base class:
  
  typedef DataBlockController<T> Controller_t;
  typedef RefCountedBlockPtr<T,BoundsChecked,Controller_t> RCBPtr_t;
  
  // Type of dynamic ID value (see comments in DataBlockController).

  typedef typename Controller_t::DynamicID_t DynamicID_t;

  // Required typedef:
  
  typedef typename RCBPtr_t::NoInitTag NoInitTag;

  //============================================================
  // Constructors
  //============================================================

  // The first several of these simply call the identical constructor
  // in the RCBPtr base class. These will invoke the appropriate 
  // DataBlockController constructors defined above. 

  DataBlockPtr() : RCBPtr_t()
  { }

  explicit DataBlockPtr(size_t size)
    : RCBPtr_t(size)
  { }

  DataBlockPtr(size_t size, const NoInitTag &tag)
    : RCBPtr_t(size,tag)
  { }
  
  // A DataBlockPtr can be initialized to a given value.

  // Note, however, that this may not give the proper memory locality
  // in a threaded execution model since initialization will occur in
  // the parse thread. For this reason, we currently don't provide the
  // obvious constructor that would specify a model and an affinity.

  DataBlockPtr(int size, const T& model)
    : RCBPtr_t(size,model)
  { }

  // A DataBlockPtr that uses foreign data:
  // Affinity is not currently specified. Need to think more
  // about handling foreign data in a NUMA environment. 

  DataBlockPtr(T* foreignData, int size)
    : RCBPtr_t(foreignData,size)
  { }

  // Initialize a block of a particular size with a DataObject having 
  // a particular affinity.
  // 
  // It would be nice to do away with the tag here, but there would
  // be an ambiguity for DataBlockPtr<int>.

  typedef typename Controller_t::WithAffinity WithAffinity_t;
  
  DataBlockPtr(int size, int affin, const WithAffinity_t&)
    : RCBPtr_t(new Controller_t(size,affin,WithAffinity_t()))
  { }

  DataBlockPtr(int size, int affin, const WithAffinity_t&,
	       const NoInitTag &tag)
    : RCBPtr_t(new Controller_t(size,affin,WithAffinity_t(),tag))
  { }

  // Constructors taking an externally supplied DataObject.
  
  DataBlockPtr(int size, DataObject_t &dobj)
    : RCBPtr_t(new Controller_t(size,dobj))
  { }
  
  DataBlockPtr(int size, const T& model, DataObject_t &dobj)
    : RCBPtr_t(new Controller_t(size,model,dobj))
  { }
  
  DataBlockPtr(int size, DataObject_t &dobj, const NoInitTag &tag)
    : RCBPtr_t(new Controller_t(size,dobj,tag))
  { }
  
  // Copy constructor.
  // This is a shallow copy, so nothing special is required.

  DataBlockPtr(const This_t& model)
    : RCBPtr_t(model)
  {  
    if (this->isValid()) this->blockControllerPtr_m->notifyOnConstruct();
  }

  // Same as above except allow the client to specify the
  // DataObject. 
  
  DataBlockPtr(const This_t& model, DataObject_t &dobj)
    : RCBPtr_t(new Controller_t(model, dobj))
  {  
    if (this->isValid()) this->blockControllerPtr_m->notifyOnConstruct();
  }

  // Opposite polarity copy constructor.

  DataBlockPtr(const That_t& model)
    : RCBPtr_t(model)
  {  
    if (this->isValid()) this->blockControllerPtr_m->notifyOnConstruct();
  }

  // Copy constructor.
  // Allow implicit conversions from base class.
  // Since the derived class's additional data is actually
  // in the underlying controller, this is safe. 

  DataBlockPtr(const RCBPtr_t& model)
    : RCBPtr_t(model)
  {  
    if (this->isValid()) this->blockControllerPtr_m->notifyOnConstruct();
  }

  // Copy constructor with an offset. 
  // Create a new DataBlockPtr that is a "view" offset
  // into the model DataBlockPtr.

  DataBlockPtr(const This_t& model, ptrdiff_t offset)
    : RCBPtr_t(model,offset)
  {  
    if (this->isValid()) this->blockControllerPtr_m->notifyOnConstruct();
  }

  //============================================================
  // Destructor
  //============================================================

  // Tell the underlying observable to notify any observers that this
  // copy is going away. Base destructor handles checking the reference
  // count and deleting the controller if necessary, so the rest is
  // trivial.

  ~DataBlockPtr()
  {  
    if (this->isValid()) this->blockControllerPtr_m->notifyOnDestruct();
  }

  //============================================================
  // Assignment operators
  //============================================================

  // Copy assignment operator. 
  // Notify the RHS that we're adding a view, the LHS that we're
  // removing a view, and perform the assignment.
  // I think this is safe. (???)

  This_t & operator=(const This_t & rhs)
  {
    if (this != &rhs)
      {
	if (rhs.isValid()) rhs.blockControllerPtr_m->notifyOnConstruct();
	if (this->isValid()) this->blockControllerPtr_m->notifyOnDestruct();
	RCBPtr_t::operator=(rhs);
      }
    return *this;
  }

  // Same as above, but for opposite polarity
  // of bounds-checking.

  This_t & operator=(const That_t & rhs)
  {
    if (this != &rhs)
      {
	if (rhs.isValid()) rhs.blockControllerPtr_m->notifyOnConstruct();
	if (this->isValid()) this->blockControllerPtr_m->notifyOnDestruct();
	RCBPtr_t::operator=(rhs);
      }
    return *this;
  }

  //============================================================
  // RCBPtr interface
  //
  // Note that we could have picked these up by inheriting
  // publicly from RCBPtr, but this isn't safe since this
  // object must manage two additional pointers using the
  // underlying RCBPtr reference count.
  //
  // Some of this can simply be exposed with using declarations.
  // However, all functions that involve a RCBPtr must be 
  // wrapped to replace the RCBPtr with a DataBlockPtr. These
  // are simple wrappers, however, as they simply call the
  // base class version.
  //
  //============================================================

  // Pointer operations
  // We override these since the base class versions
  // return the base class type. There is no danger here
  // since the base class versions actually do the right 
  // thing.

  This_t & operator++() 
  { RCBPtr_t::operator++(); return *this; }

  This_t & operator--() 
  { RCBPtr_t::operator--(); return *this; }
  
  This_t operator++(int)
  { 
    This_t tmp(*this);
    RCBPtr_t::operator++(); 
    return tmp; 
  }

  This_t operator--(int)
  { 
    This_t tmp(*this);
    RCBPtr_t::operator--(); 
    return tmp; 
  }

  // These calculate a new pointer and return it by value.
  // They do not modify the current object. 

  This_t operator+(ptrdiff_t i) const
  {
    This_t ret(*this);
    ret += i;
    return ret;
  }

  This_t operator-(ptrdiff_t i) const
  {
    This_t ret(*this);
    ret -= i;
    return ret;
  }

  // Note that this returns a pointer to the beginning of
  // the block, not to the current cursor position.

  This_t begin() const
  { return RCBPtr_t::begin(); }

  This_t end() const
  { return RCBPtr_t::end(); }


  //============================================================
  // Accessor and Mutator functions:
  //============================================================

  // Attach observers to our observable:

  void attach(SingleObserver<int> *o)
  {
    this->blockControllerPtr_m->attach(o);
  }

  // Detach the observer from our observable:

  void detach()
  {
    this->blockControllerPtr_m->detach();
  }

  // Access the smarts data object:

  inline DataObject_t* dataObject() const
  {
    return this->blockControllerPtr_m->dataObject();
  }

  // Set the data object pointer.
  // This is for internal use only!!!

  inline void dataObject(DataObject_t *obj)
  {
    this->blockControllerPtr_m->dataObject(obj);
  }

  // return the affinity for smarts:

  inline int affinity() const
  {
    return this->blockControllerPtr_m->affinity();
  }

  // set the affinity for smarts:

  inline void affinity(int affin)
  {
    this->blockControllerPtr_m->affinity(affin);
  }

  // Do two DataBlockPtr's have the same DataObject?

  bool sameDataObject(const DataBlockPtr<T>& x) const
  {
    return dataObject() == x.dataObject();
  }
    
  // Interface for locking controller's RefCounted mutex:
  
  void lockRefCount() const { this->blockControllerPtr_m->lock(); }
  void unlockRefCount() const { this->blockControllerPtr_m->unlock(); }

  //============================================================
  // DynamicID operations
  //
  // When there are several different objects using a single
  // DataBlockPtr, and those objects can perform dynamic
  // operations on the data, you must be careful to avoid
  // doing the same operation to a single DataBlockPtr more than
  // once.  We avoid this by having a "dynamic ID" value in
  // the single thing shared by all the other objects, namely
  // this DataBlockPtr.  When those objects try to do a dynamic
  // op involving this object, they first check the dynamic ID.
  // If it matches the ID of the dynamic op they are trying to
  // do, then the operation is skipped for that object.  If it
  // does NOT match, then the operation must be a new one.  In
  // that case, the other objects using this DataBlockPtr can
  // set the dynamic ID here (using setDynamicID), then perform
  // the op.  This could potentially be used for things other
  // than dynamic operations, basically any kind of operation
  // you want to do that might be done to this more than once and
  // you want to do it only once.
  //============================================================

  // Return the ID value for the most recent dynamic operation

  DynamicID_t dynamicID() const
  {
    return this->isValid() ?
      this->blockControllerPtr_m->dynamicID() :
      ObserverEvent::nullID();
  }

  // Change the ID value for the most recent dynamic operation

  void setDynamicID(DynamicID_t id)
  {
    PAssert(this->isValid());
    this->blockControllerPtr_m->setDynamicID(id);
  }

private:

  // Make That_t a friend class.

  friend class DataBlockPtr<T,!BoundsChecked>;

}; // DataBlockPtr<T,BoundsChecked>

template <class T, bool C1, bool C2>
ptrdiff_t operator-(const DataBlockPtr<T,C1> &first,
		    const DataBlockPtr<T,C2> &second)
{
  return first.currentPointer() - second.currentPointer();
}


// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_DATABLOCKPTR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DataBlockPtr.h,v $   $Author: richard $
// $Revision: 1.24 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
