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

#ifndef POOMA_ENGINE_COMPRESSIBLEBLOCK_H
#define POOMA_ENGINE_COMPRESSIBLEBLOCK_H

//-----------------------------------------------------------------------------
// Classes: 
//   CompressibleBlock<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Engine
 * @brief
 * A CompressibleBlock (CBlock) manages a block of data that can be
 * compressed to a single value.
 *
 * This data must be reference counted, and
 * so the actual data is managed by the private nested CompressibleBlock-
 * Controller (CBC) class.  CompressibleBlock has a RefCountedPtr<CBC>.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/Observable.h"
#include "Utilities/RefCounted.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/PAssert.h"
#include "Utilities/DataBlockPtr.h"
#include "Utilities/ElementProperties.h"
#include "Utilities/Statistics.h"
#include "Pooma/Pooma.h"
#include "Threads/PoomaMutex.h"

#include <stdlib.h> // for rand()

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template <class T> class SingleObserver;


/**
 * CompressibleBlock (CBlock) provides the storage for Pooma's compressible 
 * brick engine. Compressible brick engines need a data storage class that 
 * manages a block of data that has the following characteristics:
 *  -# It manages a block of data of some size.
 *  -# If all of the data in the block has the same value, then
 *     a single value can be stored. This is the "compressed" state.
 *  -# It must provide access to the uncompressed data; i.e.
 *     it is sometimes necessary to take an "uncompressed view".
 *  -# While uncompressed views exist, the data must remain 
 *     "incompressible". Once the last uncompressed view is destroyed,
 *     the block should attempt to compress itself.
 *  -# When the block compresses or uncompresses, it must notify
 *     users (compressible bricks and brickviews) so that they
 *     can sync their internal state with that of the block.
 *  -# This all needs to be thread safe. Compression/decompression
 *     usually occur when iterates are run or destroyed, so much
 *     of this happens asynchronously.
 *  -# Whether the block is compressed or uncompressed, it must appear
 *     to be the same logical object. In particular, it must maintain
 *     the same identity as viewed by the runtime system.
 *  -# It needs to have shallow copy semantics.
 *  -# If Pooma::neverCompress() returns true, then everything should
 *     always be in the uncompressed state.
 *
 * Since there can be multiple views (shallow copies) of the compressible
 * block (both of the data and the data-object), it is important that the 
 * state be maintained by a separate object. Thus, CompressibleBlock is 
 * simply an envelope that holds a ref-counted pointer to a Compressible-
 * BlockController (CBC). The CBC is a private nested class and thus does 
 * not appear directly in the user interface. However it is an important 
 * element of the design. CBlock itself is just an envelope that passes
 * messages on to the CBC.
 *
 * Given all that, here is how views of the data work: There can be 
 * compressible brick views (CompBrick and CompBrickViews) of the data. 
 * These all share the same CBlock, and they do not care whether the
 * CBlock is compressed or not. There can also be plain BrickViews
 * of the uncompressed data. When the last BrickView is destroyed, the
 * CBlock should try to compress itself. Here is a diagram, then,
 * of the possible views:
 *
 * <PRE>
 *  BrickView---->DataBlockPtr<-----CBlock<-----CompBrick
 *    BrickView------^               ^----CompBrickView
 * </PRE>
 *
 * In order for CompBrick and CBlock to have consistent compressibility
 * state, we must use the observer pattern and make it thread safe. Thus
 * we have the following relationships amongst the classes:
 *
 *  CBC observes the DataBlockPtr
 *      (~DataBlockPtr calls notify)
 *
 *  CompBrick/CompBrickView observes the CBlock (the CBC really
 *      as attach/detach simply pass the request on)
 *      (CBC::trycompress and CBC::uncompress call notify)
 */

template <class T>
class CompressibleBlock
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef CompressibleBlock<T>    This_t;
  typedef T                       Element_t;

  typedef Pooma::DataObject_t     DataObject_t;

  enum Notifier
  {
    notifyDestruct   = 0,
    notifyUncompress = 1,
    notifyCompress   = 2
  };

  //============================================================
  // Constructors
  //============================================================

  // The default constructor creates an "invalid" ref-counted
  // pointer:
  
  CompressibleBlock()
    : controller_m(0)
  { }

  // CompressibleBlockController (CBC) is a private nested class 
  // defined below...
  
  explicit CompressibleBlock(int size)
    : controller_m(new CompressibleBlockController(size))
  {  }

  CompressibleBlock(int size, int affinity)
    : controller_m(new CompressibleBlockController(size,affinity))
  {  }

  // A CompressibleBlockController can be initialized to a given value.
  // WARNING: If compressibility is turned off globally, this will result 
  // in block initialization in the parse thread, which will likely
  // lead to bad data locality.

  CompressibleBlock(int size, int affinity, const T& model)
    : controller_m(new CompressibleBlockController(size,affinity,model))
  {  }

  // Copy constructor
  // Simply copy the ref-counted pointer to the controller.
  // It does not matter if the controller gets changed asynchronously
  // during this operation. 
  
  CompressibleBlock(const CompressibleBlock &block)
    : controller_m(block.controller_m)
  {  }
    
  //============================================================
  // Destructor
  //============================================================

  // Trivial since controller_m is a smart pointer.
  
  ~CompressibleBlock()
  {
    PAssert(!controller_m.isValid() || controller_m->isValid());
  }

  //============================================================
  // Accessor and Mutator functions:
  //============================================================
  
  inline int size() const
  {
    return controller_m.size();
  }
  
  inline int capacity() const
  {
    return controller_m.capacity();
  }
  
  inline bool resize(int newsize,
		     const typename DataBlockPtr<T>::NoInitTag &)
  {
    return controller_m->resize(newsize,DataBlockPtr<T>::NoInitTag());
  }

  inline void setSize(int newsize)
  {
    controller_m->setSize(newsize);
  }


  // Access the smarts data object:

  inline DataObject_t* dataObject() const
  {
    PAssert(controller_m.isValid());
    return controller_m->dataObject();
  }

  // Return the affinity for smarts:

  inline int affinity() const
  {
    PAssert(controller_m.isValid());
    return controller_m->dataObject()->affinity();
  }

  // Get compression status:
  // Must have controller locked before doing this.
  // See comments at CompressibleBlockController::compressed().

  inline bool compressed() const
  {
    PAssert(controller_m.isValid());
    return controller_m->compressed();
  }

  // data() returns a raw pointer to the data.

  inline T* data()
  {
    PAssert(controller_m.isValid());
    return controller_m->data();
  }

  // Uncompress the data: 
  // Allocates a new block and fills it with the compressed value.

  inline void uncompress() const
  {
    PAssert(controller_m.isValid());
    controller_m->uncompress();
  }
    
  // Compress the block if we can.

  void tryCompress()
  {
    PAssert(controller_m.isValid());
    controller_m->tryCompress();
  }

  // view() uncompresses and returns a view for making BrickViews

  inline DataBlockPtr<T> view() const
  {
    PAssert(controller_m.isValid());
    return controller_m->view();
  };
  
  // Return the block pointer.
  
  DataBlockPtr<T> dataBlock() const
  {
    PAssert(controller_m.isValid());
    return controller_m->dataBlock();
  }

  // Make own copy

  void makeOwnCopy()
  {
    controller_m.makeOwnCopy();
  }

  // invalidate() invalidates the controller

  void invalidate()
  {
    controller_m.invalidate();
  }

  // Basic sanity checks

  inline bool isControllerPtrValid() const
  {
    return controller_m.isValid();
  }
  
  // NOTE: controller_m->isValid() locks and unlocks the controller, so
  // isControllerValid() MUST NOT be called while the controller is already 
  // locked. Use isControllerValidUnlocked() if you already have a lock on 
  // the controller.
  
  inline bool isControllerValid() const
  {
    return controller_m.isValid() && controller_m->isValid();
  }
  
  inline bool isControllerValidUnlocked() const
  {
    return controller_m.isValid() && controller_m->isValidUnlocked();
  }

  // Checks if the controller ptr is shared.
  // Note that this does not check if there is a BrickView of
  // the underlying block. Should it???

  inline bool isShared() const
  {
    return controller_m.isShared();
  }

  // Comparison operators:

  bool operator!=(const This_t& a) const
  {
    return controller_m != a.controller_m;
  }

  bool operator==(const This_t& a) const
  {
    return controller_m == a.controller_m;
  }

  // Attach observers to the controller:

  void attach(Observer<T*> *o)
  {
    PAssert(controller_m.isValid());
    controller_m->attach(o);
  }

  // Detach the observer from the controller:

  void detach(Observer<T*> *o)
  {
    PAssert(controller_m.isValid());
    controller_m->detach(o);
  }
  
  // Interface for external users to lock and unlock the CBC
  
  void lock() const
  {
    controller_m->lock();
  }
  
  void unlock() const
  {
    controller_m->unlock();
  }

  // Interface to set and query the number of random checks made
  // in tryCompress before the exhaustive search is started:
  
  static int randomTries() { return CompressibleBlockController::randomTries(); }
  
private:

  //---------------------------------------------------------------------------
  // class CompressibleBlockController (aka CBC)
  //
  // CBlock is just an envelope. The data is actually managed by a CBC
  // object.
  //---------------------------------------------------------------------------
   
  class CompressibleBlockController
    : public SingleObserver<int>,
      public Observable<T*>,
      public RefCounted
  {
  public:
    
    //============================================================
    // Exported typedefs and constants
    //============================================================

    typedef CompressibleBlockController      This_t;
    typedef T                                Element_t;
    typedef Pooma::DataObject_t              DataObject_t;

    //============================================================
    // Constructors and Factory Methods
    //============================================================

    // Since these are creating fresh CBCs that don't yet have any outside
    // interactions it is safe to proceed without concern for thread safety.
    
    // By default, CBCs are born compressed. 
    // If Pooma::neverCompress() is true, we must create them
    // in the uncompressed state. 
    
    // Note that Observable keeps a reference to its argument, so passing
    // the uninitialized ptr_m is OK. This is somewhat non-intuitive, though.

    CompressibleBlockController()
      : Observable<T*>(ptr_m),
        size_m(0), 
	compressible_m(true),
	ptr_m(0),           // Size zero is invalid, so don't 
	                    // point this at the compressed data.
	dataObject_m(-1), 
	ucOffset_m(-1),
	viewcount_m(0),
	countUncompressed_m(0)
    { 
      // We initialize the "compressed" value with the
      // default constructor. This is done via the ElementProperties
      // traits class so that this class can be used with objects
      // that don't have a default constructor.
      
      ElementProperties<T>::construct(&compressedData_m);
      
      if (Pooma::neverCompress())
        {
          compressible_m = false;
        }
    }

    explicit
    CompressibleBlockController(int size)
      : Observable<T*>(ptr_m),
        compressible_m(true),
        countUncompressed_m(0),
	viewcount_m(0),
	dataObject_m(-1), 
        size_m(size),
	ptr_m(&compressedData_m), 
	ucOffset_m(-1)
    { 
      ElementProperties<T>::construct(&compressedData_m);
      
      if (Pooma::neverCompress())
        {
          viewcount_m = 1;
          compressible_m = false;
          countUncompressed_m = 1;
          block_m = DataBlockPtr<T>(size_m,dataObject_m);
          ptr_m = block_m.currentPointer();
          block_m.attach(this);
        }
    }
    
    CompressibleBlockController(int size, int affinity)
      : Observable<T*>(ptr_m),
	compressible_m(true),
        countUncompressed_m(0),
	viewcount_m(0),
	dataObject_m(affinity), 
        size_m(size),
	ptr_m(&compressedData_m), 
	ucOffset_m(-1)
    { 
      ElementProperties<T>::construct(&compressedData_m);
      
      if (Pooma::neverCompress())
        {
          viewcount_m = 1;
          compressible_m = false;
          countUncompressed_m = 1;
          block_m = DataBlockPtr<T>(size_m,dataObject_m);
          ptr_m = block_m.currentPointer();
          block_m.attach(this);
        }
    }

    CompressibleBlockController(int size, int affinity, const T& value)
      : Observable<T*>(ptr_m),
	compressible_m(true),
        countUncompressed_m(0),
	viewcount_m(0),
	dataObject_m(affinity), 
        size_m(size),
	ptr_m(&compressedData_m),
	ucOffset_m(-1)
    { 
      ElementProperties<T>::construct(&compressedData_m,value);

      if (Pooma::neverCompress())
        {
          viewcount_m = 1;
          compressible_m = false;
          countUncompressed_m = 1;
          block_m = DataBlockPtr<T>(size_m,value,dataObject_m);
          ptr_m = block_m.currentPointer();
          block_m.attach(this);
        }
    }

    // The CBC is only copied if we are making a deep copy of a block,
    // via the RefCountedPtr::makeOwnCopy() method. Doing this properly 
    // requires that:
    // 1) The observable views this controller, not the model
    // 2) The new CBC is compressible since the block should have no views.
    // 3) The new CBC needs to attach the ptr to the right data.
    // 4) A new DataObject is created with the same affinity as the old one.
    //    and if uncompressed:
    // 5) The data block needs to clone itself.
    // 6) The new CBC needs to attach to observe the cloned block.
    
    CompressibleBlockController(const CompressibleBlockController& model)
      : Observable<T*>(ptr_m),
        compressible_m(!Pooma::neverCompress()),
	viewcount_m(0),
        dataObject_m(model.dataObject_m.affinity()), 
        size_m(model.size_m),
	ucOffset_m(model.ucOffset_m)
    {
      // Lock the model while we get information pertaining to it
      // being compressed or not (such data can't be initialized in
      // the initializer list).
      
      model.lock();

      bool modelCompressed = model.compressed();
      compressedData_m     = model.compressedData_m;
      block_m              = model.block_m;

      if (!modelCompressed)
	{
	  // Decrement the model's viewcount.
	  // Don't need to check for compression here. // HUH???
	  // This only gets called when makeOwnCopy is invoked.
	  // If the underlying CBC wasn't compressed prior to
	  // this, then peeling off a private copy of the 
	  // uncompressed CBC shouldn't change matters.

	  --model.viewcount_m;
	}
           
      model.unlock();
      
      if (modelCompressed)
	{
          PAssert(!Pooma::neverCompress());
          
	  ptr_m = &compressedData_m;
	  countUncompressed_m = 0;
	}
      else
	{	
	  // Copy block and fix pointer and the data object.

	  block_m.makeOwnCopy();
	  ptr_m = block_m.currentPointer();
	  block_m.dataObject(&dataObject_m);
	  countUncompressed_m = 1;

	  // Increment our view count and attach observer.

	  ++viewcount_m;
	  block_m.attach(this);	  
	}

      PAssert(isValidUnlocked());
    }

    // Destructor has to make sure to stop observing the block.
    // CBCs are never destroyed while shared, so it is safe to
    // proceed without locking.

    ~CompressibleBlockController()
    {
      PAssert(!isShared());
      PAssert(!(block_m.isValid() && block_m.isShared()));
      
      // If compressibility is not allowed, the uncompress count
      // should be one. In other words...
      
      PAssert(!Pooma::neverCompress() || countUncompressed_m == 1);

      if (!compressed())
	{
	  block_m.detach();
	}
    }

    //============================================================
    // Assignment operators
    //============================================================

    // No assignment of CBCs
    

    //============================================================
    // Accessor and Mutator functions:
    //============================================================

    // Dynamic resizing functions:
      
    inline size_t size() const
    {
      return block_m.size();
    }
    
    inline size_t capacity() const
    {
      return block_m.capacity();
    }
    
    // This seems a little dangerous. What's it for???
    // It is used to "resize" compressed arrays.
    // This should probably be handled by resize(), but later...
    
    inline void setSize(size_t size)
    {
      size_m = size;
    }

    inline bool resize(size_t newsize,
		       const typename DataBlockPtr<T>::NoInitTag &)
    {
      PAssert(block_m.isValid());
      
      if (block_m.capacity() >= newsize)
	{
	  block_m.resize(newsize, DataBlockPtr<T>::NoInitTag());
	  size_m = block_m.size();
	}
      else
	{
	  // There is not enough capacity, so make some more room ...

	  size_t nsize = newsize * sizeof(T);
	  
#if defined(POOMA_MEMORY_PAGE_SIZE)

	  // Set the size of the new area to be the requested size 
	  // rounded up to the memory page size. 
	  
	  nsize = ((nsize/POOMA_MEMORY_PAGE_SIZE)+1)*POOMA_MEMORY_PAGE_SIZE;

#endif
	  // Reserve a new DataBlockPtr with the bigger size

	  size_t n_ext = nsize/sizeof(T);
	  
	  // Make a new DataBlockPtr - remember that we need to
	  // maintain object identity!
	  
	  DataBlockPtr<T> newdata(n_ext,dataObject_m);
	  
	  newdata.resize(newsize, DataBlockPtr<T>::NoInitTag());
	  
	  // Copy the data over to the new area. 
	  // Only block_m.size() pieces of data need to be copied.
          
          T * pOld = block_m.beginPointer();
          T * pNew = newdata.beginPointer();
          const T * const pEnd = pNew + block_m.size();
          
          while (pNew != pEnd)
            {
	      ElementProperties<T>::construct(pNew++,*pOld++);
            }
            
	  block_m = newdata;
          ptr_m   = block_m.currentPointer();
	  size_m  = newsize;
	}
      return true;
    }

    // Uncompress the data: 
    // Allocates a block and fills it with the compressed value.
    // Only call uncompressUnlocked() if you already have the CBC locked.

    inline void uncompress()
    {
      lock();
      uncompressUnlocked();
      unlock();
    }
    
    void uncompressUnlocked()
    {
      if (compressed())
        {
          PAssert(compressible_m);
          PAssert(!block_m.isValid());
      
	  ++countUncompressed_m;

          // Steps in uncompressing:
          // 1) create a new block, initialized with our DataObject,
          //    and copy the compressed value into it
          // 2) notify the compressed brick views that we've uncompressed
          // 3) attach to the block so that we can observe BrickView activity
      
          block_m = DataBlockPtr<T>(size_m,dataObject_m);
          
	  ++viewcount_m;
	  PAssert(viewcount_m == 1);
	  
          ptr_m = block_m.currentPointer();
          block_m.attach(this);
          
          T *ptr = block_m.beginPointer();
          const T *const end = block_m.endPointer();
          
          while (ptr != end)
   	    {
	      ElementProperties<T>::construct(ptr++,compressedData_m);
	    }
	    
	  // Makes sure CompressibleBrick::notify doesn't try 
	  // to lock the CBC since we're already locked.
	  // A corollary is that all calls to notify should
	  // probably be made during a locked state.
	  
          Observable<T*>::notify(notifyUncompress);
          
          POOMA_INCREMENT_STATISTIC(NumUnCompresses)
        }
    }    
    
    // Compress the block if we can.
    // Only call tryCompressUnlocked() if the CBC is already locked.

    // Should this also lock viewmutex???

    inline void tryCompress()
    {
      if (!Pooma::neverCompress())
        {
          lock();
          tryCompressUnlocked();
          unlock();
        }
    }
    
    void tryCompressUnlocked()
    {
      if (!compressed() && compressible_m && !Pooma::neverCompress())
        {
          PAssert(block_m.isValid());

          int size = block_m.size();
          bool failed = false;
          
          PAssert(ucOffset_m >= -1);
          
          // If ucOffset_m > -1, then this block has failed to 
          // compress before, and ucOffset_m is the index into the 
          // block where the failure was found. We first check this 
          // location as in many applications it is likely to fail 
          // at that spot again. 
          
          if (ucOffset_m > -1 && ucOffset_m < size)
            {
              if (block_m[ucOffset_m] != block_m[0]) failed = true;
            }
          
          // Next we check a few random values. We do this since uncompressed
          // data is often "corrupted" in local spots, and a few random checks
          // may discover the uncompressibility much quicker than starting
          // at the beginning and searching the list linearly.
                    
          if (!failed)
            {
              for (int i = 0; i < randomTries(); ++i)
                {
                  int elem = rand() % size;
                  if (block_m[elem] != block_m[0])
                    {
                      failed = true;
                      ucOffset_m = elem;
                      break;
                    }
                }
            }
            
          // If the random tests didn't turn up anything, we have to do an
          // exhaustive search.
          
          if (!failed)
            {
              const T *const begin = block_m.beginPointer();
              const T *const end   = block_m.endPointer();
          
              const T * ptr;
              for (ptr = begin; ptr != end; ++ptr)
	        {
	          if (*ptr != *begin) break;
	        }
	    
              if (ptr == end) // All values are the same.
	        {
	          // Steps in compressing:
	          // 1) stop observing the block and invalidate our pointer
	          // 2) store the compressed value
	          // 3) notify compressed brick views that we've compressed
	  
	          block_m.detach();
	          block_m.invalidate();
	      
	          --viewcount_m;
	          PAssert(viewcount_m == 0);
	      
	          compressedData_m = *begin;
	          ptr_m = &compressedData_m;
	      
	          Observable<T*>::notify(notifyCompress);
	          
	          POOMA_INCREMENT_STATISTIC(NumSuccessfulTryCompresses)
	          return;
	        }
	      else
	        {
	          ucOffset_m = ptr - begin;
	        }
            }
          POOMA_INCREMENT_STATISTIC(NumUnsuccessfulTryCompresses)
        }
    }

    // view() uncompresses and returns a view for making BrickViews

    DataBlockPtr<T> view()
    {
      lock();
      if (compressed())
	{
	  uncompressUnlocked();
	}
      compressible_m = false;
      unlock();
      return block_m;
    }

    // dataBlock() just returns the block without uncompressing.
    // If data is compressed, then dataBlock().isValid() == false. 
    
    DataBlockPtr<T> dataBlock()
    {
      return block_m;
    }
      
    // data() returns a raw pointer to the data.
    // Note: Users of this function should lock the CBC since ptr_m can
    // be changed externally.

    inline T* data()
    {
      return ptr_m;
    }

    // Check if we're compressed.
    // This function checks whether the block is compressed without locking 
    // the CBC. Any function that is likely to want this information will 
    // have to lock the CBC before calling this. This is necessary since
    // locking and then calling a locking member function would result in
    // a deadlock, and calling a locking "compressed" function and then
    // taking some action (probably setting a lock) would result in a 
    // race condition. 
    
    bool compressed() const
    {
      return ptr_m == &compressedData_m;
    }
  
    // Return the smarts data object:

    DataObject_t* dataObject()
    {
      return &dataObject_m;
    }
    
    // Mutex functions

    void lock() const
    { 
      mutex_m.lock(); 
    }
    
    void unlock() const
    { 
      mutex_m.unlock(); 
    }
    
    // isValid()
    // The block should be valid if the data is uncompressed, and visa versa.
    // This is a basic sanity check that we assert on in many places.
    // A version that does not lock is provided for internal use and
    // for external use by the CompressibleBlock when it already has
    // the controller locked.
    
    bool isValidUnlocked()
    {
      return compressed() || block_m.isValid();
    }
    
    bool isValid()
    {
      lock();
      bool valid = isValidUnlocked();
      unlock();
      return valid;
    }
          
  static int randomTries() { return randomTries_s; }

  private:
      
    // Here we're notified that a brick view viewing our DataBlockPtr
    // is about to destruct.  If this will leave us with a viewcount
    // of one, then it's safe to try and compress the data.
    
    virtual void notify(const int& /* count */, const ObserverEvent &event)
    {
      
      switch (event.event())
      {
        case DataBlockController<T>::addViewEvent:
	  viewmutex_m.lock();
	  ++viewcount_m;
	  viewmutex_m.unlock();
	  break;
        case DataBlockController<T>::removeViewEvent:
	  viewmutex_m.lock();
	  --viewcount_m;
	  if (viewcount_m == 1 && !Pooma::neverCompress())
	    {
	      lock();
	      compressible_m = true;
	      tryCompressUnlocked();
	      unlock();
            }
	  viewmutex_m.unlock();
          break;
        default:
          PInsist(0,
	  "Invalid event code sent to CompressibleBlockController::notify()");
      }
    }

    //============================================================
    // Data
    //============================================================

    // compressible_m is true if it's safe to compress the data

    bool compressible_m;

    // For debugging

    int countUncompressed_m;

    // The DataBlockPtr holds the uncompressed data

    DataBlockPtr<T> block_m;
    
    // We have to keep our own reference count for the block
    // since accessing the block's refcount is difficult to
    // do in a thread-safe manner.

    mutable int viewcount_m;

    // The viewcount must be protected by its own mutex.

    mutable Pooma::Mutex_t viewmutex_m;

    // The Smarts DataObject is here. This way we'll have the
    // same DataObject regardless of whether we're compressed or
    // uncompressed. 
    
    DataObject_t dataObject_m;

    // size of the block

    int size_m;

    // The compressed data.

    T compressedData_m;

    // Pointer to the ref-counted block or to the compressed data
    // depending on compression status.

    T *ptr_m;
    
    // Last uncompressed offset - start search here.
    
    int ucOffset_m;
    
    // Mutex protection for the CBC.
    
    mutable Pooma::Mutex_t mutex_m;
    
    // The number of random compression checks to try before
    // starting the exhaustive search.
    
    enum { randomTries_s = 20 };

  }; // CompressibleBlockController
  
  // The actual controller object:

  RefCountedPtr<CompressibleBlockController> controller_m;
};

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_ENGINE_COMPRESSIBLEBLOCK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CompressibleBlock.h,v $   $Author: richard $
// $Revision: 1.30 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
