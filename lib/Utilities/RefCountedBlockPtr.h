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

#ifndef POOMA_UTILITIES_REFCOUNTEDBLOCKPTR_H
#define POOMA_UTILITIES_REFCOUNTEDBLOCKPTR_H

//-----------------------------------------------------------------------------
// Classes: 
//   RefCountedBlockPtr<T,BoundsChecked,Controller>
//   RefBlockController<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * RefCountedBlockPtr and RefBlockController classes.
 *  - RefCountedBlockPtr<T,BoundsChecked,Controller>
 *       smart pointer to reference counted block of data. 
 *       Behaves like a C array, but also provides bounds checking.
 *  - RefBlockController
 *       RefCounted class that actually manages the data and the
 *       bounds checking. 
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include <stddef.h>
#include <new>
#include <iterator>

#include "Utilities/ElementProperties.h"
#include "Utilities/RefCounted.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/PAssert.h"
#include "Pooma/Configuration.h"


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {

/**
 * This class manages the actual data allocation, reference
 * counting, and optional bounds checking for the RefCountedBlockPtr
 * class.
 * 
 * This class holds the pointer to the beginning of the allocated
 * block. It also holds a pointer to one past the end of the
 * allocated block and a pointer to one past the end of the logical
 * size of the block. These two pointers allow the block to be
 * oversized initially so that it can be resized before it has to be
 * replaced. Finally, it holds a bool, the dealloc_m flag, indicating
 * whether or not this data was actually allocated by this
 * class. The begin and logical-end pointers are used to perform
 * optional bounds checking, and dealloc_m allows RefBlockControllers 
 * to work with external data. 
 *
 * As this class inherits from RefCounted, deleting the last
 * reference to an object of this type will result in garbage
 * collection (unless dealloc_m is false).
 *
 * RefBlockController is a model for the Controller concept 
 * used by RefCountedBlockPtr defined below. 
 */

template <class T>
class RefBlockController : public RefCounted
{
public: 

  //============================================================
  // NotInitTag struct
  //============================================================

  // Nested tag class used to select no-initialization policy for
  // various constructors, etc. This can be important in multithreaded codes,
  // where the block needs to be initialized by the processor to which
  // it is assigned.

  struct NoInitTag 
  {
     NoInitTag() { }
     NoInitTag(const NoInitTag &) { }
     NoInitTag & operator=(const NoInitTag &) { return *this; }
  };


  //============================================================
  // Constructors
  //============================================================

  // This is the primary constructor - it allocates a block of
  // memory to hold an array-of-Ts of length "size". Optionally
  // takes a model used in initialization. T is not assumed to
  // possess a default constructor, so we use new to allocate raw,
  // uninitialized memory, and then initialize the memory explicitly
  // with placement operator new. This require explicit
  // destruction. The actual construction and destruction of the
  // individual elements is handled through the ElementProperties
  // traits class.

  explicit 
  RefBlockController(size_t size)
    : pBegin_m(0), pEnd_m(0), pEndOfStorage_m(0), dealloc_m(false)
  {
    // Allocate memory, and set pointers to beginning and ending.  This
    // also sets the dealloc_m flag to true.

    reallocateStorage(size, false);

    // Initialize the storage, via default constructors

    if (!ElementProperties<T>::hasTrivialDefaultConstructor)
      {
	for (T * pt = begin(); pt != end(); ++pt)
	  ElementProperties<T>::construct(pt);
      }
  }

  RefBlockController(size_t size, const T & model)
    : pBegin_m(0), pEnd_m(0), pEndOfStorage_m(0), dealloc_m(false)
  {
    // Allocate memory, and set pointers to beginning and ending.  This
    // also sets the dealloc_m flag to true.

    reallocateStorage(size, false);

    // Always use the traits class for making a copy.

    for (T * pt = begin(); pt != end(); ++pt)
      ElementProperties<T>::construct(pt, model);
  }

  // This constructor allocates raw space and sets the 
  // begin and end pointers to point to the beginning of this
  // space. This is used by the "reserve(size)" call in 
  // RefCountedBlockPtr.
  
  RefBlockController(size_t size, const NoInitTag &)
    : pBegin_m(0), pEnd_m(0), pEndOfStorage_m(0), dealloc_m(false)
  {
    // Allocate memory, and set pointers to beginning and ending.  This
    // also sets the dealloc_m flag to true.

    reallocateStorage(size, false);

    // Skip initialization in this case
  }

  // This constructor sets up a controller for storage
  // owned by somebody else. We never garbage collect 
  // such storage.

  RefBlockController(T *p, size_t size)
    : pBegin_m(p), pEnd_m(p+size), pEndOfStorage_m(p+size), dealloc_m(false)
  { }

  // Copy constructor.
  // Ordinarily, this is only used by a "makeOwnCopy" operation. 
  // However, if one wants to have a RefCountedBlockPtr<T1> where T1 itself
  // is or contains a RefCountedBlockPtr<T2>, then this
  // may occaisionally be used. When it IS used, a DEEP
  // copy is required. The RefCounted base class's copy
  // constructor properly sets the count of the new
  // class to zero.

  RefBlockController(const RefBlockController &model)
    : pBegin_m(0), pEnd_m(0), pEndOfStorage_m(0), dealloc_m(false)
  {
    // Get the size ofs the logical and allocated storage:
    
    size_t allocatedSize = model.pEndOfStorage_m - model.pBegin_m;
    size_t size          = model.end() - model.begin();
    
    // Allocate memory, and set pointers to beginning and ending.  This
    // also sets the dealloc_m flag to true.

    reallocateStorage(allocatedSize, false);
    pEnd_m = pBegin_m + size;	// adjust since it will be smaller

    // Copy over values from previous block

    T * pOld = model.begin();
    T * pNew = begin();

    while (pNew != end())
      {
	ElementProperties<T>::construct(pNew++,*pOld++);
      }
  }


  //============================================================
  // Destructor
  //============================================================

  // The destructor deletes the memory if it is owned by this class.

  ~RefBlockController()
  {
    deleteStorage();
  }


  //============================================================
  // Accessor and Mutator functions
  //============================================================

  // Resize within the limits of the allocated storage.
  // This version performs no initialization, it just adjusts the begin/end
  // pointers.  If there is not enough room, this returns false.

  bool resize(size_t newsize, const NoInitTag &)
  {
    T *pNewEnd = pBegin_m + newsize;
    
    if (pNewEnd <= pEndOfStorage_m)
      {
        pEnd_m = pNewEnd;
        return true;
      }
    else
      {
        return false;
      }
  }
  
  // Resize within the limits of the allocated storage.
  // This just adjusts the begin/end pointers, initializing storage if
  // necessary.  If there is not enough room, this returns false.

  bool resize(size_t newsize)
  {
    bool success = resize(newsize,NoInitTag());
    
    if (!ElementProperties<T>::hasTrivialDefaultConstructor)
      if (success)
	for (T * pt = begin(); pt != end(); ++pt)
          ElementProperties<T>::construct(pt);
      
    return success;
  }
  
  // Resize within the limits of the allocated storage.
  // This just adjusts the begin/end pointers, initializing storage to 'model'
  // if necessary.  If there is not enough room, this returns false.

  bool resize(size_t newsize, const T &model)
  {
    bool success = resize(newsize,NoInitTag());
    
    if (success)
      for (T * pt = begin(); pt != end(); ++pt)
        ElementProperties<T>::construct(pt, model);
      
    return success;
  }
          
  // Resize the data, copying old values over to new storage if new
  // storage must be created.

  T *resizeAndCopy(size_t newsize)
  {
    size_t oldsize = size();

    // First try to resize without mallocing

    if (!resize(newsize, NoInitTag()))
      {
	// There wasn't enough room, so remalloc first ...

	reallocateStorage(newsize, true);

	// Initialize any extra new storage
	    
	if (newsize > oldsize)
	  for (T *pt = begin() + oldsize; pt != end(); ++pt)
	    ElementProperties<T>::construct(pt);
      }

    return begin();
  }

  // Resize the data, copying old values over to new storage if new
  // storage must be created.

  T *resizeAndCopy(size_t newsize, const T &model)
  {
    size_t oldsize = size();

    // First try to resize without mallocing

    if (!resize(newsize, NoInitTag()))
      {
	// There wasn't enough room, so remalloc first ...

	reallocateStorage(newsize, true);
	
	// Initialize any extra new storage

	if (newsize > oldsize)
	  for (T *pt = begin() + oldsize; pt != end(); ++pt)
	    ElementProperties<T>::construct(pt, model);
      }

    return begin();
  }

  // Resize the data, copying old values over to new storage if new
  // storage must be created.  Return the new pointer.

  T *resizeAndCopy(size_t newsize, const NoInitTag &)
  {
    // First try to resize without mallocing

    if (!resize(newsize, NoInitTag()))
      {
	// There wasn't enough room, so remalloc first ...

	reallocateStorage(newsize, true);

	// But skip initializing any extra storage.
      }

    return begin();
  }

  // Get the pointer to the beginning of the block.

  inline T *begin() const
  {
    return pBegin_m;
  }

  inline T *end() const
  {
    return pEnd_m;
  }

  // Get info about size.
  
  inline size_t size() const
  { 
    return static_cast<size_t>(pEnd_m - pBegin_m); 
  }
  
  inline size_t capacity() const
  { 
    return static_cast<size_t>(pEndOfStorage_m - pBegin_m); 
  }
  
  inline bool empty() const
  { 
    return pEnd_m == pBegin_m; 
  }

  // Access to dealloc field:

  inline bool isMine() const
  {
    return dealloc_m;
  }

  // Check that the pointer p is a valid pointer INTO the block; 
  // i.e. that it is dereferencable. The pEnd_m pointer is a valid
  // pointer (i.e. legal), but is not dereferencable. 

  inline bool checkDeref(const T *p) const
  {
    return ((pBegin_m <= p) && (p < pEnd_m));
  }

private:
  //============================================================
  // Private utility methods
  //============================================================

  // A method to delete the existing storage.

  void deleteStorage()
  {
    if (isMine() && pBegin_m != 0)
      {
	if (!ElementProperties<T>::hasTrivialDestructor)
	  for (T *pt = begin(); pt != end(); ++pt)
	    ElementProperties<T>::destruct(pt);

	char *tmp = reinterpret_cast<char *>(pBegin_m);
	delete [] tmp; 
      }
  }

  // A method to reallocate storage, and possibly copy over the
  // old data into the new storage.

  void reallocateStorage(size_t newsize, bool copyold = false)
  {
    // pBegin, pEnd,and bEndOfStorage point to the old storage.
    // Create some new storage, and at the very end make that our
    // new storage.

    // set the size of the new area to be the requested size rounded up
    // to the memory page size, if necessary.

    T *pBeginNew = 0;
    T *pEndNew = 0;
    T *pEndOfStorageNew = 0;

    if (newsize > 0)
      {
	int nsize = newsize * sizeof(T);
#ifdef POOMA_MEMORY_PAGE_SIZE
	nsize = ((nsize/POOMA_MEMORY_PAGE_SIZE)+1)*POOMA_MEMORY_PAGE_SIZE;
#endif
	char *tmp = new char[nsize];
	pBeginNew        = reinterpret_cast<T *>(tmp);
	pEndNew          = pBeginNew + newsize;
	pEndOfStorageNew = pBeginNew + (nsize / sizeof(T));

	// Copy over old storage, if necessary

	if (copyold)
	  {
	    T * pOld = begin();
	    T * pNew = pBeginNew;
	    while (pOld != end() && pNew != pEndNew)
	      ElementProperties<T>::construct(pNew++,*pOld++);
	  }
      }

    // Deallocate the old storage, if necessary

    deleteStorage();

    // Save the new storage and indicate we own it (and must delete it).

    pBegin_m        = pBeginNew;
    pEnd_m          = pEndNew;
    pEndOfStorage_m = pEndOfStorageNew;
    dealloc_m       = true;
  }
  
  //============================================================
  // Private data
  //============================================================

  // Pointers to the actual data

  T *pBegin_m;
  T *pEnd_m;
  T *pEndOfStorage_m;

  // If this flag is true, we allocated (and must deallocate) the data.

  bool dealloc_m;
};                              // RefBlockController


/**
 *   RefCountedBlockPtr<T> is a smart-pointer class that provides
 *   reference counting for arrays of objects of type T.  As long as
 *   only RefCountedBlockPtrs are used to reference this data block,
 *   it will stay around. As soon as the last such object is deleted,
 *   the block is deleted.
 *   
 *   You create a block with: 
 *  
 *        RefCountedBlockPtr<T> p(size);
 *
 *   Then you use p like a pointer. Pointer operations are as
 *   efficient as they would be with a bare pointer as long as bounds
 *   checking is off.
 *   
 *   RefCountedBlockPtr's second template parameter, BoundsChecked, is 
 *   used to enable bounds checking.  It defaults to the value of the
 *   preprocessor symbol POOMA_BOUNDS_CHECK_DEFAULT.  This is
 *   usually set to false. However, if Pooma's configure script is
 *   passed the --bounds option, it will be set to true.
 *   
 *   The third template parameter, Controller, is the object used 
 *   to store the actual data. The concept is modeled by RefBlockController
 *   defined above, and it defaults to RefBlockController<T>.
 *   See DataBlockPtr.h for another example.
 *
 *   A RefCountedBlockPtr with bounds checking explicitly turned on is
 *   declare as:
 *   
 *       RefCountedBlockPtr<T,true> p;
 *   
 *   Similarly, a RefCountedBlockPtr with bounds checking explicitly
 *   turned off is declared as:
 *   
 *       RefCountedBlockPtr<T,false> p;
 */

template <class T, 
  bool BoundsChecked=POOMA_BOUNDS_CHECK_DEFAULT,
  class Controller=RefBlockController<T> >
class RefCountedBlockPtr
{
public:

  typedef std::random_access_iterator_tag  iterator_category;
  typedef T           value_type;
#if POOMA_NONSTANDARD_ITERATOR
  typedef ptrdiff_t   distance_type;
#else
  typedef ptrdiff_t   difference_type;
#endif
  typedef T*          pointer;
  typedef T&          reference;
  
  //============================================================
  // Exported typedefs
  //============================================================

  typedef T         Element_t;
  typedef T         Pointee_t;
  typedef ptrdiff_t Offset_t;
  
  // Convenience typedefs.

  typedef RefCountedBlockPtr<T,BoundsChecked,Controller> This_t;

  // The "sister" type, with the other bounds-checking polarity.

  typedef RefCountedBlockPtr<T,!BoundsChecked,Controller> That_t;
  
  // Nested tag class used to select no-initialization policy for
  // various constructors, etc. This can be important in multithreaded codes,
  // where the block needs to be initialized by the processor to which
  // it is assigned.
  
  struct NoInitTag 
  {
     NoInitTag() { }
     NoInitTag(const NoInitTag &) { }
     NoInitTag & operator=(const NoInitTag &) { return *this; }
  };

  //============================================================
  // Constructors
  //============================================================

  // Default constructor.

  inline RefCountedBlockPtr()
    : offset_m(0)
  { }

  // Initialize a block of a given size, optionally with a model.

  inline explicit RefCountedBlockPtr(size_t size)
    : offset_m(0),
      blockControllerPtr_m(new Controller(size))
  { }

  inline RefCountedBlockPtr(size_t size, const T & model)
    : offset_m(0),
      blockControllerPtr_m(new Controller(size,model))
  { }

  inline RefCountedBlockPtr(size_t size, const NoInitTag &)
    : offset_m(0),
      blockControllerPtr_m(new Controller(size,
#ifndef __MWERKS__
					  typename Controller::NoInitTag()))
#else
					  Controller::NoInitTag()))
#endif
  { 
#ifndef __MWERKS__
    blockControllerPtr_m->resize(0,typename Controller::NoInitTag());
#else
    blockControllerPtr_m->resize(0,Controller::NoInitTag());
#endif    
  }

  // Initialize with a user allocated pointer.
  // This turns off the deallocation and initialization, 
  // but not the bounds checking.

  inline RefCountedBlockPtr(T *p, size_t size)
    : offset_m(0),
      blockControllerPtr_m(new Controller(p, size))
  { }
  
  // Copy constructor

  inline RefCountedBlockPtr(const This_t & model)
    : offset_m(model.offset_m),
      blockControllerPtr_m(model.blockControllerPtr_m)
  { }

  // Copy constructor from a RefCountedBlockPtr of the opposite bounds
  // checking polarity (That_t).

  inline RefCountedBlockPtr(const That_t & model)
    : offset_m(model.offset_m),
      blockControllerPtr_m(model.blockControllerPtr_m)
  { }

  // Copy constructor with offset
  // This was added so that BrickView Engines could
  // initialize their RefCountedBlockPtr without using
  // the operator+ defined below.  (It creates a temporary
  // which adds to compile time.)

  inline RefCountedBlockPtr(const This_t & model, Offset_t offset)
    : offset_m(model.offset_m + offset),
      blockControllerPtr_m(model.blockControllerPtr_m)
  { }

  // Inline destructor for better performance

  inline ~RefCountedBlockPtr() {}

  //============================================================
  // Assignment
  //============================================================

  inline This_t & operator=(const This_t & rhs)
  {
    blockControllerPtr_m = rhs.blockControllerPtr_m;
    offset_m = rhs.offset_m;
    return *this;
  }

  inline This_t & operator=(const That_t & rhs)
  {
    blockControllerPtr_m = rhs.blockControllerPtr_m;
    offset_m = rhs.offset_m;
    return *this;
  }

  //============================================================
  // Accessors
  //============================================================

  // Comparison functions

  template<class T1, bool BoundsChecked1, class Controller1>
  inline bool operator==(
    const RefCountedBlockPtr<T1, BoundsChecked1, Controller1>&) const 
  {
    return false;
  }

  template<class T1, bool BoundsChecked1, class Controller1>
  inline bool operator!=(
    const RefCountedBlockPtr<T1, BoundsChecked1, Controller1>&) const 
  {
    return true;
  }

  inline bool operator==(const This_t& a) const 
  {
    return (beginPointer() == a.beginPointer() && offset() == a.offset());
  }

  inline bool operator!=(const This_t& a) const 
  {
    return (beginPointer() != a.beginPointer() || offset() != a.offset());
  }

  inline bool operator<(const This_t& a) const 
  {
    PAssert(beginPointer() == a.beginPointer());
    return offset() < a.offset();
  }

  inline bool operator>(const This_t& a) const 
  {
    PAssert(beginPointer() == a.beginPointer());
    return offset() > a.offset();
  }

  inline bool operator<=(const This_t& a) const 
  {
    PAssert(beginPointer() == a.beginPointer());
    return offset() <= a.offset();
  }

  inline bool operator>=(const This_t& a) const 
  {
    PAssert(beginPointer() == a.beginPointer());
    return offset() >= a.offset();
  }

  // One for the opposite bounds-checking polarity

  inline bool operator==(const That_t& a) const 
  {
    return (beginPointer() == a.beginPointer() && offset() == a.offset());
  }

  inline bool operator!=(const That_t& a) const 
  {
    return (beginPointer() != a.beginPointer() || offset() != a.offset());
  }

  // Dereferencing operations

  // All operations that attempt to dereference the pointer do
  // bounds-checking first if BoundsChecked is true.  Because
  // BoundsChecked is a compile time parameter, the bounds-checking
  // code should be completely removed by the compiler if
  // BoundsChecked is false.

  // These are const member functions since "const Array" must be able 
  // to modify its data. RefCountedBlockPtr is thus really like a 
  // non-const iterator. We don't currently need a ConstRefCountedBlockPtr,
  // but one could easily be made, and would probably serve as a base
  // class for this class.

  inline T& operator*() const
  {
    T *p = currentPointer();
    boundsAssert(p);
    return *p;
  }

  // Array Indexing.
  
  inline T& operator[](Offset_t i) const
  {
    T *p = currentPointer() + i;
    boundsAssert(p);
    return *p;
  }

  // Member selection

  inline T *operator->() const
  {
    T *p = currentPointer();
    boundsAssert(p);
    return p;
  }

  //============================================================
  // Mutators
  //============================================================

  // RefCountedBlockPtr provides all the usual pointer manipulation
  // functions.  Bounds are not checked here - only on dereferencing.

  inline This_t & operator++()
  {
    ++offset_m;
    return *this;
  }

  inline This_t & operator--()
  {
    --offset_m;
    return *this;
  }

  inline This_t operator++(int)
  {
    This_t save(*this);
    ++offset_m;
    return save;
  }

  inline This_t operator--(int)
  {
    This_t save(*this);
    --offset_m;
    return save;
  }

  inline void operator+=(Offset_t i)
  {
    offset_m += i;
  }

  inline void operator-=(Offset_t i)
  {
    offset_m -= i;
  }

  //============================================================
  // Other iterator functionality.
  //============================================================

  // These calculate a new pointer and return it by value.
  // They do not modify the current object. 

  inline This_t operator+(Offset_t i) const
  {
    This_t ret(*this);
    ret += i;
    return ret;
  }

  inline This_t operator-(Offset_t i) const
  {
    This_t ret(*this);
    ret -= i;
    return ret;
  }

  // Note that this returns a pointer to the beginning of
  // the block, not to the current cursor position.

  inline This_t begin() const
  {
    return This_t(*this, -offset_m);
  }

  inline This_t end() const
  {
    return This_t(*this, size() - offset_m);
  }

  //============================================================
  // Utility Mutators
  //============================================================
  
  // Reserve space. Allows setting aside some uninitialized space
  // for future resizing. This is only valid for uninitialized
  // (i.e. "invalid") objects.

  void reserve(size_t size)
  {
    PAssert(!isValid());    
    typedef typename Controller::NoInitTag NoInit_t;
    blockControllerPtr_m = new Controller(size, NoInit_t());
    blockControllerPtr_m->resize(0,NoInit_t());
    offset_m = 0;
  }

  // The next versions of resize only resize within the currently
  // allocated storage, they do not ever cause a reallocate of memory.
  // The offset is not affected.

  bool resize(size_t size, const NoInitTag &)
  {
    PAssert(isValid());
    typedef typename Controller::NoInitTag NoInit_t;
    return blockControllerPtr_m->resize(size, NoInit_t());
  }
  
  bool resize(size_t size)
  {
    PAssert(isValid());
    return blockControllerPtr_m->resize(size);
  }
  
  bool resize(size_t size, const T &model)
  {
    PAssert(isValid());
    return blockControllerPtr_m->resize(size, model);
  }

  // The next versions of resize WILL copy over old data into
  // newly allocated storage if a realloc is necessary.
  // The offset is not affected.

  void resizeAndCopy(size_t size, const NoInitTag &)
  {
    PAssert(isValid());
    typedef typename Controller::NoInitTag NoInit_t;
    blockControllerPtr_m->resizeAndCopy(size, NoInit_t());
  }

  void resizeAndCopy(size_t size)
  {
    PAssert(isValid());
    blockControllerPtr_m->resizeAndCopy(size);
  }

  void resizeAndCopy(size_t size, const T &model)
  {
    PAssert(isValid());
    blockControllerPtr_m->resizeAndCopy(size, model);
  }

  // Used to do actual destruction, when needed.

  void invalidate()
  {
    blockControllerPtr_m.invalidate();
    offset_m = 0;
  }

  // Make private copy of data pointed to by this reference.
  // makeOwnCopy() returns *this for use in chained computations.

  This_t & makeOwnCopy() 
  {
    blockControllerPtr_m.makeOwnCopy(); 
    return *this;
  }

  //============================================================
  // Utility Accessors
  //============================================================

  // Return our current offset from the beginning of the allocated block.
  // If the offset is non-zero, it means that:
  //   1. This is for a view of the data, looking past the beginning.
  //   2. Or we're using a non-zero-based view, and we have a negative
  //      offset to compensate that our domain*strides will come out
  //      to start at a value greater than zero.

  inline Offset_t offset() const
  {
    return offset_m;
  }

  // Note that isValid queries the state of the controller.
  // It does NOT do bounds checking! 

  inline bool isValid() const
  {
    return blockControllerPtr_m.isValid();
  }

  // Check to see if the block is shared.

  inline bool isShared() const
  {
    return blockControllerPtr_m.isShared();
  }

  // Return the current value of the reference count.

  inline int count() const
  {
    return isValid() ? blockControllerPtr_m.count() : 0;
  }

  // Get info about controller's size.
  
  inline size_t size() const
  { 
    return isValid() ? blockControllerPtr_m->size() : 0; 
  }
  
  inline size_t capacity() const 
  { 
    return isValid() ? blockControllerPtr_m->capacity() : 0; 
  }
  
  inline bool empty() const 
  { 
    return isValid() ? blockControllerPtr_m->empty() : true; 
  }

  // Check to see if pointer is at beginning of block.

  inline bool isAtBeginning() const 
  { 
    return (offset() == 0);
  }

  // Return true if this data is owned by the controller, and will
  // be deallocated by the controller.

  inline bool isMine() const
  {
    return isValid() ? blockControllerPtr_m->isMine() : true;
  }

  // Here are the evil member function that grant access to the raw
  // underlying pointers. Use with care as these pointers are NOT
  // reference counted!
  
  // Note that these are const member functions even though they return 
  // non-const T*. Like the dereference operations above, this is necessary
  // since "const" RCBPtr implies that the pointer object itself is
  // not changed, not that the block itself is const.

  inline T *beginPointer() const
  {
    PAssert(isValid());
    return blockControllerPtr_m->begin();
  }

  inline T *endPointer() const
  {
    PAssert(isValid());
    return blockControllerPtr_m->end();
  }

  inline T *currentPointer() const
  {
    return beginPointer() + offset();
  }

protected:

  // Make That_t a friend class.

  friend class RefCountedBlockPtr<T,!BoundsChecked,Controller>;

  //============================================================
  // Private utility methods
  //============================================================

  // Initialize with an already constructed controller. This allows
  // for derived classes with special controllers to initialize their
  // controllers, which may take additional constructor arguments.
  
  RefCountedBlockPtr(Controller *con)
    : offset_m(0), blockControllerPtr_m(con)
  { }

  // Utility function: Check bounds and throw assertion if
  // checkDeref fails. BoundsChecked is a compile time constant. This
  // must expand to nothing if BoundsChecked == false.

  inline void boundsAssert(T *p) const
  {
    if (BoundsChecked)
      {
        PInsist(isValid() && blockControllerPtr_m->checkDeref(p),
                "RefCountedBlockPtr: Bounds Violation.");
      }
  }

  //============================================================
  // Private data
  //============================================================

  // RefCountedBlockPtr has two data items: an offset into the block,
  // and a reference counted pointer to a Controller object
  // that manages the storage for the block.

  Offset_t offset_m;

  RefCountedPtr<Controller> blockControllerPtr_m;
};


template <class T, bool C1, bool C2, class Controller>
inline ptrdiff_t
operator-(const RefCountedBlockPtr<T,C1,Controller> &first,
          const RefCountedBlockPtr<T,C2,Controller> &second)
{
  return first.currentPointer() - second.currentPointer();
}

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_REFCOUNTEDBLOCKPTR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RefCountedBlockPtr.h,v $   $Author: richard $
// $Revision: 1.64 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
