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

#ifndef POOMA_UTILITIES_REFCOUNTEDPTR_H
#define POOMA_UTILITIES_REFCOUNTEDPTR_H

//-----------------------------------------------------------------------------
// Classes: 
//   RefCountedPtr<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 *  RefCountedPtr<T> - reference counted pointer-to-T.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/ElementProperties.h"
#include "Pooma/Configuration.h"


///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {


/**
 * RefCountedPtr<T> is a smart-pointer class that provides reference
 * counting for objects of type T. T must provide the same interface
 * and semantics as RefCounted, which is usually accomplished by
 * inheriting from RefCounted.
 */  

template <class T>
class RefCountedPtr
{
public:
  
  //============================================================
  // Exported typedefs
  //============================================================

  typedef RefCountedPtr<T> This_t;
  typedef T Pointee_t;
  
  //============================================================
  // Constructors
  //============================================================

  // Default constructor initializes pointer to NULL

  RefCountedPtr() : ptr_m(0) { }

  // Main constructor - takes a pointer to an already created
  // RefCounted object.

  RefCountedPtr(T * const pT)
    : ptr_m(pT) 
    { if (isValid()) ptr_m->addReference(); }

  // Copy constructors copy the raw pointer and increment 
  // the reference count. 

  RefCountedPtr(const This_t &model)
    : ptr_m(model.ptr_m)
    { if (isValid()) ptr_m->addReference(); }


  //============================================================
  // Destructor
  //============================================================

  ~RefCountedPtr();
  
  //============================================================
  // Assignment operators
  //============================================================

  // Assignment operators increment the reference count.

  RefCountedPtr & operator=(const RefCountedPtr &);
  RefCountedPtr & operator=(T *);

  //============================================================
  // Accessors and Mutators
  //============================================================

  // Two ways to use these pointers: 
  // member selection and dereferencing.

  inline T * operator->() const { return ptr_m;  }
  inline T & operator*()  const { return *ptr_m; }

  // Comparison functions

  bool operator==(const This_t& a) const 
  { return ptr_m == a.ptr_m; }

  bool operator!=(const This_t& a) const 
  { return ptr_m != a.ptr_m; }

  // Removes reference and sets pointer field to NULL.

  void invalidate(); 

  // Check to see if a pointer is valid:

  inline bool isValid() const { return ptr_m != 0; }

  // Check to see if the pointer is shared.

  inline bool isShared() const { return ptr_m->isShared(); }

  // Return the current value of the reference count.

  inline int count() const { return ptr_m->count(); }

  // Make private copy of data pointed to by this reference.

  RefCountedPtr<T> & makeOwnCopy();

  // Interoperability with non-POOMA code may require access to the
  // raw data pointer. Thus we provide the following accessor
  // functions. These should be used with care as the returned pointer
  // is not reference counted.

  inline T * rawPointer() { return ptr_m; }
  inline const T * rawPointer() const { return ptr_m; }


#if ! POOMA_NO_TEMPLATE_FRIENDS

private:

  // Make RefCountedBlockPtr a friend.

  template <class T2, bool val, class Controller>
  friend class RefCountedBlockPtr;

#endif

  // The pointer itself.

  T * ptr_m;

};


//////////////////////////////////////////////////////////////////////
//
// Inline implementation of the functions for RefCountedPtr<T>
//
//////////////////////////////////////////////////////////////////////

/// void invalidate()
/// Delete reference and set pointer field to NULL

template <class T>
inline void RefCountedPtr<T>::invalidate()
{ 
  if ( isValid() && ptr_m->removeRefAndCheckGarbage() )
    delete ptr_m;
  ptr_m = 0;
}

/// ~RefCountedPtr
/// Destructor. Work done by invalidate().

template <class T>
inline RefCountedPtr<T>::~RefCountedPtr()
{
  invalidate();
}

/// RefCountedPtr<T>& operator=(const RefCountedPtr<T>& rhs)
/// Assignment of RefCountedPtr<T>

template <class T>
inline RefCountedPtr<T> & 
RefCountedPtr<T>::operator=(const RefCountedPtr<T>& rhs)
{
  // Check self-assignment.

  if (ptr_m != rhs.ptr_m) 
    {
      // First unlink from the one we're pointing to now, and collect
      // garbage if that was the last reference.  

      // (This can potentially throw, but only if *ptr_m is already in
      // an inconsistent state [refcount == 0, but still in
      // existence], so the fact that *this will be left in an
      // inconsistent state isn't a big deal.)

      if ( isValid() && ptr_m->removeRefAndCheckGarbage() )
	delete ptr_m;

      // Now assign the new one.
      
      ptr_m = rhs.ptr_m;
      if ( isValid() ) ptr_m->addReference();
    }

  return *this;
}

/// RefCountedPtr<T>& operator=(T *pp)
/// Assignment of pointer to T to a RefCountedPtr<T>.

template <class T>
inline RefCountedPtr<T> & 
RefCountedPtr<T>::operator=(T *pp)
{
  // If we already manage this pointer, just return.

  if (ptr_m != pp)
    {
      // First unlink from the one we're pointing to now, and collect
      // garbage if that was the last reference.

      // (This can potentially throw, but only if *ptr_m is already in
      // an inconsistent state [refcount == 0, but still in
      // existence], so the fact that *this will be left in an
      // inconsistent state isn't a big deal.)

      if ( isValid() && ptr_m->removeRefAndCheckGarbage() )
	delete ptr_m;

      // Now assign the new one.
 
      ptr_m = pp;
      if ( isValid() ) ptr_m->addReference();
    }

  return *this;
}

/// RefCountedPtr<T>& makeOwnCopy()
/// If we aren't the sole owner of the data, make a private copy.
/// Returns itself for use in chained expressions.

template <class T>
inline RefCountedPtr<T> & 
RefCountedPtr<T>::makeOwnCopy()
{
  // If more than one thing is referring to this one

  if ( isValid() && ptr_m->isShared() )
    {
      // First try to allocate new memory and assign to a temporary
      // pointer. If this throws, *this will still be in a consistent
      // state. 

      // ElementProperties<T>::clone() is used to allow specialization 
      // for T's that have shallow copy semantics, or that provide a
      // virtual clone method.

      T * temp = ElementProperties<T>::clone(*ptr_m);

      // Remove reference from copy-ee. It was shared so there is no
      // garbage to collect. (This can throw, but only if *ptr_m was
      // already corrupt.)

      ptr_m->removeReference();

      // Assign pointer to new object.

      ptr_m = temp;

      // Increment reference for new copy.

      ptr_m->addReference();
    }

  return *this;
}

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_REFCOUNTEDPTR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RefCountedPtr.h,v $   $Author: richard $
// $Revision: 1.18 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
