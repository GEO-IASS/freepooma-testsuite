// -*- C++ -*-
//
// Copyright (C) 2004  Richard Guenther
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

#ifndef CHEETAH_MATCHINGHANDLER_SERIALIZE_H
#define CHEETAH_MATCHINGHANDLER_SERIALIZE_H

//-----------------------------------------------------------------------------
// Classes:
//   Cheetah
//   Serialize<Tag, T>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Serialize is a simple class that serializes/unserializes items to/from
// a buffer.  It can be partially specialized for different types T,
// or for different general tags Tag.  Provided tags are:
//
// 1. 'CHEETAH' is a simple tag type for the default case used by other parts
//    of Cheetah.  Objects are instantiated in place in the provided buffer.
// 3. 'ARRAY' serializes arrays.  API changes a little from other
//    serialize tags as array length must be provided in serialize methods.
//    Objects are instantiated in place in the provided buffer.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files:
//-----------------------------------------------------------------------------

#include <new>
#include <string>


namespace Cheetah {

//----------------------------------------------------------------------
//
// class Serialize
//
// Serialize is a class that can be specialized to pack and unpack
// items of type T to/from a provided buffer of bytes.  It is used by
// the MatchingHandler to prepare and use data sent between MatchingHandler
// send and request calls.  It has two template parameters: a tag, and a data
// type.  The tag can be used to specialize to different categories of
// serialize operations; the data type indicates the type of data that
// will be packed or unpacked.
//
// Serialize specializations should define the following four static
// functions:
//
//   // Return the storage needed to pack the item of type T
//   static int size(const T &item);
//
//   // Pack an item of type T into the given buffer.  Return space used.
//   static int pack(const T &item, char *buffer);
//
//   // Unpack an item of type T from the given buffer.  Set the given
//   // pointer to point at this item.  Return bytes unpacked.
//   static int unpack(T* &p, char *buffer);
//
//   // Delete the item pointed to by the given pointer, that was
//   // unpacked with a previous call to unpack().
//   static void cleanup(T *p);
//
// There is a general template for this class that does nothing, 
// one specialization for a tag 'CHEETAH'.
//
//----------------------------------------------------------------------


//----------------------------------------------------------------------
// Returns padding necessary for word alignment.
//----------------------------------------------------------------------
static inline int padding(int size)
{
  int extra = size % sizeof(void*);
  return (extra == 0) ? 0 : sizeof(void*) - extra;
}


//----------------------------------------------------------------------
// CHEETAH serialize specialization
//----------------------------------------------------------------------

// The general tag type used to specialize Serialize later.

struct CHEETAH
{
  inline CHEETAH() { }
  inline ~CHEETAH() { }
};


// The general template, that does nothing.

template<class Tag, class T>
class Serialize { };


// A specialization for the CHEETAH tag that provides some default ability
// to pack items.

template<class T>
class Serialize< ::Cheetah::CHEETAH, T>
{
public:
  // Return the storage needed to pack the item of type T.
  // For the default case, this is just sizeof(T), but perhaps rounded
  // up to be pointer-word-size aligned.

  static inline int size(const T &)
  {
    return sizeof(double) * ((sizeof(T) + sizeof(double) - 1) / sizeof(double));
  }

  // Pack an item of type T into the given buffer.  Return space used.
  // By default, this just does a placement-new into the buffer,
  // assuming the storage required is sizeof(T).

  static inline int pack(const T &item, char *buffer)
  {
    new ((void*)buffer) T(item);
    return size(item);
  }

  // Unpack an item of type T from the given buffer.  Set the given
  // pointer to point at this item.  Return bytes unpacked.
  // By default, this just recasts the current buffer pointer.

  static inline int unpack(T* &p, char *buffer)
  {
    p = reinterpret_cast<T *>(buffer);
    return size(*p);
  }

  // Delete the item pointed to by the given pointer, that was
  // unpacked with a previous call to unpack().
  // By default, this just runs the destructor on the data, which for
  // many things will do nothing.

  static inline void cleanup(T *p)
  {
    p->~T();
  }
};


//----------------------------------------------------------------------
// ARRAY serialize specialization
//----------------------------------------------------------------------

struct ARRAY
{
  inline ARRAY() { }
  inline ~ARRAY() { }
};


// A specialization for the POINTER tag that provides marshaling of
// arrays.

template<class T>
class Serialize< ::Cheetah::ARRAY, T>
{
public:

  // Return the storage needed to pack count items of type T,
  // This includes the bytes needed to store the size of the array.

  static inline int size(const T* items, const int& count)
  {
    int arraySize = count*sizeof(T);
    return ( Serialize<CHEETAH, int>::size(count)
            + arraySize + padding(arraySize) );
  }

  // Pack an item of type T into the given buffer.  Return space used.
  // By default, this just does a placement-new into the buffer,
  // assuming the storage required is sizeof(T).

  static inline int pack(const T* items, char* buffer, const int& count)
  {
     int n = Serialize<CHEETAH, int>::pack(count, buffer);
     memcpy(n+buffer, items, count*sizeof(T));
     return size(items, count);
  }

  // Unpack an item of type T from the given buffer.  Set the given
  // pointer to point at this item.  Return bytes unpacked.

  static inline int unpack(T* &p, char *buffer, int& count)
  {
     int* iPtr;
     int n = Serialize<CHEETAH, int>::unpack(iPtr, buffer);
     count = *iPtr;
     p = reinterpret_cast<T *>(n+buffer);
     return size(p, count);
  }

  // Delete the item pointed to by the given pointer, that was unpacked with a
  //  previous call to unpack(). By default, this just runs the destructor on
  // the data, which for many things will do nothing.  Memory has been
  // allocated from the provided buffer so no freeing of memory need be done
  // here.

  static inline void cleanup(T *p)
  {
    p->~T();
  }
};


//
// This class is used so that serialization routines can be specialized
// for either delegation (WrappedBool<true>) or CHEETAH 
// (WrappedBool<false>).
//

template<bool flag> class WrappedBool
{
public:
  WrappedBool()  {}
  ~WrappedBool() {}
};

} // namespace Cheetah

#endif // CHEETAH_MATCHINGHANDLER_SERIALIZE_H
