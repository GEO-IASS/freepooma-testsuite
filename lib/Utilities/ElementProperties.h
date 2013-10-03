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

#ifndef POOMA_UTILITIES_ELEMENTPROPERTIES_H
#define POOMA_UTILITIES_ELEMENTPROPERTIES_H

//-----------------------------------------------------------------------------
// Classes:
//   ElementProperties<T> and specializations
//   TrivialElementProperties<T>
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/** @file
 * @ingroup Utilities
 * @brief
 * Traits class for determining, and possibly modifying, the
 * construction and destruction properties of elements of type T.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <new>

/**
 * Traits class for determining, and possibly modifying, the
 * construction and destruction properties of elements of type T.
 *
 * In detail, this class serves several purposes:
 *  - First it allows RefCountedBlockPtr to optimize away the
 *   constructor and destructor calls for classes with "trivial"
 *   default constructors and destructors (e.g. the native C data
 *   types, but see below for details).
 *  - Second, certain types can be safely copied with memcpy. We 
 *   refer to these as "concrete" types. Such types should set
 *   the concrete trait to true.
 *  - Third, it allows specializations to provide special construct,
 *   clone, and destruct methods that override the default behavior.
 *   The primary reason for this capability is to allow RefCountedPtr
 *   and RefCountedBlockPtr to store deep copies of (and potentially
 *   make further deep copies of) objects that have shallow copy
 *   semantics. This can be done by specializing ElementProperties to
 *   provide construct and clone methods that make deep copies of the
 *   model object.
 *  - Finally, one might want RefCountedPtr<T> to point to an object
 *   that inherits from T. In such situations, asking the
 *   RefCountedPtr to make a deep copy of its pointee would, with the
 *   default behavior, cause the object to be sliced (as T's copy
 *   constructor would be used). If T has a "virtual constructor" (a
 *   virtual clone method), then one can specialize ElementProperties'
 *   clone() method to call the virtual constructor and make the
 *   proper copy.
 *
 * The first capability is provided by defining the two bool fields:
 *  - enum { hasTrivialDefaultConstructor };
 *  - enum { hasTrivialDestructor };
 *
 * hasTrivialDefaultConstructor is true for data types whose default
 * constructors have the same semantics as the native C data types;
 * i.e. they do nothing: no allocations, no default values, etc.
 * Normally RefCountedBlockPtr calls placement operator new to
 * initizlize objects in the space that it allocates and manages.
 * However, this is unnecessary overhead for types whose default
 * constructor does nothing. Thus if hasTrivialDefaultConstructor is
 * true, RefCountedBlockPtr will leave memory uninitialized in the
 * default case.
 *
 * Versions of ElementProperties for the most common C data types are
 * defined below. Similar specializations might also be useful for
 * other statically sized data types, such as TinyArrays. (Note that
 * one could optionally define the specializations for native data
 * types to initialize memory with some obviously corrupt value in
 * order to help track initialization problems during debugging,
 * although purify is probably a better tool for such investigations.)
 *
 * Similarly, hasTrivialDestructor == true causes RefCountedBlockPtr
 * to skip the explicit destructor calls that are normally necessary
 * when destroying an object created with placement new.  
 * This will almost always have the same value as
 * hasTrivialDefaultConstructor, but the additional flexibility
 * carries no additional cost so it was included.
 *
 * The class must also define the following static functions:
 *  - inline static void construct(T * addr, const T & model)
 *  - inline static void T * clone(const T & model);
 *  - inline static void construct(T * addr)
 *  - inline static void destruct(T *addr)
 *
 * If the "trivial" flags are true, then the last two functions will
 * never be called, but they must be defined or the compiler will
 * complain. In these cases it is best to define the functions to
 * throw an exception (see below).
 *    
 * The non-specialized ElementProperties<T> class defines both flags
 * to be false.  It defines the construct methods to use the default
 * and copy constructors with placement new, respectively, under the
 * assumption that these will make deep copies. Finally, it defines
 * the destruct method to explicitly invoke the destructor on the
 * object.
 */

template <class T>
struct ElementProperties
{ 
  //---------------------------------------------------------------------------
  // Convenience typedef (not needed here, but useful for
  // specializations, which may want to copy parts of this definition)

  typedef T This_t;

  //---------------------------------------------------------------------------
  // By default, we assume that the default constructor and the
  // destructor do something.

  enum { hasTrivialDefaultConstructor = false };

  enum { hasTrivialDestructor = false };

  //---------------------------------------------------------------------------
  // We specialize this struct for concrete types. These are types that
  // have no pointers, etc., so that their data can be copied with routines
  // such as std::copy or std::memcpy.

  enum { concrete = false };

  //---------------------------------------------------------------------------
  // Sometimes it is necessary to know if a type is one of the C basic
  // types. The following trait answers this question.

  enum { basicType = false };

  //---------------------------------------------------------------------------
  // Since the above are false, the code will use placement new for
  // initialization, and thus must explicitly call the destructor. The
  // following are the default methods for doing these things.

  static void construct(This_t * addr)
  {
    new (addr) This_t();
  }

  static void construct(This_t * addr, const This_t & model)
  {
    new (addr) This_t(model);
  }

  static This_t * clone(const This_t &model)
  {
    return new This_t(model);
  }

  static void destruct(This_t * addr)
  {
    addr->~This_t();
  }
};


/**
 * Concrete types that have trivial default construction and destruction
 * semantics can just inherit from this:
 */

template <class T>
struct TrivialElementPropertiesBase
{
  typedef T This_t;

  enum { hasTrivialDefaultConstructor = true };

  enum { hasTrivialDestructor = true };

  enum { concrete = true };
  
  static void construct(This_t * addr, const This_t & model)
  {
    new (addr) This_t(model);
  }

  static This_t * clone(const This_t &model)
  {
    return new This_t(model);
  }

  static void construct(This_t *addr)
  {
    new (addr) This_t();
  }

  static void destruct(This_t *)
  {
    PInsist(0,"TrivialElementProperties<T>::destruct(addr) not allowed!");
  }
};

template <class T>
struct TrivialElementProperties : public TrivialElementPropertiesBase<T>
{
  enum { basicType = false };
};


/**
 * Basic types are the C basic types. This is the same as
 * TrivialElementProperties, but with basicType == true;
 */

template <class T>
struct BasicTypeProperties : public TrivialElementPropertiesBase<T>
{
  enum { basicType = true };
};
  

/**
 * Classes that have shallow copy semantics and "makeOwnCopy" methods
 * can specialize ElementProperties<T> by simply inheriting from this.
 */

template <class T>
struct MakeOwnCopyProperties
{
  typedef T This_t;

  enum { hasTrivialDefaultConstructor = false };

  enum { hasTrivialDestructor = false };

  enum { concrete = false };

  enum { basicType = false };

  static void construct(This_t * addr)
  {
    new (addr) This_t;
    addr->makeOwnCopy();
  }

  static void construct(This_t * addr, const This_t & model)
  {
    new (addr) This_t(model);
    addr->makeOwnCopy();
  }

  static This_t * clone(const This_t &model)
  {
    This_t * temp = new This_t(model);
    temp->makeOwnCopy();
    return temp;
  }

  static void destruct(This_t * addr)
  {
    addr->~This_t();
  }

};


//-----------------------------------------------------------------------------
// Specializations for standard C++ types, which do not do initialization
// in their "default constructors".

template <>
struct ElementProperties<bool>   : public BasicTypeProperties<bool>
{ };

template <>
struct ElementProperties<char>   : public BasicTypeProperties<char>
{ };

template <>
struct ElementProperties<unsigned char>   
  : public BasicTypeProperties<unsigned char>
{ };

template <>
struct ElementProperties<short>  : public BasicTypeProperties<short>
{ };

template <>
struct ElementProperties<unsigned short>
  : public BasicTypeProperties<unsigned short>
{ };

template <>
struct ElementProperties<int>    : public BasicTypeProperties<int>
{ };

template <>
struct ElementProperties<unsigned int>
  : public BasicTypeProperties<unsigned int>
{ };

template <>
struct ElementProperties<long>   : public BasicTypeProperties<long>
{ };

#if POOMA_HAS_LONG_LONG
template <>
struct ElementProperties<long long>   : public BasicTypeProperties<long long>
{ };
#endif

template <>
struct ElementProperties<unsigned long>
  : public BasicTypeProperties<unsigned long>
{ };

template <>
struct ElementProperties<float>  : public BasicTypeProperties<float>
{ };

template <>
struct ElementProperties<double> : public BasicTypeProperties<double>
{ };

namespace std {
  template <class Float> class complex;
}

template <class FloatType>
struct ElementProperties<std::complex<FloatType> > 
  : public TrivialElementProperties<std::complex<FloatType> >
{ };


// } // namespace POOMA
///////////////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_ELEMENTPROPERTIES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ElementProperties.h,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
