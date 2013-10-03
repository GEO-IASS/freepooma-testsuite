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

/** @file
 * @ingroup Utilities
 * @brief
 * A mixin class for providing fast new and delete.
 */

#ifndef POOMA_UTILITIES_POOLED_H
#define POOMA_UTILITIES_POOLED_H

// Classes:
// Pooled
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/Pool.h"
#include "Pooma/Configuration.h"


/**
 * This mixin class provides two functions: new and delete.
 * It has a static Pool from which it gets and returns memory.
 *
 * You use this class by inheriting from it like so:
 *
 * class A : public Pooled<A>
 * {
 *   ...
 * };
 *
 * Pooled is templated on the class that inherits from it so that it
 * will know the size of the blocks to request from the pool. 
 *
 * This technique will not be correct for a class B which inherits
 * from A, so Pooled can only be used for classes which are leaves in
 * the inheritance heirarchy.
 */

template<class T>
class Pooled
{
#if POOMA_POOLED
public:

  // Allocate memory by getting it from the pool.
  inline void* operator new(size_t) { return pool_s.alloc(); }

  // Return memory to the pool.
  inline void operator delete(void *p, size_t) { if (p) pool_s.free(p); }

  // Provide a placement version that works.
  inline void* operator new(size_t, void* ptr) { return ptr; }

  // Need a matching delete
#if ! POOMA_NO_PLACEMENT_DELETE
  inline void operator delete(void *, void *) { }
#endif

private:

  // The pool which all objects of type T use.
  static Pool pool_s;
#endif
};


// Declare the pool for objects of type T.
#if POOMA_POOLED
template<class T>
Pool Pooled<T>::pool_s(sizeof(T));
#endif

//////////////////////////////////////////////////////////////////////

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Pooled.h,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
