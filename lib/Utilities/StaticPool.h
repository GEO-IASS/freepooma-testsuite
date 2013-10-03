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
 * A class that manages a static Pool in which the block size is
 * a template parameter.
 *
 * If you just create a Pool as a static object in each of many different
 * pooled classes, you end up with potentially a large number of different 
 * pools. In particular, if you pool of expression objects, you will have a
 * different pool for each kind of expression object, which is inefficient
 * because many different expression types will have the same size, and
 * could therefore share a pool.
 *
 * class StaticPool<S> has a static Pool of size S.  Strictly speaking,
 * it has a pool of size S', where S' is rounded up to a multiple of 8 bytes.
 * All the StaticPools that round up to size S' share the same pool.
 *
 * This is done by having the Pool be static data in a base class 
 * RoundedStaticPool<SP>, where SP is S rounded up.
 *
 * Usage: When you want a chunk of memory for an object of type T, you say:
 *
 *     T* p = StaticPool<sizeof(T)>::alloc() 
 *
 * To free that memory you say:
 *
 *     StaticPool<sizeof(T)>::free(p);
 */


#ifndef POOMA_UTILITIES_STATIC_POOL_H
#define POOMA_UTILITIES_STATIC_POOL_H

//-----------------------------------------------------------------------------
// Classes:
// StaticPool
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Utilities/Pool.h"

/**
 * All this does is define alloc and free as static functions,
 * and the static pool itself.
 */

template<int SP>
class RoundedStaticPool
{
public:

  // Get a block of memory.
  static void *alloc() { return pool_s.alloc(); }

  // Return a block of memory.
  static void free(void *p) { pool_s.free(p); }

private:

  // This class stores the pool.
  static Pool pool_s;

  // Forbid construction by making this private.
  RoundedStaticPool() {}

};

//
// Declare the storage for the static pool.
//

template<int SP> Pool RoundedStaticPool<SP>::pool_s(SP);

/**
 * This is a wrapper class on RoundedStaticPool, which just rounds up
 * its input block size and inherits from RoundedStaticPool.
 *
 * It doesn't need to do anything else since it inherits the alloc
 * and free functions.
 */

template<class T>
class StaticPool 
: public RoundedStaticPool<(sizeof(T)%8 ? sizeof(T)+8-(sizeof(T)%8) : sizeof(T)) >
{
private:

  // Forbid construction by making this private:
  StaticPool() {}

};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: StaticPool.h,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
