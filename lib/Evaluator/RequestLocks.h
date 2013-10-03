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

//-----------------------------------------------------------------------------
// Class:
// WriteRequest
// DataObjectRequest<WriteRequest>
// ReadRequest
// DataObjectRequest<ReadRequest>
// WriteRelease
// DataObjectRequest<WriteRelease>
// ReadRelease
// DataObjectRequest<ReadRelease>
// CountBlocks
// DataObjectRequest<CountBlocks>
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_REQUESTLOCKS_H
#define POOMA_EVALUATOR_REQUESTLOCKS_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * Classes and functors that are necessary to request locks on an expression. 
 *
 * DataObjectRequest<RequestType> is defined for 4 request types here:
 * - WriteRequest: request a write lock
 * - ReadRequest: request a read lock
 * - WriteRelease: release a write lock
 * - ReadRelease: release a read lock
 * - CountBlocks: count the number of data objects in an array
 *
 * DataObjectRequest is defined in Engine/DataObject.h.  It acts as a PETE
 * functor tag and as a tag that can be used by array's message function.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Threads/PoomaSmarts.h"
#include "Engine/DataObject.h"
#include "Engine/EngineFunctor.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Request type tags:
//-----------------------------------------------------------------------------

struct WriteRequest {};
struct ReadRequest {};
struct WriteRelease {};
struct ReadRelease {};
struct CountBlocks {};

//-----------------------------------------------------------------------------
// DataObjectRequest<WriteRequest>
// Used to request write locks.
// Also stores a pointer to the data object so that we can compare to it
// while requesting read locks.  If the same object appears on the left and
// right hand sides, then requesting both read and write locks would lead to
// a deadlock.  We just make the write lock.
//-----------------------------------------------------------------------------

template<>
class DataObjectRequest<WriteRequest>
{
public:
  // return an int (not used, but we need a dummy value in ForEach)
  typedef int Type_t;

  // We don't use the return so NullCombine it.
  typedef NullCombine Combine_t;

  // Constructor takes the iterate which is used in making the request.

  DataObjectRequest(Pooma::Iterate_t& iterate)
    : iterate_m(iterate),
      lhs1_m(NULL), lhs2_m(NULL)
  { }

  // If there's a data object then request the lock.
  inline Type_t operator()(Pooma::DataObject_t* obj) const
  {
    if ((obj != lhs1_m) && (obj != lhs2_m))
    {
      if (lhs1_m == NULL)
      {
	lhs1_m = obj;
      }
      else
      {
	if (lhs2_m == NULL)
	{
	  lhs2_m = obj;
	}
	else
	{
	  PAssert(false);
	}
      }
      obj->request(iterate_m,Pooma::SmartsTag_t::Write);
    }
    return 0;
  }

  // If no object then just return.
  inline Type_t defaultValue() const
  {
    return 0;
  }

  // access the data object
  Pooma::DataObject_t* dataObject1() const { return lhs1_m; }
  Pooma::DataObject_t* dataObject2() const { return lhs2_m; }

  // access the iterate
  Pooma::Iterate_t& iterate() const { return iterate_m; }

private:
  // Iterate that we're requesting locks for
  mutable Pooma::Iterate_t& iterate_m;

  // cached pointer to the data object(s) written to
  // (used when requesting reads).

  mutable Pooma::DataObject_t* lhs1_m;
  mutable Pooma::DataObject_t* lhs2_m;
};

//-----------------------------------------------------------------------------
// DataObjectRequest<ReadRequest>
// Used to request read locks.
// Checks to see if the object was on the lhs, and if so doesn't request lock.
//-----------------------------------------------------------------------------

template<>
class DataObjectRequest<ReadRequest>
{
public:
  // return an int (not used, but we need a dummy value in ForEach)
  typedef int Type_t;

  // We don't use the return so NullCombine it.
  typedef NullCombine Combine_t;

  // Constructor takes the iterate and a pointer to the lhs data object
  DataObjectRequest(const DataObjectRequest<WriteRequest>& write)
    : iterate_m(write.iterate()),
      lhs1_m(write.dataObject1()),
      lhs2_m(write.dataObject2())
  { }

  // Constructor takes the iterate which is used in making the request.

  DataObjectRequest(Pooma::Iterate_t& iterate)
    : iterate_m(iterate),
      lhs1_m(NULL), lhs2_m(NULL)
  { }

  // If there's a data object, then we first compare it to the lhs
  // data object.  If it's the same object, the we
  // don't need to request this lock.
  inline Type_t operator()(Pooma::DataObject_t* obj) const
  {
    if ((lhs1_m != obj) && (lhs2_m != obj))
    {
      obj->request(iterate_m, Pooma::SmartsTag_t::Read);
    }
    return 0;
  }

  inline Type_t defaultValue() const
  {
    return 0;
  }

private:
  Pooma::Iterate_t& iterate_m;
  Pooma::DataObject_t* lhs1_m;
  Pooma::DataObject_t* lhs2_m;
};

//-----------------------------------------------------------------------------
// DataObjectRequest<WriteRelease>
// Used to release write locks.
// As with read locks, we cache the lhs data object pointer.
//-----------------------------------------------------------------------------

template<>
class DataObjectRequest<WriteRelease>
{
public:
  typedef int Type_t;
  typedef NullCombine Combine_t;

  DataObjectRequest()
    : lhs1_m(NULL), lhs2_m(NULL)
  { }

  inline Type_t operator()(Pooma::DataObject_t* obj) const
  {
    if ((obj != lhs1_m) && (obj != lhs2_m))
    {
      if (lhs1_m == NULL)
      {
	lhs1_m = obj;
      }
      else
      {
	if (lhs2_m == NULL)
	{
	  lhs2_m = obj;
	}
	else
	{
	  PAssert(false);
	}
      }
      obj->release(Pooma::SmartsTag_t::Write);
    }
    return 0;
  }

  inline Type_t defaultValue() const
  {
    return 0;
  }

  Pooma::DataObject_t* dataObject1() const { return lhs1_m; }
  Pooma::DataObject_t* dataObject2() const { return lhs2_m; }

private:
  mutable Pooma::DataObject_t* lhs1_m;
  mutable Pooma::DataObject_t* lhs2_m;
};

//-----------------------------------------------------------------------------
// DataObjectRequest<ReadRelease>
// Used to release read locks.
// As with requesting locks, we check if the data object was on the lhs.
//-----------------------------------------------------------------------------

template<>
class DataObjectRequest<ReadRelease>
{
public:
  typedef int Type_t;
  typedef NullCombine Combine_t;

  DataObjectRequest()
    : lhs1_m(NULL), lhs2_m(NULL)
  { }

  DataObjectRequest(const DataObjectRequest<WriteRelease>& write)
    : lhs1_m(write.dataObject1()), lhs2_m(write.dataObject2())
  { }

  inline Type_t operator()(Pooma::DataObject_t* obj) const
  {
    if ((lhs1_m != obj) && (lhs2_m != obj))
    {
      obj->release(Pooma::SmartsTag_t::Read);
    }
    return 0;
  }

  inline Type_t defaultValue() const
  {
    return 0;
  }

private:
  Pooma::DataObject_t* lhs1_m;
  Pooma::DataObject_t* lhs2_m;
};

//-----------------------------------------------------------------------------
// DataObjectRequest<CountBlocks>
// Used to count the total number of locks we're going to make.
//-----------------------------------------------------------------------------

template<>
class DataObjectRequest<CountBlocks>
{
public:
  // The result of the count is an int.
  typedef int Type_t;

  // Use SumCombine to get the total # of blocks
  typedef SumCombine Combine_t;

  DataObjectRequest() { }

  // if there's a data object then count 1
  inline Type_t operator()(Pooma::DataObject_t*) const
  {
    return 1;
  }

  // otherwise count 0
  inline Type_t defaultValue() const
  {
    return 0;
  }
};


//-----------------------------------------------------------------------------
// Code for ExpressionApply is provided here, since the lock requests can
// be implemented using expressionApply (the affinity access still needs to
// go through engineFunctor).
//-----------------------------------------------------------------------------

template<class Eng, class Tag>
struct DefaultExpressionApply<Eng, DataObjectRequest<Tag> >
{
  enum { hasDataObject = Eng::hasDataObject };

  inline static
  int apply(const Eng &e,
	    const ExpressionApply<DataObjectRequest<Tag> > &request)
  {
    DataObjectApply<hasDataObject>::apply(e, request.tag());
    return 0;
  }
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_REQUESTLOCKS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RequestLocks.h,v $   $Author: richard $
// $Revision: 1.26 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
