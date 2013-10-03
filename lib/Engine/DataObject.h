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
// DataObjectApply
// DataObjectRequest<RequestTag>
// LeafFunctor<DataObjectRequest,Leaf>
// GetDataObject
// DataObjectRequest<GetDataObject>
// BlockAffinity
// DataObjectRequest<BlockAffinity>
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_DATAOBJECT_H
#define POOMA_ENGINE_DATAOBJECT_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Engine
 * @brief
 * Traits and tags necessary for dealing with the Smarts data object
 * inside some engines.
 *
 * Engine_t::dataObject is true if the engine has
 * a smarts data object, false otherwise.  Anything that requires information
 * about smarts data objects should query this trait first.
 *
 * The tag needs to be of type
 * DataObjectTag<RequestType>, which satisfy the interface of array message
 * tags and ForEach LeafFunctor tags.  In this file we define two request types,
 * GetDataObject and BlockAffinity to return the data object and the affinity:
 *
 * <PRE>
 * Pooma::DataObject_t obj = forEach(array,DataObjectRequest<GetDataObject>())
 * int affinity = forEach(array,DataObjectRequest<BlockAffinity>())
 * </PRE>
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/PoomaSmarts.h"
#include "Engine/EngineFunctor.h"
#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * DataObjectApply is used to wrap the two cases of whether an object has a
 * smarts data object or not.
 */

template<bool Block>
struct DataObjectApply
{ };

template<>
struct DataObjectApply<false>
{
  // engine has no data object, so return default value
  template<class Engine,class Functor>
  inline static
  typename Functor::Type_t
  apply(const Engine&, const Functor& functor)
  {
    return functor.defaultValue();
  }
};

template<>
struct DataObjectApply<true>
{
  // engine has a data object, so pass request via message
  template<class Engine, class Functor>

  inline static
  typename Functor::Type_t
  apply(const Engine& engine, const Functor& functor)
  {
    return functor(engine.dataObject());
  }
};

/**
 *  This class has two functions.
 *   -# It's a message functor that can be handed to arrays.  Arrays that have
 *      data objects will support the function:
 *      array.dataBlockRequest(datablocktag);
 *      To implement the message function, the engine should call operator() on
 *      the DataObjectRequest with a pointer to the data object.
 *   -# It's a LeafFunctor tag class that allows us to apply smarts data object
 *      to expressions.
 *
 *  In cases where the engine contains an expression, the request will be
 *  passed on to engines in the expression that have data objects.  For example,
 *  if we use DataObjectRequest to request a lock on a stencil engine, it uses
 *  forEach to request the same lock on all the engines in the expression
 *  contained inside the stencil.
 *
 *  To implement a version of DataObjectRequest<Request>, the following
 *  interface is required:
 *   - typedef ... Type_t;
 *   - typedef ... Combine_t;
 *   - Type_t operator()(Pooma::DataObject_t*) const;
 *   - Type_t defaultValue() const;
 *
 *  Type_t is the return type of the functor.  Combine_t is a ForEach combine
 *  tag that will be used if the engine contains an expression.  operator()
 *  lets you compute the return value given a pointer to the actual data
 *  object. defaultValue() returns the value that should be returned when
 *  there is no data object.
 */

template<class RequestType>
class DataObjectRequest
{};

template<class Eng, class RequestType>
struct EngineFunctorDefault<Eng, DataObjectRequest<RequestType> >
{
  enum { hasDataObject = Eng::hasDataObject };

  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;

  static inline
  Type_t apply(const Eng &e, const DataObjectRequest<RequestType> &request)
  {
    return DataObjectApply<hasDataObject>::apply(e, request);
  }
};

/**
 * LeafFunctors for DataObjectRequest.
 * For scalars, we return the default value provided by the functor.
 */

template<class RequestType,class T>
struct EngineFunctorScalar<T, DataObjectRequest<RequestType> >
{
  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;
  inline static
  Type_t apply(const T &, const DataObjectRequest<RequestType> &tag)
  {
    return tag.defaultValue();
  }
};

/**
 * DataObjectRequest<BlockAffinity>
 * Used to get the affinity for an array.
 */

struct BlockAffinity { };


/**
 *   This trivial combiner returns the left-most object in an expression.
 */

struct AffinityCombine
{
  AffinityCombine() { }
  AffinityCombine(const AffinityCombine &) { }
};

template<class Op>
struct Combine2<int, int, Op, AffinityCombine>
{
  typedef int Type_t;
  inline static
  Type_t combine(int a, int b, AffinityCombine)
  {
    return a;
  }
};

template<>
class DataObjectRequest<BlockAffinity>
{
public:
  // This functor returns an affinity
  typedef int Type_t;

  // Affinities combine to return the left-most affinity.  It might make
  // more sense to perform a more intelligent combination.
  // (Currently the affinity for an iterate comes from the LHS. If the LHS
  // has multiple parts, then we're only interested in the leftmost array.)

  typedef AffinityCombine Combine_t;

  DataObjectRequest() { }

  // Just return the DataObject pointer if there is one, otherwise
  // return a null pointer.
  inline Type_t operator()(Pooma::DataObject_t* obj) const
  {
    return obj->affinity();
  }

  inline Type_t defaultValue() const
  {
    return (-1);
  }
};

#endif // POOMA_ENGINE_DATAOBJECT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DataObject.h,v $   $Author: richard $
// $Revision: 1.32 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
