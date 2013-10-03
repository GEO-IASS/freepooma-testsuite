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
//   CollectFromContexts<T, Op>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Tulip
 * @brief
 * CollectFromContext encapsulates functionality like MPI_Gather.
 */

#ifndef POOMA_MESSAGING_COLLECTFROMCONTEXTS_H
#define POOMA_MESSAGING_COLLECTFROMCONTEXTS_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tulip/Messaging.h"
#include "Utilities/PAssert.h"

#include <vector>


#if !POOMA_MESSAGING

template<class T>
class CollectFromContexts
{
public:

  CollectFromContexts(const T &val, int context = 0, bool valid = true)
    {
      PAssert(valid);
      PAssert(context == 0);
      value_m = val;
    }
    
  T &operator[](int i)
    {
      PAssert(i == 0);
      return value_m;
    }
    
  T operator[](int i) const
    {
      PAssert(i == 0);
      return value_m;
    }

private:
  T value_m;
};
 

#else // POOMA_MESSAGING


/**
 * This class associates a value with a flag that indicates whether or not
 * it is valid. It takes special care to not read the value if it is invalid.
 * It also encodes the context the value was created on.
 */

template<class T>
class CollectionValue
{
public:

  CollectionValue(bool valid, const T &val)
    : valid_m(valid), context_m(Pooma::context())
    {
      if (valid_m)
	val_m = val;
    }

  CollectionValue(bool valid, const T &val, int con)
    : valid_m(valid), context_m(con)
    {
      if (valid_m)
	val_m = val;
    }

  CollectionValue(const CollectionValue<T> &model)
    {
      valid_m = model.valid();
      context_m = model.context();
      if (valid_m)
	val_m = model.value();
    }

  CollectionValue<T> &operator=(const CollectionValue<T> &rhs)
    {
      if (&rhs != this)
	{
	  valid_m = rhs.valid();
	  context_m = rhs.context();
	  if (valid_m)
	    val_m = rhs.value();
	}

      return *this;
    }

  bool valid() const { return valid_m; }
  int context() const { return context_m; }
  const T &value() const { PAssert(valid()); return val_m; }

private:

  bool valid_m;
  int context_m;
  T val_m;
};


namespace Cheetah {

/**
 * This class is used to serialize CollectionValue<T> objects, taking care
 * not to send invalid values.
 */

template<class T>
class Serialize<CHEETAH, CollectionValue<T> >
{
public:

  static inline int size(const CollectionValue<T> &v)
  {
    int nBytes = Serialize<CHEETAH, bool>::size(v.valid());
    nBytes += Serialize<CHEETAH, int>::size(v.context());
    if (v.valid())
      nBytes += Serialize<CHEETAH, T>::size(v.value());

    return nBytes;
  }

  static inline int pack(const CollectionValue<T> &v, char *buffer)
  {
    int nBytes = Serialize<CHEETAH, bool>::pack(v.valid(), buffer); 
    nBytes += Serialize<CHEETAH, int>::pack(v.context(), buffer + nBytes); 

    if (v.valid())
      {
	nBytes += Serialize<CHEETAH, T>::pack(v.value(), buffer + nBytes); 
      }

    return nBytes;
  }

  static inline int unpack(CollectionValue<T>* &vp, char *buffer)
  {
    bool *pvalid;
    int con, *pcon = &con;
    T val, *pval = &val;

    int nBytes = Serialize<CHEETAH, bool>::unpack(pvalid, buffer);

    nBytes += Serialize<CHEETAH, int>::unpack(pcon, buffer + nBytes); 

    if (*pvalid)
      {
	nBytes += Serialize<CHEETAH, T>::unpack(pval, buffer + nBytes); 
      }

    vp = new CollectionValue<T>(*pvalid, *pval, *pcon);

    if (*pvalid)
      Serialize<CHEETAH, T>::cleanup(pval);

    return nBytes;
  }
  
  static inline void cleanup(CollectionValue<T> *vp)
  {
    delete vp;
  }
};

} // namespace Cheetah


#if POOMA_CHEETAH
/**
 * This struct holds a few static quantities that are shared by all
 * instantiations of CollectFromContexts<T>. In particular, we want to 
 * maintain a running sequence of tags across all instantiations. 
 */

struct CollectFromContextsBase
{
  // We use this as a counter to generate tags.

  static int tagBase_m;
};
#endif

/**
 * This class is used to collect all valid values from all contexts.
 */

template<class T>
class CollectFromContexts
{
  typedef CollectFromContexts<T> This_t;

public:

  // All the work happens in the constructor. If we're on the "to" context,
  // we set up to receive messages from all of the other contexts. The 
  // receive() handler performs the collection from contexts incrementally as
  // we get the messages in. We poll until everything shows up. If we're
  // not on the "to" context, we send our value to there.

  // Things are slightly more complicated by the fact that we don't want
  // to read 'val' unless it is a valid value. Values don't have to
  // valid because not all contexts necessarily contribute to the collection.

#if POOMA_MPI

  CollectFromContexts(const T &val, int toContext = 0, bool valid = true)
    : toContext_m(toContext), data_m(Pooma::contexts())
    {
      typedef Cheetah::Serialize<Cheetah::CHEETAH, CollectionValue<T> > Serialize_t;
      CollectionValue<T> v(valid, val);
      // We need to get at the maximum size we need to transfer per context.
      // With the valid/invalid mechanism we can't use size(v) for this, and
      // for dynamic types like Grid<> we can't use CV<T>(true, T()) either.
      // So for these cases we need to communicate the maximum size needed
      // (but we might be able to optimize this with appropriate type tags).
      int thislength = Serialize_t::size(v);
      thislength = (thislength+7)&~7; // round to qword size
      int length;
      MPI_Allreduce(&thislength, &length, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      char *buffer = new char[length]; 
      char *recvbuffer = NULL;
      if (Pooma::context() == toContext)
        recvbuffer = new char[length*Pooma::contexts()];
      Serialize_t::pack(v, buffer);
      MPI_Gather(buffer, length, MPI_CHAR,
                 recvbuffer, length, MPI_CHAR,
                 toContext, MPI_COMM_WORLD);
      delete[] buffer;
      if (Pooma::context() == toContext) {
        for (int i=0; i<Pooma::contexts(); ++i) {
          CollectionValue<T> *v2;
          Serialize_t::unpack(v2, recvbuffer+i*length);
          if (v2->valid())
	    data_m[i] = v2->value();
          Serialize_t::cleanup(v2);
        }
        delete[] recvbuffer;
      }
    }

#elif POOMA_CHEETAH

  CollectFromContexts(const T &val, int toContext = 0, bool valid = true)
    : toContext_m(toContext), data_m(Pooma::contexts())
    {
      int tagBase = CollectFromContextsBase::tagBase_m;
      CollectFromContextsBase::tagBase_m += Pooma::contexts();

      if (Pooma::context() == toContext)
	{
	  toReceive_m = Pooma::contexts();
	  
	  int fromContext;
	  for (fromContext = 0; fromContext < Pooma::contexts(); ++fromContext)
	    {
	      if (fromContext != toContext)
		{
		  Pooma::collectionHandler()->
		    request(fromContext, tagBase + fromContext, receive, this);
		}
	      else
		{
		  CollectionValue<T> v(valid, val);
		  receive(this, v);
		}
	    }
	  
	  while (toReceive_m != 0)
	    {
	      Pooma::poll();
	    }
	}
      else
	{
	  CollectionValue<T> v(valid, val);

	  Pooma::collectionHandler()->
	    send(toContext, tagBase + Pooma::context(), v);
	}
    }

private:

  // Handler function for cheetah.

  static void receive(This_t *me, CollectionValue<T> &v)
  {
    if (v.valid())
      {
        me->data_m[v.context()] = v.value();
      }

    me->toReceive_m--;
  }

  // The number of messages we're receiving.

  int toReceive_m;
  
#endif

public:

  T &operator[](int i)
    {
      PAssert(Pooma::context() == toContext_m);
      PAssert(i >= 0 and i < Pooma::contexts());
      return data_m[i];
    }
    
  T operator[](int i) const
    {
      PAssert(Pooma::context() == toContext_m);
      PAssert(i >= 0 and i < Pooma::contexts());
      return data_m[i];
    }

private:

  // The actual value we're reducing.

  std::vector<T> data_m;

  // The context we're reducing on.
  
  int toContext_m;

};

#endif // POOMA_MESSAGING

#endif     // POOMA_MESSAGING_COLLECTFROMCONTEXTS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CollectFromContexts.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
