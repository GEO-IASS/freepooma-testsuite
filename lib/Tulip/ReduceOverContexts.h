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
//   ReduceOverContexts<T, Op>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Tulip
 * @brief
 * ReduceOverContexts encapsulates functionality like MPI_Reduce
 * and MPI_Allreduce by means of the ReduceOverContexts::broadcast() method.
 */

#ifndef POOMA_CHEETAH_REDUCEOVERCONTEXTS_H
#define POOMA_CHEETAH_REDUCEOVERCONTEXTS_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tulip/Messaging.h"
#include "Evaluator/OpMask.h"
#include "Tulip/RemoteProxy.h"
#include "Evaluator/OpMask.h"


/**
 * This class associates a value with a flag that indicates whether or not
 * it is valid. It takes special care to not read the value if it is invalid.
 */

template<class T>
class ReductionValue
{
public:

  ReductionValue(bool valid, const T &val)
    : valid_m(valid)
    {
      if (valid_m)
	val_m = val;
    }

  ReductionValue(const ReductionValue<T> &model)
    {
      valid_m = model.valid();
      if (valid_m)
	val_m = model.value();
    }

  ReductionValue<T> &operator=(const ReductionValue<T> &rhs)
    {
      if (&rhs != this)
	{
	  valid_m = rhs.valid();
	  if (valid_m)
	    val_m = rhs.value();
	}

      return *this;
    }

  bool valid() const { return valid_m; }
  const T &value() const { PAssert(valid()); return val_m; }
  T &value() { PAssert(valid()); return val_m; }

private:

  bool valid_m;
  T val_m;
};


#if POOMA_MESSAGING

namespace Cheetah {

/**
 * This class is used to serialize ReductionValue<T> objects, taking care
 * not to send invalid values.
 */

template<class T>
class Serialize<CHEETAH, ReductionValue<T> >
{
public:

  static inline int size(const ReductionValue<T> &v)
  {
    int nBytes = Serialize<CHEETAH, bool>::size(v.valid());
    if (v.valid())
      nBytes += Serialize<CHEETAH, T>::size(v.value());

    return nBytes;
  }

  static inline int pack(const ReductionValue<T> &v, char *buffer)
  {
    int nBytes = Serialize<CHEETAH, bool>::pack(v.valid(), buffer); 

    if (v.valid())
      {
	nBytes += Serialize<CHEETAH, T>::pack(v.value(), buffer + nBytes); 
      }

    return nBytes;
  }

  static inline int unpack(ReductionValue<T>* &vp, char *buffer)
  {
    bool *pvalid;
    T val, *pval = &val;

    int nBytes = Serialize<CHEETAH, bool>::unpack(pvalid, buffer);

    if (*pvalid)
      {
	nBytes += Serialize<CHEETAH, T>::unpack(pval, buffer + nBytes); 
      }

    vp = new ReductionValue<T>(*pvalid, *pval);

    if (*pvalid)
      Serialize<CHEETAH, T>::cleanup(pval);

    return nBytes;
  }
  
  static inline void cleanup(ReductionValue<T> *vp)
  {
    delete vp;
  }
};

} // namespace Cheetah

#endif


#if POOMA_CHEETAH
/**
 * This struct holds a few static quantities that are shared by all
 * instantiations of ReduceOverContexts<T>. In particular, we want to 
 * maintain a running sequence of tags across all instantiations. 
 */

struct ReduceOverContextsBase
{
  // We use this as a counter to generate tags.

  static int tagBase_m;
};
#endif

/**
 * This class is used to implement the final reduction over contexts used
 * in Reduction<RemoteMultiPatchTag>::evaluate().
 */

template<class T, class ReductionOp>
class ReduceOverContexts
{
  typedef ReduceOverContexts<T, ReductionOp> This_t;

public:

  // All the work happens in the constructor. If we're on the "to" context,
  // we set up to receive messages from all of the other contexts. The 
  // receive() handler performs the reduction over contexts incrementally as
  // we get the messages in. We poll until everything shows up. If we're
  // not on the "to" context, we send our value to there.

  // Things are slightly more complicated by the fact that we don't want
  // to read 'val' unless it is a valid value. Values don't have to
  // valid because not all contexts necessarily contribute to the reduction.

#if POOMA_MESSAGING

#if POOMA_CHEETAH

  ReduceOverContexts(const T &val, int toContext = 0, bool valid = true)
    : valid_m(false), toContext_m(toContext)
    {
      int tagBase = ReduceOverContextsBase::tagBase_m;
      ReduceOverContextsBase::tagBase_m += Pooma::contexts();

      if (Pooma::context() == toContext)
	{
	  toReceive_m = Pooma::contexts();
	  int fromContext;
	  for (fromContext = 0; fromContext < Pooma::contexts(); ++fromContext)
	    {
	      if (fromContext != toContext)
		{
		  Pooma::reductionHandler()->
		    request(fromContext, tagBase + fromContext, receive, this);
		}
	      else
		{
		  ReductionValue<T> v(valid, val);
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
	  ReductionValue<T> v(valid, val);

	  Pooma::reductionHandler()->
	    send(toContext, tagBase + Pooma::context(), v);
	}
    }

#elif POOMA_MPI

  ReduceOverContexts(const T &val, int toContext = 0, bool valid = true)
    : toContext_m(toContext)
    {
      typedef Cheetah::Serialize<Cheetah::CHEETAH, ReductionValue<T> > Serialize_t;
      ReductionValue<T> v(valid, val);
      // invalid size is different (doh!), so use some default for size
      // strictly speaking this is incorrect, too (see CollectOverContexts),
      // but we might not have reduction ops over dynamic sized objects...
      int length = Serialize_t::size(ReductionValue<T>(true, T()));
      length = (length+7)&~7; // round to qword size
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
          if (i == toContext) // this we already have in v
            continue;
          ReductionValue<T> *v2;
          Serialize_t::unpack(v2, recvbuffer+i*length);
          if (v2->valid()) {
	    if (!v.valid())
	      v = *v2;
	    else
	      ReductionOp()(v.value(), v2->value());
          }
          Serialize_t::cleanup(v2);
        }
        delete[] recvbuffer;
        if (v.valid())
          value_m = v.value();
      }
   }

#endif

  // FIXME: with a different API we could use MPI_AllGather here...
  void broadcast(T &val)
    {
      RemoteProxy<T> broadcast(value_m, toContext_m);
      val = broadcast;
    }
      
#else

  ReduceOverContexts(const T &val, int = 0, bool valid = true)
    {
      PAssert(valid);
      value_m = val;
    }
    
  void broadcast(T &val)
    {
      val = value_m;
    }

#endif // POOMA_MESSAGING

  inline operator T() const { return value_m; }

private:

  // Handler function for cheetah.

  static void receive(This_t *me, ReductionValue<T> &v)
  {
    if (v.valid())
      {
	if (!me->valid_m)
	  {
	    me->value_m = v.value();
	    me->valid_m = true;
	  }
	else
	  {
	    ReductionOp()(me->value_m, v.value());
	  }
      }

    me->toReceive_m--;
  }
 
  // The actual value we're reducing.

  T value_m;

  // If it's valid.

  bool valid_m;

  // The number of messages we're receiving.

  int toReceive_m;
  
  // The context we're reducing on.
  
  int toContext_m;
};

#endif     // POOMA_CHEETAH_REDUCEOVERCONTEXTS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ReduceOverContexts.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
