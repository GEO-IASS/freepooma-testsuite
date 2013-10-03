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
// TagGenerator
//-----------------------------------------------------------------------------

#ifndef POOMA_TULIP_MESSAGING_H
#define POOMA_TULIP_MESSAGING_H

/** @file
 * @ingroup Tulip
 * @brief
 * Functions and classes needed to support interaction with the cheetah
 * messaging library.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"

#if POOMA_MPI
# include "Tulip/CheetahSerialize.h"
# include <mpi.h>
#endif

#if POOMA_CHEETAH
# include "Cheetah/Cheetah.h"
#endif

#include <vector>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------


/**
 * There exist situations where pooma needs to generate identifying tags for
 * data being transferred from one context to another.
 * In order to generate the correct data-flow, pooma needs to manage a set of
 * tags that identify bricks being transferred from context to context.
 */

class TagGenerator
{
public:

  TagGenerator()
    : send_m(1), receive_m(1)
  {
    send_m[0]    = 0;
    receive_m[0] = 0;
  }

  TagGenerator(int n)
    : send_m(n), receive_m(n)
  {
    int i;
    for (i = 0; i < n; ++i)
      {
        send_m[i]    = 0;
        receive_m[i] = 0;
      }
  }

  int send(int otherContext)
  {
    int tag = send_m[otherContext];
    send_m[otherContext]++;
    return tag;
  }

  int receive(int otherContext)
  {
    int tag = receive_m[otherContext];
    receive_m[otherContext]++;
    return tag;
  }

private:

  std::vector<int> send_m;
  std::vector<int> receive_m;
};


#if POOMA_MESSAGING

namespace Cheetah {

  /**
   * This class is used to serialize std::vector<T> objects. First the size
   * and then the elements are sent.
   */

  template<class T>
  class Serialize<CHEETAH, std::vector<T> >
  {
  public:

    static int size(const std::vector<T> &v)
    {
      int nBytes = Serialize<CHEETAH, Size_t>::size(v.size());
      for (int i = 0; i < v.size(); i++)
	nBytes += Serialize<CHEETAH, T>::size(v[i]);
      
      return nBytes;
    }
    
    static int pack(const std::vector<T> &v, char *buffer)
    {
      int nBytes = Serialize<CHEETAH, Size_t>::pack(v.size(), buffer); 
      for (int i = 0; i < v.size(); i++)
	{
	  T vi = v[i];
	  nBytes += Serialize<CHEETAH, T>::pack(vi, buffer + nBytes); 
	}
      
      return nBytes;
    }
    
    static int unpack(std::vector<T>* &vp, char *buffer)
    {
      Size_t *psize;
      T val, *pval = &val;
      
      int nBytes = Serialize<CHEETAH, Size_t>::unpack(psize, buffer);
      vp = new std::vector<T>;
      vp->reserve(*psize);
      for (int i = 0; i < *psize; i++)
	{
	  nBytes += Serialize<CHEETAH, T>::unpack(pval, buffer + nBytes); 
	  vp->push_back(*pval);
	}

      return nBytes;
    }
    
    static inline void cleanup(std::vector<T> *vp)
    {
      delete vp;
    }
    
  private:
    
    typedef typename std::vector<T>::size_type Size_t;
    
  };
  
} // namespace Cheetah

#endif // #if POOMA_MESSAGING

namespace Pooma {

extern int expectedMessages_g;

#if POOMA_CHEETAH

extern Cheetah::MatchingHandler *collectionHandler_g;
extern Cheetah::MatchingHandler *indexHandler_g;
extern Cheetah::MatchingHandler *reductionHandler_g;
extern Cheetah::MatchingHandler *remoteEngineHandler_g;
extern Cheetah::MatchingHandler *particleSwapHandler_g;


inline Cheetah::MatchingHandler *collectionHandler()
{
  return collectionHandler_g;
}

inline Cheetah::MatchingHandler *indexHandler()
{
  return indexHandler_g;
}

inline Cheetah::MatchingHandler *reductionHandler()
{
  return reductionHandler_g;
}

inline Cheetah::MatchingHandler *remoteEngineHandler()
{
  return remoteEngineHandler_g;
}

inline Cheetah::MatchingHandler *particleSwapHandler()
{
  return particleSwapHandler_g;
}

#endif // #if POOMA_CHEETAH

void initializeCheetahHelpers(int contexts);
void finalizeCheetahHelpers();

int sendTag(int context);

int receiveTag(int context);

inline void addIncomingMessage()
{
  expectedMessages_g++;
}

inline void gotIncomingMessage()
{
  expectedMessages_g--;
}

inline bool incomingMessages()
{
  return (expectedMessages_g > 0);
}

}

#endif     // POOMA_TULIP_MESSAGING_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Messaging.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
