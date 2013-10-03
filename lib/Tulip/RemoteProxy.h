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
// RemoteProxy<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Tulip
 * @brief
 * This is like MPI_Bcast.
 *
 * It moves a value from one context to all others.
 * Special about this is that assigns to a RemoteProxy object
 * on the owning context is performed to the underlying data,
 * while on the remote contexts it is just done to the proxy
 * object.
 */

#ifndef POOMA_CHEETAH_REMOTE_PROXY_H
#define POOMA_CHEETAH_REMOTE_PROXY_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tulip/Messaging.h"
#include "Domain/Loc.h"
#include "Tiny/Vector.h"
#include "Functions/ComponentAccess.h"


// For Cheetah support we need to mark more types not delegate.

#if POOMA_MESSAGING
namespace Cheetah {

  /**
   * CHEETAH specializations for STL strings
   */

  template<>
  class Serialize< ::Cheetah::CHEETAH, std::string>
  { 
  public:
    static inline int size(const std::string& str) 
    {
      return Serialize<ARRAY, char>::size(0, str.length());
    }
    
    static int pack(const std::string &str, char* buffer)
    {
      return Serialize<ARRAY, char>::pack(str.data(), buffer, str.length());
    }

    static int unpack(std::string* &str, char* buffer)
    {
      char* ptr;
      int size;

      int n = Serialize<ARRAY, char>::unpack(ptr, buffer, size);
      str = new std::string(ptr, size);
      return n;
    }

    static void cleanup(std::string* str) { delete str; }
  };

} // namespace Cheetah
#endif


#if POOMA_CHEETAH
struct RemoteProxyBase
{
  /// If we need a remote value, then this flag lets us know when it's
  /// ready.  This value is static because it is used to block the parse
  /// thread until the data is received.

  static bool ready_m;

  /// We only need one tag for all the remote proxies.  Perhaps this could
  /// be packaged with the handler for remote proxies.

  static int tag_m;
};
#endif

/**
 * This class is the return type of the remote brick engine operator().
 * We need an object that lets us assign to data on this context, but that
 * can also contain the data that came from another context, and that prevents
 * you from writing to that data.
 *
 * Outstanding:
 *  -# maybe need a copy constructor etc.  right now it's ok as long as proxies
 *     don't live too long.
 *  -# other operator functions.
 *
 * Usage:
 *
 * A RemoteProxy must be constructed with a value and the context that the
 * value belongs to.
 */

template<class T>
class RemoteProxy
{
public:
  typedef RemoteProxy<T> This_t;

  /// All the work happens in the remote proxy constructor.
  /// If we're on the right context, we store a pointer to the
  /// value and broadcast the value to the other contexts.
  /// Otherwise we receive the value from the owning context.

#if POOMA_CHEETAH

  RemoteProxy(T &val, int owningContext = 0)
  {
    int tag = RemoteProxyBase::tag_m++;
    if (Pooma::context() == owningContext)
    {
      value_m = &val;

      int toContext;
      for (toContext = 0; toContext < Pooma::contexts(); ++toContext)
      {
	if (toContext != Pooma::context())
	{
	  Pooma::indexHandler()->sendWith(Cheetah::CHEETAH(), toContext, tag, val);
	}
      }
    }
    else
    {
      storedValue_m = val;
      value_m = &storedValue_m;

      RemoteProxyBase::ready_m = false;

      Pooma::indexHandler()->requestWith(Cheetah::CHEETAH(), owningContext, tag,
				         This_t::receive, this);

      while (!RemoteProxyBase::ready_m)
      {
	Pooma::poll();
      }
    }
  }

private:
  // Handler function for Cheetah.

  static void receive(This_t *me, T &value)
  {
    me->storedValue_m = value;
    RemoteProxyBase::ready_m = true;
  }

public:
#elif POOMA_MPI

  RemoteProxy(T &val, int owningContext = 0)
  {
    typedef Cheetah::Serialize<Cheetah::CHEETAH, T> Serialize_t;
    int length = Serialize_t::size(val);
    // Only the owningContext can possibly know the actual length for
    // types like std::vector<>. Maybe we can conditionalize this extra
    // communication on a tag field in the Cheetah::Serialize type.
    MPI_Bcast(&length, 1, MPI_INT, owningContext, MPI_COMM_WORLD);
    char *buffer = new char[length];
    if (Pooma::context() == owningContext)
      Serialize_t::pack(val, buffer);
    MPI_Bcast(buffer, length, MPI_CHAR, owningContext, MPI_COMM_WORLD);
    if (Pooma::context() == owningContext) {
      value_m = &val;
    } else {
      T *nval;
      Serialize_t::unpack(nval, buffer);
      storedValue_m = *nval;
      value_m = &storedValue_m;
      Serialize_t::cleanup(nval);
    }
    delete[] buffer;
  }

#else

  RemoteProxy(T &val, int owningContext = 0)
  {
    if (Pooma::context() == owningContext)
    {
      value_m = &val;
    }
  }

#endif

  RemoteProxy(const RemoteProxy<T> &s)
  {
    if (s.value_m != &s.storedValue_m)
    {
      value_m = s.value_m;
    }
    else
    {
      storedValue_m = s.value();
      value_m = &storedValue_m;
    }
  }

  inline
  operator T() const { return *value_m; }

  inline
  T &value() { return *value_m; }

  inline
  const T &value() const { return *value_m; }

  template<class S>
  inline
  RemoteProxy<T> &operator=(const S &s)
  {
    *value_m = s;
    return *this;
  }

  template<class S>
  inline
  RemoteProxy<T> &operator=(const RemoteProxy<S> &s)
  {
    *value_m = s.value();
    return *this;
  }

  inline
  RemoteProxy<T> &operator=(const RemoteProxy<T> &s)
  {
    *value_m = s.value();
    return *this;
  }

  inline
  typename ComponentAccess<T, Loc<1> >::Element_t
  operator()(int i) const
  {
    return ComponentAccess<T, Loc<1> >::index(value(), Loc<1>(i));
  }

  inline
  typename ComponentAccess<T, Loc<1> >::ElementRef_t
  operator()(int i)
  {
    return ComponentAccess<T, Loc<1> >::indexRef(value(), Loc<1>(i));
  }

private:
 
 // Pointer to the actual value represented by this proxy.

  T *value_m;

  // If we get a remote value, we store it here.

  T storedValue_m;

};

/// Some mathematical operations.
/// These probably need to be improved to promote the types.

template<class T, class S>
inline T
operator*(const RemoteProxy<T> &t, const S &s)
{
  return t.value() * s;
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CHEETAH_REMOTE_PROXY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RemoteProxy.h,v $   $Author: richard $
// $Revision: 1.22 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
