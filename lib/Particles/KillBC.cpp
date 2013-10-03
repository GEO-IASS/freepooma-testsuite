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
// ParticleBC template definitions for KillBC specialization.
//-----------------------------------------------------------------------------

// include files

#include "Domain/Interval.h"
#include "DynamicArray/DynamicArray.h"
#include "Domain/IndirectionList.h"
#include "Evaluator/PatchFunction.h"
#include "Utilities/PAssert.h"
#include "Tiny/Vector.h"


// KillBC functors for general elements and specialized for Vectors.
// This is the general case in which we assume the Object to be a
// Particles object.

//-----------------------------------------------------------------------------
// KillBCFunc functor
// Apply the kill BC in a patchwise fashion.
// The object of the KillBCFunc functor *must* be a Particles object
// or a DynamicArray (i.e., something that can perform destroys).
//-----------------------------------------------------------------------------

template <class T, class Object>
struct KillBCFunc
{
  // Typedefs
  typedef Object Object_t;

  // Constructor
  KillBCFunc(const T& min, const T& max, Object_t& object)
    : min_m(min), max_m(max), object_m(object)
  {
  }

  // apply method implements BC
  template <class ArrayPatch>
  void apply(const ArrayPatch& sub, int node)
  {
    // store kill list as a DynamicBrick
    Array<1,int,Brick> killlist(sub.domain().size());
    
    // loop over elements in this patch and apply BC
    int from = sub.domain()[0].first();
    int to = sub.domain()[0].last();
    int killed = 0;
    for (int i = from; i <= to; ++i)
      {
        // check both boundaries
        if (sub.read(i) < min_m || sub.read(i) > max_m)
          {
            // add this particle's local index to the kill list
            killlist(killed++) = i;
          }
      }

    // convert the kill list into an IndirectionList and cache it
    IndirectionList<int> klist(killlist(Interval<1>(killed)));
    object_m.deferredDestroy(klist, node);
  }

  // bounds for BC
  T min_m, max_m;

  // store reference to object so we can invoke its destroy method
  Object_t& object_m;
};

//-----------------------------------------------------------------------------
// KillBCFunc functor specialized for Vectors
// Apply the kill BC in a patchwise fashion.
// Each component is treated separately.
//-----------------------------------------------------------------------------

template <int Dim, class T, class E, class Object>
struct KillBCFunc< Vector<Dim,T,E>, Object >
{
  // Typedefs
  typedef Vector<Dim,T,E> Type_t;
  typedef Object Object_t;

  // Constructor
  KillBCFunc(const Type_t& min, const Type_t& max, Object_t& object)
    : min_m(min), max_m(max), object_m(object)
  {
  }

  // apply method implements BC
  template <class ArrayPatch>
  void apply(const ArrayPatch& sub, int node)
  {
    // store kill list as a DynamicBrick
    Array<1,int,Brick> killlist(sub.domain().size());
    
    // loop over elements in this patch and apply BC
    int from = sub.domain()[0].first();
    int to = sub.domain()[0].last();
    int killed = 0;
    for (int i = from; i <= to; ++i)
      {
        bool wasKilled = false;
        int d = 0;
        while (d < Dim && !wasKilled)
          {
            // check both boundaries
            if (sub.read(i)(d) < min_m(d) || sub.read(i)(d) > max_m(d))
              {
                // add this particle's local index to the kill list
                killlist(killed++) = i;
                wasKilled = true;
              }
            ++d;
          }
      }

    // convert the kill list into an IndirectionList and cache it
    IndirectionList<int> klist(killlist(Interval<1>(killed)));
    object_m.deferredDestroy(klist, node);
  }

  // bounds for BC
  Type_t min_m, max_m;

  // store reference to object so we can invoke its destroy method
  Object_t& object_m;
};


// KillBC functors for general elements and specialized for Vectors.
// This is the special case in which the Object is a DynamicArray.

//-----------------------------------------------------------------------------
// KillBCFunc functor
// Apply the kill BC in a patchwise fashion.
// The object of the KillBCFunc functor *must* be a Particles object
// or a DynamicArray (i.e., something that can perform destroys).
//-----------------------------------------------------------------------------

template <class T1, class T2, class E>
struct KillBCFunc< T1, DynamicArray<T2,E> >
{
  // Typedefs
  typedef DynamicArray<T2,E> Object_t;

  // Constructor
  KillBCFunc(const T1& min, const T1& max, Object_t& object)
    : min_m(min), max_m(max), object_m(object)
  {
  }

  // apply method implements BC
  template <class ArrayPatch>
  void apply(const ArrayPatch& sub, int node)
  {
    // store kill list as a DynamicBrick
    Array<1,int,Brick> killlist(sub.domain().size());
    
    // loop over elements in this patch and apply BC
    int from = sub.domain()[0].first();
    int to = sub.domain()[0].last();
    int killed = 0;
    for (int i = from; i <= to; ++i)
      {
        // check both boundaries
        if (sub.read(i) < min_m || sub.read(i) > max_m)
          {
            // add this particle's local index to the kill list
            killlist(killed++) = i;
          }
      }

    // convert the kill list into an IndirectionList and cache it
    IndirectionList<int> klist(killlist(Interval<1>(killed)));
    object_m.destroy(klist, node);
  }

  // bounds for BC
  T1 min_m, max_m;

  // store reference to object so we can invoke its destroy method
  Object_t& object_m;
};

//-----------------------------------------------------------------------------
// KillBCFunc functor specialized for Vectors
// Apply the kill BC in a patchwise fashion.
// Each component is treated separately.
//-----------------------------------------------------------------------------

template <int Dim, class T1, class E1, class T2, class E2>
struct KillBCFunc< Vector<Dim,T1,E1>, DynamicArray<T2,E2> >
{
  // Typedefs
  typedef Vector<Dim,T1,E1> Type_t;
  typedef DynamicArray<T2,E2> Object_t;

  // Constructor
  KillBCFunc(const Type_t& min, const Type_t& max, Object_t& object)
    : min_m(min), max_m(max), object_m(object)
  {
  }

  // apply method implements BC
  template <class ArrayPatch>
  void apply(const ArrayPatch& sub, int node)
  {
    // store kill list as a DynamicBrick
    Array<1,int,Brick> killlist(sub.domain().size());
    
    // loop over elements in this patch and apply BC
    int from = sub.domain()[0].first();
    int to = sub.domain()[0].last();
    int killed = 0;
    for (int i = from; i <= to; ++i)
      {
        bool wasKilled = false;
        int d = 0;
        while (d < Dim && !wasKilled)
          {
            // check both boundaries
            if (sub(i)(d) < min_m(d) || sub(i)(d) > max_m(d))
              {
                // add this particle's local index to the kill list
                killlist(killed++) = i;
                wasKilled = true;
              }
            ++d;
          }
      }
    // convert the kill list into an IndirectionList and cache it
    IndirectionList<int> klist(killlist(Interval<1>(killed)));
    object_m.destroy(klist, node);
  }

  // bounds for BC
  Type_t min_m, max_m;

  // store reference to object so we can invoke its destroy method
  Object_t& object_m;
};


//-----------------------------------------------------------------------------
// void applyBoundaryCondition()
// Apply the "kill" boundary condition to the subject and object.
//-----------------------------------------------------------------------------

template <class Subject, class Object, class T>
void
ParticleBC< Subject, Object, KillBC<T> >::
applyBoundaryCondition(int pid)
{
  // get limits of KillBC range
  T min = bc_m.min(), max = bc_m.max();

  // loop over local patches and apply the BC using functor...
  KillBCFunc<T,Object> bcfun(min,max,object_m);
  if (pid < 0)
    {
      PatchFunction< KillBCFunc<T,Object>,
                     PatchParticle1<false> > patchfun(bcfun);
      patchfun.block(subject_m);
    }
  else
    {
      bcfun.apply(subject_m.patchLocal(pid)(), pid);
    }
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: KillBC.cpp,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
