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
// ParticleBC template definitions for PeriodicBC specialization.
//-----------------------------------------------------------------------------

// include files

#include "Domain/Interval.h"
#include "Array/Array.h"
#include "Evaluator/PatchFunction.h"
#include "Tiny/Vector.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// PeriodicBCFunc functor
// Apply the periodic BC in a patchwise fashion.
//-----------------------------------------------------------------------------

template <class T>
struct PeriodicBCFunc
{
  // Constructor
  PeriodicBCFunc(const T& min, const T& max)
    : min_m(min), max_m(max) { }

  // apply method implements BC
  template <class ArrayPatch1, class ArrayPatch2>
  void apply(const ArrayPatch1& obj, const ArrayPatch2& sub, int)
  {
    PAssert(sub.domain() == obj.domain());
    T diff = max_m - min_m;
    int from = sub.domain()[0].first();
    int to = sub.domain()[0].last();
    for (int i = from; i <= to; ++i)
      {
        // check lower boundary
        if (sub(i) < min_m)
          obj(i) = sub(i) + diff;

        // check upper boundary
        else if (sub(i) > max_m)
          obj(i) = sub(i) - diff;
      }
  }

  // bounds for BC
  T min_m, max_m;
};


//-----------------------------------------------------------------------------
// PeriodicBCFunc functor, specialized for Vector's
// Apply the periodic BC in a patchwise fashion.  Each component is
// treated separately.
//-----------------------------------------------------------------------------

template <int Dim, class T, class E>
struct PeriodicBCFunc< Vector<Dim, T, E> >
{
  typedef Vector<Dim, T, E> Type_t;

  // Constructor
  PeriodicBCFunc(const Type_t& min, const Type_t& max)
    : min_m(min), max_m(max) { }

  // apply method implements BC
  template <class ArrayPatch1, class ArrayPatch2>
  void apply(const ArrayPatch1& obj, const ArrayPatch2& sub, int)
  {
    PAssert(sub.domain() == obj.domain());
    Type_t diff = max_m - min_m;
    int from = sub.domain()[0].first();
    int to = sub.domain()[0].last();
    for (int i = from; i <= to; ++i)
      {
        for (int d = 0; d < Dim; ++d)
          {
            if (sub(i)(d) < min_m(d))
              {
                // check lower boundary
                obj(i)(d) = sub(i)(d) + diff(d);
              }
            else if (sub(i)(d) > max_m(d))
              {
                // check upper boundary
                obj(i)(d) = sub(i)(d) - diff(d);
              }
          }
      }
  }

  // bounds for BC
  Type_t min_m, max_m;
};


//-----------------------------------------------------------------------------
// void applyBoundaryCondition()
// Apply the periodic boundary condition to the subject.
//-----------------------------------------------------------------------------

template <class Subject, class Object, class T>
void
ParticleBC< Subject, Object, PeriodicBC<T> >::
applyBoundaryCondition(int pid)
{
  // get limits of periodic range
  T min = bc_m.min(), max = bc_m.max();

  // loop over local patches and apply BC using functor
  PeriodicBCFunc<T> bcfun(min,max);
  if (pid < 0)
    {
      PatchFunction< PeriodicBCFunc<T>,
                     PatchParticle2<true,false> > patchfun(bcfun);
      patchfun.block(object_m,subject_m);
    }
  else
    {
      bcfun.apply(object_m.patchLocal(pid)(), subject_m.patchLocal(pid)(), pid);
    }
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PeriodicBC.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
