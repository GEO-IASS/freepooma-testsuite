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
// Iota
//
// function:
// iota
//-----------------------------------------------------------------------------

#ifndef POOMA_POOMA_INDICES_H
#define POOMA_POOMA_INDICES_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Pooma
 * @brief
 * iota(domain) is a handy function that returns an Array that contains
 * an array of vectors whose elements correspond to index values.
 *
 * For example, iota(Interval<2>(10,10))(3,4) is Vector<2,int>(3,4).
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tiny/Vector.h"
#include "Engine/IndexFunctionEngine.h"
#include "Array/Array.h"
#include "Pooma/View.h"

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

// Ack!!! Pooma names aren't supposed to be all caps.

class IotaFunctor
{
public:
  inline IotaFunctor() { }
  inline IotaFunctor(const IotaFunctor &) { }
  inline IotaFunctor &operator=(const IotaFunctor &) { return *this; }

  inline
  Vector<1, int> operator()(int i1) const
  {
    return Vector<1, int>(i1);
  }

  inline
  Vector<2, int> operator()(int i1, int i2) const
  {
    return Vector<2, int>(i1, i2);
  }

  inline
  Vector<3, int> operator()(int i1, int i2, int i3) const
  {
    return Vector<3, int>(i1, i2, i3);
  }

};

template<int Dim>
struct Iota
{
  typedef Array<Dim, Vector<Dim, int>, IndexFunction<IotaFunctor> > 
    Iota_t;
  typedef typename ComponentView<Loc<1>, Iota_t>::Type_t Index_t;
};

template<int Dim>
inline
typename Iota<Dim>::Iota_t
iota(const Interval<Dim> &domain)
{
  typedef typename Iota<Dim>::Iota_t Iota_t;
  return Iota_t(domain);
}

template<int Dim>
inline
typename Iota<Dim>::Index_t
iotaIndex(const Interval<Dim> &domain, int i)
{
  typedef typename Iota<Dim>::Iota_t Iota_t;
  return Iota_t(domain).comp(i);
}

inline
Array<1, Vector<1, int>, IndexFunction<IotaFunctor> >
iota(int i1);

inline
Array<1, Vector<1, int>, IndexFunction<IotaFunctor> >
iota(int i1)
{
  return Array<1, Vector<1, int>,
    IndexFunction<IotaFunctor> >(Interval<1>(i1));
}

inline
Array<2, Vector<2, int>, IndexFunction<IotaFunctor> >
iota(int i1, int i2);

inline
Array<2, Vector<2, int>, IndexFunction<IotaFunctor> >
iota(int i1, int i2)
{
  return Array<2, Vector<2, int>,
    IndexFunction<IotaFunctor> >(Interval<2>(i1, i2));
}

inline
Array<3, Vector<3, int>, IndexFunction<IotaFunctor> >
iota(int i1, int i2, int i3);

inline
Array<3, Vector<3, int>, IndexFunction<IotaFunctor> >
iota(int i1, int i2, int i3)
{
  return Array<3, Vector<3, int>,
    IndexFunction<IotaFunctor> >(Interval<3>(i1, i2, i3));
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_POOMA_INDICES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Indices.h,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:17:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
