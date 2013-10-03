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
// RNGValue
// RNGSeed
// RNGAdvance
// RNGAdvanceProxy<RNG>
//-----------------------------------------------------------------------------

#ifndef POOMA_FUNCTIONS_RNG_COMPONENT_H
#define POOMA_FUNCTIONS_RNG_COMPONENT_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Utilities
 * @brief
 * UNSUPPORTED.
 *
 * These component accessors allow you to access values from arrays of random
 * numbers and to advance them.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class T, class Components> struct ComponentAccess;

//-----------------------------------------------------------------------------
//
// RNGValue:
//
// RNGValue lets you access the value of random numbers stored in an array.
//
// Array<2, RNG> a;
//
// b = a.comp(RNGValue());
//
// While a.comp(RNGValue()) can be an Array (not just a Array), it is an
// error to attempt to assign to the values.
//-----------------------------------------------------------------------------

struct RNGValue { };

template<class RNG>
struct ComponentAccess< RNG, RNGValue >
{
  typedef typename RNG::Type_t Element_t;
  typedef Element_t ElementRef_t;

  static inline ElementRef_t indexRef(RNG &rng, const RNGValue &)
  {
    return rng.value();
  }
  
  static inline Element_t index(const RNG &rng, const RNGValue &)
  {
    return rng.value();
  }
};

//-----------------------------------------------------------------------------
//
// RNGAdvance:
//
// RNGAdvance lets you advance the value of random numbers stored in an array.
//
// Array<2, RNG> a;
//
// a.comp(RNGAdvance()) = 5;
//
// (advance the random numbers 5 times at every location in a)
//
//-----------------------------------------------------------------------------

struct RNGAdvance { };

template<class RNG>
struct RNGAdvanceProxy
{ 
  RNGAdvanceProxy(RNG &rng)
    : rng_m(rng)
  { }

  inline RNGAdvanceProxy &operator=(int i)
  {
    rng_m.advance(i);
    return *this;
  }

  inline const RNGAdvanceProxy &operator=(int i) const
  {
    rng_m.advance(i);
    return *this;
  }

  inline operator int() const
  {
    return 0;
  }

  mutable RNG &rng_m;
};

template<class RNG>
struct ComponentAccess< RNG, RNGAdvance >
{
  typedef int Element_t;
  typedef RNGAdvanceProxy<RNG> ElementRef_t;
  
  static inline ElementRef_t indexRef(RNG &rng, const RNGAdvance &)
  {
    return ElementRef_t(rng);
  }
  
  static inline Element_t index(const RNG &rng, const RNGAdvance &)
  {
    return 0;
  }
};

//-----------------------------------------------------------------------------
//
// RNGSeed:
//
// RNGSeed lets you set and get the seeds of random numbers stored in an array.
//
// Array<2, RNG> a;
//
// oldSeed = a.comp(RNGSeed());
// a.comp(RNGSeed()) = newSeed;
//
// (Some RNGs may mess with the seed value when you set it to make it more
// random.  In that case, the values you read can be different from the values
// you assign.)
//-----------------------------------------------------------------------------

struct RNGSeed { };

template<class RNG>
struct RNGSeedProxy
{
  RNGSeedProxy(RNG &rng)
    : rng_m(rng)
  { }

  inline RNGSeedProxy &operator=(long i)
  {
    rng_m.seed(i);
    return *this;
  }

  inline const RNGSeedProxy &operator=(long i) const
  {
    rng_m.seed(i);
    return *this;
  }

  inline operator long() const
  {
    return rng_m.seed();
  }

  RNG &rng_m;
};

template<class RNG>
struct ComponentAccess< RNG, RNGSeed >
{
  typedef int Element_t;
  typedef RNGSeedProxy<RNG> ElementRef_t;
  
  static inline ElementRef_t indexRef(RNG &rng, const RNGSeed &)
  {
    return ElementRef_t(rng);
  }
  
  static inline Element_t index(const RNG &rng, const RNGSeed &)
  {
    return rng.seed();
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_FUNCTIONS_RNG_COMPONENT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RNGComponent.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:49 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
