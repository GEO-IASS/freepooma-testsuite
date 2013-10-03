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
// RNGSimple
//-----------------------------------------------------------------------------

#ifndef POOMA_FUNCTIONS_RNG_SIMPLE_H
#define POOMA_FUNCTIONS_RNG_SIMPLE_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Utilities
 * @brief
 * UNSUPPORTED.
 *
 * Simple class that implements random number generator a la
 * Numerical Recipes, in the range [0...1].
 */

class RNGSimple
{
public:

  // return type
  typedef double Type_t;

public:

  // default constructor

  RNGSimple(int adv = 0)
    : currentRand_m(randShift + 1)
  {
    advance(adv);
  }

  // copy constructor

  RNGSimple(const RNGSimple& rng)
    : currentRand_m(rng.currentRand_m)
  { }

  // destructor

  ~RNGSimple()
  { }

  //   advance indicates number of times to advance random number source

  inline void advance(int advance = 1)
  {
    int iadv;
    for (iadv = 0; iadv < advance; iadv++)
    {
      currentRand_m = (currentRand_m * randMultiplier + randShift) % randModulus;
    }
  }

  // set seed to user-specified value, plus shift to ensure it is large

  inline void seed(long seed)
  {
    currentRand_m = (seed + randShift) % randModulus;
  }

  inline long seed() const
  {
    return currentRand_m;
  }

  // return the next pseudo-random number (from 0 ... 1)

  inline Type_t value() const
  {
    return ( Type_t(currentRand_m) / Type_t(randModulus) );
  }

  // return the period of the RNG

  static Type_t period()
  {
    return Type_t(randModulus);
  }

private:

  mutable long currentRand_m;

  enum { randModulus    = 714025 };
  enum { randMultiplier = 1366 };
  enum { randShift      = 150889 };
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_FUNCTIONS_RNG_SIMPLE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RNGSimple.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:49 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
