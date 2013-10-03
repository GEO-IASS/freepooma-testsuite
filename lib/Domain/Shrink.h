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

#ifndef POOMA_DOMAIN_SHRINK_H
#define POOMA_DOMAIN_SHRINK_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Interval<Dim> shrinking and growing domains.
 *
 * Shrinking and growing domains can be done asymmetrically by one of
 * the shrinkLeft, shrinkRight, growLeft or growRight variants.  Symmetric
 * shrinking and growing can be done using the overloaded shrink function.
 *
 * Examples:
 * - shrinkRight(Interval<1>(0, 4), 1) == Interval<1>(0, 3)
 * - growLeft(Interval<1>(0, 4), 1) == Interval<1>(-1, 4)
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Loc.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/// Deprecated. Use shrinkRight().

template<int Dim>
Interval<Dim> &
shrinkRightInPlace(Interval<Dim> &dom, const Loc<Dim> &s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() - s[d].first();
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use shrinkRight().

template<int Dim>
Interval<Dim> &
shrinkRightInPlace(Interval<Dim> &dom, int s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() - s;
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use growRight().

template<int Dim>
Interval<Dim> &
growRightInPlace(Interval<Dim> &dom, const Loc<Dim> &s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() + s[d].first();
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use growRight().

template<int Dim>
Interval<Dim> &
growRightInPlace(Interval<Dim> &dom, int s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() + s;
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use shrinkLeft().

template<int Dim>
Interval<Dim> &
shrinkLeftInPlace(Interval<Dim> &dom, const Loc<Dim> &s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + s[d].first();
      int b = dom[d].last();
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use shrinkLeft().

template<int Dim>
Interval<Dim> &
shrinkLeftInPlace(Interval<Dim> &dom, int s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + s;
      int b = dom[d].last();
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use growLeft().

template<int Dim>
Interval<Dim> &
growLeftInPlace(Interval<Dim> &dom, const Loc<Dim> &s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - s[d].first();
      int b = dom[d].last();
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}

/// Deprecated. Use growLeft().

template<int Dim>
Interval<Dim> &
growLeftInPlace(Interval<Dim> &dom, int s)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - s;
      int b = dom[d].last();
      dom[d] = Interval<1>(a, b);
    }
  return dom;
}



/// Shrinks the Interval dom from the right by s[i] in direction i.

template<int Dim>
inline Interval<Dim> 
shrinkRight(const Interval<Dim> &dom, const Loc<Dim> &s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() - s[d].first();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Shrinks the Interval dom from the right by s in every direction.

template<int Dim>
inline Interval<Dim> 
shrinkRight(const Interval<Dim> &dom, int s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() - s;
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Grows the Interval dom to the right by s[i] in direction i.

template<int Dim>
inline Interval<Dim> 
growRight(const Interval<Dim> &dom, const Loc<Dim> &s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() + s[d].first();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Grows the Interval dom to the right by s in every direction.

template<int Dim>
inline Interval<Dim> 
growRight(const Interval<Dim> &dom, int s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first();
      int b = dom[d].last() + s;
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}


/// Shrinks the Interval dom from the left by s[i] in direction i.

template<int Dim>
inline Interval<Dim> 
shrinkLeft(const Interval<Dim> &dom, const Loc<Dim> &s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + s[d].first();
      int b = dom[d].last();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Shrinks the Interval dom from the left by s in every direction.
template<int Dim>
inline Interval<Dim> 
shrinkLeft(const Interval<Dim> &dom, int s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + s;
      int b = dom[d].last();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Grows the Interval dom to the left by s[i] in direction i.

template<int Dim>
inline Interval<Dim> 
growLeft(const Interval<Dim> &dom, const Loc<Dim> &s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - s[d].first();
      int b = dom[d].last();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Grows the Interval dom to the left by s in every direction.

template<int Dim>
inline Interval<Dim> 
growLeft(const Interval<Dim> &dom, int s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - s;
      int b = dom[d].last();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}


/// Grows the interval dom by s on each side/dim.

template<int Dim>
inline Interval<Dim>
grow(const Interval<Dim> &dom, int s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - s;
      int b = dom[d].last() + s;
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Grows the interval dom by s[i] on each side in direction i.

template<int Dim>
inline Interval<Dim>
grow(const Interval<Dim> &dom, const Loc<Dim> &s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - s[d].first();
      int b = dom[d].last() + s[d].first();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}


/// Shrinks the interval dom by s on each side/dim.

template<int Dim>
inline Interval<Dim>
shrink(const Interval<Dim> &dom, int s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + s;
      int b = dom[d].last() - s;
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}

/// Shrinks the interval dom by s[i] on each side in direction i.

template<int Dim>
inline Interval<Dim>
shrink(const Interval<Dim> &dom, const Loc<Dim> &s)
{
  Interval<Dim> ret = Pooma::NoInit();
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + s[d].first();
      int b = dom[d].last() - s[d].first();
      ret[d] = Interval<1>(a, b);
    }
  return ret;
}



//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_SHRINK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Shrink.h,v $   $Author: richi $
// $Revision: 1.11 $   $Date: 2004/11/29 12:23:51 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
