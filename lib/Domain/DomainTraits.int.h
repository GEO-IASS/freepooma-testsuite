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

#ifndef POOMA_DOMAIN_DOMAIN_TRAITS_INT_H
#define POOMA_DOMAIN_DOMAIN_TRAITS_INT_H

//-----------------------------------------------------------------------------
// Class:
//   DomainTraits<T>, for T = char, short, int, long, and unsigned versions
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainTraits<int> is a specialization of the general DomainTraits
 * class, for the case of integers.
 *
 * Integers can often be used in
 * combination with other domain types, and the traits class defines how
 * they interact with those other domains.  Generally, an int gets promoted
 * to a Loc or Interval, based on the context, but in some cases the
 * int is used directly.
 *
 * This also defines the same kind of traits for long, char and short.
 * Both are treated just as if they were int's for Domain purposes.  The
 * same is true for the unsigned versions of these integral types.
 *
 * DomainTraits<int> stores the characteristics of integers when they are
 * used to specify domains in expressions and object constructors.  An int
 * behaves differently in different situations; for example, when used to
 * construct a new, separate domain object, an int acts just like a Loc<1>
 * object.  When used to construct an Array, an int acts just like an
 * Interval<1> object, for an interval [0 ... intval - 1].
 *
 * Since we do not ever have int as a domain "type", just in constructors
 * and subset objects, then we need only define the bare minimum of info
 * in this traits class.
 *
 * Identical traits are defined for long, char and short.  Both integral types
 * are treated just as if they were ints; they get converted to Interval
 * or Loc objects in almost all cases.  The same is true to unsigned versions.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template <int Dim> class Interval;
template <> class Interval<1>;
template <int Dim> class Loc;
template <> class Loc<1>;


template<>
struct DomainTraits<char>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<unsigned char>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<short>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<unsigned short>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<int>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<unsigned int>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<long>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };

template<>
struct DomainTraits<unsigned long>
  : public DomainTraitsScalar<Interval<1>, int, Loc<1> > { };


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_TRAITS_INT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainTraits.int.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
