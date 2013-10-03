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

/** @file
 * @ingroup Utilities
 * @brief
 * A tag for checking whether the terms in an expression have
 * conforming domains.
 */

#ifndef POOMA_UTILITIES_CONFORM_H
#define POOMA_UTILITIES_CONFORM_H

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class T, class A> struct LeafFunctor;
template<class T> class Scalar;

/**
 * When evaluating expressions, we need to check to make sure that the
 * left hand side has the same domain as the right.  To do that, we
 * build a ConformTag functor with the domain from the lhs, and
 * drop it down the rhs tree.  It makes sure that the lengths of the
 * domains are all the same.
 */

template<int D>
class ConformTag
{
public:
  template<class Domain>
  ConformTag(const Domain& domain) 
  {
    for (int i=0; i<D; ++i)
      lengths_m[i] = domain[i].length();
  }
  int length(int i) const { return lengths_m[i]; }
private:
  int lengths_m[D];
};

/// @name conforms(Domain,ConformTag)
/// Check to see whether a given domain conforms with the given
/// ConformTag.  
/// 
/// We just check to see if the length of each dimension is the same.
/// We specialize this for ranks 1 thru 7 for efficiency.
//@{

template<class Domain>
bool conforms(const Domain &d, const ConformTag<1> &ct)
{
  return d.length() == ct.length(0);
}

template<class Domain>
bool conforms(const Domain &d, const ConformTag<2> &ct)
{
  return (d[0].length() == ct.length(0))
      && (d[1].length() == ct.length(1));
}

template<class Domain>
bool conforms(const Domain &d, const ConformTag<3> &ct)
{
  return (d[0].length() == ct.length(0)) 
      && (d[1].length() == ct.length(1)) 
      && (d[2].length() == ct.length(2));
}

template<class Domain>
bool conforms(const Domain &d, const ConformTag<4> &ct)
{
  return (d[0].length() == ct.length(0)) 
      && (d[1].length() == ct.length(1)) 
      && (d[2].length() == ct.length(2)) 
      && (d[3].length() == ct.length(3));
}

template<class Domain>
bool conforms(const Domain &d, const ConformTag<5> &ct)
{
  return (d[0].length() == ct.length(0)) 
      && (d[1].length() == ct.length(1))
      && (d[2].length() == ct.length(2))
      && (d[3].length() == ct.length(3))
      && (d[4].length() == ct.length(4));
}

template<class Domain>
bool conforms(const Domain &d, const ConformTag<6> &ct)
{
  return (d[0].length() == ct.length(0)) 
      && (d[1].length() == ct.length(1))
      && (d[2].length() == ct.length(2))
      && (d[3].length() == ct.length(3))
      && (d[3].length() == ct.length(4))
      && (d[5].length() == ct.length(5));
}

template<class Domain>
bool conforms(const Domain &d, const ConformTag<7> &ct)
{
  return (d[0].length() == ct.length(0)) 
      && (d[1].length() == ct.length(1))
      && (d[2].length() == ct.length(2))
      && (d[3].length() == ct.length(3))
      && (d[3].length() == ct.length(4))
      && (d[5].length() == ct.length(5))
      && (d[6].length() == ct.length(6));
}

//@}

/// Scalars conform with anything, so always return true.

template<int D, class T>
struct LeafFunctor<Scalar<T>, ConformTag<D> >
{
  typedef bool Type_t;
  static Type_t apply(const Scalar<T> &, const ConformTag<D> &)
  {
    return true;
  }
};

#endif // POOMA_UTILITIES_CONFORM_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Conform.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
