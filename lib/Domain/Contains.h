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

#ifndef POOMA_DOMAIN_CONTAINS_H
#define POOMA_DOMAIN_CONTAINS_H

//-----------------------------------------------------------------------------
// Function:
//   bool contains(domain, domain);
// Class:
//   ContainsDomain<T1,T2,Dim>
//   ContainsDomainSingle<T1,T2,bool>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * bool contains(domain,domain) is a global function which determines if
 * the points in the second domain are all points which are in the first
 * domain.
 *
 * If there is even just one point in the second not in the first,
 * then this returns false.  Note that the order is important: if
 * contains(a,b) is true, then the only way that contains(b,a) can be true is
 * if a == b.  The order for the query is: is b contained within a?
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * ContainsDomainSingle<T1,T2,bool strided>::contains(a,b) compares two
 * domains a and b of type T1 and T2, and returns true if a contains b.
 * a and b are assumed to be 1D domains, and this struct is used by
 * ContainsDomain for each dimension in a multidimensional contains operation.
 * The final boolean template parameter is used to specialize the calculation
 * to the following two cases:
 *   strided == false: one or both of the domains has unit stride.  In
 *                     this case, the computation is quite simple: check if
 *                     the endpoints of b lie within the endpoints of a.
 *   strided == true: neither domain has unit stride.  This is more complicated
 *                    since it is possible that even if the endpoints of b
 *                    are contained in a, that all the points in b are not
 *                    found in a.  The striding of a may lead to it not
 *                    referring to points in b.  Only do this calculation when
 *                    absolutely necessary.
 *
 *
 * The default (unit-stride) version of ContainsDomainSingle, which assumes
 * that both arguments to 'contains' are 1D domains with unit stride
 */

template<class T1, class T2, bool strided>
struct ContainsDomainSingle {
  static bool contains(const T1 &a, const T2 &b) {
    return (a.min() <= b.min() && a.max() >= b.max());
  }
};

/**
 * The non-unit-stride version of ContainsDomainSingle, which does extra
 * work for the case where a and b do not have unit stride.
 */

template<class T1, class T2>
struct ContainsDomainSingle<T1,T2,true> {
  static bool contains(const T1 &a, const T2 &b) {
    // types for elements in these domains
    typedef typename DomainTraits<T1>::Element_t E1_t;
    typedef typename DomainTraits<T2>::Element_t E2_t;

    // Find min and max values of type domains
    E1_t a0 = a.min();
    E1_t a1 = a.max();
    E1_t  s = a.stride();
    E2_t b0 = b.min();
    E2_t b1 = b.max();
    E2_t  t = b.stride();
    if (s < 0)
      s = -s;
    if (t < 0)
      t = -t;

    // We can do a quick short-circuit check to make sure they do not overlap
    // at all just from their endpoints.  If they don't even do this, we can
    // quit and say they do not touch.
    bool quicktest = (a0 <= b0 && a1 >= b1);
    if (!quicktest || s == 1)
      return quicktest;

    // OK, the endpoints of a contain those of b, and we must find out if
    // all the points in b are found in a.  This will be true if:
    //   1. The stride of b is a multiple of the stride of a
    //   2. The endpoints of b are found in a
    // If either of these conditions are false, a does not contain b 
    return (t % s == 0) && ((b0-a0) % s == 0) && ((a1-b1) % s == 0);
  }
};


/**
 * ContainsDomain implements a basic template meta-program to
 * compare each dimension separately of the multidimensional domains for
 * whether a contains b.  It uses ContainsDomainSingle to do the single-domain
 * comparison, telling that struct whether the domains have unit stride
 * or not.  A general version of ContainsDomain is defined, to compare the
 * domains in the 'Dim' dimension, and then a specialization is provided
 * for Dim==1 that stops the metaprogram recursion.
 */

template<class T1, class T2, int Dim>
struct ContainsDomain {
  // domain has non unit stride
  enum { strided = !DomainTraits<T1>::unitStride };

  // which dimension we're working with
  enum { n = Dim - 1 };

  static bool contains(const T1 &a, const T2 &b) {
    // the types for the two domains used in ContainsDomainSingle may be a
    // little different than T1 and T2, since ContainsDomainSingle works
    // with 1D domains, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // compare the results for the 'Dim' dimension, and then for lower dims
    return ContainsDomainSingle<Dom1_t,Dom2_t,strided>::contains(
      DomainTraits<T1>::getDomain(a,n), DomainTraits<T2>::getDomain(b,n)) &&
      ContainsDomain<T1,T2,n>::contains(a,b);
  }
};

template<class T1, class T2>
struct ContainsDomain<T1,T2,1> {
  // domain has non unit stride
  enum { strided = !DomainTraits<T1>::unitStride };

  static bool contains(const T1 &a, const T2 &b) {
    // the types for the two domains used in ContainsDomainSingle may be a
    // little different than T1 and T2, since ContainsDomainSingle works
    // with 1D domains, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // compare the results for the 'Dim' dimension, and then for lower dims
    return ContainsDomainSingle<Dom1_t,Dom2_t,strided>::contains(
      DomainTraits<T1>::getDomain(a,0), DomainTraits<T2>::getDomain(b,0));
  }
};


/**
 * bool contains(domain1, domain2) is one of the domain calculus routines
 * used to analyze domains to determine their relative characteristics.  It
 * returns true if ALL the points in domain2 are found in the set of points
 * which form domain1.  If so, this returns true, else it returns false.
 *
 * The implementation of contains is deferred to the ContainsDomain
 * struct, which performs the contains comparison for each dimension and and's
 * the results together.
 */

template<class T1, class T2>
inline bool contains(const T1 &a, const T2 &b)
{
  CTAssert((int)DomainTraits<T1>::dimensions == DomainTraits<T2>::dimensions);
  return ContainsDomain<T1,T2,DomainTraits<T1>::dimensions>::contains(a, b);
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_CONTAINS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Contains.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
