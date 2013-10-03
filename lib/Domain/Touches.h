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

#ifndef POOMA_DOMAIN_TOUCHES_H
#define POOMA_DOMAIN_TOUCHES_H

//-----------------------------------------------------------------------------
// Function:
//   bool touches(domain, domain);
// Class:
//   TouchesDomain<T1,T2,Dim>
//   TouchesDomainSingle<T1,T2,bool>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * bool touches(domain,domain) is a global function which determines if
 * two domains d1 and d2 overlap in any way.
 *
 * 'Overlap' means that there is
 * at least one point which resides in both domains.  If this is the case,
 * it returns true, otherwise false.
 *
 * touches uses a partially-specialized struct 'TouchesDomain' to do the work;
 * TouchesDomain implements a template metaprogram to check that each dimension
 * touches and and's the results together.  TouchesDomain, in turn, uses
 * TouchesDomainSingle to perform the touches calculation for two 1D domains,
 * with specialization for the case where one or both of the domains have
 * unit stride, and for the case where neither have unit stride (which is
 * a much more complicated computation).
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Domain/DomainCalculus.h"
#include "Utilities/PAssert.h"
#include "Utilities/WrappedInt.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * TouchesDomainSingle<T1,T2,bool strided>::touches(a,b) compares two
 * domains a and b of type T1 and T2, and returns true if they touch.
 * a and b are assumed to be 1D domains, and this struct is used by
 * TouchesDomain for each dimension in a multidimensional touches calculation.
 * The final boolean template parameter is used to specialize the calculation
 * to the following two cases:
 *  - strided == false: one or both of the domains has unit stride.  In
 *                     this case, the computation is quite simple: check if
 *                     the endpoints overlap in any way.
 *  - strided == true: neither domain has unit stride.  This is more complicated
 *                    since it is possible that even if the endpoints overlap,
 *                    they will not have any points in common due to the
 *                    striding.  We only do this long calculation when
 *                    absolutely necessary.
 *
 * The default (unit-stride) version of TouchesDomainSingle, which assumes
 * that both arguments to 'touches' are 1D domains with unit stride
 */

template<class T1, class T2, bool strided>
struct TouchesDomainSingle {
  static bool touches(const T1 &a, const T2 &b) {
    return (a.min() <= b.max() && a.max() >= b.min());
  }
};

/**
 * The non-unit-stride version of TouchesDomainSingle, which does extra
 * work for the case where a and b do not have unit stride.
 */

template<class T1, class T2>
struct TouchesDomainSingle<T1,T2,true> {
  static bool touches(const T1 &a, const T2 &b) {
    // We can do a quick short-circuit check to make sure they do not overlap
    // at all just from their endpoints.  If they don't even do this, we can
    // quit and say they do not touch.  If they do have overlapping endpoints,
    // though, and at least one has a unity stride (+1 or -1), we can just
    // report the results of the quick test.
    bool quicktest = TouchesDomainSingle<T1,T2,false>::touches(a, b);
    if (!quicktest || a.stride() == 1 || a.stride() == (-1) ||
	b.stride() == 1 || b.stride() == (-1)) {
      return quicktest;
    }

    // OK, the endpoints overlap in some way, and we must find out if there
    // is at least one point in a that is also in b.  Use the general
    // routine 'findLeftCommonEndpoint' to do this; if it returns true,
    // then there is a possible left endpoint for a domain that contains
    // points from both a and b.  If it returns false, then they don't touch.
    int endpoint = 0;
    return findLeftCommonEndpoint(a.min(), a.max(), a.stride(),
				  b.min(), b.max(), b.stride(), endpoint);
  }
};


/**
 * TouchesDomain implements a basic template meta-program to
 * compare each dimension separately of the multidimensional domains for
 * whether they touch.  It uses TouchesDomainSingle to do the single-domain
 * comparison, telling that struct whether the domains have unit stride
 * or not.  A general version of TouchesDomain is defined, to compare the
 * domains in the 'Dim' dimension, and then a specialization is provided
 * for Dim==1 that stops the metaprogram recursion.
 */

template<class T1, class T2, int Dim>
struct TouchesDomain {
  // both domains have non unit stride
  enum { strided = 
    !DomainTraits<T1>::unitStride && !DomainTraits<T2>::unitStride };

  static bool touches(const T1 &a, const T2 &b) {
    // the types for the two domains used in TouchesDomainSingle may be a
    // little different than T1 and T2, since TouchesDomainSingle works
    // with 1D domains, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // compare the results for the 'Dim' dimension, and then for lower dims
    return TouchesDomainSingle<Dom1_t,Dom2_t,strided>::touches(
      DomainTraits<T1>::getDomain(a,Dim-1), 
      DomainTraits<T2>::getDomain(b,Dim-1)) &&
      TouchesDomain<T1,T2,Dim-1>::touches(a,b);
  }
};

template<class T1, class T2>
struct TouchesDomain<T1,T2,1> {
  // both domains have non unit stride
  enum { strided = 
    !DomainTraits<T1>::unitStride && !DomainTraits<T2>::unitStride };

  static bool touches(const T1 &a, const T2 &b) {
    // the types for the two domains used in TouchesDomainSingle may be a
    // little different than T1 and T2, since TouchesDomainSingle works
    // with 1D domains, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // compare the results for the 'Dim' dimension, and then for lower dims
    return TouchesDomainSingle<Dom1_t,Dom2_t,strided>::touches(
      DomainTraits<T1>::getDomain(a,0), DomainTraits<T2>::getDomain(b,0));
  }
};


template<class T1, class T2, int Dim1, int Dim2>
inline bool touches2(const T1 &, const T2 &, const WrappedInt<Dim1> &,
  const WrappedInt<Dim2> &)
{
  return false;
}

template<class T1, class T2, int Dim>
inline bool touches2(const T1 &a, const T2 &b, const WrappedInt<Dim> &,
  const WrappedInt<Dim> &)
{
  return TouchesDomain<T1,T2,Dim>::touches(a, b);
}

/// bool touches(domain1, domain2) is one of the domain calculus routines
/// used to analyze domains to determine their relative characteristics.  It
/// returns true if there is at least one point in domain1 which is also in
/// domain2.  Otherwise, if they are completely disjoint, it returns false.
/// domain1 and domain2 must have the same number of dimensions; if they don't,
/// it results in an assertion failure.  The comparison is done for each
/// dimension; if any dimension fails to have overlapping domains for the
/// two domain objects, then the whole domains do not touch.
///
/// The implementation of touches is deferred to the TouchesDomain
/// struct, which performs the touches comparison for each dimension and and's
/// the results together.

template<class T1, class T2>
inline bool touches(const T1 &a, const T2 &b)
{
  if (a.empty() || b.empty())
    return false;
  else 
    return touches2(a, b, WrappedInt<DomainTraits<T1>::dimensions>(),
		    WrappedInt<DomainTraits<T2>::dimensions>());
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_TOUCHES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Touches.h,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
