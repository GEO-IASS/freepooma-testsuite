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

#ifndef POOMA_DOMAIN_INTERSECT_H
#define POOMA_DOMAIN_INTERSECT_H

//-----------------------------------------------------------------------------
// Function:
//   domain intersect(domain, domain);
// Class:
//   IntersectDomain<T1,T2,T3,Dim>
//   IntersectDomainSingle<T1,T2,T3,bool>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////


/** @file
 * @ingroup Domain
 * @brief
 * domain3 intersect(domain1,domain2) is a global function which determines if
 * two domains intersect with each other, and if so, what their
 * intersection is.
 *
 * By 'intersect', we mean the set of points which are
 * found in BOTH domains, expressed as a new domain.  If they have no
 * points in common, this returns an empty domain.  The type of domain
 * returned is the most general domain type which can store the information
 * from domain1 and domain2.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Domain/DomainCalculus.h"
#include "Domain/NewDomain.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * IntersectDomainSingle<T1,T2,T3,int Dim,bool strided>::intersect(a,b,c)
 * finds the intersection of two domains a and b of type T1 and T2, and
 * puts the domain intersection in c[Dim-1]. If there are no points which these
 * domains have in common, this returns an empty domain.  The domains a and b
 * are assumed to be 1D domains.
 * The final boolean template parameter is used to specialize the calculation
 * to the following two cases:
 *  - strided == false: both of the domains have unit stride.  In
 *                     this case, the computation is quite simple: check if
 *                     the endpoints overlap in any way, and return it.
 *  - strided == true: one or both domains have non-unit stride.
 *                    This is more complicated
 *                    since it is possible that even if the endpoints overlap,
 *                    they will not have any points in common due to the
 *                    striding.  We only do this long calculation when
 *                    absolutely necessary.
 *
 * The default (unit-stride) version of IntersectDomainSingle, which assumes
 * that both arguments to 'intersect' are 1D domains with unit stride
 */

template<class T1, class T2, class T3, int Dim, bool strided>
struct IntersectDomainSingle {
  static void intersect(const T1 &a, const T2 &b, T3 &c) {
    // types for elements in these domains
    typedef typename DomainTraits<T1>::Element_t E1_t;
    typedef typename DomainTraits<T2>::Element_t E2_t;

    // find min and max of the domains
    E1_t a0 = a.min();
    E1_t a1 = a.max();
    E2_t b0 = b.min();
    E2_t b1 = b.max();

    // make sure they touch in some way
    if (a1 < b0 || a0 > b1)
      return;

    // they overlap; find where
    if (b0 > a0)
      a0 = b0;
    if (b1 < a1)
      a1 = b1;

    // fill in the domain with the correct value
    typedef typename DomainTraits<T3>::OneDomain_t Dom_t;
    c[Dim-1] = Dom_t(a0, a1);
  }
};

/**
 * The non-unit-stride version of IntersectDomainSingle, which does extra
 * work for the case where a or b do not have unit stride.
 */

template<class T1, class T2, class T3, int Dim>
struct IntersectDomainSingle<T1,T2,T3,Dim,true> {
  static void intersect(const T1 &a, const T2 &b, T3 &c) {
    // types for elements in these domains
    typedef typename DomainTraits<T1>::Element_t E1_t;
    typedef typename DomainTraits<T2>::Element_t E2_t;

    // if both strides are +1, we can do the quick operation
    E1_t s = a.stride();
    E2_t t = b.stride();
    if (s == 1 && t == 1) {
      IntersectDomainSingle<T1,T2,T3,Dim,false>::intersect(a,b,c);
      return;
    }

    // make sure they touch in some way
    E1_t a0 = a.min();
    E1_t a1 = a.max();
    E2_t b0 = b.min();
    E2_t b1 = b.max();
    if (a1 < b0 || a0 > b1)
      return;

    // OK, the endpoints overlap in some way, and we must find out if there
    // are any points in both a and b.  If so, the result will have a
    // stride which is the least-common-multiple of the strides of a and b,
    // and endpoints which match this stride and the unit-stride intersection
    // endpoints.  We use the general routine 'findIntersectionEndpoints'
    // to do this; if it returns true, then the endpoints and stride
    // parameters are set to new values, otherwise there are no points
    // in common.
    int i1, i2, is;
    if (findIntersectionEndpoints(a0, a1, s, b0, b1, t, i1, i2, is)) {
      // there was an intersection; so, adjust i1, i2, and is to have
      // the right direction.   When i1, i2, and s come back from this
      // routine, it is always the case that i1 <= i2, and is > 0.  This
      // might need to be reversed.
      typedef typename DomainTraits<T3>::OneDomain_t Dom_t;
      if (s < 0)
        c[Dim-1] = Dom_t(i2, i1, -is);
      else
        c[Dim-1] = Dom_t(i1, i2, is);
    }
  }
};



/**
 * IntersectDomain implements a basic template meta-program to
 * intersect each dimension separately of the multidimensional domains.
 * It uses IntersectDomainSingle to do the single-domain intersection,
 * telling that struct whether the domains have unit stride
 * or not.  A general version of IntersectDomain is defined, to intersect the
 * domains in the 'Dim' dimension, and then a specialization is provided
 * for Dim==1 that stops the metaprogram recursion.
 */

template<class T1, class T2, class T3, int Dim>
struct IntersectDomain {
  // either domain has non unit stride
  enum { strided = 
    !DomainTraits<T1>::unitStride || !DomainTraits<T2>::unitStride };

  static void intersect(const T1 &a, const T2 &b, T3 &c) {
    // the types for the two domains used in IntersectDomainSingle may be a
    // little different than T1 and T2, since IntersectDomainSingle works
    // with 1D domains for a and b, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // intersect the 'Dim' dimension, and then the lower dims
    IntersectDomainSingle<Dom1_t,Dom2_t,T3,Dim,strided>::intersect(
      DomainTraits<T1>::getDomain(a,Dim-1), 
      DomainTraits<T2>::getDomain(b,Dim-1), c);
    IntersectDomain<T1,T2,T3,Dim-1>::intersect(a,b,c);
  }
};

template<class T1, class T2, class T3>
struct IntersectDomain<T1,T2,T3,1> {
  // either domain has non unit stride
  enum { strided = 
    !DomainTraits<T1>::unitStride || !DomainTraits<T2>::unitStride };

  static void intersect(const T1 &a, const T2 &b, T3 &c) {
    // the types for the two domains used in IntersectDomainSingle may be a
    // little different than T1 and T2, since IntersectDomainSingle works
    // with 1D domains, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // intersect the lowest dimension
    IntersectDomainSingle<Dom1_t,Dom2_t,T3,1,strided>::intersect(
      DomainTraits<T1>::getDomain(a,0), DomainTraits<T2>::getDomain(b,0), c);
  }
};



/**
 * a simple struct used to figure out the return type when intersecting
 * types T1 and T2.  It defines a typedef 'Type_t' for what the return type is.
 * Note that we use the 'DomainChangeDim' mechanism after we find out the
 * type when combining T1 and T2, since the combined type will have a
 * dimension of dim(T1) + dim(T2), and we want the dim to be the same as T1.
 */
template<class T1, class T2>
struct IntersectReturnType {
  typedef typename NewDomain2<T1,T2>::Type_t Combine_t;
  typedef typename 
    DomainChangeDim<Combine_t,DomainTraits<T1>::dimensions>::NewType_t Type_t;
};

/**
 * domain3 intersect(domain1,domain2) is a global function which determines if
 * two domains intersect with each other, and if so, what their
 * intersection is.  By 'intersect', we mean the set of points which are
 * found in BOTH domains, expressed as a new domain.  If they have no
 * points in common, this returns an empty domain.  The type of domain
 * returned is the most general domain type which can store the information
 * from domain1 and domain2.  For example, if you intersect an Interval
 * with a Range, the most general type is Range.  The same rules that are
 * used to combine domain's together are used to determine the return type.
 *
 * The implementation of intersect is deferred to the IntersectDomain
 * struct, which performs the intersection for each dimension and and's
 * the results together.
 */

template<class T1, class T2>
inline typename IntersectReturnType<T1,T2>::Type_t
intersect(const T1 &a, const T2 &b)
{
  typedef typename IntersectReturnType<T1,T2>::Type_t T3;
  CTAssert((int)DomainTraits<T1>::dimensions == DomainTraits<T2>::dimensions);
  T3 c;
  IntersectDomain<T1,T2,T3,DomainTraits<T1>::dimensions>::intersect(a, b, c);
  return c;
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_INTERSECT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Intersect.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
