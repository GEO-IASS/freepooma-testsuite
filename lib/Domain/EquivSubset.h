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

#ifndef POOMA_DOMAIN_EQUIV_SUBSET_H
#define POOMA_DOMAIN_EQUIV_SUBSET_H

//-----------------------------------------------------------------------------
// Function:
//   domain equivSubset(domain, domain, domain);
// Class:
//   EquivSubsetDomain<T1,T2,T3,Dim>
//   EquivSubsetDomainSingle<T1,T2,T3,bool>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * domain4 equivSubset(domain1,domain2,domain3) is a global function
 * which finds the 'equivalent subset' for domain3 given a linear relationship
 * between domain1 and domain2.
 *
 * For example, given the relationship
 *    I --> 2I - 1
 * and a domain  3J, then the equivalent subset is
 *    3J --> 2(3J) - 1 = 6J - 1
 * The returned domain type is the most general type which could hold
 * the data in the domains 1, 2 and 3.
 *
 * EquivSubsetDomainSingle<T1,T2,T3,T4,int Dim,bool strided>::equiv(a,b,c,d)
 * finds the equivalent subset for c given the relationship between a and b,
 * for just the Dim'th dimension.  a, b, and c are assumed to be 1D domains.
 * The final boolean template parameter is used to specialize the calculation
 * to the following two cases:
 *   - strided == false: all the domains have unit stride.
 *   - strided == true: one or more domains have non-unit stride.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Domain/NewDomain.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * The default (unit-stride) version of EquivSubsetDomainSingle, which assumes
 * that all arguments to 'equiv' are domains with unit stride.  The fourth
 * argument, the returned equivalent subset domain, is assumed to have been
 * set equal to the third domain beforehand.  It will be offset and
 * scaled in the same way that the second is from the first.
 */

template<class T1, class T2, class T3, int Dim, bool strided>
struct EquivSubsetDomainSingle {
  static void equiv(const T1 &a, const T2 &b, T3 &d) {
    // Since this has unit stride, that is the only possible difference.
    // So, just offset by the difference in a and b.
    d[Dim-1] += (b.first() - a.first());
  }
};

/**
 * The non-unit-stride version of EquivSubsetDomainSingle, which does extra
 * work for the case where a, b, or c do not have unit stride.  The fourth
 * argument, the returned equivalent subset domain, is assumed to have been
 * set equal to the third domain beforehand.  It will be offset and
 * scaled in the same way that the second is from the first.
 *
 * To calculate the non-unit-stride equiv. subset, we do this:
 *   -# Look for the relationship from a --> b, in the form:
 *         b = m a + k
 *      Calculate the change in stride, m, by
 *         m = stride(b) / stride(a)
 *      which must be an integer.  Calculate the offset, k, by
 *         k = b - m a
 *      using the first point in the domains a and b.
 *   -# Apply the linear transformation to the third domain to get the fourth:
 *         d = m c + k
 */

template<class T1, class T2, class T3, int Dim>
struct EquivSubsetDomainSingle<T1,T2,T3,Dim,true> {
  static void equiv(const T1 &a, const T2 &b, T3 &d) {
    // types for elements in these domains
    typedef typename DomainTraits<T3>::Element_t E3_t;

    // Find the ratio of the stride of b to the stride of a.  This must
    // be an integral value.
    E3_t m = b.stride() / a.stride();
    PAssert(m * a.stride() == b.stride());

    // Find the offset from a to b.
    E3_t k = b.first() - m * a.first();

    // Apply this change to the final argument.
    d[Dim-1] *= m;
    d[Dim-1] += k;
  }
};


/**
 * EquivSubsetDomain implements a basic template meta-program to find the
 * equiv subset of each dimension separately of the multidimensional domains.
 * It uses EquivSubsetDomainSingle to do the single-domain calculations,
 * telling that struct whether the domains have unit stride
 * or not.  A general version of EquivSubsetDomain is defined, to calculate the
 * subset in the 'Dim' dimension, and then a specialization is provided
 * for Dim==1 that stops the metaprogram recursion.
 *
 * Since this is calculating a new domain d given domains c , b, and a so that
 * c --> d in the same way that a --> b, we can set things here so that we
 * first set d == c, and then modify d accordingly.  So, this struct assumes
 * that d == c already, and does not need to have c provided in another var.
 */

template<class T1, class T2, class T3, int Dim>
struct EquivSubsetDomain {
  // determine if the third domain type has non unit stride
  enum { strided = !DomainTraits<T3>::unitStride };

  static void equiv(const T1 &a, const T2 &b, T3 &c) {
    // the types for the two domains used in EquivSubsetDomainSingle may be a
    // little different than T1 and T2, since EquivSubsetDomainSingle works
    // with 1D domains for a and b, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // calculate the 'Dim' dimension, and then the lower dims
    EquivSubsetDomainSingle<Dom1_t,Dom2_t,T3,Dim,strided>::equiv(
      DomainTraits<T1>::getDomain(a,Dim-1), 
      DomainTraits<T2>::getDomain(b,Dim-1), c);
    EquivSubsetDomain<T1,T2,T3,Dim-1>::equiv(a,b,c);
  }
};

template<class T1, class T2, class T3>
struct EquivSubsetDomain<T1,T2,T3,1> {
  // determine if the third domain type has non unit stride
  enum { strided = !DomainTraits<T3>::unitStride };

  static void equiv(const T1 &a, const T2 &b, T3 &c) {
    // the types for the two domains used in EquivSubsetDomainSingle may be a
    // little different than T1 and T2, since EquivSubsetDomainSingle works
    // with 1D domains, not N-D domains.
    typedef typename DomainTraits<T1>::OneDomain_t Dom1_t;
    typedef typename DomainTraits<T2>::OneDomain_t Dom2_t;

    // calculate the lowest dimension
    EquivSubsetDomainSingle<Dom1_t,Dom2_t,T3,1,strided>::equiv(
      DomainTraits<T1>::getDomain(a,0), DomainTraits<T2>::getDomain(b,0), c);
  }
};


/**
 * A simple struct used to figure out the return type when examining
 * types T1,T2,T3.  It defines a typedef 'Type_t' for what the return type is.
 * Note that we use the 'DomainChangeDim' mechanism after we find out the
 * type when combining T1,T2,T3, since the combined type will have a
 * dimension of dim(T1) + dim(T2) + dim(T3), and we want the dim to be
 * the same as T1.
 */
template<class T1, class T2, class T3>
struct EquivSubsetReturnType {
  typedef typename NewDomain3<T1,T2,T3>::Type_t Combine_t;
  typedef typename 
    DomainChangeDim<Combine_t,DomainTraits<T1>::dimensions>::NewType_t Type_t;
};

/**
 * domain4 equivSubset(domain1,domain2,domain3) is a global function
 * which finds the 'equivalent subset' for domain3 given a linear relationship
 * between domain1 and domain2.  For example, given the relationship
 *    I --> 2I - 1
 * and a domain  3J, then the equivalent subset is
 *    3J --> 2(3J) - 1 = 6J - 1
 * The returned domain type is the most general type which could hold
 * the data in the domains 1, 2 and 3.  The NewDomain3 struct is used
 * to find the type of domain which would result if the three types
 * were combined, but with the proper number of dimensions.
 *
 * The implementation of equivSubset is deferred to the EquivSubsetDomain
 * struct, which performs the intersection for each dimension and and's
 * the results together.
 */

template<class T1, class T2, class T3>
inline typename EquivSubsetReturnType<T1,T2,T3>::Type_t
equivSubset(const T1 &a, const T2 &b, const T3 &c)
{
  typedef typename EquivSubsetReturnType<T1,T2,T3>::Type_t T4;
  CTAssert((int)DomainTraits<T1>::dimensions == DomainTraits<T2>::dimensions);
  CTAssert((int)DomainTraits<T1>::dimensions == DomainTraits<T3>::dimensions);
  T4 d = c;
  EquivSubsetDomain<T1,T2,T4,DomainTraits<T1>::dimensions>::equiv(a, b, d);
  return d;
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_EQUIV_SUBSET_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EquivSubset.h,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
