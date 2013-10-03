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

#ifndef POOMA_DOMAIN_SPLIT_H
#define POOMA_DOMAIN_SPLIT_H

//-----------------------------------------------------------------------------
// Function:
//   bool split(domain, domain, domain);
// Class:
//   SplitDomain<T,Dim>
//   SplitDomainSingle<T,Dim,bool>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * void split(domain,domain,domain) is a global function which splits the
 * first argument into two separate domains, roughly in the moddle.
 *
 * If the first argument has zero length, this does nothing.  If the first
 * argument has a length of one, the second argument is a copy of the
 * first, and the second is set to be empty.
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
 * SplitDomainSingle<T,Dim,bool strided>::split(a,b,c) splits just the
 * Dim dimension of the first argument into the second and third argument.
 * It is specialized on the third parameter indicating whether the
 * domain has unit stride or not, and whether the type is int or not.
 *
 * The default (unit-stride) version of SplitDomainSingle, which assumes
 * that the domains have unit stride.
 */

template<class T, int Dim, bool strided>
struct SplitDomainSingle {
  static void split(const T &a, int axis, T &b, T &c) {
    // types for elements in these domains
    typedef typename DomainTraits<T>::Element_t E1_t;

    // a typedef for a single-dimension domain
    typedef typename DomainTraits<T>::OneDomain_t OneDomain_t;

    // split operation depends on length
    if (axis != (Dim - 1)) {
      b[Dim-1] = a[Dim-1];
      c[Dim-1] = a[Dim-1];
    } else if (a[Dim-1].length() < 2) {
      b[Dim-1] = a[Dim-1];
    } else {
      E1_t  a0 = a[Dim-1].first();
      E1_t  a1 = a[Dim-1].last();
      E1_t mid = a0 + a[Dim-1].length()/2;
      b[Dim-1] = OneDomain_t(a0, mid-1);
      c[Dim-1] = OneDomain_t(mid, a1);
    }
  }

  static void split(const T &a, int axis, int leftLength, T &b, T &c) {
    // types for elements in these domains
    typedef typename DomainTraits<T>::Element_t E1_t;

    // a typedef for a single-dimension domain
    typedef typename DomainTraits<T>::OneDomain_t OneDomain_t;

    // split operation depends on length
    if (axis != (Dim - 1)) {
      b[Dim-1] = a[Dim-1];
      c[Dim-1] = a[Dim-1];
    } else if (a[Dim-1].length() < 2) {
      b[Dim-1] = a[Dim-1];
    } else {
      E1_t  a0 = a[Dim-1].first();
      E1_t  a1 = a[Dim-1].last();
      E1_t mid = a0 + leftLength;
      b[Dim-1] = OneDomain_t(a0, mid-1);
      c[Dim-1] = OneDomain_t(mid, a1);
    }
  }

  static void split(const T &a, T &b, T &c) { split(a, Dim-1, b, c); }
};

/**
 * The non-unit-stride version of SplitDomainSingle.
 */

template<class T, int Dim>
struct SplitDomainSingle<T,Dim,true> {
  static void split(const T &a, int axis, T &b, T &c) {
    // types for elements in these domains
    typedef typename DomainTraits<T>::Element_t E1_t;

    // a typedef for a single-dimension domain
    typedef typename DomainTraits<T>::OneDomain_t OneDomain_t;

    // split operation depends on length
    if (axis != (Dim - 1)) {
      b[Dim-1] = a[Dim-1];
      c[Dim-1] = a[Dim-1];
    } else if (a[Dim-1].length() < 2) {
      b[Dim-1] = a[Dim-1];
    } else {
      E1_t  a0 = a[Dim-1].first();
      E1_t  a1 = a[Dim-1].last();
      E1_t   s = a[Dim-1].stride();
      E1_t mid = a0 + (a[Dim-1].length()/2 * s);
      b[Dim-1] = OneDomain_t(a0, mid - s, s);
      c[Dim-1] = OneDomain_t(mid, a1, s);
    }
  }

  static void split(const T &a, int axis, int leftLength, T &b, T &c) {
    // types for elements in these domains
    typedef typename DomainTraits<T>::Element_t E1_t;

    // a typedef for a single-dimension domain
    typedef typename DomainTraits<T>::OneDomain_t OneDomain_t;

    // split operation depends on length
    if (axis != (Dim - 1)) {
      b[Dim-1] = a[Dim-1];
      c[Dim-1] = a[Dim-1];
    } else if (a[Dim-1].length() < 2) {
      b[Dim-1] = a[Dim-1];
    } else {
      E1_t  a0 = a[Dim-1].first();
      E1_t  a1 = a[Dim-1].last();
      E1_t   s = a[Dim-1].stride();
      E1_t mid = a0 + (leftLength * s);
      b[Dim-1] = OneDomain_t(a0, mid - s, s);
      c[Dim-1] = OneDomain_t(mid, a1, s);
    }
  }

  static void split(const T &a, T &b, T &c) { split(a, Dim-1, b, c); }
};

/**
 * Special version of SplitDomainSingle for int's, which must be
 * handled uniquely.
 */

template<int Dim, bool strided>
struct SplitDomainSingle<int,Dim,strided> {
  static void split(int a, int, int &b, int &c) {
    // splitting an integer just means to copy a to b, and set c to zero
    b = a;
    c = 0;
  }
  static void split(int a, int, int, int &b, int &c) {
    b = a;
    c = 0;
  }
  static void split(int a, int &b, int &c) {
    b = a;
    c = 0;
  }
};


/**
 * SplitDomain implements a basic template meta-program to split
 * each dimension separately of the multidimensional domain.
 * It uses SplitDomainSingle to do the single-domain splits,
 * telling that struct whether the domain has unit stride or not.
 * A general version of SplitDomain is defined, to split the
 * domain in the 'Dim' dimension, and then a specialization is provided
 * for Dim==1 that stops the metaprogram recursion.
 */

template<class T, int Dim>
struct SplitDomain {
  // domain has non unit stride
  enum { strided = !DomainTraits<T>::unitStride };

  static void split(const T &a, T &b, T &c) {
    // Split along the current dimension
    SplitDomainSingle<T,Dim,strided>::split(a, b, c);

    // Split the remaining dimensions
    SplitDomain<T,Dim-1>::split(a, b, c);
  }

  static void split(const T &a, int axis, T &b, T &c) {
    // Split along the current dimension
    SplitDomainSingle<T,Dim,strided>::split(a, axis, b, c);

    // Split the remaining dimensions
    SplitDomain<T,Dim-1>::split(a, axis, b, c);
  }

  static void split(const T &a, int axis, int leftLength, T &b, T &c) {
    // Split along the current dimension
    SplitDomainSingle<T,Dim,strided>::split(a, axis, leftLength, b, c);

    // Split the remaining dimensions
    SplitDomain<T,Dim-1>::split(a, axis, leftLength, b, c);
  }
};

template<class T>
struct SplitDomain<T,1> {
  // domain has non unit stride
  enum { strided = !DomainTraits<T>::unitStride };

  static void split(const T &a, T &b, T &c) {
    // Split along the current dimension
    SplitDomainSingle<T,1,strided>::split(a, b, c);
  }

  static void split(const T &a, int axis, T &b, T &c) {
    // Split along the current dimension
    SplitDomainSingle<T,1,strided>::split(a, axis, b, c);
  }

  static void split(const T &a, int axis, int leftLength, T &b, T &c) {
    // Split along the current dimension
    SplitDomainSingle<T,1,strided>::split(a, axis, leftLength, b, c);
  }
};


/// void split(domain,domain,domain) is a global function which splits the
/// first argument into two separate domains, roughly in the middle.  If
/// the first argument has zero length, this does nothing.  If the first
/// argument has a length of one, the second argument is a copy of the
/// first, and the second is set to be empty.
///
/// The implementation of split is deferred to the SplitDomain
/// struct, which performs the split for each dimension.

template<class T>
inline void split(const T &a, T &b, T &c)
{
  SplitDomain<T,DomainTraits<T>::dimensions>::split(a, b, c);
}

/// void split(domain,axis,domain,domain) is a global function which splits
/// the first argument into two separate domains just along the Nth axis instead
/// of along all axes.  Otherwise it is the same as the other global split.

template<class T>
inline void split(const T &a, int axis, T &b, T &c)
{
  SplitDomain<T,DomainTraits<T>::dimensions>::split(a, axis, b, c);
}

/// void split(domain,axis,leftLength,domain,domain) is a global function which
/// splits the first argument into two separate domains just along the Nth axis
/// with specified size of the left part of the domain.
/// Otherwise it is the same as the other global split.

template<class T>
inline void split(const T &a, int axis, int leftLength, T &b, T &c)
{
  SplitDomain<T,DomainTraits<T>::dimensions>::split(a, axis, leftLength, b, c);
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_SPLIT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Split.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
