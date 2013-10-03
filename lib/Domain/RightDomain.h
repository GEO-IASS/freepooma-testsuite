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

#ifndef POOMA_DOMAIN_RIGHT_DOMAIN_H
#define POOMA_DOMAIN_RIGHT_DOMAIN_H

//-----------------------------------------------------------------------------
// Class:
// RightDomain<int>
// DomainTraits<RightDomain<int> >
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * RightDomain is one of the domain wildcards, which are used when constructing
 * other domains using specific combination rules.
 *
 * RightDomain means to use
 * the ending endpoint of the domain of a second 'reference' domain, with
 * a new user-provided left endpoint, when constructing a new
 * domain.  It is also used when constructing new domains with no other
 * arguments to mean that the domain should not be initialized, which can
 * save considerable time in some circumstances.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Domain/Loc.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * RightDomain is a special domain class which is used as a 'wildcard'.
 * Wildcards are useful when constructing new domains based on some other
 * 'reference' domain, which is done when doing things like making a new
 * view of an Array.  Wildcard domains use the reference domain to determine
 * what the 'final' domain should be.  RightDomain refers to 'use the right
 * endpoint of the reference domain, with newly provided left endpoint,
 * as the new domain values'.
 *
 * RightDomain can be used as one of the arguments to the 'combineSlice' or
 * 'fillSlice' routines in the NewDomain combiners, in which case the user-
 * supplied reference domain is used with the 'setWildcardDomain' method
 * of the domain being filled to get the final domain settings.
 *
 * Wildcard domains in general can also be used in the constructors for
 * regular domain objects.  If they are given, they indicate that those
 * dimensions should not be initialized, which can be helpful to avoid
 * extra unneeded work when the domain will be filled with new values very
 * soon.
 */

template<int Dim>
class RightDomain
{
public:
  //
  // Typedefs and static data
  //

  typedef RightDomain<Dim> Domain_t;
  typedef RightDomain<1>   OneDomain_t;
  typedef int              Element_t;

  enum { dimensions = Dim };

  //
  // Constructors.
  //

  // default constructor
  RightDomain() : endpoints_m(Pooma::NoInit()) { CTAssert(Dim > 0); }

  // copy constructor
  RightDomain(const RightDomain<Dim> &d) : endpoints_m(d.endpoints_m) {
    CTAssert(Dim > 0);
  }

  // templated constructors for 1 ... 7 arguments, used to fill up the
  // values for the endpoints
  template<class T1>
  explicit RightDomain(const T1 &a)
    : endpoints_m(a) { CTAssert(Dim > 0); }

  template<class T1, class T2>
  RightDomain(const T1 &a, const T2 &b)
    : endpoints_m(a,b) { CTAssert(Dim > 0); }

  template<class T1, class T2, class T3>
  RightDomain(const T1 &a, const T2 &b, const T3 &c)
    : endpoints_m(a,b,c) { CTAssert(Dim > 0); }

  template<class T1, class T2, class T3, class T4>
  RightDomain(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
    : endpoints_m(a,b,c,d) { CTAssert(Dim > 0); }

  template<class T1, class T2, class T3, class T4, class T5>
  RightDomain(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e)
    : endpoints_m(a,b,c,d,e) { CTAssert(Dim > 0); }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6>
  RightDomain(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
      const T6 &f)
    : endpoints_m(a,b,c,d,e,f) { CTAssert(Dim > 0); }

  template<class T1, class T2, class T3, class T4, class T5,
           class T6, class T7>
  RightDomain(const T1 &a, const T2 &b, const T3 &c, const T4 &d, const T5 &e,
      const T6 &f, const T7 &g)
    : endpoints_m(a,b,c,d,e,f,g) { CTAssert(Dim > 0); }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~RightDomain() { }

  //
  // Domain-like methods
  //

  // Get the Nth element of our domain.  This only returns a copy; there
  // is no way to modify a RightDomain after it has been constructed, except
  // through the setDomain method and the operator=.
  OneDomain_t operator[](int n) const { return OneDomain_t(endpoints_m[n]); }

  // setDomain: change our RightDomain to be the newly provided one.
  void setDomain(const RightDomain<Dim> &d) { endpoints_m = d.endpoints_m; }

  // given another reference domain, return the proper values for first,
  // length, and stride
  template<class T>
  typename DomainTraits<T>::Element_t first(const T &u) const {
    // for a RightDomain, first is taken from our own first
    return endpoints_m[0].first();
  }
  int first(int u) const { return u; }

  template<class T>
  typename DomainTraits<T>::Element_t length(const T &u) const {
    // for a RightDomain, length is determined by our left endpoint and
    // the given argument's last.  We do this by making a temporary of type
    // T with the proper endpoints and calling length on that.  If the values
    // for the endpoints are inconsistent with the domain type T, it will
    // be an error.
    CTAssert(Dim == 1);
    T dom(endpoints_m[0].first(), u.last(), u.stride());
    return dom.length();
  }
  int length(int u) const {
    CTAssert(Dim == 1);
    Interval<1> dom(endpoints_m[0].first(), u);
    return dom.length();
  }

  template<class T>
  typename DomainTraits<T>::Element_t stride(const T &u) const {
    // for a RightDomain, the stride is the same as the given domain's stride
    return u.stride();
  }
  int stride(int) const { return 1; }


  //
  // operator=
  //

  RightDomain<Dim> &operator=(const RightDomain<Dim> &d) {
    endpoints_m = d.endpoints_m;
    return *this;
  }

protected:

private:
  // our list of right endpoints
  Loc<Dim> endpoints_m;
};


/**
 * DomainTraits<RightDomain<Dim>> provides traits info about RightDomain,
 * which is one of the domain wildcards.  It has a quite stripped-down
 * selection of traits, the basic ones needed to allow wildcards to be used
 * in the construction of regular and strided domains.  This includes the
 * dimension and the type of the wildcard, and an enum indicating that it is
 * a wildcard.  Also, getDomain returns a 1D element of the N-dimensional
 * list of wildcards.
 */

template<int Dim>
struct DomainTraits< RightDomain<Dim> >
{
  // necessary typedefs: the domain type when this is copied
  typedef RightDomain<Dim> Domain_t;
  typedef RightDomain<1>   OneDomain_t;
  typedef RightDomain<1>   PointDomain_t;
  typedef RightDomain<Dim> AskDomain_t;

  // necessary static data
  enum { domain          = true };
  enum { dimensions      = Dim,
	 sliceDimensions = Dim };
  enum { wildcard        = true };
  enum { singleValued    = false };

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).  Since RightDomain does not store any
  // data, we can just return a new copy.
  static OneDomain_t getDomain(const Domain_t &d, int n) {
    return d[n];
  }

  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a new copy
  // of PointDomain_t, since this object does not have any data.
  static PointDomain_t getPointDomain(const Domain_t &d, int n) {
    return d[n];
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_RIGHT_DOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RightDomain.h,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
