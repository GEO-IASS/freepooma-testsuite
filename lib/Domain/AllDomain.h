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

#ifndef POOMA_DOMAIN_ALL_DOMAIN_H
#define POOMA_DOMAIN_ALL_DOMAIN_H

//-----------------------------------------------------------------------------
// Class:
// AllDomain<int>
// DomainTraits<AllDomain<int> >
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 *  @ingroup Domain
 *
 * AllDomain is one of the domain wildcards, which are used when constructing
 * other domains using specific combination rules.  AllDomain means to use
 * the entire domain of a second 'reference' domain when constructing a new
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
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * AllDomain is a special domain class which is used as a 'wildcard'.
 * Wildcards are useful when constructing new domains based on some other
 * 'reference' domain, which is done when doing things like making a new
 * view of an Array.  Wildcard domains use the reference domain to determine
 * what the 'final' domain should be.  AllDomain refers to 'use the same
 * exact domain values as the reference domain'.
 *
 * AllDomain can be used as one of the arguments to the 'combineSlice' or
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
class AllDomain
{
public:
  //
  // Typedefs and static data
  //

  typedef AllDomain<Dim> Domain_t;
  typedef AllDomain<Dim> NewDomain1_t;
  typedef AllDomain<1>   OneDomain_t;
  typedef int            Element_t;

  enum { dimensions = Dim };

  //
  // Constructors.
  //

  // default constructor
  AllDomain() {
    CTAssert(Dim > 0);
  }

  // copy constructor
  AllDomain(const AllDomain<Dim> &) {
    CTAssert(Dim > 0);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~AllDomain() { }

  //
  // Domain-like methods
  //

  // Get the Nth element of our domain
  OneDomain_t operator[](int) const { return OneDomain_t(); }

  // setDomain: for AllDomain, this does nothing, since there is nothing
  // to set.  There is only one this we can set this with, and that is
  // another AllDomain
  void setDomain(const AllDomain<Dim> &) { }

  // given another reference domain, return the proper values for first,
  // length, and stride
  template<class T>
  typename DomainTraits<T>::Element_t first(const T &u) const {
    return u.first();
  }
  int first(int u) const { return u; }

  template<class T>
  typename DomainTraits<T>::Element_t length(const T &u) const {
    return u.length();
  }
  int length(int) const { return 1; }

  template<class T>
  typename DomainTraits<T>::Element_t stride(const T &u) const {
    return u.stride();
  }
  int stride(int) const { return 1; }


  //
  // operator=, which does nothing here.
  //

  AllDomain<Dim> &operator=(const AllDomain<Dim> &) { return *this; }

protected:

private:
};


/**
 * DomainTraits<AllDomain<Dim>> provides traits information about AllDomain,
 * which is one of the domain wildcards.  It has a quite stripped-down
 * selection of traits, the basic ones needed to allow wildcards to be used
 * in the construction of regular and strided domains.  This includes the
 * dimension and the type of the wildcard, and an enum indicating that it is
 * a wildcard.  Also, getDomain returns a 1D element of the N-dimensional
 * list of wildcards.
 */

template<int Dim>
struct DomainTraits< AllDomain<Dim> >
{
  // necessary typedefs: the domain type when this is copied
  typedef AllDomain<Dim> Domain_t;
  typedef AllDomain<Dim> NewDomain1_t;
  typedef AllDomain<1>   OneDomain_t;
  typedef AllDomain<1>   PointDomain_t;
  typedef AllDomain<Dim> AskDomain_t;

  // necessary static data
  enum { domain = true };
  enum { dimensions = Dim,
	 sliceDimensions = Dim };
  enum { wildcard = true };
  enum { singleValued = false };

  // get the Nth element of the domain, and return a OneDomain_t
  // object with it (here, as a copy).  Since AllDomain does not store any
  // data, we can just return a new copy.
  static OneDomain_t getDomain(const Domain_t &, int) {
    return OneDomain_t();
  }

  // convert from the Nth element of the domain to a single point, if
  // possible, and return a PointDomain_t.  Here, we just return a new copy
  // of PointDomain_t, since this object does not have any data.
  static PointDomain_t getPointDomain(const Domain_t &, int) {
    return PointDomain_t();
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_ALL_DOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: AllDomain.h,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:16:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
