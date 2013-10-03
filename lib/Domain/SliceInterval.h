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

#ifndef POOMA_DOMAIN_SLICE_INTERVAL_H
#define POOMA_DOMAIN_SLICE_INTERVAL_H

//-----------------------------------------------------------------------------
// Class:
// SliceInterval<int,int>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Slice domain objects are domains which have N dimensions worth of 1D
 * domain data, but really represent the result of taking an M-dimensional
 * slice (M < N) of another N dimensional domain.
 *
 * SliceInterval<N,M> is is basically an array of N Interval<1> objects,
 * but it also knows that only M of these are full domains, and that N-M
 * domains are actually referring to single points.  You can retrieve all
 * N 1D domains as a normal Interval<N> object, or the smaller slice domain
 * as an Interval<M> object.
 *
 * SliceInterval defers most of its implementation to the
 * SliceDomain<DomainTraits<SliceInterval>> base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.SliceInterval.h"
#include "Domain/SliceDomain.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * SliceInterval<N,M> is a special form of domain object which stores a
 * total domain and a slice domain, both as Interval<N> or <M> objects.
 * It is a subclass of SliceDomain, which is a general base class for all
 * sliced domain objects.  A sliced domain stores a total domain and a
 * subset of the total domain, called the slice domain, which has a number
 * of dimensions M < less than the total number of dimensions N.  Other
 * than being a specialization to Interval domains, SliceInterval has
 * the exact same interface as described in SliceDomain.
 *
 * SliceInterval has only a default constructor and a copy constructor.
 * The default constructor initializes the two internal domain objects with
 * their default constructors.  It is generally initialized with the
 * 'fillSlice' or 'combineSlice' methods in NewDomain.h, which combine
 * other non-slice domains into a SliceDomain subclass.
 *
 * A traits class DomainTraits<SliceInterval<N,M>> is used to provide
 * information to the rest of the code about SliceInterval.  The traits
 * for SliceDomain subclasses are given in DomainTraits.SliceInterval.h and
 * similar files.
 *
 * There is no specialization to 1D for SliceInterval, since it is not
 * needed.  Typically, after a SliceInterval is created by whatever
 * means, you call either 'getTotalDomain()' or 'getSliceDomain()' to
 * retrieve a reference to the relevant domain, and then use that domain
 * as normal.  getTotalDomain() returns an Interval<N> ref, while
 * getSliceDomain() returns an Interval<M> ref.
 */

template<int Dim, int SliceDim>
class SliceInterval
  : public SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > >
{
public:
  //
  // Constructors.
  //

  // default constructor : initialize to an empty slice domain
  SliceInterval() { }

  // skip initialization constructor ... if a Pooma::NoInit object
  // is provided, skip initialization
  SliceInterval(const Pooma::NoInit &d)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > >(d) { }

  // copy constructor ... base class does all the work
  SliceInterval(const SliceInterval<Dim,SliceDim> &nd)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > >(nd) {
  }

  // constructors to create slices directly.
  
  template <class Base, class D1, class D2>
  SliceInterval(const Base &baseDomain, const D1 &d1, const D2 &d2)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain2<D1,D2> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2);
  }

  template <class Base, class D1, class D2, class D3>
  SliceInterval(const Base &baseDomain, const D1 &d1, const D2 &d2, 
                const D3 &d3)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain3<D1,D2,D3> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2, d3);
  }

  template <class Base, class D1, class D2, class D3, 
            class D4>
  SliceInterval(const Base &baseDomain, const D1 &d1, const D2 &d2, 
                const D3 &d3, const D4 &d4)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain4<D1,D2,D3,D4> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2, d3, d4);
  }

  template <class Base, class D1, class D2, class D3, 
            class D4, class D5>
  SliceInterval(const Base &baseDomain, const D1 &d1, const D2 &d2, 
                const D3 &d3, const D4 &d4, const D5 &d5)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain5<D1,D2,D3,D4,D5> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2, d3, d4, d5);
  }

  template <class Base, class D1, class D2, class D3, 
            class D4, class D5, class D6>
  SliceInterval(const Base &baseDomain, const D1 &d1, const D2 &d2, 
                const D3 &d3, const D4 &d4, const D5 &d5, const D6 &d6)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain6<D1,D2,D3,D4,D5,D6> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2, d3, d4, d5, d6);
  }

  template <class Base, class D1, class D2, class D3, 
            class D4, class D5, class D6, class D7>
  SliceInterval(const Base &baseDomain, const D1 &d1, const D2 &d2, 
                const D3 &d3, const D4 &d4, const D5 &d5, const D6 &d6, 
                const D7 &d7)
    : SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain7<D1,D2,D3,D4,D5,D6,D7> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2, d3, d4, d5, d6, d7);
  }

  //
  // Destructor.  For this class there is nothing to do.
  //

  ~SliceInterval() { }

  //
  // operator= ... just use the op= in the base class
  //

  SliceInterval<Dim,SliceDim> &
    operator=(const SliceInterval<Dim,SliceDim> &nd) {
      SliceDomain<DomainTraits<SliceInterval<Dim,SliceDim> > >::operator=(nd);
      return *this;
  }

protected:

private:

};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_SLICE_INTERVAL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SliceInterval.h,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
