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

#ifndef POOMA_DOMAIN_SLICE_RANGE_H
#define POOMA_DOMAIN_SLICE_RANGE_H

//-----------------------------------------------------------------------------
// Class:
// SliceRange<int,int>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * Slice domain objects are domains which have N dimensions worth of 1D
 * domain data, but really represent the result of taking an M-dimensional
 * slice (M < N) of another N dimensional domain.
 *
 * SliceRange<N,M> is is basically an array of N Range<1> objects,
 * but it also knows that only M of these are full domains, and that N-M
 * domains are actually referring to single points.  You can retrieve all
 * N 1D domains as a normal Range<N> object, or the smaller slice domain
 * as a Range<M> object.
 *
 * SliceRange defers most of its implementation to the
 * SliceDomain<DomainTraits<SliceRange>> base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.SliceRange.h"
#include "Domain/SliceDomain.h"
#include "Utilities/NoInit.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * SliceRange<N,M> is a special form of domain object which stores a
 * total domain and a slice domain, both as Range<N> or <M> objects.
 * It is a subclass of SliceDomain, which is a general base class for all
 * sliced domain objects.  A sliced domain stores a total domain and a
 * subset of the total domain, called the slice domain, which has a number
 * of dimensions M < less than the total number of dimensions N.  Other
 * than being a specialization to Range domains, SliceRange has
 * the exact same interface as described in SliceDomain.
 *
 * SliceRange has only a default constructor and a copy constructor.
 * The default constructor initializes the two internal domain objects with
 * their default constructors.  It is generally initialized with the
 * 'fillSlice' or 'combineSlice' methods in NewDomain.h, which combine
 * other non-slice domains into a SliceDomain subclass.
 *
 * A traits class DomainTraits<SliceRange<N,M>> is used to provide
 * information to the rest of the code about SliceRange.  The traits
 * for SliceDomain subclasses are given in DomainTraits.SliceRange.h and
 * similar files.
 *
 * There is no specialization to 1D for SliceRange, since it is not
 * needed.  Typically, after a SliceRange is created by whatever
 * means, you call either 'getTotalDomain()' or 'getSliceDomain()' to
 * retrieve a reference to the relevant domain, and then use that domain
 * as normal.  getTotalDomain() returns a Range<N> ref, while
 * getSliceDomain() returns a Range<M> ref.
 */

template<int Dim, int SliceDim>
class SliceRange
  : public SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > >
{
public:
  //
  // Constructors.
  //

  // default constructor : initialize to an empty slice domain
  SliceRange() { }

  // skip initialization constructor ... if any version of AllDomain
  // is provided, skip initialization
  SliceRange(const Pooma::NoInit &e)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > >(e) { }

  // copy constructor ... base class does all the work
  SliceRange(const SliceRange<Dim,SliceDim> &nd)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > >(nd) {
  }
  
  // copy constructor from a SliceInterval
  SliceRange(const SliceInterval<Dim,SliceDim> &nd)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > >(nd) {
  }
  
  // constructors to create slices directly.
  
  template <class Base, class D1>
  SliceRange(const Base &baseDomain, const D1 &d1)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain1<D1> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      // Why is this CTAssert commented out here?  -- JCC
      // CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1);
  }
  
  template <class Base, class D1, class D2>
  SliceRange(const Base &baseDomain, const D1 &d1, const D2 &d2)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
        Pooma::NoInit())
  {
      typedef NewDomain2<D1,D2> NewDomain_t;
      typedef typename NewDomain_t::SliceType_t SliceType_t;
      CTAssert(DomainTraits<SliceType_t>::dimensions == Dim);
      CTAssert(DomainTraits<SliceType_t>::sliceDimensions == SliceDim);
      NewDomain_t::fillSlice(*this, baseDomain, d1, d2);
  }

  template <class Base, class D1, class D2, class D3>
  SliceRange(const Base &baseDomain, const D1 &d1, const D2 &d2, const D3 &d3)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
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
  SliceRange(const Base &baseDomain, const D1 &d1, const D2 &d2, const D3 &d3, 
             const D4 &d4)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
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
  SliceRange(const Base &baseDomain, const D1 &d1, const D2 &d2, const D3 &d3, 
             const D4 &d4, const D5 &d5)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
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
  SliceRange(const Base &baseDomain, const D1 &d1, const D2 &d2, const D3 &d3, 
             const D4 &d4, const D5 &d5, const D6 &d6)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
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
  SliceRange(const Base &baseDomain, const D1 &d1, const D2 &d2, const D3 &d3, 
             const D4 &d4, const D5 &d5, const D6 &d6, const D7 &d7)
    : SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > > (
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

  ~SliceRange() { }

  //
  // operator= ... just use the op= in the base class
  //

  SliceRange<Dim,SliceDim> &
    operator=(const SliceRange<Dim,SliceDim> &nd) {
      SliceDomain<DomainTraits<SliceRange<Dim,SliceDim> > >::operator=(nd);
      return *this;
  }

protected:

private:

};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_SLICE_RANGE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SliceRange.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
