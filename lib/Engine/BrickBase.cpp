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
//-----------------------------------------------------------------------------
// BrickBase and BrickViewBase non-inline template definitions.
//-----------------------------------------------------------------------------

#include "Engine/BrickBase.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

#if !defined(NOPassert)
#include <algorithm>
#endif

///////////////////////////////////////////////////////////////////////////////

namespace Pooma {

///////////////////////////////////////////////////////////////////////////////
//
// BrickBase Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// BrickBase(const Interval<Dim> &domain)
//
// Constructs a BrickBase for a multidimensional domain described by 
// the Interval<Dim> domain.
//
//-----------------------------------------------------------------------------

template <int Dim> 
BrickBase<Dim>::BrickBase(const Domain_t &dom, bool compressible)
  : layout_m(dom), compressibleBase_m(compressible)
{
  // Compute the strides and offset.

  strides_m[0] = 1;
  off_m        = -domain()[0].first();

  for (int d = 1; d < Dim; ++d)
    {
      strides_m[d] = strides_m[d-1]*domain()[d-1].length();
      off_m       -= domain()[d].first()*strides_m[d];
    }
  
  for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
}

//-----------------------------------------------------------------------------
//
// BrickBase(const Node<Domain_t> &node)
//
// Constructs a BrickBase using the specified Node.
//
//-----------------------------------------------------------------------------

template <int Dim> 
BrickBase<Dim>::BrickBase(const Node<Domain_t> &node, bool compressible)
  : layout_m(node), compressibleBase_m(compressible)
{
  // Compute the strides and offset.

  strides_m[0] = 1;
  off_m        = -domain()[0].first();

  for (int d = 1; d < Dim; ++d)
    {
      strides_m[d] = strides_m[d-1]*domain()[d-1].length();
      off_m       -= domain()[d].first()*strides_m[d];
    }
  
  for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
}

//-----------------------------------------------------------------------------
//
// BrickBase(const domainLayout<Dim> &layout)
//
// Constructs a BrickBase using the specified DomainLayout.
//
//-----------------------------------------------------------------------------

template <int Dim> 
BrickBase<Dim>::BrickBase(const Layout_t &layout, bool compressible)
  : layout_m(layout), compressibleBase_m(compressible)
{
  // Compute the strides and offset.

  strides_m[0] = 1;
  off_m        = -domain()[0].first();

  for (int d = 1; d < Dim; ++d)
    {
      strides_m[d] = strides_m[d-1]*domain()[d-1].length();
      off_m       -= domain()[d].first()*strides_m[d];
    }
  
  for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
}


///////////////////////////////////////////////////////////////////////////////
//
// BrickViewBase Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// BrickViewBase Constructors
//
// This is where the real complexity lives. There are two categories
// of views:
// 
// Slices (Dim < BaseDim):
//   These can be created either by subsetting a BrickBase with a 
//   SliceInterval or a SliceRange, by subsetting a BrickViewBase
//   with any rectangular domain (Interval, Range, SliceInterval,
//   SliceRange).
//
// Non-slices (Dim == BaseDim):
//   These are created by subsetting a BrickBase with an Interval
//   or Range, or by subsetting another non-sliced view with an
//   Interval or Range. These are handled by partial specializations
//   that fully specialize on the stride since that allows us to
//   only define the useful constructors.
// 
// We do the more general cases (the slices) first.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// BrickViewBase constructors (for sliced BrickViewBase objects)
//
// These just call the various sliceInit functions defined below.
// (Interval types are converted to Range types. This costs a few
// extra multiplies since the Intervals had unit stride in each
// direction, but that savings is not worth the code duplication);
// 
// Note that slice-of-slice constructors are inline member templates. They
// delegate their work to the DoubleSliceHelper functions defined below.
//
//-----------------------------------------------------------------------------

template <int Dim>
BrickViewBase<Dim>::
BrickViewBase(const BrickBase<Dim> &bbase,
	      const Interval<Dim> &dom)
: domain_m(Pooma::NoInit()),
  baseOffset_m(bbase.offset()), compressibleBase_m(bbase.compressibleBase())
{    
  viewInit(bbase, Range<Dim>(dom));
}

template <int Dim>
BrickViewBase<Dim>::
BrickViewBase(const BrickBase<Dim> &bbase,
	      const Range<Dim> &dom)
: domain_m(Pooma::NoInit()),
  baseOffset_m(bbase.offset()), compressibleBase_m(bbase.compressibleBase())
{    
  viewInit(bbase, dom);
}

template <int Dim>
BrickViewBase<Dim>::
BrickViewBase(const This_t &bvbase, const Range<Dim> &domain)
: domain_m(Pooma::NoInit()),  
  baseOffset_m(bvbase.baseOffset()), 
  compressibleBase_m(bvbase.compressibleBase())
{
  sliceInit(bvbase, domain);
}

template <int Dim>
BrickViewBase<Dim>::
BrickViewBase(const This_t &bvbase, const Interval<Dim> &domain)
: domain_m(Pooma::NoInit()),  
  baseOffset_m(bvbase.baseOffset()), 
  compressibleBase_m(bvbase.compressibleBase())
{
  sliceInit(bvbase, Range<Dim>(domain));
}

template <int Dim>
BrickViewBase<Dim>::
BrickViewBase(const This_t &bvbase, bool compressible)
{
  *this = bvbase;
  compressibleBase_m = compressible;
  if (!compressible) restoreStrides();
}   

//-----------------------------------------------------------------------------
//
// sliceInit(const BrickViewBase<Dim,Dim,UStrd> &, 
//           const Range<Dim> &)
//
// Helper function used in taking a non-sliced view of a sliced view.
//
//-----------------------------------------------------------------------------

template <int Dim>
void
BrickViewBase<Dim>::
sliceInit(const This_t &bvbase, const Range<Dim> &domain)
{
  // Compute the strides and domain.
    
  for (int d = 0; d < Dim; ++d)
    {
      domain_m[d]   = Interval<1>(domain[d].length());
      strides_m[d]  = bvbase.ostrides_m[d] * domain[d].stride();
      baseOffset_m += domain[d].first() * bvbase.ostrides_m[d];
    }

  for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
}

//-----------------------------------------------------------------------------
//
// DoubleSliceHelper::init
//
// Helper functions used in initializing a slice of a slice.
//
//-----------------------------------------------------------------------------

template <int Dim, int Dim2>
void DoubleSliceHelper<Dim,Dim2>::
init(Interval<Dim> &domain, 
     int *strides,
     int &baseOffset,
     const BrickViewBase<Dim2> &bvbase,
     const SliceInterval<Dim2,Dim> &dom)
{  
  SliceRange<Dim2,Dim> tmp = dom;
  init(domain, strides, baseOffset, bvbase, tmp);
}

template <int Dim, int Dim2>
void DoubleSliceHelper<Dim,Dim2>::
init(Interval<Dim> &domain, 
     int *strides,
     int &baseOffset, 
     const BrickViewBase<Dim2> &bvbase,
     const SliceRange<Dim2,Dim> &dom)
{  
  // Compute the domain and strides.
  // The domain is an Interval with the length of each component
  // equal to the length of the corrsponding domain in the SliceRange.
  // The strides are calculated by multiplying the strides in 
  // the non-ignorable directions of the engine being viewed, by the
  // strides in the SliceRange that is doing the viewing. Since we must
  // skip over the ignorable dimensions in the engine being viewed,
  // we write this as a loop over the viewed engine's dimensions and
  // just do nothing for the ignorable ones.
    
  typedef typename SliceRange<Dim2,Dim>::TotalDomain_t TotalDomain_t;
  const TotalDomain_t &totDomain = dom.totalDomain();
    
  int d, dt;
  for (dt = 0, d = 0; dt < Dim2; ++dt)
    {
      if (!dom.ignorable(dt))
        {
          PAssert(d < Dim);
          domain[d]  = Interval<1>(totDomain[dt].length());
          strides[d] = bvbase.originalStrides()[dt] * totDomain[dt].stride();
          ++d;
        }
      baseOffset += totDomain[dt].first() * bvbase.originalStrides()[dt];
    }
  PAssert(d == Dim);
   
}

} // namespace Pooma

///////////////////////////////////////////////////////////////////////////////

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BrickBase.cpp,v $   $Author: richi $
// $Revision: 1.12 $   $Date: 2004/11/29 16:21:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
