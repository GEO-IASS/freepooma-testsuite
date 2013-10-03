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
// BrickViewBase test code
// This code tests the functionality of BrickViewBases that are created
// with SliceIntervals. 
//-----------------------------------------------------------------------------

#include "Engine/BrickBase.h"

#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Domain/SliceInterval.h"

//-----------------------------------------------------------------------------
// sliceTester functions
// 
// These functions test various properties of a slice of a brick. 
// I do not currently test the offset function as I don't know
// an easy way to do that in a dimension independent fashion.
//
// The sliceTest function tests a particular slice of a particular
// BrickBase.
//
// The sliceTester functions take a reference domain and a set
// of generic "domains", construct the appropriate BrickBase
// and construct the appropriate SliceInterval using NewDomainN,
// and then call sliceTest.
//-----------------------------------------------------------------------------

template <int Dim, int BaseDim>
void sliceTest(Pooma::Tester &t, 
               const Pooma::BrickBase<BaseDim> &A,
               const SliceInterval<BaseDim,Dim> &slice)
{
  Pooma::BrickViewBase<Dim> AV(A,slice);

  Interval<Dim> domain = Pooma::NoInit();
  
  int d;

  // First check the domain information: domain,
  // and ignorableDomain. These can be calculated directly 
  // from the slice domain.  

  for (d = 0; d < Dim; ++d)
    domain[d] = Interval<1>(slice.sliceDomain()[d].length());
    
  t.check(AV.domain() == domain);
  
  // Views are always zero-based...
  
  for (d = 0; d < Dim; ++d)
    {
      t.check(AV.first(d) == 0);
      t.check(AV.domain()[d].first() == 0);
    }

  // Finally we test the offset calculation. We've already 
  // confirmed that the strides array is right, so this is
  // simply a matter of summing up the offset in a particular
  // direction multiplied by the stride in that direction. 
  
  typedef typename Interval<Dim>::iterator Iterator_t;
  
  Iterator_t ploc = AV.domain().begin();
  while (ploc != AV.domain().end())
    {
      int off = 0;
      const Loc<Dim> &loc = *ploc;
      for (d = 0; d < Dim; ++d)
        off += loc[d].first() * AV.strides()[d];
      t.check(AV.offset(loc) == off);
      ++ploc;
    }
}

template <int BaseDim, class D1, class D2>
void sliceTester(Pooma::Tester &t, 
                 const Interval<BaseDim> &domain,
                 const D1 &d1, const D2 &d2)
{
  Pooma::BrickBase<BaseDim> A(domain);
  typename NewDomain2<D1,D2>::SliceType_t slice(domain,d1,d2);
  sliceTest(t,A,slice);
}
  
template <int BaseDim, class D1, class D2, class D3>
void sliceTester(Pooma::Tester &t, 
                 const Interval<BaseDim> &domain,
                 const D1 &d1, const D2 &d2, const D3 &d3)
{
  Pooma::BrickBase<BaseDim> A(domain);
  typename NewDomain3<D1,D2,D3>::SliceType_t slice(domain,d1,d2,d3);
  sliceTest(t,A,slice);
}
  
template <int BaseDim, class D1, class D2, class D3, class D4>
void sliceTester(Pooma::Tester &t, 
                 const Interval<BaseDim> &domain,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4)
{
  Pooma::BrickBase<BaseDim> A(domain);
  typename NewDomain4<D1,D2,D3,D4>::SliceType_t slice(domain,d1,d2,d3,d4);
  sliceTest(t,A,slice);
}
  
template <int BaseDim, class D1, class D2, class D3, class D4,
                       class D5>
void sliceTester(Pooma::Tester &t, 
                 const Interval<BaseDim> &domain,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4,
                 const D5 &d5)
{
  Pooma::BrickBase<BaseDim> A(domain);
  typename NewDomain5<D1,D2,D3,D4,D5>::SliceType_t 
    slice(domain,d1,d2,d3,d4,d5);
  sliceTest(t,A,slice);
}
  
template <int BaseDim, class D1, class D2, class D3, class D4,
                       class D5, class D6>
void sliceTester(Pooma::Tester &t, 
                 const Interval<BaseDim> &domain,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4,
                 const D5 &d5, const D6 &d6)
{
  Pooma::BrickBase<BaseDim> A(domain);
  typename NewDomain6<D1,D2,D3,D4,D5,D6>::SliceType_t 
    slice(domain,d1,d2,d3,d4,d5,d6);
  sliceTest(t,A,slice);
}
  
template <int BaseDim, class D1, class D2, class D3, class D4,
                       class D5, class D6, class D7>
void sliceTester(Pooma::Tester &t, 
                 const Interval<BaseDim> &domain,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4,
                 const D5 &d5, const D6 &d6, const D7 &d7)
{
  Pooma::BrickBase<BaseDim> A(domain);
  typename NewDomain7<D1,D2,D3,D4,D5,D6,D7>::SliceType_t 
    slice(domain,d1,d2,d3,d4,d5,d6,d7);
  sliceTest(t,A,slice);
}



int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif

    using Pooma::BrickBase;
    using Pooma::BrickViewBase;
    
    tester.out() << "\nTesting sliced BrickViewBase (single-slice only)." 
                 << std::endl;

// Check all combinations of views of BrickBase<D> for D > 1:

// First create BrickBase objects for D = 2, ..., 7.
// We'll take 1-6 dimensional views of these.

    Interval<1> L1(-5,5);
    Interval<2> L2(L1,L1);
    Interval<3> L3(L1,L1,L1);
    Interval<4> L4(L1,L1,L1,L1);
    Interval<5> L5(L1,L1,L1,L1,L1);
    Interval<6> L6(L1,L1,L1,L1,L1,L1);
    Interval<7> L7(L1,L1,L1,L1,L1,L1,L1);
    
    BrickBase<2> A2(L2);
    BrickBase<3> A3(L3);
    BrickBase<4> A4(L4);
    BrickBase<5> A5(L5);
    BrickBase<6> A6(L6);
    BrickBase<7> A7(L7);
    
// Now invoke the tester function for the domain and the
// various combinations of 0 and the domain I1, which is a 
// subset of the domain L1.

    Interval<1> I1(-1,1);
    
    sliceTester(tester,A2.domain(), 0,I1);
    sliceTester(tester,A2.domain(),I1, 0);

    sliceTester(tester,A3.domain(), 0, 0,I1);
    sliceTester(tester,A3.domain(), 0,I1, 0);
    sliceTester(tester,A3.domain(), 0,I1,I1);
    sliceTester(tester,A3.domain(),I1, 0, 0);
    sliceTester(tester,A3.domain(),I1, 0,I1);
    sliceTester(tester,A3.domain(),I1,I1, 0);

    sliceTester(tester,A4.domain(), 0, 0, 0,I1);
    sliceTester(tester,A4.domain(), 0, 0,I1, 0);
    sliceTester(tester,A4.domain(), 0, 0,I1,I1);
    sliceTester(tester,A4.domain(), 0,I1, 0, 0);
    sliceTester(tester,A4.domain(), 0,I1, 0,I1);
    sliceTester(tester,A4.domain(), 0,I1,I1, 0);
    sliceTester(tester,A4.domain(), 0,I1,I1,I1);
    sliceTester(tester,A4.domain(),I1, 0, 0, 0);
    sliceTester(tester,A4.domain(),I1, 0, 0,I1);
    sliceTester(tester,A4.domain(),I1, 0,I1, 0);
    sliceTester(tester,A4.domain(),I1, 0,I1,I1);
    sliceTester(tester,A4.domain(),I1,I1, 0, 0);
    sliceTester(tester,A4.domain(),I1,I1, 0,I1);
    sliceTester(tester,A4.domain(),I1,I1,I1, 0);

    sliceTester(tester,A5.domain(), 0, 0, 0, 0,I1);
    sliceTester(tester,A5.domain(), 0, 0, 0,I1, 0);
    sliceTester(tester,A5.domain(), 0, 0, 0,I1,I1);
    sliceTester(tester,A5.domain(), 0, 0,I1, 0, 0);
    sliceTester(tester,A5.domain(), 0, 0,I1, 0,I1);
    sliceTester(tester,A5.domain(), 0, 0,I1,I1, 0);
    sliceTester(tester,A5.domain(), 0, 0,I1,I1,I1);
    sliceTester(tester,A5.domain(), 0,I1, 0, 0, 0);
    sliceTester(tester,A5.domain(), 0,I1, 0, 0,I1);
    sliceTester(tester,A5.domain(), 0,I1, 0,I1, 0);
    sliceTester(tester,A5.domain(), 0,I1, 0,I1,I1);
    sliceTester(tester,A5.domain(), 0,I1,I1, 0, 0);
    sliceTester(tester,A5.domain(), 0,I1,I1, 0,I1);
    sliceTester(tester,A5.domain(), 0,I1,I1,I1, 0);
    sliceTester(tester,A5.domain(), 0,I1,I1,I1,I1);
    sliceTester(tester,A5.domain(),I1, 0, 0, 0, 0);
    sliceTester(tester,A5.domain(),I1, 0, 0, 0,I1);
    sliceTester(tester,A5.domain(),I1, 0, 0,I1, 0);
    sliceTester(tester,A5.domain(),I1, 0, 0,I1,I1);
    sliceTester(tester,A5.domain(),I1, 0,I1, 0, 0);
    sliceTester(tester,A5.domain(),I1, 0,I1, 0,I1);
    sliceTester(tester,A5.domain(),I1, 0,I1,I1, 0);
    sliceTester(tester,A5.domain(),I1, 0,I1,I1,I1);
    sliceTester(tester,A5.domain(),I1,I1, 0, 0, 0);
    sliceTester(tester,A5.domain(),I1,I1, 0, 0,I1);
    sliceTester(tester,A5.domain(),I1,I1, 0,I1, 0);
    sliceTester(tester,A5.domain(),I1,I1, 0,I1,I1);
    sliceTester(tester,A5.domain(),I1,I1,I1, 0, 0);
    sliceTester(tester,A5.domain(),I1,I1,I1, 0,I1);
    sliceTester(tester,A5.domain(),I1,I1,I1,I1, 0);

    sliceTester(tester,A6.domain(), 0, 0, 0, 0, 0,I1);
    sliceTester(tester,A6.domain(), 0, 0, 0, 0,I1, 0);
    sliceTester(tester,A6.domain(), 0, 0, 0, 0,I1,I1);
    sliceTester(tester,A6.domain(), 0, 0, 0,I1, 0, 0);
    sliceTester(tester,A6.domain(), 0, 0, 0,I1, 0,I1);
    sliceTester(tester,A6.domain(), 0, 0, 0,I1,I1, 0);
    sliceTester(tester,A6.domain(), 0, 0, 0,I1,I1,I1);
    sliceTester(tester,A6.domain(), 0, 0,I1, 0, 0, 0);
    sliceTester(tester,A6.domain(), 0, 0,I1, 0, 0,I1);
    sliceTester(tester,A6.domain(), 0, 0,I1, 0,I1, 0);
    sliceTester(tester,A6.domain(), 0, 0,I1, 0,I1,I1);
    sliceTester(tester,A6.domain(), 0, 0,I1,I1, 0, 0);
    sliceTester(tester,A6.domain(), 0, 0,I1,I1, 0,I1);
    sliceTester(tester,A6.domain(), 0, 0,I1,I1,I1, 0);
    sliceTester(tester,A6.domain(), 0, 0,I1,I1,I1,I1);
    sliceTester(tester,A6.domain(), 0,I1, 0, 0, 0, 0);
    sliceTester(tester,A6.domain(), 0,I1, 0, 0, 0,I1);
    sliceTester(tester,A6.domain(), 0,I1, 0, 0,I1, 0);
    sliceTester(tester,A6.domain(), 0,I1, 0, 0,I1,I1);
    sliceTester(tester,A6.domain(), 0,I1, 0,I1, 0, 0);
    sliceTester(tester,A6.domain(), 0,I1, 0,I1, 0,I1);
    sliceTester(tester,A6.domain(), 0,I1, 0,I1,I1, 0);
    sliceTester(tester,A6.domain(), 0,I1, 0,I1,I1,I1);
    sliceTester(tester,A6.domain(), 0,I1,I1, 0, 0, 0);
    sliceTester(tester,A6.domain(), 0,I1,I1, 0, 0,I1);
    sliceTester(tester,A6.domain(), 0,I1,I1, 0,I1, 0);
    sliceTester(tester,A6.domain(), 0,I1,I1, 0,I1,I1);
    sliceTester(tester,A6.domain(), 0,I1,I1,I1, 0, 0);
    sliceTester(tester,A6.domain(), 0,I1,I1,I1, 0,I1);
    sliceTester(tester,A6.domain(), 0,I1,I1,I1,I1, 0);
    sliceTester(tester,A6.domain(), 0,I1,I1,I1,I1,I1);
    sliceTester(tester,A6.domain(),I1, 0, 0, 0, 0, 0);
    sliceTester(tester,A6.domain(),I1, 0, 0, 0, 0,I1);
    sliceTester(tester,A6.domain(),I1, 0, 0, 0,I1, 0);
    sliceTester(tester,A6.domain(),I1, 0, 0, 0,I1,I1);
    sliceTester(tester,A6.domain(),I1, 0, 0,I1, 0, 0);
    sliceTester(tester,A6.domain(),I1, 0, 0,I1, 0,I1);
    sliceTester(tester,A6.domain(),I1, 0, 0,I1,I1, 0);
    sliceTester(tester,A6.domain(),I1, 0, 0,I1,I1,I1);
    sliceTester(tester,A6.domain(),I1, 0,I1, 0, 0, 0);
    sliceTester(tester,A6.domain(),I1, 0,I1, 0, 0,I1);
    sliceTester(tester,A6.domain(),I1, 0,I1, 0,I1, 0);
    sliceTester(tester,A6.domain(),I1, 0,I1, 0,I1,I1);
    sliceTester(tester,A6.domain(),I1, 0,I1,I1, 0, 0);
    sliceTester(tester,A6.domain(),I1, 0,I1,I1, 0,I1);
    sliceTester(tester,A6.domain(),I1, 0,I1,I1,I1, 0);
    sliceTester(tester,A6.domain(),I1, 0,I1,I1,I1,I1);
    sliceTester(tester,A6.domain(),I1,I1, 0, 0, 0, 0);
    sliceTester(tester,A6.domain(),I1,I1, 0, 0, 0,I1);
    sliceTester(tester,A6.domain(),I1,I1, 0, 0,I1, 0);
    sliceTester(tester,A6.domain(),I1,I1, 0, 0,I1,I1);
    sliceTester(tester,A6.domain(),I1,I1, 0,I1, 0, 0);
    sliceTester(tester,A6.domain(),I1,I1, 0,I1, 0,I1);
    sliceTester(tester,A6.domain(),I1,I1, 0,I1,I1, 0);
    sliceTester(tester,A6.domain(),I1,I1, 0,I1,I1,I1);
    sliceTester(tester,A6.domain(),I1,I1,I1, 0, 0, 0);
    sliceTester(tester,A6.domain(),I1,I1,I1, 0, 0,I1);
    sliceTester(tester,A6.domain(),I1,I1,I1, 0,I1, 0);
    sliceTester(tester,A6.domain(),I1,I1,I1, 0,I1,I1);
    sliceTester(tester,A6.domain(),I1,I1,I1,I1, 0, 0);
    sliceTester(tester,A6.domain(),I1,I1,I1,I1, 0,I1);
    sliceTester(tester,A6.domain(),I1,I1,I1,I1,I1, 0);

    sliceTester(tester,A7.domain(), 0, 0, 0, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0, 0,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0, 0,I1,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1, 0,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(), 0,I1,I1,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0, 0,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1, 0,I1,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1, 0,I1,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0,I1,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1, 0,I1,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1, 0, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1, 0, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1, 0,I1, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1, 0,I1,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1,I1, 0, 0);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1,I1, 0,I1);
    sliceTester(tester,A7.domain(),I1,I1,I1,I1,I1,I1, 0);

#if POOMA_EXCEPTIONS
  }
  catch(const char *err) 
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
  catch(const Pooma::Assertion &err)
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
#endif    
  int ret = tester.results("brickviewbase_test2");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brickviewbase_test2.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
