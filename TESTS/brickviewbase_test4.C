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
// This version tests slices of slices using SliceInterval.
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

template <int Dim, int Dim2>
void sliceTest(Pooma::Tester &t, 
               const Pooma::BrickViewBase<Dim2> &A,
               const SliceInterval<Dim2,Dim> &slice)
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

template <int Dim, class D1, class D2>
void sliceTester(Pooma::Tester &t, 
                 const Pooma::BrickViewBase<Dim> &BV,
                 const D1 &d1, const D2 &d2)
{
  typename NewDomain2<D1,D2>::SliceType_t slice(BV.domain(),d1,d2);
  sliceTest(t,BV,slice);
}
  
template <int Dim, class D1, class D2, class D3>
void sliceTester(Pooma::Tester &t, 
                 const Pooma::BrickViewBase<Dim> &BV,
                 const D1 &d1, const D2 &d2, const D3 &d3)
{
  typename NewDomain3<D1,D2,D3>::SliceType_t slice(BV.domain(),d1,d2,d3);
  sliceTest(t,BV,slice);
}
  
template <int Dim, class D1, class D2, class D3, class D4>
void sliceTester(Pooma::Tester &t, 
                 const Pooma::BrickViewBase<Dim> &BV,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4)
{
  typename NewDomain4<D1,D2,D3,D4>::SliceType_t slice(BV.domain(),d1,d2,d3,d4);
  sliceTest(t,BV,slice);
}
  
template <int Dim, class D1, class D2, class D3, class D4,
                   class D5>
void sliceTester(Pooma::Tester &t, 
                 const Pooma::BrickViewBase<Dim> &BV,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4,
                 const D5 &d5)
{
  typename NewDomain5<D1,D2,D3,D4,D5>::SliceType_t 
    slice(BV.domain(),d1,d2,d3,d4,d5);
  sliceTest(t,BV,slice);
}
  
template <int Dim, class D1, class D2, class D3, class D4,
                   class D5, class D6>
void sliceTester(Pooma::Tester &t, 
                 const Pooma::BrickViewBase<Dim> &BV,
                 const D1 &d1, const D2 &d2, const D3 &d3, const D4 &d4,
                 const D5 &d5, const D6 &d6)
{
  typename NewDomain6<D1,D2,D3,D4,D5,D6>::SliceType_t 
    slice(BV.domain(),d1,d2,d3,d4,d5,d6);
  sliceTest(t,BV,slice);
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
    
    tester.out() 
      << "\nTesting double-sliced BrickViewBase with Interval slices." 
      << std::endl;

// First create BrickBase objects for D = 2, ..., 7.
// We'll take 1-6 dimensional views of these and then
// run slice tests for lots of permuations of slices of
// these views.

    Interval<1> L1(-5,5);
    Interval<2> L2(L1,L1);
    Interval<3> L3(L1,L1,L1);
    Interval<4> L4(L1,L1,L1,L1);
    Interval<5> L5(L1,L1,L1,L1,L1);
    Interval<6> L6(L1,L1,L1,L1,L1,L1);
    Interval<7> L7(L1,L1,L1,L1,L1,L1,L1);
    
    BrickBase<3> A3(L3);
    BrickBase<4> A4(L4);
    BrickBase<5> A5(L5);
    BrickBase<6> A6(L6);
    BrickBase<7> A7(L7);
    
    SliceInterval<3,2> SI32(A3.domain(), 0,L1,L1);
    SliceInterval<4,3> SI43(A4.domain(), 0,L1,L1,L1);
    SliceInterval<4,2> SI42(A4.domain(), 0, 0,L1,L1);
    SliceInterval<5,4> SI54(A5.domain(), 0,L1,L1,L1,L1);
    SliceInterval<5,3> SI53(A5.domain(), 0, 0,L1,L1,L1);
    SliceInterval<5,2> SI52(A5.domain(), 0, 0, 0,L1,L1);
    SliceInterval<6,5> SI65(A6.domain(), 0,L1,L1,L1,L1,L1);
    SliceInterval<6,4> SI64(A6.domain(), 0, 0,L1,L1,L1,L1);
    SliceInterval<6,3> SI63(A6.domain(), 0, 0, 0,L1,L1,L1);
    SliceInterval<6,2> SI62(A6.domain(), 0, 0, 0, 0,L1,L1);
    SliceInterval<7,6> SI76(A7.domain(), 0,L1,L1,L1,L1,L1,L1);
    SliceInterval<7,5> SI75(A7.domain(), 0, 0,L1,L1,L1,L1,L1);
    SliceInterval<7,4> SI74(A7.domain(), 0, 0, 0,L1,L1,L1,L1);
    SliceInterval<7,3> SI73(A7.domain(), 0, 0, 0, 0,L1,L1,L1);
    SliceInterval<7,2> SI72(A7.domain(), 0, 0, 0, 0, 0,L1,L1);
    
    BrickViewBase<2> AV32(A3,SI32);
    BrickViewBase<3> AV43(A4,SI43);
    BrickViewBase<2> AV42(A4,SI42);
    BrickViewBase<4> AV54(A5,SI54);
    BrickViewBase<3> AV53(A5,SI53);
    BrickViewBase<2> AV52(A5,SI52);
    BrickViewBase<5> AV65(A6,SI65);
    BrickViewBase<4> AV64(A6,SI64);
    BrickViewBase<3> AV63(A6,SI63);
    BrickViewBase<2> AV62(A6,SI62);
    BrickViewBase<6> AV76(A7,SI76);
    BrickViewBase<5> AV75(A7,SI75);
    BrickViewBase<4> AV74(A7,SI74);
    BrickViewBase<3> AV73(A7,SI73);
    BrickViewBase<2> AV72(A7,SI72);

// Now invoke the tester function for the first view and the
// various combinations of 0 and the domain I1, which is a 
// subset of the view's domain (Interval<1>(L1.length())).

    Interval<1> I1(4,8);
    
    sliceTester(tester,AV32, 0,I1);
    sliceTester(tester,AV32,I1, 0);
    sliceTester(tester,AV42, 0,I1);
    sliceTester(tester,AV42,I1, 0);
    sliceTester(tester,AV52, 0,I1);
    sliceTester(tester,AV52,I1, 0);
    sliceTester(tester,AV62, 0,I1);
    sliceTester(tester,AV62,I1, 0);
    sliceTester(tester,AV72, 0,I1);
    sliceTester(tester,AV72,I1, 0);

    sliceTester(tester,AV43, 0, 0,I1);
    sliceTester(tester,AV43, 0,I1, 0);
    sliceTester(tester,AV43, 0,I1,I1);
    sliceTester(tester,AV43,I1, 0, 0);
    sliceTester(tester,AV43,I1, 0,I1);
    sliceTester(tester,AV43,I1,I1, 0);
    sliceTester(tester,AV53, 0, 0,I1);
    sliceTester(tester,AV53, 0,I1, 0);
    sliceTester(tester,AV53, 0,I1,I1);
    sliceTester(tester,AV53,I1, 0, 0);
    sliceTester(tester,AV53,I1, 0,I1);
    sliceTester(tester,AV53,I1,I1, 0);
    sliceTester(tester,AV63, 0, 0,I1);
    sliceTester(tester,AV63, 0,I1, 0);
    sliceTester(tester,AV63, 0,I1,I1);
    sliceTester(tester,AV63,I1, 0, 0);
    sliceTester(tester,AV63,I1, 0,I1);
    sliceTester(tester,AV63,I1,I1, 0);
    sliceTester(tester,AV73, 0, 0,I1);
    sliceTester(tester,AV73, 0,I1, 0);
    sliceTester(tester,AV73, 0,I1,I1);
    sliceTester(tester,AV73,I1, 0, 0);
    sliceTester(tester,AV73,I1, 0,I1);
    sliceTester(tester,AV73,I1,I1, 0);

    sliceTester(tester,AV54, 0, 0, 0,I1);
    sliceTester(tester,AV54, 0, 0,I1, 0);
    sliceTester(tester,AV54, 0, 0,I1,I1);
    sliceTester(tester,AV54, 0,I1, 0, 0);
    sliceTester(tester,AV54, 0,I1, 0,I1);
    sliceTester(tester,AV54, 0,I1,I1, 0);
    sliceTester(tester,AV54, 0,I1,I1,I1);
    sliceTester(tester,AV54,I1, 0, 0, 0);
    sliceTester(tester,AV54,I1, 0, 0,I1);
    sliceTester(tester,AV54,I1, 0,I1, 0);
    sliceTester(tester,AV54,I1, 0,I1,I1);
    sliceTester(tester,AV54,I1,I1, 0, 0);
    sliceTester(tester,AV54,I1,I1, 0,I1);
    sliceTester(tester,AV54,I1,I1,I1, 0);
    sliceTester(tester,AV64, 0, 0, 0,I1);
    sliceTester(tester,AV64, 0, 0,I1, 0);
    sliceTester(tester,AV64, 0, 0,I1,I1);
    sliceTester(tester,AV64, 0,I1, 0, 0);
    sliceTester(tester,AV64, 0,I1, 0,I1);
    sliceTester(tester,AV64, 0,I1,I1, 0);
    sliceTester(tester,AV64, 0,I1,I1,I1);
    sliceTester(tester,AV64,I1, 0, 0, 0);
    sliceTester(tester,AV64,I1, 0, 0,I1);
    sliceTester(tester,AV64,I1, 0,I1, 0);
    sliceTester(tester,AV64,I1, 0,I1,I1);
    sliceTester(tester,AV64,I1,I1, 0, 0);
    sliceTester(tester,AV64,I1,I1, 0,I1);
    sliceTester(tester,AV64,I1,I1,I1, 0);
    sliceTester(tester,AV74, 0, 0, 0,I1);
    sliceTester(tester,AV74, 0, 0,I1, 0);
    sliceTester(tester,AV74, 0, 0,I1,I1);
    sliceTester(tester,AV74, 0,I1, 0, 0);
    sliceTester(tester,AV74, 0,I1, 0,I1);
    sliceTester(tester,AV74, 0,I1,I1, 0);
    sliceTester(tester,AV74, 0,I1,I1,I1);
    sliceTester(tester,AV74,I1, 0, 0, 0);
    sliceTester(tester,AV74,I1, 0, 0,I1);
    sliceTester(tester,AV74,I1, 0,I1, 0);
    sliceTester(tester,AV74,I1, 0,I1,I1);
    sliceTester(tester,AV74,I1,I1, 0, 0);
    sliceTester(tester,AV74,I1,I1, 0,I1);
    sliceTester(tester,AV74,I1,I1,I1, 0);

    sliceTester(tester,AV65, 0, 0, 0, 0,I1);
    sliceTester(tester,AV65, 0, 0, 0,I1, 0);
    sliceTester(tester,AV65, 0, 0, 0,I1,I1);
    sliceTester(tester,AV65, 0, 0,I1, 0, 0);
    sliceTester(tester,AV65, 0, 0,I1, 0,I1);
    sliceTester(tester,AV65, 0, 0,I1,I1, 0);
    sliceTester(tester,AV65, 0, 0,I1,I1,I1);
    sliceTester(tester,AV65, 0,I1, 0, 0, 0);
    sliceTester(tester,AV65, 0,I1, 0, 0,I1);
    sliceTester(tester,AV65, 0,I1, 0,I1, 0);
    sliceTester(tester,AV65, 0,I1, 0,I1,I1);
    sliceTester(tester,AV65, 0,I1,I1, 0, 0);
    sliceTester(tester,AV65, 0,I1,I1, 0,I1);
    sliceTester(tester,AV65, 0,I1,I1,I1, 0);
    sliceTester(tester,AV65, 0,I1,I1,I1,I1);
    sliceTester(tester,AV65,I1, 0, 0, 0, 0);
    sliceTester(tester,AV65,I1, 0, 0, 0,I1);
    sliceTester(tester,AV65,I1, 0, 0,I1, 0);
    sliceTester(tester,AV65,I1, 0, 0,I1,I1);
    sliceTester(tester,AV65,I1, 0,I1, 0, 0);
    sliceTester(tester,AV65,I1, 0,I1, 0,I1);
    sliceTester(tester,AV65,I1, 0,I1,I1, 0);
    sliceTester(tester,AV65,I1, 0,I1,I1,I1);
    sliceTester(tester,AV65,I1,I1, 0, 0, 0);
    sliceTester(tester,AV65,I1,I1, 0, 0,I1);
    sliceTester(tester,AV65,I1,I1, 0,I1, 0);
    sliceTester(tester,AV65,I1,I1, 0,I1,I1);
    sliceTester(tester,AV65,I1,I1,I1, 0, 0);
    sliceTester(tester,AV65,I1,I1,I1, 0,I1);
    sliceTester(tester,AV65,I1,I1,I1,I1, 0);
    sliceTester(tester,AV75, 0, 0, 0, 0,I1);
    sliceTester(tester,AV75, 0, 0, 0,I1, 0);
    sliceTester(tester,AV75, 0, 0, 0,I1,I1);
    sliceTester(tester,AV75, 0, 0,I1, 0, 0);
    sliceTester(tester,AV75, 0, 0,I1, 0,I1);
    sliceTester(tester,AV75, 0, 0,I1,I1, 0);
    sliceTester(tester,AV75, 0, 0,I1,I1,I1);
    sliceTester(tester,AV75, 0,I1, 0, 0, 0);
    sliceTester(tester,AV75, 0,I1, 0, 0,I1);
    sliceTester(tester,AV75, 0,I1, 0,I1, 0);
    sliceTester(tester,AV75, 0,I1, 0,I1,I1);
    sliceTester(tester,AV75, 0,I1,I1, 0, 0);
    sliceTester(tester,AV75, 0,I1,I1, 0,I1);
    sliceTester(tester,AV75, 0,I1,I1,I1, 0);
    sliceTester(tester,AV75, 0,I1,I1,I1,I1);
    sliceTester(tester,AV75,I1, 0, 0, 0, 0);
    sliceTester(tester,AV75,I1, 0, 0, 0,I1);
    sliceTester(tester,AV75,I1, 0, 0,I1, 0);
    sliceTester(tester,AV75,I1, 0, 0,I1,I1);
    sliceTester(tester,AV75,I1, 0,I1, 0, 0);
    sliceTester(tester,AV75,I1, 0,I1, 0,I1);
    sliceTester(tester,AV75,I1, 0,I1,I1, 0);
    sliceTester(tester,AV75,I1, 0,I1,I1,I1);
    sliceTester(tester,AV75,I1,I1, 0, 0, 0);
    sliceTester(tester,AV75,I1,I1, 0, 0,I1);
    sliceTester(tester,AV75,I1,I1, 0,I1, 0);
    sliceTester(tester,AV75,I1,I1, 0,I1,I1);
    sliceTester(tester,AV75,I1,I1,I1, 0, 0);
    sliceTester(tester,AV75,I1,I1,I1, 0,I1);
    sliceTester(tester,AV75,I1,I1,I1,I1, 0);

    sliceTester(tester,AV76, 0, 0, 0, 0, 0,I1);
    sliceTester(tester,AV76, 0, 0, 0, 0,I1, 0);
    sliceTester(tester,AV76, 0, 0, 0, 0,I1,I1);
    sliceTester(tester,AV76, 0, 0, 0,I1, 0, 0);
    sliceTester(tester,AV76, 0, 0, 0,I1, 0,I1);
    sliceTester(tester,AV76, 0, 0, 0,I1,I1, 0);
    sliceTester(tester,AV76, 0, 0, 0,I1,I1,I1);
    sliceTester(tester,AV76, 0, 0,I1, 0, 0, 0);
    sliceTester(tester,AV76, 0, 0,I1, 0, 0,I1);
    sliceTester(tester,AV76, 0, 0,I1, 0,I1, 0);
    sliceTester(tester,AV76, 0, 0,I1, 0,I1,I1);
    sliceTester(tester,AV76, 0, 0,I1,I1, 0, 0);
    sliceTester(tester,AV76, 0, 0,I1,I1, 0,I1);
    sliceTester(tester,AV76, 0, 0,I1,I1,I1, 0);
    sliceTester(tester,AV76, 0, 0,I1,I1,I1,I1);
    sliceTester(tester,AV76, 0,I1, 0, 0, 0, 0);
    sliceTester(tester,AV76, 0,I1, 0, 0, 0,I1);
    sliceTester(tester,AV76, 0,I1, 0, 0,I1, 0);
    sliceTester(tester,AV76, 0,I1, 0, 0,I1,I1);
    sliceTester(tester,AV76, 0,I1, 0,I1, 0, 0);
    sliceTester(tester,AV76, 0,I1, 0,I1, 0,I1);
    sliceTester(tester,AV76, 0,I1, 0,I1,I1, 0);
    sliceTester(tester,AV76, 0,I1, 0,I1,I1,I1);
    sliceTester(tester,AV76, 0,I1,I1, 0, 0, 0);
    sliceTester(tester,AV76, 0,I1,I1, 0, 0,I1);
    sliceTester(tester,AV76, 0,I1,I1, 0,I1, 0);
    sliceTester(tester,AV76, 0,I1,I1, 0,I1,I1);
    sliceTester(tester,AV76, 0,I1,I1,I1, 0, 0);
    sliceTester(tester,AV76, 0,I1,I1,I1, 0,I1);
    sliceTester(tester,AV76, 0,I1,I1,I1,I1, 0);
    sliceTester(tester,AV76, 0,I1,I1,I1,I1,I1);
    sliceTester(tester,AV76,I1, 0, 0, 0, 0, 0);
    sliceTester(tester,AV76,I1, 0, 0, 0, 0,I1);
    sliceTester(tester,AV76,I1, 0, 0, 0,I1, 0);
    sliceTester(tester,AV76,I1, 0, 0, 0,I1,I1);
    sliceTester(tester,AV76,I1, 0, 0,I1, 0, 0);
    sliceTester(tester,AV76,I1, 0, 0,I1, 0,I1);
    sliceTester(tester,AV76,I1, 0, 0,I1,I1, 0);
    sliceTester(tester,AV76,I1, 0, 0,I1,I1,I1);
    sliceTester(tester,AV76,I1, 0,I1, 0, 0, 0);
    sliceTester(tester,AV76,I1, 0,I1, 0, 0,I1);
    sliceTester(tester,AV76,I1, 0,I1, 0,I1, 0);
    sliceTester(tester,AV76,I1, 0,I1, 0,I1,I1);
    sliceTester(tester,AV76,I1, 0,I1,I1, 0, 0);
    sliceTester(tester,AV76,I1, 0,I1,I1, 0,I1);
    sliceTester(tester,AV76,I1, 0,I1,I1,I1, 0);
    sliceTester(tester,AV76,I1, 0,I1,I1,I1,I1);
    sliceTester(tester,AV76,I1,I1, 0, 0, 0, 0);
    sliceTester(tester,AV76,I1,I1, 0, 0, 0,I1);
    sliceTester(tester,AV76,I1,I1, 0, 0,I1, 0);
    sliceTester(tester,AV76,I1,I1, 0, 0,I1,I1);
    sliceTester(tester,AV76,I1,I1, 0,I1, 0, 0);
    sliceTester(tester,AV76,I1,I1, 0,I1, 0,I1);
    sliceTester(tester,AV76,I1,I1, 0,I1,I1, 0);
    sliceTester(tester,AV76,I1,I1, 0,I1,I1,I1);
    sliceTester(tester,AV76,I1,I1,I1, 0, 0, 0);
    sliceTester(tester,AV76,I1,I1,I1, 0, 0,I1);
    sliceTester(tester,AV76,I1,I1,I1, 0,I1, 0);
    sliceTester(tester,AV76,I1,I1,I1, 0,I1,I1);
    sliceTester(tester,AV76,I1,I1,I1,I1, 0, 0);
    sliceTester(tester,AV76,I1,I1,I1,I1, 0,I1);
    sliceTester(tester,AV76,I1,I1,I1,I1,I1, 0);

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
  int ret = tester.results("brickviewbase_test4");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brickviewbase_test4.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
