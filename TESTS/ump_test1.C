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
// ump test 1
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Layout/UniformGridLayout.h"
#include "Pooma/UMPArrays.h"
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  // Create the total domain.
  
  Interval<1> x(12);
  Interval<1> domain(x);
  
  // Create the block sizes.
  
  Loc<1> blocks(3), blocks2(4);

  // Create the partitioners.
  
  UniformGridPartition<1> partition(blocks), partition2(blocks2);
  
  // Create the layouts.
  
  UniformGridLayout<1> layout(domain, partition, ReplicatedTag());
  UniformGridLayout<1> layout2(domain, partition2, ReplicatedTag());

  Array<1, double, MultiPatch<UniformTag,Brick> > a(layout);
  Array<1, double, MultiPatch<UniformTag,Brick> > a2(layout2);
  Array<1, double, MultiPatch<UniformTag,CompressibleBrick> > ac(layout);  
  
  // Store some stuff.
  
  typedef Interval<1>::iterator Iterator;
  Iterator b0 = domain[0].begin();
  Iterator e0 = domain[0].end();
  for (Iterator i0=b0; i0!=e0; ++i0)
    a(i0->first()) = ac(i0->first()) = double(i0->first());

  tester.out() << a << std::endl;
  tester.out() << ac << std::endl;
    
  Array<1, double, BrickView> b = a(*(layout.beginGlobal() + 2));
  tester.out() << b << std::endl;
    
  // Create a view of a multipatch.
  
  Range<1> vdom(3,11,2), xdom(1,3,2), zdom; 
  Array<1, double, MultiPatchView<UniformTag,Brick,1> > av = a(vdom);
  tester.out() << av << std::endl;
  
  // Create an Intersector object and use it.
  
  Intersector<1> inter;
  inter.intersect(a.engine());
  inter.intersect(a2.engine());
  
  Intersector<1>::const_iterator p = inter.begin();
  while (p != inter.end())
    {
      tester.out() << a(*p) << ac(*p) << std::endl;
      ++p;
    }
    
  // Play with view layouts.
  
  UniformGridLayoutView<1,1> vlayout(layout, vdom);

#if POOMA_NO_OSTREAM_ITERATOR_1ARG
  std::ostream_iterator<Node<Range<1>,Interval<1> >, char> os(tester.out().stream(), "\n");
#else
  std::ostream_iterator<Node<Range<1>,Interval<1> > > os(tester.out().stream(), "\n");
#endif

  vlayout.touches(xdom, os);

  tester.out() << vlayout << std::endl;

  UniformGridLayoutView<1,1> vvlayout(vlayout, Interval<1>(1,2));
  tester.out() << vvlayout << std::endl;
  
  Interval<1> I(6);
  Interval<5> I5(I,I,I,I,I);
  Loc<5> blocks5(2,2,2,2,2);
  UniformGridPartition<5> partition5(blocks5);   
  UniformGridLayout<5> layout5(I5, partition5, ReplicatedTag());
  
  tester.out() << layout5 << std::endl;
  
  typedef NewDomain5<int, Range<1>, int, AllDomain<1>, Interval<1> > 
  NewDomain_t;
  typedef NewDomain_t::SliceType_t SliceType_t;
  SliceType_t VD;
  AllDomain<1> A;
  Interval<1> I1(1,3);
  Range<1> R(0,4,2);
  NewDomain_t::fillSlice(VD, layout5.domain(), 2, R, 1, A, I1);
  UniformGridLayoutView<3,5> vlayout3(layout5, VD);

  // vlayout3.domain() should be:     [0:2:1,0:4:1,0:2:1]
  // vlayout3.baseDomain() should be: [2:2:1,0:4:2,1:1:1,0:5:1,1:3:1]
  
  tester.out() << vlayout3 << std::endl;

  typedef NewDomain3<int, Range<1>, Interval<1> > NewDomain2_t;
  typedef NewDomain2_t::SliceType_t SliceType2_t;
  SliceType2_t VD2;
  Interval<1> I2(1,2);
  Range<1> R2(0,2,2);
  NewDomain2_t::fillSlice(VD2, vlayout3.domain(), 0, R2, I2);
  UniformGridLayoutView<2,5> vvlayout2(vlayout3, VD2);

  // vvlayout2.domain() should be:     [0:1:1,0:1:1]
  // vvlayout2.baseDomain() should be: [2:2:1,0:0:1,1:1:1,0:2:2,2:3:1]
  
  tester.out() << vvlayout2 << std::endl;

  int ret = tester.results("ump_test1");
  Pooma::finalize();
  return ret;
}
	     
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ump_test1.cpp,v $   $Author: richard $
// $Revision: 1.30 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
