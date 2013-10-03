// -*- C++  -*-
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
// ----------------------------------------------------------------------
// Grid-based Multi-Patch Array's test 1
// ----------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Layout/GridLayout.h"
#include "Pooma/BrickArrays.h"
#include "Engine/MultiPatchEngine.h"
//#include "Evaluator/MultiPatchEval.h"
#include "iostream"
#include "Domain/Grid.h"
#include "Utilities/Tester.h"

#include <vector>
using std::vector;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  // Create the total domain.

  Interval<1> x(12);
  Interval<1> domain(x);
  
  // Create the block sizes.
  
  Loc<1> blocks(3), blocks2(4), blocks_5(5);

  tester.out() << "Initial domain = " << domain << ", blocks = " << blocks << std::endl;

  // Create the partitioners.
  
  GridPartition<1> partition(blocks), partition2(blocks2),partition_5(blocks_5);

  tester.out() << "Created GridPartition<1> p1 = " << partition << std::endl;
  tester.out() << "Created GridPartition<1> p2 = " << partition2 << std::endl;
  tester.out() << "Created GridPartition<1> p_5 = " << partition_5 << std::endl;

  // Create the layouts.
  
  GridLayout<1> layout(domain, partition, ReplicatedTag());
  GridLayout<1> layout2(domain, partition2, ReplicatedTag());
  GridLayout<1> layout_5(domain,partition_5, ReplicatedTag());

  tester.out() << "Created GridLayout<1> l1 = " << layout << std::endl;
  tester.out() << "Created GridLayout<1> l2 = " << layout2 << std::endl;
  tester.out() << "Created GridLayout<1> l_5 = " << layout_5 << std::endl;
 
  Range<1> range(domain[0].first(), domain[0].last()+1,
		 domain[0].size()/blocks[0].first());
  Grid<1> grid(range);
  tester.out() << "Created Grid<1> = " << grid << std::endl;
  
  Range<1> range2(domain[0].first(), domain[0].last()+1,
	          domain[0].size()/blocks2[0].first());
  Grid<1> grid2(range2);

  Array<1, double, MultiPatch<GridTag,Brick> > a(layout);
  Array<1, double, MultiPatch<GridTag,Brick> > a2(layout2);
  Array<1, double, MultiPatch<GridTag,CompressibleBrick> > ac(layout);  
  
  Array<1, double, MultiPatch<GridTag,Brick> > g(layout);
  Array<1, double, MultiPatch<GridTag,Brick> > g2(layout2);
  Array<1, double, MultiPatch<GridTag,CompressibleBrick> > gc(layout);  
 
  // Store some stuff.
  
  typedef Interval<1>::iterator Iterator;
  Iterator b0 = domain[0].begin();
  Iterator e0 = domain[0].end();
  for (Iterator i0=b0; i0!=e0; ++i0)
  {
    int j=i0->first();
    a(j) = ac(j) = g(j) = gc(j) = double(j);
  }

  tester.out()  << a << std::endl;
  tester.out()  << ac << std::endl;

  tester.out()  << g << std::endl;
  tester.out()  << gc << std::endl;

  Array<1, double, BrickView> b = a(*(layout.beginGlobal() + 2));
  tester.out()  << " view b "<<std::endl;
  tester.out()  << b << std::endl;

  Array<1, double, BrickView> gb = g(*(layout.beginGlobal() + 2));
  tester.out()  << " view gb "<<std::endl;
  tester.out()  << gb << std::endl;
    
  // Create a view of a multipatch.
  
  Range<1> vdom(3,11,2), xdom(1,3,2), zdom; 
  Array<1, double, MultiPatchView<GridTag,Brick,1> > av = a(vdom);
  tester.out()  << " view av " <<std::endl;
  tester.out()  << av << std::endl;

  Array<1, double, MultiPatchView<GridTag,Brick,1> > gv = g(vdom);
  tester.out()  << " view gv " <<std::endl;
  tester.out()  << gv << std::endl;

  // Create an Intersector object and use it.

  Intersector<1> inter;
  inter.intersect(a.engine());
  inter.intersect(a2.engine());

  Intersector<1>::const_iterator p = inter.begin();
  tester.out()  << " intersect " <<std::endl;
  while (p != inter.end())
    {
      tester.out()  << a(*p) << ac(*p) << std::endl;
      ++p;
    }
    
  // Play with view layouts.  
  typedef Node<Interval<1> >         Node_t1;
  vector<Node_t1> domains1;
 
  GridLayoutView<1,1> vlayout(layout, vdom);
  
#if POOMA_NO_OSTREAM_ITERATOR_1ARG
  std::ostream_iterator<Node<Range<1>,Interval<1> >, char> os(tester.out().stream(), "\n");
#else
  std::ostream_iterator<Node<Range<1>,Interval<1> > > os(tester.out().stream(), "\n");
#endif

  vlayout.touches(xdom, os );

  tester.out()  << " before write of vlayout " <<std::endl;

  tester.out()  << vlayout << std::endl;

  GridLayoutView<1,1> vvlayout(vlayout, Interval<1>(1,2));

  tester.out()  << " before write of vvlayout " <<std::endl;

  tester.out()  << vvlayout << std::endl;
  

  Interval<1> I(6);
  Array<1,int> al(3);
  al(0)=0;al(1)=4,al(2)=6;
  IndirectionList<int> il(al);


  Interval<4> I4(I,I,I,I);
  Loc<4> blocks4(2,2,2,2);
  Grid<4> grid4(il,il,il,il);
  
  Interval<5> I5(I,I,I,I,I);
  Loc<5> blocks5(2,2,2,2,2);
  Grid<5> grid5(il,il,il,il,il);

  Interval<6> I6(I,I,I,I,I,I);
  Loc<6> blocks6(2,2,2,2,2,2);
  Grid<6> grid6(il,il,il,il,il,il);

  Interval<7> I7(I,I,I,I,I,I,I);
  Loc<7> blocks7(2,2,2,2,2,2,2);
  Grid<7> grid7(il,il,il,il,il,il,il);

  GridPartition<4> partition4(grid4);
  GridPartition<5> partition5(grid5);
  GridPartition<6> partition6(grid6);
  GridPartition<7> partition7(grid7);
  
  GridLayout<4> layout4(I4, partition4, ReplicatedTag());
  GridLayout<5> layout5(I5, partition5, ReplicatedTag());
  GridLayout<6> layout6(I6, partition6, ReplicatedTag());
  GridLayout<7> layout7(I7, partition7, ReplicatedTag());
  
  tester.out()  << layout4 << std::endl;
  tester.out()  << layout5 << std::endl;
  tester.out()  << layout6 << std::endl;
  tester.out()  << layout7 << std::endl;

  Interval<4> t4(5,5,5,5);
  Interval<5> t5(5,5,5,5,5);
  Interval<6> t6(5,5,5,5,5,5);
  Interval<7> t7(5,5,5,5,5,5,5);
  
#if POOMA_NO_OSTREAM_ITERATOR_1ARG
  std::ostream_iterator<Node<Interval<4> >, char> osd(tester.out().stream(), "\n");
  std::ostream_iterator<Node<Interval<5> >, char> ose(tester.out().stream(), "\n");
  std::ostream_iterator<Node<Interval<6> >, char> osf(tester.out().stream(), "\n");
  std::ostream_iterator<Node<Interval<7> >, char> osg(tester.out().stream(), "\n");
#else
  std::ostream_iterator<Node<Interval<4> > > osd(tester.out().stream(), "\n");
  std::ostream_iterator<Node<Interval<5> > > ose(tester.out().stream(), "\n");
  std::ostream_iterator<Node<Interval<6> > > osf(tester.out().stream(), "\n");
  std::ostream_iterator<Node<Interval<7> > > osg(tester.out().stream(), "\n");
#endif

  layout4.touches(t4, osd);

  layout5.touches(t5, ose);

  layout6.touches(t6, osf);

  layout7.touches(t7, osg);



  typedef NewDomain5<int, Range<1>, int, AllDomain<1>, Interval<1> > 
  NewDomain_t;
  typedef NewDomain_t::SliceType_t SliceType_t;
  SliceType_t VD;
  AllDomain<1> A;
  Interval<1> I1(1,3);
  Range<1> R(0,4,2);
  NewDomain_t::fillSlice(VD, layout5.domain(), 2, R, 1, A, I1);
  GridLayoutView<3,5> vlayout3(layout5, VD);

  // vlayout3.domain() should be:     [0:2:1,0:4:1,0:2:1]
  // vlayout3.baseDomain() should be: [2:2:1,0:4:2,1:1:1,0:5:1,1:3:1]
  
  tester.out()  << vlayout3 << std::endl;

  typedef NewDomain3<int, Range<1>, Interval<1> > NewDomain2_t;
  typedef NewDomain2_t::SliceType_t SliceType2_t;
  SliceType2_t VD2;
  Interval<1> I2(1,2);
  Range<1> R2(0,2,2);
  NewDomain2_t::fillSlice(VD2, vlayout3.domain(), 0, R2, I2);
  GridLayoutView<2,5> vvlayout2(vlayout3, VD2);

  // vvlayout2.domain() should be:     [0:1:1,0:1:1]
  // vvlayout2.baseDomain() should be: [2:2:1,0:0:1,1:1:1,0:2:2,2:3:1]
  
  tester.out()  << vvlayout2 << std::endl;
  tester.out()  << "\nAll Done!" << std::endl;
  
  int ret = tester.results("gmp_test1");
  Pooma::finalize();
  return ret;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: gmp_test1.cpp,v $   $Author: richard $
// $Revision: 1.20 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
