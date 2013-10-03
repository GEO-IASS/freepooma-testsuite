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
#include "Pooma/Pooma.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Touches.h"
#include "Domain/Split.h"
#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/EquivSubset.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include <iostream>
#include "Utilities/Tester.h"
int main(int argc, char *argv[]) 
{

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> n1(1,5);
  Interval<1> n2(4,8);
  Interval<1> n3(10,20);
  Interval<2> a(n1,n2);
  Interval<3> b(n1,n2,n3);

  Range<1> r1(1,5);
  Range<1> r2(4,8,2);
  Range<1> r3(5,9,2);
  Range<1> r4(10,20,5);
  Range<2> ra(r1,r2);
  Range<2> rb(r1,r3);
  Range<3> rc(r1,r2,r3);

  tester.out() << "1: touches(" << a[0] << "," << a[1] << ") ? ";
  tester.out() << touches(a[0], a[1]) << std::endl;
  tester.check( touches(a[0], a[1]) );
  tester.out() << "0: touches(" << a[0] << "," << b[2] << ") ? ";
  tester.out() << touches(a[0], b[2]) << std::endl;
  tester.check( touches(a[0], b[2])==0);
  tester.out() << "1: touches(" << a[0] << "," << ra[0] << ") ? ";
  tester.out() << touches(a[0], ra[0]) << std::endl;
  tester.check(touches(a[0], ra[0]));
  tester.out() << "1: touches(" << ra[0] << "," << ra[1] << ") ? ";
  tester.out() << touches(ra[0], ra[1]) << std::endl;
  tester.check( touches(ra[0], ra[1]));
  tester.out() << "0: touches(" << r2 << "," << r3 << ") ? ";
  tester.out() << touches(r2, r3) << std::endl;
  tester.check( touches(r2, r3)==0);
  tester.out() << "0: touches(" << ra << "," << rb << ") ? ";
  tester.out() << touches(ra, rb) << std::endl;
  tester.check(  touches(ra, rb) ==0);
  tester.out() << "1: touches(" << rc << "," << rc << ") ? ";
  tester.out() << touches(rc, rc) << std::endl;
  tester.check( touches(rc, rc) );
  tester.out() << "------------------------------------" << std::endl;

  tester.check(" touches ", true);

  Interval<1> c1(1,10);
  Interval<1> c2(3,8);
  Interval<1> c3(5,15);
  Interval<2> ca(c1, c1);
  Interval<2> cb(c1, c2);
  Range<1>    cr1(2,20,2);
  Range<1>    cr2(4,16,4);
  Range<1>    cr3(3,15,2);
  Range<1>    cr4(5,15,5);

  tester.out() << "1: contains(" << c1 << "," << c2 << ") ? ";
  tester.out() << contains(c1,c2) << std::endl;
  tester.check(contains(c1,c2));
  tester.out() << "0: contains(" << c2 << "," << c1 << ") ? ";
  tester.out() << contains(c2,c1) << std::endl;
  tester.check(contains(c2,c1)==0);
  tester.out() << "0: contains(" << c1 << "," << c3 << ") ? ";
  tester.out() << contains(c1,c3) << std::endl;
  tester.check(contains(c1,c3)==0);
  tester.out() << "1: contains(" << ca << "," << cb << ") ? ";
  tester.out() << contains(ca,cb) << std::endl;
  tester.check(contains(ca,cb));
  tester.out() << "0: contains(" << cb << "," << ca << ") ? ";
  tester.out() << contains(cb,ca) << std::endl;
  tester.check(contains(cb,ca)==0);
  tester.out() << "1: contains(" << cr1 << "," << cr2 << ") ? ";
  tester.out() << contains(cr1,cr2) << std::endl;
  tester.check( contains(cr1,cr2));
  tester.out() << "0: contains(" << cr1 << "," << cr3 << ") ? ";
  tester.out() << contains(cr1,cr3) << std::endl;
  tester.check(contains(cr1,cr3)==0);
  tester.out() << "1: contains(" << c3 << "," << cr4 << ") ? ";
  tester.out() << contains(c3,cr4) << std::endl;
  tester.check(contains(c3,cr4));
  tester.out() << "0: contains(" << cr4 << "," << c3 << ") ? ";
  tester.out() << contains(cr4,c3) << std::endl;
  tester.check(contains(cr4,c3)==0);
  tester.out() << "------------------------------------" << std::endl;

  Interval<2> s1, s2;
  Range<2>    sr1, sr2;

  split(cb, s1, s2);
  tester.out() << "split(" << cb << ") = " << s1 << " and " << s2 << std::endl;
  tester.check(s1==Interval<2>(Interval<1>(1,5),Interval<1>(3,5)));
  tester.check(s2==Interval<2>(Interval<1>(6,10),Interval<1>(6,8)));

 

  split(rb, sr1, sr2);
  tester.out() << "split(" << rb << ") = " << sr1 << " and " << sr2 << std::endl;
  tester.check(sr1==Range<2>(Range<1>(1,2),Range<1>(5,5,2)));
  tester.check(sr2==Range<2>(Range<1>(3,5),Range<1>(7,9,2)));

  tester.out() << "------------------------------------" << std::endl;

  tester.out() << "intersect(" << cb << "," << ca << ") = ";
  tester.out() << intersect(cb,ca) << std::endl;
  tester.check(intersect(cb,ca)==Interval<2>(Interval<1>(1,10),
					     Interval<1>(3,8)));

  tester.out() << "intersect(" << rb << "," << ra << ") = ";
  tester.out() << intersect(rb,ra) << std::endl;


  Range<1> i1(1,16,3);
  Range<1> i2(17,3,-2);
  tester.out() << "intersect(" << i1 << "," << i2 << ") = ";
  tester.out() << intersect(i1,i2) << std::endl;
  tester.check( intersect(i1,i2) == Range<1>(7,14,6));

  tester.out() << "intersect(" << i2 << "," << i1 << ") = ";
  tester.out() << intersect(i2,i1) << std::endl;
  tester.check( intersect(i2,i1) == Range<1>(13,7,-6));

  tester.out() << "------------------------------------" << std::endl;

  Interval<1> eq1(1,5);
  Range<1> eq2 = -2 * eq1 + 3;
  Range<1> eq3(-8,8,4);
  Range<1> eq4 = 3 * eq1;
  Range<1> eq5 = 2 * eq1;
  Range<1> eq6 = 6 * eq1 + 1;

  tester.out() << "For " << eq1 << " --> " << eq4 << ", then " << eq3 << " --> ";
  tester.out() << equivSubset(eq1,eq4,eq3) << std::endl;

  tester.out() << "For " << eq4 << " --> " << eq6 << ", then " << eq3 << " --> ";
  tester.out() << equivSubset(eq4,eq6,eq3) << std::endl;

  tester.out() << "For " << eq1 << " --> " << eq2 << ", then " << eq3 << " --> ";
  tester.out() << equivSubset(eq1,eq2,eq3) << std::endl;

  tester.out() << "------------------------------------" << std::endl;

  NewDomain3<Interval<1>, Interval<1>, int>::SliceType_t  ba;
  //  tester.out() << "Created initial slice domain ba = " << ba << std::endl;
  ba = NewDomain3<Interval<1>, Interval<1>, int>::combineSlice(ba,eq1,eq1,7);
  tester.out() << "After taking slice, ba = " << ba << std::endl;

  int retval = tester.results("Domain Calc");
  Pooma::finalize();
  return retval;


}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: domaincalc.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
