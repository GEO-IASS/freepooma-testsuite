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

#include <iostream>

#include "Pooma/Pooma.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/IndirectionList.h"
#include "Domain/Grid.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include "Domain/Region.h"
#include "Domain/Touches.h"
#include "Domain/Contains.h"
#include "Domain/Split.h"
#include "Domain/Intersect.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[]) 
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << "Starting domain test." << std::endl << std::endl;

  Interval<5> i5(90,90,90,90,90);
  tester.out() << " size of 90^5 is" << i5.size()<<std::endl;
  Interval<5> ii5(101,101,101,101,101);
  tester.out() << " size of 101^5 is" << ii5.size()<<std::endl;


  Loc<1> loc1(9);
  Loc<5> loc5(9);
  Interval<1> interval1(1,10);
  Interval<5> interval5(interval1,interval1,interval1,interval1,interval1);
  Range<1> r1(1,9,2);
  Range<5> r5(r1,r1,r1,r1,r1);
  
  Loc<1> foo = loc1 + loc1;

  tester.check("Loc<1> addition ",foo==Loc<1>(18));

  foo = loc1 - loc1;

  tester.check("Loc<1> subtraction ",foo==Loc<1>(0));

  Loc<5> goo = loc1 + loc5;
 
  tester.check("Loc<1> + Loc<5> ",goo==Loc<5>(18));

  goo = loc1 - loc5;
  tester.check("Loc<1> - Loc<5> ",goo==Loc<5>(0));

  goo = loc5 - loc1;

  tester.check("Loc<5> - Loc<1> ",goo==Loc<5>(0));

  Interval<1> ir1 = interval1 + loc1;

  tester.check("Interval<1> + Loc<1> ",ir1 == Interval<1>(10,19));

  ir1 = loc1 + interval1;

  tester.check("Loc<1> + Interval<1> ",ir1 == Interval<1>(10,19));

  Range<1> range1(1,13,3);


  Interval<5> zero = interval5 - loc5;
  Interval<1> t(-8,1);
  tester.check("Interval<5> - Loc<1> ",zero == Interval<5>(t,t,t,t,t) );


  Range<5> rzero = loc5 - interval5;
  Range<1> tr(8,-1,-1);
  tester.check("Loc<1> - Interval<1> ",rzero == Range<5>(tr,tr,tr,tr,tr));


  Range<5> rprod = loc5 * interval5;
  tr=Range<1>(9,90,9);
  tester.check(" Loc<5> * Interval<5> ",rprod == Range<5>(tr,tr,tr,tr,tr));

  rprod = interval5 * loc5;
  tester.check(" Interval<5> * loc<5> ",rprod == Range<5>(tr,tr,tr,tr,tr));

  rprod = interval5 * loc1;
  tester.check(" Interval<5> * Loc<1> ",rprod == Range<5>(tr,tr,tr,tr,tr)); 

  rprod = loc1 * interval5;
  tester.check(" Loc<1> * Interval<5> ",rprod == Range<5>(tr,tr,tr,tr,tr));

  Range<1> rprod1 = range1 * 5;
  tester.check(" Range<1> * (int) 5",rprod1 == Range<1>(5,65,15));

  rprod1 = range1 * loc1;
  tester.check(" Range<1> * Loc<1> ",rprod1 == Range<1>(9,117,27));

  Range<5> sprod = interval5 * 9 ;
  sprod = 9 * interval5;
  tester.check(" Interval<5> * (int) 9 ", sprod == Range<5>(tr,tr,tr,tr,tr));


  Loc<1> slprod = loc1 * 9;
  slprod = 9 * loc1;
  tester.check(" Loc<1> * (int) 9 ", slprod == Loc<1>(81));

  slprod = loc1 / 9;
  tester.check("Loc<1> / (int) 9 " ,slprod == Loc<1>(1));

  Loc<5> sl5prod = loc5 * 9;
  tester.check("Loc<5> * (int) 9 ", sl5prod == Loc<5>(81,81,81,81,81));

  sl5prod = 9 * loc5;
  tester.check(" (int) 9 * Loc<5>  ", sl5prod == Loc<5>(81,81,81,81,81));
 
  sl5prod = loc5 / 9;
  tester.check("Loc<5> / (int) 9 ", sl5prod == Loc<5>(1,1,1,1,1));

  Range<5> sdiv = interval5 / .5;
  Range<1> ttt(2,20,2);
  tester.check(" Interval<5> / 0.5 ",sdiv == Range<5>(ttt,ttt,ttt,ttt,ttt));
  
 

  sdiv = interval5 / loc1 ;
  tester.out()<<sdiv<<std::endl;
  
  sdiv = interval5 / loc5;
  tester.out()<<sdiv<<std::endl;

  Range<1> s1div = range1 / 3;
  tester.out()<< s1div <<std::endl;

  s1div = range1 / loc1;
  tester.out()<< s1div <<std::endl;
 
  Range<1> rr1 = r1 + loc1;
  tester.check("Range<1> + Loc<1> ",rr1 == Range<1>(10,19,2) );

  rr1 = loc1 + r1;
  tester.check(" Loc<1> + Range<1> ",rr1 == Range<1>(10,19,2) );
 

  tester.out() << "Testing NewDomain<*> combine methods:" << std::endl;
  tester.out() << "-------------------------------------" << std::endl;
  {
    Interval<3> t1(Interval<1>(0,0),Interval<1>(0,2),Interval<1>(0,4));
     
    tester.check("  NewDomain3<int,int,int>::combine(1,3,5)",
		 NewDomain3<int,int,int>::combine(1,3,5) == t1);

    Interval<3> t2(Interval<1>(2,2),Interval<1>(4,4),Interval<1>(0,5));

    tester.check( "  NewDomain2<Loc<2>,int>::combine(Loc<2>(2,4),6) ",
		  NewDomain2<Loc<2>,int>::combine(Loc<2>(2,4),6)==t2);

    Range<7> t3(Interval<1>(12),Range<1>(5),Range<1>(10),Range<1>(15),
		Interval<1>(0,0),Interval<1>(0,1),Interval<1>(20));


    tester.out() <<t3<<std::endl;
    tester.out() <<  NewDomain4<int,Range<3>,Interval<2>,Interval<1> >::
      combine(12,Range<3>(Range<1>(5),Range<1>(10),Range<1>(15)),
	      Interval<2>(1,2), Interval<1>(20)) <<std::endl;

    tester.check(
		 "  NewDomain4<int,Range<3>,Interval<2>,Interval<1> >::combine",
		 NewDomain4<int,Range<3>,Interval<2>,Interval<1> >::
                 combine(12,Range<3>(Range<1>(5),Range<1>(10),Range<1>(15)),
                         Interval<2>(1,2), Interval<1>(20)) == t3);
  }

  tester.out() << std::endl;
  tester.out() << "Testing Loc<N>:" << std::endl;
  tester.out() << "---------------" << std::endl;
  {
    Loc<1> a;

    a = 3;
    tester.check("  after a = 3 : a = " , a == Loc<1>(3));

    Loc<2> b(a,a);
    tester.check("  2D Loc<2> b(a,a) = " , b == Loc<2>(Loc<1>(3),Loc<1>(3)) );

    b[0] = 2;
    tester.check("  after b[0] = 2 : b = ", b ==Loc<2>(Loc<1>(2),Loc<1>(3)) );

    Loc<2> bb = b;

    b += b;
    tester.check("  after b += b : b = ", b == bb*2);

    b += a;
    tester.check("  after b += a : b = " , b == Loc<2>(Loc<1>(7),Loc<1>(9)) );

    Loc<2> foo =  2 + b*a - 3;
    tester.check("  result of 2 + b * a - 3 = ",Loc<2>(Loc<1>(20),Loc<1>(26)) == foo );

    tester.check("  result of (b == b) = ", (b == b) );
    tester.check("  result of (b != b) = ", !(b != b) );

    Loc<3> c(b,10);
    tester.check("  3D Loc<3> c(b,10) = " , c == Loc<3>(Loc<1>(7),Loc<1>(9),Loc<1>(10)) );

    tester.check("  c[1].length() = " , c[1].length() == 1 );
    tester.check("  -c = ",-c == Loc<3>(Loc<1>(-7),Loc<1>(-9),Loc<1>(-10)) );
    tester.check("  results of ++c = ", (++c) == Loc<3>(Loc<1>(8),Loc<1>(10),Loc<1>(11)) );

    Loc<1>::iterator locitera = c[2].begin();
    Loc<1>::iterator lociterb = c[2].end();
    tester.out()<<"  c[2].begin = " << *locitera <<std::endl;

    tester.out()<<"  Iterating over c[2]: values = ";
    while (locitera != lociterb)
      tester.out()<<*locitera++ << " ";
    tester.out()<<std::endl;

  
    long val1 = 3;
    char val2 = 7;
    Loc<2> typesloc(val1, val2);
    tester.check("  Creating Loc from long and char: Loc<2>(3L, 7c) " ,typesloc == Loc<2>(Loc<1>(3),Loc<1>(7)) );

  
    typesloc = Loc<1>(4);
    tester.check("  Setting the above 2D Loc to Loc<1>(4) ",typesloc == Loc<2>(Loc<1>(4),Loc<1>(4))  );

    short val3 = 8;
    typesloc = val3;
    tester.check("  Setting the above 2D Loc to (short)8 ", typesloc == Loc<2>(Loc<1>(8),Loc<1>(8)) );
  }

  tester.out() << std::endl;
  tester.out() << "Testing Interval<N>:" << std::endl;
  tester.out() << "--------------------" << std::endl;
  {
    Interval<1> a;
   
    a = 3;
    tester.out() << "  after a = 3 : a = " << a << std::endl;

    Interval<2> b(a,a);
    tester.out() << "  2D Interval<2> b(a,a) = " << b << std::endl;

    b[0] = Interval<1>(2,5);
    tester.out() << "  after b[0] = (2,5) : b = " << b << std::endl;

    b += Loc<2>(1,2);
    tester.out() << "  after b += Loc<2>(1,2) : b = " << b << std::endl;

    //    b *= 2;
    //    tester.out() << "  after b *= 2 : b = " << b << std::endl;

    tester.out() << "  result of 2 + b - 3 = " << 2 + b - 3 << std::endl;

    tester.out() << "  result of (b == b) = " << (b == b) << std::endl;
    tester.out() << "  result of (b != b) = " << (b != b) << std::endl;
    tester.out() << "  result of (b  < b) = " << (b  < b) << std::endl;
    tester.out() << "  result of (b >= b) = " << (b >= b) << std::endl;

    Interval<3> c(10,b);
    tester.out() << "  3D Interval<3> c(10,b) = " << c << std::endl;

    tester.out() << "  c[1].length() = " << c[1].length() << std::endl;
    tester.out() << "  -c = " << -c << std::endl;
    tester.out() << "  results of ++c = " << ++c << std::endl;

    Interval<1>::iterator locitera = c[1].begin();
    Interval<1>::iterator lociterb = c[1].end();
    tester.out() << "  c[1].begin = " << *locitera << std::endl;

    tester.out() << "  Iterating over c[1]: values = ";
    while (locitera != lociterb)
      tester.out() << *locitera++ << " ";
    tester.out() << std::endl;

    tester.out()
      << "  Creating Interval from long and char: Interval<2>(3L, 7c) = ";
    long val1 = 3;
    char val2 = 7;
    Interval<2> typesloc(val1, val2);
    tester.out() << typesloc << std::endl;

    tester.out() << "  Setting the above 2D Interval to Loc<1>(4) = ";
    typesloc = Loc<1>(4);
    tester.out() << typesloc << std::endl;

    tester.out() << "  Setting the above 2D Interval to (short)8 = ";
    short val3 = 8;
    typesloc = val3;
    tester.out() << typesloc << std::endl;

    tester.out() << "  firsts for this domain  = " << typesloc.firsts()
		 << std::endl;
    tester.out() << "  lasts for this domain   = " << typesloc.lasts()
		 << std::endl;
    tester.out() << "  strides for this domain = " << typesloc.strides()
		 << std::endl;
    tester.out() << "  lengths for this domain = " << typesloc.lengths()
		 << std::endl;
    tester.out() << "  mins for this domain    = " << typesloc.mins()
		 << std::endl;
    tester.out() << "  maxes for this domain   = " << typesloc.maxes()
		 << std::endl;
  }

  tester.out() << std::endl;
  tester.out() << "Testing Range<N>:" << std::endl;
  tester.out() << "--------------------" << std::endl;
  {
    Range<1> a;
    
    a = 3;
    tester.out() << "  after a = 3 : a = " << a << std::endl;

    Range<2> b(a,a);
    tester.out() << "  2D Range<2> b(a,a) = " << b << std::endl;

    b[0] = Range<1>(2,5);
    tester.out() << "  after b[0] = (2,5) : b = " << b << std::endl;

    b += Loc<2>(1,2);
    tester.out() << "  after b += Loc<2>(1,2) : b = " << b << std::endl;

    b *= 2;
    tester.out() << "  after b *= 2 : b = " << b << std::endl;

    tester.out() << "  result of 2 + b - 3 = " << 2 + b - 3 << std::endl;

    tester.out() << "  result of (b == b) = " << (b == b) << std::endl;
    tester.out() << "  result of (b != b) = " << (b != b) << std::endl;
    tester.out() << "  result of (b  < b) = " << (b  < b) << std::endl;
    tester.out() << "  result of (b >= b) = " << (b >= b) << std::endl;

    Range<3> c(10,b);
    tester.out() << "  3D Range<3> c(10,b) = " << c << std::endl;

    tester.out() << "  c[1].length() = " << c[1].length() << std::endl;
    tester.out() << "  -c = " << -c << std::endl;
    tester.out() << "  results of ++c = " << ++c << std::endl;

    Range<1>::iterator locitera = c[1].begin();
    Range<1>::iterator lociterb = c[1].end();
    tester.out() << "  c[1].begin = " << *locitera << std::endl;

    tester.out() << "  Iterating over c[1]: values = ";
    while (locitera != lociterb)
      tester.out() << *locitera++ << " ";
    tester.out() << std::endl;

    tester.out() << "  checking b = " << b << std::endl;
    Range<4> d(Interval<2>(a,Interval<1>(5,10)), b);
    tester.out() << "  4D Range<4> d(a,(5,10),b) = " << d << std::endl;

    Range<6> e(Interval<2>(a,Interval<1>(5,10)), b, b);
    tester.out() << "  6D Range<6> e(a,(5,10),b,b) = " << e << std::endl;
    tester.out() << "  6D Range<6> f(b,d) = " << Range<6>(b,d) << std::endl;

    tester.out() << "  Creating Range from long and char: Range<2>(3L, 7c) = ";
    long val1 = 3;
    char val2 = 7;
    Range<2> typesloc(val1, val2);
    tester.out() << typesloc << std::endl;

    tester.out() << "  Setting the above 2D Range to Loc<1>(4) = ";
    typesloc = Loc<1>(4);
    tester.out() << typesloc << std::endl;

    tester.out() << "  Setting the above 2D Range to (short)8 = ";
    short val3 = 8;
    typesloc = val3;
    tester.out() << typesloc << std::endl;
  }

  tester.out() << std::endl;
  tester.out() << "Testing SliceInterval<N>:" << std::endl;
  tester.out() << "-------------------------" << std::endl;
  {
    SliceInterval<2,1> a;
  
    Interval<1> b1(1,5);
    Interval<1> b2(8,9);
    NewDomain2<Interval<1>,Interval<1> >::SliceType_t b;
    b = NewDomain2<Interval<1>,Interval<1> >::combineSlice(b, b1, b2);
    tester.out() << "  combineSlice b(1:5,8:9) = " << b << std::endl;

    NewDomain2<int,Interval<1> >::SliceType_t bs;
    bs = NewDomain2<int,Interval<1> >::combineSlice(bs,7,b2);
    tester.out() << "  combineSlice bs(7,8:9) = " << bs << std::endl;

    NewDomain4<int,Interval<2>,int,Interval<1> >::SliceType_t b3s;
    b3s = NewDomain4<int,Interval<2>,int,Interval<1> >::
      combineSlice(b3s,7,b,2,b1);
    tester.out() << "  combineSlice bs(7,1:5,8:9,2,1:5) = " << b3s << std::endl;
  }

  tester.out() << std::endl;
  tester.out() << "Testing SliceRange<N>:" << std::endl;
  tester.out() << "-------------------------" << std::endl;
  {
    SliceRange<3,1> a;
  
    Interval<1> b1(1,5);
    Range<1>    b2(2,8,2);
    NewDomain2<Interval<1>,Range<1> >::SliceType_t b;
    b = NewDomain2<Interval<1>,Range<1> >::combineSlice(b,b1,b2);
    tester.out() << "  combineSlice b(1:5,2:8:2) = " << b << std::endl;

    NewDomain3<int,Range<1>,int>::SliceType_t bs;
    bs = NewDomain3<int,Range<1>,int>::combineSlice(bs,7,b2,3);
    tester.out() << "  combineSlice bs(7,2:8:2,3) = " << bs << std::endl;

    NewDomain4<Loc<1>,Range<2>,int,Interval<1> >::SliceType_t b3s;
    b3s = NewDomain4<Loc<1>,Range<2>,int,Interval<1> >::
      combineSlice(b3s,Loc<1>(7),b,2,b1);
    tester.out() << "  combineSlice bs(7,1:5,2:8:2,2,1:5) = " << b3s
		 << std::endl;
  }

  tester.out() << std::endl;
  tester.out() << "Testing Region<N,double>:" << std::endl;
  tester.out() << "-------------------------" << std::endl;
  {
    Region<3> a;
  
    Region<1> a1(3);
    tester.out() << "  Region<1>(3) a1 = " << a1 << std::endl;

    Region<2> b(2.0, Region<1,double>(1.0, 1.5));
    tester.out() << "  Region<2>(2.0, Region<1>(1.0, 1.5)) b = " << b
		 << std::endl;

    Interval<1> b1(1,5);
    Range<1>    b2(2,8,2);
    tester.out() << "  combine(b, 2:8:2) = ";
    tester.out() << NewDomain2<Region<2>,Range<1> >::combine(b,b2) << std::endl;

    NewDomain2<Interval<1>,Range<1> >::fill(b, b1, b2);
    tester.out() << "  fill(1:5, 2:8:2) = " << b << std::endl;

    b *= 2;
    tester.out() << "  b *= 2 ==> b = " << b << std::endl;

    b += Loc<2>(3,4);

    tester.out() << "  b += Loc<2>(3,4) ==> b = " << b << std::endl;

    tester.out() << "  result of (b == b) = " << (b == b) << std::endl;
    tester.out() << "  result of (b != b) = " << (b != b) << std::endl;
    tester.out() << "  result of (b  < b) = " << (b  < b) << std::endl;
    tester.out() << "  result of (b >= b) = " << (b >= b) << std::endl;

    Region<1,double> a2(3,5);
    Region<1,double> a3(3.5,4);
    Range<1> r2(2,10,2);

    tester.out() << "  touches([3:5], [2:10:2]) = " << touches(a2,r2)
		 << std::endl;
    tester.out() << "  contains([3:5], [2:10:2]) = " << contains(a2,r2)
		 << std::endl;
    tester.out() << "  contains([3:5], [3.5,4]) = " << contains(a2,a3)
		 << std::endl;
    tester.out() << "  intersect([3:5], [2:10:2]) = " << intersect(a2,r2)
		 << std::endl;
    tester.out() << "  intersect([3:5], [3.5,4]) = " << intersect(a2,a3)
		 << std::endl;

    Region<1,double> a4, a5;
    split(a3,a4,a5);
    tester.out() << "  split([3.5,4]) ==> " << a4 << ", " << a5 << std::endl;
  }

  tester.out() << std::endl;
  tester.out() << "Testing Region<N,float>:" << std::endl;
  tester.out() << "-------------------------" << std::endl;
  {
    Region<3,float> a;
  
    Region<1,float> a1(3);
    tester.out() << "  Region<1>(3) a1 = " << a1 << std::endl;

    Region<2,float> b(2.0, Region<1,float>(1.0, 1.5));
    tester.out() << "  Region<2>(2.0, Region<1>(1.0, 1.5)) b = " << b
		 << std::endl;

    Interval<1> b1(1,5);
    Range<1>    b2(2,8,2);
    tester.out() << "  combine(b, 2:8:2) = ";
    tester.out() << NewDomain2<Region<2,float>,Range<1> >::combine(b,b2)
		 << std::endl;

    NewDomain2<Interval<1>,Range<1> >::fill(b, b1, b2);
    tester.out() << "  fill(1:5, 2:8:2) = " << b << std::endl;

    b *= 2;
    tester.out() << "  b *= 2 ==> b = " << b << std::endl;

    b += Loc<2>(3,4);

    tester.out() << "  b += Loc<2>(3,4) ==> b = " << b << std::endl;

    tester.out() << "  result of (b == b) = " << (b == b) << std::endl;
    tester.out() << "  result of (b != b) = " << (b != b) << std::endl;
    tester.out() << "  result of (b  < b) = " << (b  < b) << std::endl;
    tester.out() << "  result of (b >= b) = " << (b >= b) << std::endl;

    Region<1,float> a2(3,5);
    Region<1,float> a3(3.5,4);
    Region<1,double> a3d(3.5,7);
    Range<1> r2(2,10,2);

    tester.out() << "  touches([3:5], [2:10:2]) = " << touches(a2,r2)
		 << std::endl;
    tester.out() << "  contains([3:5], [2:10:2]) = " << contains(a2,r2)
		 << std::endl;
    tester.out() << "  contains([3:5], [3.5,4]) = " << contains(a2,a3)
		 << std::endl;
    tester.out() << "  intersect([3:5], [2:10:2]) = " << intersect(a2,r2)
		 << std::endl;
    tester.out() << "  intersect([3:5], [3.5,7]) = " << intersect(a2,a3d)
		 << std::endl;

    Region<1,float> a4, a5;
    split(a3,a4,a5);
    tester.out() << "  split([3.5,4]) ==> " << a4 << ", " << a5 << std::endl;
  }
 
  int ret = tester.results("domaintest");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: domaintest.cpp,v $   $Author: richard $
// $Revision: 1.19 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
