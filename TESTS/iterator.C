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
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"
#include "Array/Array.h"
#include "Engine/BrickEngine.h"
#include "Domain/IndirectionList.h"


int main(int argc, char *argv[]) {
  int i, j;

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  {
    Grid<5> foo;
  }

  { 
    Grid<5> goo;
    goo[1]=Grid<1>(8);
  }


  // Test iterating over 1D domains.

  tester.out() << "Single-dimension domain iterators:" << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  Loc<1> a(3);
  Interval<1> b(1, 4);
  Range<1> c(4, 8, 2);

  tester.out() << "Iterating over Loc<1> = " << a << ":" << std::endl;

  Loc<1>::iterator ai = a.begin();
  for (i = 3; ai != a.end(); ++ai, ++i)
    {
      tester.out() << "  " << *ai;
      tester.check(*ai == i);
    }
  tester.out() << std::endl;

  tester.out() << "Iterating over Interval<1> = " << b << ":" << std::endl;

  Interval<1>::iterator bi = b.begin();
  for (i = 1; bi != b.end(); ++bi, ++i)
    {
      tester.out() << "  " << *bi;
      tester.check(*bi == i);
    }
  tester.out() << std::endl;

  tester.out() << "Iterating over Range<1> = " << c << ":" << std::endl;

  Range<1>::iterator ci = c.begin();
  for (i = 4; ci != c.end(); ++ci, i += 2)
    {
      tester.out() << "  " << *ci;
      tester.check(*ci == i);
    }
  tester.out() << std::endl;

  tester.out() << "Testing operator-> on domain iterator:" << std::endl;
  ci = c.begin();
  tester.out() << "  ci->first() == " << ci->first() << " (should be ";
  tester.out() << (*ci).first() << ")" << std::endl;
  tester.check(ci->first() == (*ci).first());

  // Test iterating over 2D domains.

  tester.out() << "\nTwo-dimensional domain iterators:" << std::endl;
  tester.out() << "---------------------------------" << std::endl;

  Interval<2> b2(b, b);
  Range<2> c2(b, c);

  tester.out() << "Iterating over Interval<2> = " << b2 << ":" << std::endl;

  Interval<2>::iterator b2i = b2.begin();
  for (j = 1; b2i != b2.end(); ++j)
    {
      for (i = 1; i <= 4 && b2i != b2.end(); ++b2i, ++i)
	{
	  tester.out() << "  " << *b2i;
	  tester.check(*b2i == Loc<2>(i, j));
	}
      tester.out() << std::endl;
    }
  tester.check(b2i == b2.end());

  tester.out() << "Iterating over Range<2> = " << c2 << ":" << std::endl;

  Range<2>::iterator c2i = c2.begin();
  for (j = 4; c2i != c2.end(); j += 2)
    {
      for (i = 1; i <= 4 && c2i != c2.end(); ++c2i, ++i)
	{
	  tester.out() << "  " << *c2i;
	  tester.check(*c2i == Loc<2>(i, j));
	}
      tester.out() << std::endl;
    }
  tester.check(c2i == c2.end());
  
  tester.out()<< " Testing blockIterator on Grid<2>(Range<2>) " <<std::endl;
  Grid<2> d2(c2);

  Grid<2>::blockIterator d2bi(d2);
  Grid<2>::blockIterator d2biend;
  while(d2bi!=d2biend)
    {
      tester.out()<<" "<<*d2bi;
      tester.out()<<" "<<d2bi.index();
      tester.out()<<" "<<d2bi.point();
      tester.out()<<std::endl;
      ++d2bi;
    }



  Array<1,int> ar1(Interval<1>(0,5));
  ar1(0)=0;ar1(1)=3;ar1(2)=4;ar1(3)=7;ar1(4)=8;ar1(5)=10;
  Array<1,int> ar2(Interval<1>(0,5));
  ar2(0)=0;ar2(1)=1;ar2(2)=2;ar2(3)=6;ar2(4)=8;ar2(5)=10;

  IndirectionList<int> il1(ar1),il2(ar2);

  il1+=0;
  il2+=0;

  tester.out()<<" indirections lists used to make Grid<2> " <<std::endl;
  for(int i=0;i<il1.size();i++)
    tester.out() << il1(i);
  tester.out() << std::endl;
  for(int i=0;i<il1.size();i++)
    tester.out() << il2(i);
  tester.out() << std::endl;

  Grid<2> g2(il1,il2);
  Grid<2>::blockIterator g2i(g2);
  Grid<2>::blockIterator g2end;


 tester.out()<< " Testing blockIterator on Grid<2> " <<std::endl;
 
  while(g2i!=g2end)
    {
      tester.out()<< " " << *g2i;
      tester.out()<< " " << g2i.index();
      tester.out()<< " " << g2i.point();
      tester.out()<<std::endl;
      ++g2i;
    }



  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Domain Iterators");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: iterator.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
