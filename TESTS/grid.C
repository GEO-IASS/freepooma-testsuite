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


int main(int argc, char *argv[]) {
  int i, j;

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Test constructing and iterating over Grid objects

  tester.out() << "Grid domain tests:" << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  // Construct an IndirectionList and a Range

  IndirectionList<int> list(4);
  list(0) = 2;
  list(1) = 5;
  list(2) = 6;
  list(3) = 9;

  tester.out() << "Created IndirectionList =";
  for (i=0; i < list.length(); ++i)
    tester.out() << "  " << list(i);
  tester.out() << std::endl;
  tester.check(list.length() == 4);

  Range<1> range(8, 4, -2);

  tester.out() << "Created Range = " << range << std::endl;

  // Construct a 1D Grid from the IndirectionList and from the range.

  Grid<1> g1(list);
  Grid<1> g2(range);

  tester.out() << "Created Grid<1> from list = " << g1 << std::endl;
  Grid<1>::iterator g1i = g1.begin();
  for (i=0 ; g1i != g1.end(); ++g1i, ++i)
    {
      tester.check(*g1i == list(i));
    }

  tester.out() << "Created Grid<1> from range = " << g2 << std::endl;
  Grid<1>::iterator  g2i = g2.begin();
  Range<1>::iterator ri = range.begin();
  for ( ; g2i != g2.end(); ++g2i, ++ri)
    {
      tester.check(*g2i == *ri);
    }

  // Construct a 2D Grid from the IndirectionList and the range, and
  // from a list and a Grid<1>

  Grid<2> g3(range, list);
  Grid<2> g4(list, g1);

  tester.out() << "Created Grid<2> from range and list:" << std::endl;
  Grid<2>::iterator g3i = g3.begin();
  for (j = 0; j < g3[1].length(); ++j)
    {
      Range<1>::iterator ri = range.begin();
      for (i = 0; i < g3[0].length(); ++g3i, ++i, ++ri)
	{
	  tester.out() << "  " << *g3i;
	  tester.check(*g3i == Loc<2>(*ri, list(j)));
	}
      tester.out() << std::endl;
    }
  tester.check(g3i == g3.end());

  tester.out() << "Created Grid<2> from list and Grid<1>:" << std::endl;
  g3i = g4.begin();
  for (j = 0; j < g4[1].length(); ++j)
    {
      for (i = 0; i < g4[0].length(); ++g3i, ++i)
	{
	  tester.out() << "  " << *g3i;
	  tester.check(*g3i == Loc<2>(list(i), list(j)));
	}
      tester.out() << std::endl;
    }
  tester.check(g3i == g4.end());

  // Test +=, -=, etc. operations

  tester.out() << "\nArithmetic operations:" << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  Grid<1> a1;

  tester.out() << "Testing g1 += 4:" << std::endl;
  a1 = g1;
  a1 += 4;
  tester.out() << "  original = " << g1 << std::endl;
  tester.out() << "       new = " << a1 << std::endl;
  for (i=0; i < a1.length(); ++i)
    tester.check(a1(i) == (g1(i) + 4));

  tester.out() << "Testing g1 -= 2:" << std::endl;
  a1 = g1;
  a1 -= 2;
  tester.out() << "  original = " << g1 << std::endl;
  tester.out() << "       new = " << a1 << std::endl;
  for (i=0; i < a1.length(); ++i)
    tester.check(a1(i) == (g1(i) - 2));

  tester.out() << "Testing g1 *= -3:" << std::endl;
  a1 = g1;
  a1 *= -3;
  tester.out() << "  original = " << g1 << std::endl;
  tester.out() << "       new = " << a1 << std::endl;
  for (i=0; i < a1.length(); ++i)
    tester.check(a1(i) == (g1(i) * -3));

  tester.out() << "Testing g1 /= 2:" << std::endl;
  a1 = g1;
  a1 /= 2;
  tester.out() << "  original = " << g1 << std::endl;
  tester.out() << "       new = " << a1 << std::endl;
  for (i=0; i < a1.length(); ++i)
    tester.check(a1(i) == (g1(i) / 2));

  tester.out() << "Testing g2 += (5,10):" << std::endl;
  Grid<2> g5(g1, g1);
  Grid<2> g6(g5);
  Loc<2> val(5,10);
  g6 += val;
  tester.out() << "  original = " << g5 << std::endl;
  tester.out() << "       new = " << g6 << std::endl;
  for (j=0; j < 2; ++j)
    for (i=0; i < g5[j].length(); ++i)
      tester.check(g6[j](i) == (g5[j](i) + val[j](i)));


  tester.out() << "\nBlock iterator:" << std::endl;
  tester.out() << "-------------------------------------------" << std::endl;

  tester.out() << "Blocks in Grid<2> = " << g4 << ":" << std::endl;
  Grid<2>::blockIterator bi = g4.beginBlock();
  for (j = 0; j < (g4[1].length() - 1); ++j)
    {
      for (i = 0; i < (g4[0].length() - 1); ++bi, ++i)
	{
	  tester.out() << "  " << *bi;
	  tester.out() << " (point = " << bi.point() << ", index = ";
	  tester.out() << bi.index() << ")" << std::endl;
	  tester.check(*bi==Interval<2>(Interval<1>(g4[0](i), g4[0](i+1)-1),
					Interval<1>(g4[1](j), g4[1](j+1)-1)));
	}
    }
  tester.check(bi == g4.endBlock());

  tester.out() << "\nBlocks in Grid<1> = " << g1 << ":" << std::endl;
  Grid<1>::blockIterator bbi = g1.beginBlock();
  for (i=0; i < (g1.length() - 1); ++i, ++bbi)
    {
      tester.out() << "  " << *bbi;
      tester.out() << " (point = " << bbi.point() << ", index = ";
      tester.out() << bbi.index() << ")" << std::endl;
      tester.check(*bbi==Interval<1>(g1(i), g1(i+1) - 1));
    }
  tester.check(bbi == g1.endBlock());

  Range<1> decrange(9,3,-3);
  tester.out() << "\nBlocks in Range<1> = " << decrange << ":" << std::endl;
  Range<1>::blockIterator dri = decrange.beginBlock();
  for (i=0; i < (decrange.length() - 1); ++i, ++dri)
    {
      tester.out() << "  " << *dri;
      tester.out() << " (point = " << dri.point() << ", index = ";
      tester.out() << dri.index() << ")" << std::endl;
      tester.check(*dri==Interval<1>(decrange(i+1) + 1, decrange(i)));
    }
  tester.check(dri == decrange.endBlock());
  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Grid domain tests.");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: grid.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
