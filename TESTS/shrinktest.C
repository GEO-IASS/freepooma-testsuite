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
// DynamicArray test 1: Create/destroy operations, using SharedBrick engines.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Domain/Shrink.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Domain shrink and grow functions." << std::endl;
  tester.out() << "----------------------------------------------" << std::endl;

  int i;
  Interval<2> a1(10,10),a2(Interval<1>(5,20),Interval<1>(5,20));
  Interval<2> b,c;
  b = shrinkRight(a1,Loc<2>(1,1));
  c = shrinkRight(a2,Loc<2>(1,1));
  tester.out() << b << "," << c << std::endl;
  tester.check(b == Interval<2>(9,9));
  tester.check(c == Interval<2>(Interval<1>(5,19),Interval<1>(5,19)));

  Interval<2> d,e;
  d = growRight(b,Loc<2>(2,2));
  e = growRight(c,Loc<2>(2,2));
  tester.out() << d << "," << e << std::endl;
  tester.check(d == Interval<2>(11,11));
  tester.check(e == Interval<2>(Interval<1>(5,21),Interval<1>(5,21)));


  Interval<2> I(Interval<1>(2, 4), Interval<1>(1, 5));

  tester.check("shrinkRight(D, Loc<>)",
	       shrinkRight(I, Loc<2>(1, 0)),
	       Interval<2>(Interval<1>(2, 3), Interval<1>(1, 5)));
  tester.check("shrinkRight(D, int)",
	       shrinkRight(I, 2),
	       Interval<2>(Interval<1>(2, 2), Interval<1>(1, 3)));
  tester.check("growRight(D, Loc<>)",
	       growRight(I, Loc<2>(1, 0)),
	       Interval<2>(Interval<1>(2, 5), Interval<1>(1, 5)));
  tester.check("growRight(D, int)",
	       growRight(I, 2),
	       Interval<2>(Interval<1>(2, 6), Interval<1>(1, 7)));
  tester.check("shrinkLeft(D, Loc<>)",
	       shrinkLeft(I, Loc<2>(1, 0)),
	       Interval<2>(Interval<1>(3, 4), Interval<1>(1, 5)));
  tester.check("shrinkLeft(D, int)",
	       shrinkLeft(I, 2),
	       Interval<2>(Interval<1>(4, 4), Interval<1>(3, 5)));
  tester.check("growLeft(D, Loc<>)",
	       growLeft(I, Loc<2>(1, 0)),
	       Interval<2>(Interval<1>(1, 4), Interval<1>(1, 5)));
  tester.check("growLeft(D, int)",
	       growLeft(I, 2),
	       Interval<2>(Interval<1>(0, 4), Interval<1>(-1, 5)));
  tester.check("grow(D, Loc<>)",
	       grow(I, Loc<2>(1, 2)),
	       Interval<2>(Interval<1>(1, 5), Interval<1>(-1, 7)));
  tester.check("grow(D, int)",
	       grow(I, 1),
	       Interval<2>(Interval<1>(1, 5), Interval<1>(0, 6)));
  tester.check("shrink(D, Loc<>)",
	       shrink(I, Loc<2>(1, 2)),
	       Interval<2>(Interval<1>(3, 3), Interval<1>(3, 3)));
  tester.check("shrink(D, int)",
	       shrink(I, 1),
	       Interval<2>(Interval<1>(3, 3), Interval<1>(2, 4)));

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "----------------------------------------------" << std::endl;
  int retval = tester.results("Domain shrink");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: shrinktest.cpp,v $   $Author: richi $
// $Revision: 1.6 $   $Date: 2004/11/29 12:23:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
