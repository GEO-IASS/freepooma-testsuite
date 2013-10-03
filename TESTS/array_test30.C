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
// array_test30: verify correctness of igc updates
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Arrays.h"
#include "Utilities/Tester.h"
#include <iostream>


template <class A1, class A2>
bool test(Pooma::Tester& tester,
	  const A1& a_mp, const A1& b_mp,
	  const A2& a_sp, const A2& b_sp,
	  const Loc<2>& delta1, const Loc<2>& delta2,
	  bool initial_f, const Loc<2>& initial)
{
  static int sequence = 0;
  Interval<2> I;

  // initialize rhs arrays, ensure wrong igc values
  // via sequence number.
  I = b_sp.totalDomain();
  b_sp(I) = sequence + iota(I).comp(0) + I[0].size()*iota(I).comp(1);
  b_mp.engine().setGuards(0);
  b_mp(I) = b_sp(I);

  // if requested, force initial update of a set of igcs
  if (initial_f) {
    b_sp(b_sp.physicalDomain()) = b_mp(b_sp.physicalDomain()+initial);
    b_sp(I) = sequence + iota(I).comp(0) + I[0].size()*iota(I).comp(1);
    Pooma::blockAndEvaluate();
  }

  // do calculation both sp and mp
  I = a_sp.physicalDomain();
  a_sp(I) = b_sp(I+delta1) - b_sp(I+delta2);
  a_mp(I) = b_mp(I+delta1) - b_mp(I+delta2);

  // check the results are the same everywhere
  bool res = all(a_sp(I) == a_mp(I));
  tester.out() << "For deltas " << delta1 << " and " << delta2 << " ";
  if (initial_f)
    tester.out() << "with initial " << initial << " ";
  tester.check("result is", res);
  if (!res) {
    int n = b_mp.layout().sizeGlobal();
    for (int i=0; i<n; ++i) {
      Array<2, int, Remote<Brick> > b(b_mp.engine().globalPatch(i));
      tester.out() << "Brick " << i << " " << intersect(b.domain(), b_mp.physicalDomain())
	      	   << " on context " << b.engine().owningContext()
		   << " is\n" << b(intersect(b.totalDomain(), b_mp.physicalDomain()))
		   << std::endl;
    }
    tester.out() << "Aborting." << std::endl;
    return false;
  }

  sequence++;

  return true;
}


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<2> domain(12, 12);
  UniformGridLayout<2> layout_mp(domain, Loc<2>(3, 3),
				 GuardLayers<2>(2), DistributedTag());
  DomainLayout<2> layout_sp(domain, GuardLayers<2>(2));

  Array<2, int, MultiPatch<UniformTag, Remote<Brick> > >
    a_mp(layout_mp), b_mp(layout_mp);
  Array<2, int, Brick>
    a_sp(layout_sp), b_sp(layout_sp);

  // all 5^4 == 625 uninitialized cases
  for (int d1i = -2; d1i <= 2; ++d1i)
    for (int d1j = -2; d1j <= 2; ++d1j)
      for (int d2i = -2; d2i <= 2; ++d2i)
	for (int d2j = -2; d2j <= 2; ++d2j)
	  if (!test(tester, a_mp, b_mp, a_sp, b_sp,
		    Loc<2>(d1i, d1j), Loc<2>(d2i, d2j),
		    false, Loc<2>(0)))
	    goto out;

  // all 5^4 == 625 initialized cases with simplified expression
  for (int ii = -2; ii <= 2; ++ii)
    for (int ij = -2; ij <= 2; ++ij)
      for (int d1i = -2; d1i <= 2; ++d1i)
        for (int d1j = -2; d1j <= 2; ++d1j)
	  if (!test(tester, a_mp, b_mp, a_sp, b_sp,
		    Loc<2>(d1i, d1j), Loc<2>(d1i, d1j),
		    true, Loc<2>(ii, ij)))
	    goto out;

 out:
  tester.out() << "Best testing is done with all 1 to 9 processes" << std::endl;

  int retval = tester.results("array_test30");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test30.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
