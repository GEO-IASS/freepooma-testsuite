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
// Test of PatchSizeSyncer
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Tulip/Messaging.h"
#include "Tulip/CollectFromContexts.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const int numContexts = Pooma::contexts();
  const int myContext   = Pooma::context();

  tester.out() << "Running with " << numContexts << " contexts." << std::endl;

  CollectFromContexts<int> ranks(2*(Pooma::context()+1));
  if (Pooma::context() == 0) {
    bool check = true;
    for (int i=0; i<Pooma::contexts(); ++i)
      if (ranks[i] != 2*(i+1)) {
	tester.out() << "[" << i << "] should be "
		     << 2*(i+1) << ", but is " << ranks[i] << "\n";
        check = false;
      }
    tester.check("Collecting ranks", check);
  }

  // We can't do the following test on !MESSAGING, as invalid data on
  // context 0 is not supported in this case.
#if POOMA_MESSAGING
  CollectFromContexts<int> ranks2(Pooma::context()+1, 0,
			          Pooma::context() > 0
				  && Pooma::context() < Pooma::contexts()-1);
  if (Pooma::context() == 0) {
    bool check = true;
    for (int i=1; i<Pooma::contexts()-1; ++i)
      if (ranks2[i] != i+1) {
	tester.out() << "[" << i << "] should be "
		     << (i+1) << ", but is " << ranks[i] << "\n";
        check = false;
      }
    tester.check("Collecting ranks, but not first and last", check);
  }
#endif
 
  int ret = tester.results("CollectFromContextsTest");
  Pooma::finalize();

  return ret;
}

