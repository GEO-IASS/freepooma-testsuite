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
// Paws test 8: Send and Receive an int and double set of scalars,
//              plus a dynamic array, in conjunction with test 7.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"

#if POOMA_PAWS
#include "Pooma/Paws.h"
#endif // POOMA_PAWS

#include "Pooma/DynamicArrays.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Paws DynmaicArray send/receive test B\n";
  tester.out() << "----------------------------------------------------";
  tester.out() << std::endl;

#if POOMA_PAWS

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  int iters = 10;

  // DynamicArrays for reference ...

  Interval<1> refdomain(100);
  Loc<1> refblocks(2);
  GridPartition<1> refgpar(refblocks);
  LocalMapper<1> refcmap(refgpar);
  DynamicLayout reflayout(refdomain, refgpar, refcmap);
  DynamicArray<float, MultiPatch<DynamicTag,Dynamic> > refa1(reflayout);
  DynamicArray<int, MultiPatch<DynamicTag,Dynamic> > refa2(reflayout);
  DynamicArray<double, Dynamic> refa3(30);

  // DynamicArrays to receive ...

  Interval<1> domain(3);
  Loc<1> blocks(3);
  GridPartition<1> gpar(blocks);
  LocalMapper<1> cmap(gpar);
  DynamicLayout layout(domain, gpar, cmap);
  DynamicLayout layout2(domain, gpar, cmap);
  DynamicArray<float, MultiPatch<DynamicTag,Dynamic> >  a1(layout);
  DynamicArray<int, MultiPatch<DynamicTag,Dynamic> >    a2(layout);
  DynamicArray<double, MultiPatch<DynamicTag,Dynamic> > a3(layout2);

  // Initialize the arrays

  refa1 = 1 + iota(refa1.domain()).comp(0);
  refa2 = 1000 + refa1;
  refa3 = 4.5;
  a1 = 0;
  a2 = 0;
  a3 = 0;
  Pooma::blockAndEvaluate();

  // Create a Paws connection

  tester.out() << "Creating PawsConnection object ..." << std::endl;
  Connection<Paws> *paws = new Connection<Paws>("test8", argc, argv);
  tester.out() << "Finished creating PawsConnection object." << std::endl;

  // Establish connections for the two scalars

  tester.out() << "Connecting s1 = " << s1 << " for input ..." << std::endl;
  ConnectorBase *s1p = paws->connectScalar("s1", s1, ConnectionBase::in);
  tester.out() << "Connecting s2 = " << s2 << " for output ..." << std::endl;
  ConnectorBase *s2p = paws->connectScalar("s2", s2, ConnectionBase::out);
  tester.out() << "Connecting iters = " << iters << " for input ...";
  tester.out() << std::endl;
  ConnectorBase *iterp = paws->connectScalar("iters",iters,ConnectionBase::in);

  // Establish connections for the arrays

  tester.out() << "Connecting a1 = " << a1 << " for input ..." << std::endl;
  paws->connect("a1", a1, ConnectionBase::in);
  tester.out() << "Connecting a2 = " << a2 << " for input ..." << std::endl;
  paws->connect("a2", a2, ConnectionBase::in);
  tester.out() << "Connecting a3 = " << a3 << " for input ..." << std::endl;
  paws->connect("a3", a3, ConnectionBase::in);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  paws->ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Modify s1, and update

  s1 *= 2;
  tester.out() << "Updating current s1 = " << s1 << " and s2 = " << s2;
  tester.out() << ", plus arrays ..." << std::endl;
  paws->update();

  // Report the results

  tester.out() << "Received update.  New values:" << std::endl;
  tester.out() << "  s1 = " << s1 << " (should be " << origs1 << ")\n";
  tester.out() << "  s2 = " << s2 << " (should be " << origs2 << ")\n";
  tester.out() << std::endl;
  tester.check("s1 OK", s1 == origs1);
  tester.check("s2 OK", s2 == origs2);

  // Disconnect the scalars

  int connections = paws->size();
  tester.out() << "Disconnecting scalars ..." << std::endl;
  delete s1p;
  delete s2p;
  delete iterp;
  tester.check("3 less connections", paws->size() == (connections - 3));

  // Do, in a loop, updates of the receiver. Add one to the arrays each time,
  // plus delete the second element.

  int myiters = iters;
  while (myiters-- > 0)
    {
      refa1 += 1;
      refa2 += 1;
      refa3 += 1;
      Pooma::blockAndEvaluate();
      refa1.destroy(Interval<1>(1,1), ShiftUp());
      refa1.sync();

      tester.out() << "Receiving for iters = " << myiters << std::endl;
      paws->update();

      // Compare to reference

      tester.check("a1 size", a1.domain().size() == refa1.domain().size());
      tester.check("a2 size", a2.domain().size() == refa2.domain().size());
      tester.check("a3 size", a3.domain().size() == refa3.domain().size());
      int a1msd = sum((a1 - refa1)*(a1 - refa1));
      int a2msd = sum((a2 - refa2)*(a2 - refa2));
      int a3msd = sum((a3 - refa3)*(a3 - refa3));
      tester.check("a1 MSD", a1msd == 0);
      tester.check("a2 MSD", a2msd == 0);
      tester.check("a3 MSD", a3msd == 0);
    }

  // Delete PAWS connection, disconnecting us from the other code.

  tester.out() << "Deleting Connection<Paws> object ..." << std::endl;
  delete paws;

#else // POOMA_PAWS

  tester.out() << "Please configure with --paws to use this test code!"
	       << std::endl;

#endif // POOMA_PAWS

  // Finish up and report results

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Paws DynamicArray send/receive test A");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: paws_test8.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:24 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
