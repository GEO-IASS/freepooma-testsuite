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
// Paws test 7: Send and Receive an int and double set of scalars,
//              plus a dynamic array, in conjunction with test 8.
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
  tester.out() << argv[0] << ": Paws DynmaicArray send/receive test A\n";
  tester.out() << "----------------------------------------------------";
  tester.out() << std::endl;

#if POOMA_PAWS

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  int iters = 10;

  // DynamicArrays to send and receive ...

  Interval<1> domain(100);
  Loc<1> blocks(4);
  GridPartition<1> gpar(blocks);
  LocalMapper<1> cmap(gpar);
  DynamicLayout layout(domain, gpar, cmap);
  DynamicArray<float, MultiPatch<DynamicTag,Dynamic> > a1(layout);
  DynamicArray<int, MultiPatch<DynamicTag,Dynamic> > a2(layout);
  DynamicArray<double, Dynamic> a3(30);

  // Initialize the arrays

  a1 = 1 + iota(a1.domain()).comp(0);
  a2 = 1000 + a1;
  a3 = 4.5;
  Pooma::blockAndEvaluate();

  // Create a Paws connection

  tester.out() << "Creating PawsConnection object ..." << std::endl;
  Connection<Paws> *paws = new Connection<Paws>("test7", argc, argv);
  tester.out() << "Finished creating PawsConnection object." << std::endl;

  // Establish connections for the two scalars

  tester.out() << "Connecting s1 = " << s1 << " for output ..." << std::endl;
  ConnectorBase *s1p = paws->connectScalar("s1", s1, ConnectionBase::out);
  tester.out() << "Connecting s2 = " << s2 << " for input ..." << std::endl;
  ConnectorBase *s2p = paws->connectScalar("s2", s2, ConnectionBase::in);
  tester.out() << "Connecting iters = " << iters << " for output ...";
  tester.out() << std::endl;
  ConnectorBase *iterp=paws->connectScalar("iters",iters,ConnectionBase::out);

  // Establish connections for the arrays

  tester.out() << "Connecting a1 = " << a1 << " for output ..." << std::endl;
  paws->connect("a1", a1, ConnectionBase::out);
  tester.out() << "Connecting a2 = " << a2 << " for output ..." << std::endl;
  paws->connect("a2", a2, ConnectionBase::out);
  tester.out() << "Connecting a3 = " << a3 << " for output ..." << std::endl;
  paws->connect("a3", a3, ConnectionBase::out);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  paws->ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Modify s2, and update

  s2 *= 2;
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
      a1 += 1;
      a2 += 1;
      a3 += 1;
      Pooma::blockAndEvaluate();
      a1.destroy(Interval<1>(1,1), ShiftUp());
      a1.sync();

      tester.out() << "Sending for iters = " << myiters << std::endl;
      paws->update();
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
// $RCSfile: paws_test7.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:24 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
