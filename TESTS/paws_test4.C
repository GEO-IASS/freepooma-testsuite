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
// Paws test 4: Send and Receive an int and double set of scalars,
//              plus a fixed-size 3D Array, in conjunction with test 3.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"

#if POOMA_PAWS
#include "Pooma/Paws.h"
#endif // POOMA_PAWS

#include "Pooma/Arrays.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Paws Array send/receive test B" << std::endl;
  tester.out() << "--------------------------------------------" << std::endl;

#if POOMA_PAWS

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  int iters = 0, citers = 10;

  // Arrays to send and receive ... use different layouts in the two
  // test codes.

  Loc<3> blocks(2,2,1);
  Interval<3> domain(2,4,8);
  Interval<3> subdomain(1, 2, 2);
  GridLayout<3> layout(domain, blocks, ReplicatedTag());
  Array<3, float, Brick> a1(domain);
  Array<3, int, MultiPatch<GridTag,Brick> > a2(layout);
  Array<3, float, Brick> a3(subdomain);
  Array<3, float, Brick> ca1(domain);
  Array<3, int, Brick> ca2(domain);
  Array<3, float, Brick> ca3(subdomain);

  // Initialize the arrays to be all zero, and the compare arrays to be
  // what we expect from the sender.

  a1 = 0;
  a2 = 0;
  a3 = 0;
  ca1 = 100 * (iota(domain).comp(2) + 1) + 10 * (iota(domain).comp(1) + 1) +
    iota(domain).comp(0) + 1;
  ca2 = ca1 + 1000;
  ca3 = ca1(subdomain);

  // Create a Paws connection

  tester.out() << "Creating PawsConnection object ..." << std::endl;
  Connection<Paws> *paws = new Connection<Paws>("test4", argc, argv);
  tester.out() << "Finished creating PawsConnection object." << std::endl;

  // Establish connections for the two scalars

  tester.out() << "Connecting s1 = " << s1 << " for input ..." << std::endl;
  ConnectorBase *s1p = paws->connectScalar("s1", s1, ConnectionBase::in);
  tester.out() << "Connecting s2 = " << s2 << " for output ..." << std::endl;
  ConnectorBase *s2p = paws->connectScalar("s2", s2, ConnectionBase::out);
  tester.out() << "Connecting iters = " << iters << " for input ...";
  tester.out() << std::endl;
  ConnectorBase *ip = paws->connectScalar("iters", iters, ConnectionBase::in);

  // Establish connections for the two arrays; also connect up a view of
  // the first array

  tester.out() << "Connecting a1 = " << a1 << " for input ..." << std::endl;
  paws->connect("a1", a1, ConnectionBase::in);
  tester.out() << "Connecting a2 = " << a2 << " for input ..." << std::endl;
  paws->connect("a2", a2, ConnectionBase::in);
  tester.out() << "Connecting a3 = " << a3 << " for input ..." << std::endl;
  paws->connect("a1view", a3, ConnectionBase::in);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  paws->ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Modify s1, and update

  s1 *= 2;
  tester.out() << "Updating current s1 = " << s1 << " and s2 = " << s2;
  tester.out() << " ..." << std::endl;
  paws->update();

  // Report the results

  tester.out() << "Received update.  New values:" << std::endl;
  tester.out() << "  s1 = " << s1 << " (should be " << origs1 << ")\n";
  tester.out() << "  s2 = " << s2 << " (should be " << origs2 << ")\n";
  tester.out() << "  iters = " << iters << " (should be " << citers << ")\n";
  tester.out() << std::endl;
  tester.check("s1 OK", s1 == origs1);
  tester.check("s2 OK", s2 == origs2);
  tester.check("iters OK", iters == citers);

  // Also report Array results

  tester.out() << "Received Arrays as well.  New values:" << std::endl;
  tester.out() << "  a1 = " << a1 << std::endl;
  tester.out() << "  a2 = " << a2 << std::endl;
  tester.out() << "  a3 = " << a3 << std::endl;
  tester.check("a1 OK", all(a1 == ca1));
  tester.check("a2 OK", all(a2 == ca2));
  tester.check("a3 OK", all(a3 == ca3));

  // Disconnect the scalars

  int connections = paws->size();
  tester.out() << "Disconnecting scalars ..." << std::endl;
  delete s1p;
  delete s2p;
  delete ip;
  tester.check("3 less connections", paws->size() == (connections - 3));

  // Do, in a loop, updates of the receiver. Add one to the arrays each time.

  int myiters = iters;
  while (myiters-- > 0)
    {
      ca1 += 1;
      ca2 += 1;
      ca3 += 1;
      tester.out() << "Receiving for iters = " << myiters << std::endl;
      paws->update();
      tester.out() << "Receive complete." << std::endl;
      tester.check("a1 OK", all(a1 == ca1));
      tester.check("a2 OK", all(a2 == ca2));
      tester.check("a3 OK", all(a3 == ca3));
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
  int retval = tester.results("Paws Array send/receive test B");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: paws_test4.cpp,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:23 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
