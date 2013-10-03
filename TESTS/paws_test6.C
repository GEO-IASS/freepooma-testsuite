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
// Paws test 6: Send and Receive an int and double set of scalars,
//              plus a fixed-size 2D Field, in conjunction with test 5.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"

#if POOMA_PAWS
#include "Pooma/Paws.h"
#endif // POOMA_PAWS

#include "Pooma/Fields.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Paws Field send/receive test B" << std::endl;
  tester.out() << "--------------------------------------------" << std::endl;

#if POOMA_PAWS

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  int iters = 0, citers = 10;

  // Fields to send and receive ... use different layouts in the two
  // test codes.

  Loc<2> blocks(2,1);
  Interval<2> domain(6,2);
  Interval<2> subdomain(3, 2);
  Vector<2,double> origin(2.0);
  Vector<2,double> spacings(0.2);
  RectilinearMesh<2> mesh(domain, origin, spacings);

  // Create the geometry and layout.

  typedef DiscreteGeometry<Vert, RectilinearMesh<2> > Geometry_t;
  Geometry_t geom(mesh);
  GridLayout<2> layout(domain, blocks, ReplicatedTag());

  // Now create the Fields

  Field<Geometry_t, float, Brick> a1(geom);
  Field<Geometry_t, int, MultiPatch<GridTag,Brick> > a2(geom, layout);
  Array<2, float, Brick> a3(subdomain);
  Array<2, float, Brick> ca1(domain);
  Array<2, int, Brick> ca2(domain);
  Array<2, float, Brick> ca3(subdomain);

  // Initialize the arrays to be all zero, and the compare arrays to be
  // what we expect from the sender.

  a1 = 0;
  a2 = 0;
  a3 = 0;
  ca1 = 10 * (iota(domain).comp(1) + 1) + iota(domain).comp(0) + 1;
  ca2 = ca1 + 1000;
  ca3 = ca1(subdomain);

  // Create a Paws connection

  tester.out() << "Creating PawsConnection object ..." << std::endl;
  Connection<Paws> *paws = new Connection<Paws>("test6", argc, argv);
  tester.out() << "Finished creating PawsConnection object." << std::endl;

  // Establish connections for the two scalars

  tester.out() << "Connecting s1 = " << s1 << " for input ..." << std::endl;
  ConnectorBase *s1p = paws->connectScalar("s1", s1, ConnectionBase::in);
  tester.out() << "Connecting s2 = " << s2 << " for output ..." << std::endl;
  ConnectorBase *s2p = paws->connectScalar("s2", s2, ConnectionBase::out);
  tester.out() << "Connecting iters = " << iters << " for input ...";
  tester.out() << std::endl;
  ConnectorBase *iterp=paws->connectScalar("iters", iters, ConnectionBase::in);

  // Establish connections for the two fields; also connect up a view of
  // the first fields

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

  // Also report Field results

  tester.out() << "Received Fields as well.  New values:" << std::endl;
  tester.out() << "  a1 = " << a1 << std::endl;
  tester.out() << "  a2 = " << a2 << std::endl;
  tester.out() << "  a3 = " << a3 << std::endl;
  tester.check("a1 OK", all(a1.array() == ca1));
  tester.check("a2 OK", all(a2.array() == ca2));
  tester.check("a3 OK", all(a3 == ca3));

  // Disconnect the scalars

  int connections = paws->size();
  tester.out() << "Disconnecting scalars ..." << std::endl;
  paws->disconnect(s1p);
  paws->disconnect(s2p);
  paws->disconnect(iterp);
  tester.check("3 less connections", paws->size() == (connections - 3));

  // Do, in a loop, updates of the receiver. Add one to the arrays each time.

  int runiters =  iters;
  while (runiters-- > 0)
    {
      ca1 += 1;
      ca2 += 1;
      ca3 += 1;
      tester.out() << "Receiving for iters = " << iters << std::endl;
      paws->update();
      tester.out() << "Receive complete." << std::endl;
      tester.check("a1 OK", all(a1.array() == ca1));
      tester.check("a2 OK", all(a2.array() == ca2));
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
  int retval = tester.results("Paws Field send/receive test B");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: paws_test6.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:23 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
