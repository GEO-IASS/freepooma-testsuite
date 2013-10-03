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
// Paws test 5: Send and Receive an int and double set of scalars,
//              plus a fixed-size 2D Field, in conjunction with test 6.
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
  tester.out() << argv[0] << ": Paws Field send/receive test A" << std::endl;
  tester.out() << "--------------------------------------------" << std::endl;

#if POOMA_PAWS

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  int iters = 10;

  // Fields to send and receive ... use different layouts in the two
  // test codes.

  Loc<2> blocks(3,2);
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

  Field<Geometry_t, float, MultiPatch<GridTag,Brick> > a1(geom, layout);
  Field<Geometry_t, int, Brick> a2(geom);

  // Initialize the arrays

  a1 = 10 * (iota(domain).comp(1) + 1) + iota(domain).comp(0) + 1;
  a2 = a1 + 1000;

  // Create a Paws connection

  tester.out() << "Creating PawsConnection object ..." << std::endl;
  Connection<Paws> *paws = new Connection<Paws>("test5", argc, argv);
  tester.out() << "Finished creating PawsConnection object." << std::endl;

  // Establish connections for the two scalars

  tester.out() << "Connecting s1 = " << s1 << " for output ..." << std::endl;
  ConnectorBase *s1p = paws->connectScalar("s1", s1, ConnectionBase::out);
  tester.out() << "Connecting s2 = " << s2 << " for input ..." << std::endl;
  ConnectorBase *s2p = paws->connectScalar("s2", s2, ConnectionBase::in);
  tester.out() << "Connecting iters = " << iters << " for output ...";
  tester.out() << std::endl;
  ConnectorBase *iterp=paws->connectScalar("iters",iters,ConnectionBase::out);

  // Establish connections for the two fields; also connect up a view of
  // the first fields

  tester.out() << "Connecting a1 = " << a1 << " for output ..." << std::endl;
  paws->connect("a1", a1, ConnectionBase::out);
  tester.out() << "Connecting a2 = " << a2 << " for output ..." << std::endl;
  paws->connect("a2", a2, ConnectionBase::out);
  tester.out() << "Connecting a1(" << subdomain << ") = " << a1(subdomain);
  tester.out() << " for output ..." << std::endl;
  paws->connect("a1view", a1(subdomain), ConnectionBase::out);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  paws->ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Modify s2, and update

  s2 *= 2;
  tester.out() << "Updating current s1 = " << s1 << " and s2 = " << s2;
  tester.out() << ", plus fields ..." << std::endl;
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
  paws->disconnect(s1p);
  paws->disconnect(s2p);
  paws->disconnect(iterp);
  tester.check("3 less connections", paws->size() == (connections - 3));

  // Do, in a loop, updates of the receiver. Add one to the arrays each time.

  while (iters-- > 0)
    {
      a1 += 1;
      a2 += 1;
      tester.out() << "Sending for iters = " << iters << std::endl;
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
  int retval = tester.results("Paws Field send/receive test A");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: paws_test5.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:23 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
