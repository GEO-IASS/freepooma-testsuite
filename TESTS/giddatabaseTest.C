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
// Test of global ID database
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Layout/UniformGridLayout.h"
#include "Layout/GlobalIDDataBase.h"
#include "Utilities/Tester.h"

#include <iostream>

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Tests of global ID database." << std::endl;
  tester.out() << "------------------------------------------------" << std::endl;

  int size = 120;

  Interval<1> domain(size);

  UniformGridPartition<1> partition1(Loc<1>(10));
  UniformGridLayout<1> layout1(domain,partition1, ReplicatedTag());

  UniformGridPartition<1> partition2(Loc<1>(6));
  UniformGridLayout<1> layout2(domain,partition2, ReplicatedTag());

  std::vector<INode<1> > inodes;
  GlobalIDDataBase gidStore;

  UniformGridLayout<1>::const_iterator p = layout1.beginGlobal();
  while (p != layout1.endGlobal())
  {
    inodes.push_back(INode<1>(*p, layout1.ID(), &gidStore));
    ++p;
  }
  int i;
  int ni = inodes.size();
  for (i = 0; i < ni; ++i)
  {
    layout2.touches(inodes[i].domain(), std::back_inserter(inodes),
		    TouchesConstructINode<1>(layout2.ID(),inodes[i].key(),
					     &gidStore)
		    );
  }
  inodes.erase(inodes.begin(), inodes.begin() + ni);

  int gid1,gid2;
  int gidStore1,gidStore2;
  for (i = 0; i < inodes.size(); ++i)
  {
    gid1 = layout1.globalID(inodes[i].domain().firsts());
    gid2 = layout2.globalID(inodes[i].domain().firsts());
    gidStore1 = inodes[i].globalID(layout1.ID());
    gidStore2 = inodes[i].globalID(layout2.ID());

    tester.check(gid1 == gidStore1);
    tester.check(gid2 == gidStore2);

    tester.out() << "domain " << inodes[i].domain()
		 << ", key " << inodes[i].key()
		 << ", gid #1 - (" << gid1 << " == " << gidStore1 << ")"
		 << ", gid #2 - (" << gid2 << " == " << gidStore2 << ")"
		 << std::endl;
  }

  gidStore.print(tester.out());

  tester.out() << "------------------------------------------------"
	       << std::endl;

  int retval = tester.results("giddatabaseTest");
  Pooma::finalize();
  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: giddatabaseTest.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
