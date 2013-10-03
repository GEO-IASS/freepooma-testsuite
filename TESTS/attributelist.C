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
// Particles test: AttributeList and Attribute
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Particles/AttributeList.h"
#include "Domain/Interval.h"
#include "Layout/DynamicLayout.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "DynamicArray/DynamicArray.h"

#include <iostream>

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": AttributeList operations" << std::endl;
  tester.out() << "------------------------------------------------"
               << std::endl;

  // Create some Attributes

  tester.out() << "Creating DynamicArray objects for attributes ..."
               << std::endl;
  Interval<1> D(10);
  int blocks = 4;
  DynamicLayout layout(D,blocks);
  tester.out() << "DynamicLayout object:\n" << layout << std::endl;
#if POOMA_MESSAGING
  typedef MultiPatch< DynamicTag, Remote<Dynamic> > EngineTag_t;
#else
  typedef MultiPatch<DynamicTag,Dynamic> EngineTag_t;
#endif
  DynamicArray<int,EngineTag_t>   a1(layout);
  DynamicArray<long,EngineTag_t>  a2(layout);
  DynamicArray<float,EngineTag_t> a3(layout);

  // Initialize the arrays with scalars.
  // Block since we're starting scalar code.

  Pooma::blockAndEvaluate();
  
  tester.out() << "Initializing DynamicArray objects ..."
               << std::endl;
  int i;
  for (i=0; i < D.size(); ++i)
    {
      a1(i) = 10 + i;
      a2(i) = 100 + i;
      a3(i) = 0.1 * i;
    }
  tester.out() << "Initialization complete:" << std::endl;
  tester.out() << "  a1 = " << a1 << std::endl;
  tester.out() << "  a2 = " << a2 << std::endl;
  tester.out() << "  a3 = " << a3 << std::endl;

  // Add these to an AttributeList

  tester.out() << "Adding DynamicArray's to the AttributeList ..."
               << std::endl;
  AttributeList attriblist;
  attriblist.add(a1);
  attriblist.add(a2);
  attriblist.add(a3);
  tester.out() << "Added " << attriblist.size() << " attributes." << std::endl;
  tester.check(attriblist.size() == 3);

  // Delete some of the elements in the attributes

  tester.out() << "Deleting even-numbered elements ..." << std::endl;
  tester.out() << "Domain size before destroy = " << layout.domain().size()
               << std::endl;
  Range<1> killlist(D.first(), D.last(), 2);
  layout.destroy(killlist, BackFill());
  layout.sync();
  tester.out() << "Domain size after destroy = " << layout.domain().size()
               << std::endl;
  tester.check(layout.domain().size() == (D.size() - killlist.size()));

  // Loop through the attributes, printing them out

  tester.out() << "Current contents of attributes:" << std::endl;
  for (i=0; i < attriblist.size(); ++i)
    {
      tester.out() << "  attrib[" << i << "] = ";
      tester.out() << *(attriblist.attribute(i)) << std::endl;
    }

  // Multiply values together for some attributes

  tester.out() << "Multiplying a2 *= (a1 + a3) ..." << std::endl;
  a2 *= (a1 + a3);
  tester.out() << "Results:" << std::endl;
  tester.out() << attriblist << std::endl;

  // Return resulting error code and exit; Tester will shut down POOMA.

  tester.out() << "------------------------------------------------"
               << std::endl;
  int retval = tester.results("AttributeList operations");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: attributelist.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:00 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
