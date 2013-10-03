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
// DomainLayout test: Create and use DomainLayout objects
//-----------------------------------------------------------------------------


#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"
#include "Layout/DomainLayout.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": DomainLayout operations." << std::endl;
  tester.out() << "----------------------------------------" << std::endl;

  Interval<1> I1(0,9);
  Interval<2> I2(I1,I1);
  
  DomainLayout<2> layout1(I2, GuardLayers<2>(2));
  DomainLayout<2> layout2(layout1);

  tester.out() << layout1 << std::endl;
  tester.out() << layout2 << std::endl;
 
  tester.out() << layout1.externalGuards() << std::endl;
  tester.out() << layout2.externalGuards() << std::endl;

  tester.check("correct external guards",
	       layout1.externalGuards() == layout2.externalGuards());
 
  tester.check("correct domains",
	       layout1.domain() == layout2.domain());
 
  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("DomainLayout operations");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: domainLayout.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:55 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
