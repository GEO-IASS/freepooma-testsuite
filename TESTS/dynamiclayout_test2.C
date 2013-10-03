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
// Test of DynamicLayout
//-----------------------------------------------------------------------------

// Include files

#include "Layout/DynamicLayout.h"
#include "Partition/GridPartition.h"
#include "Utilities/Tester.h"
#include "Pooma/Pooma.h"
#include "Domain/Range.h"
#include "Domain/Grid.h"
#if POOMA_NO_STRINGSTREAM
	#include <strstream>
#else
        #include <sstream>
#endif

#define BARRIER

#ifndef BARRIER
#if POOMA_CHEETAH
# define BARRIER Pooma::controller()->barrier()
#else
# define BARRIER
#endif
#endif

void printLayout(Inform &out, const std::string msg, const DynamicLayout &l);

void printLayout(Inform &out, const std::string msg, const DynamicLayout &l)
{
  DynamicLayout::const_iterator pos;
  const int numContexts = Pooma::contexts();
  const int myContext   = Pooma::context();

  BARRIER;

  out.setOutputContext(0);
  out << msg << std::endl;
  out.setOutputContext(-1);

  BARRIER;

// This looks like a silly loop, but with the BARRIER it causes the contexts
// to print things in order.

  for (int i = 0; i < numContexts; ++i)
    {
      if (myContext == i)
	for (pos = l.beginLocal(); pos != l.endLocal(); ++pos)
	  out << pos->domain() << std::endl;
      BARRIER;
    }

  out.setOutputContext(0);
  out << "Total Domain = " << l.domain() << std::endl;

  BARRIER;
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  const int numContexts = Pooma::contexts();
  const int myContext   = Pooma::context();

  tester.out() << "Testing dynamic ops in DynamicLayout class . . .\n";
  tester.out() << "Running with " << numContexts << " contexts." << std::endl;

  Interval<1> domain;
  int numBlocks = numContexts * 5;

  tester.out() << "Creating DynamicLayout with domain " << domain
               << " and " << numBlocks << " blocks." << std::endl;
               
  BARRIER;

  Loc<1> foo(numBlocks);
  GridPartition<1> gp(foo);
  DistributedMapper<1> cmap(gp);
  DynamicLayout layout(domain, gp, cmap);

  DynamicLayout::iterator pos;

  std::string msg("Here are the patch domains for the initial partitioning:");
  printLayout(tester.out(), msg, layout);

  // Creates 35 elements in the first patch of each subdomain.

  for (int i = 0; i < numContexts; ++i)
    {
      if (i == myContext)
        {
          layout.create(35,0);
	  layout.create(10,1);
#if POOMA_NO_STRINGSTREAM
	  std::strstream s;
#else
          std::ostringstream s;
#endif
          s << "Here are the patch domains after adding elements\n"
            << "to the first two patches on context " << i 
	    << ", before syncing.";
          printLayout(tester.out(), s.str(), layout);
	  BARRIER;

          layout.sync();
#if POOMA_NO_STRINGSTREAM
	  std::strstream ss;
#else
          std::ostringstream ss;
#endif       
	  
          ss << "Here are the patch domains on context " << i 
	    << ", after syncing.";
          printLayout(tester.out(), ss.str(), layout);
          BARRIER;
        }
    }

  BARRIER;

  int ret = tester.results("DynamicLayout Test2");
  Pooma::finalize();

  return ret;
}

