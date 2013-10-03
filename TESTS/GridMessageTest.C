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
// Test of Grid Serialize specialization and JAC's intro to sending
// messages with MatchingHandler 8-)
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Grid.h"
#include "Domain/Range.h"

#if POOMA_MESSAGING
#include "Tulip/Messaging.h"
#endif

#define BARRIER

#ifndef BARRIER
# if POOMA_CHEETAH
#  define BARRIER Pooma::controller()->barrier()
# else
#  define BARRIER
# endif
#endif

bool gotIt = false;

typedef Grid<1> Send_t;

void receiveGrid(Send_t *lg, Send_t &rg)
{
  *lg = rg;
  gotIt = true;
}
 
int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

#if POOMA_CHEETAH

  typedef Cheetah::MatchingHandler Handler_t;
  
  const int numContexts = Pooma::contexts();
  const int myContext   = Pooma::context();

  Handler_t *handler = new Cheetah::MatchingHandler(*Pooma::controller());

  tester.out() << "Testing Grid messages . . .\n";
  tester.out() << "Running with " << numContexts << " contexts." << std::endl;

  int start = myContext * 10;
  int end   = (myContext + 1) * 10;

  Range<1> r(start, end, 2);

  Grid<1> foo(r);

  tester.out() << "Here are our Grids..." << std::endl;

  BARRIER;

  tester.out().setOutputContext(-1);
  tester.out() << foo << std::endl;

  // Here's the message pattern - we're just sending in a ring:

  int toContext   = (myContext + 1)               % numContexts;
  int fromContext = (myContext + numContexts - 1) % numContexts;

  BARRIER;

  tester.out() << "Node "             << myContext 
	       << ";   Sending to "     << toContext 
	       << ";   Receiving from " << fromContext << std::endl;

  int msgtag = 0;

  handler->send(toContext, msgtag, foo);

  Send_t bar;

  handler->request(fromContext, msgtag, receiveGrid, &bar);

  while (!gotIt) Pooma::poll();

  BARRIER;

  tester.out().setOutputContext(0);
  tester.out() << "Here are the Grids we received:" << std::endl;

  BARRIER;

  tester.out().setOutputContext(-1);
  tester.out() << bar << std::endl;

  start = fromContext * 10;
  end   = (fromContext + 1) * 10;

  Range<1> ans(start, end, 2);

  BARRIER;

  tester.check(bar == ans);

  delete handler;

#else

  tester.out() << "This test requires Cheetah" << std::endl;

#endif

  int ret = tester.results("GridMessage Test");
  Pooma::finalize();

  return ret;
}

