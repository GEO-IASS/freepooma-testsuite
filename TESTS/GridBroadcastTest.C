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

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Grid.h"
#include "Domain/Range.h"
#include "Tulip/RemoteProxy.h"

#define BARRIER

#ifndef BARRIER
# if POOMA_CHEETAH
#  define BARRIER Pooma::controller()->barrier()
# else
#  define BARRIER
# endif
#endif

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << "Testing Grid broadcasting . . ." << std::endl;

  // Hmmm... If I comment out the following output statement,
  // then the Inform object does not correctly limit itself to
  // context 0 only, nor does the barrier appear to work.
  // What is going on here???

  // JCC: I think I figured out what was wrong.  The first output 
  // statement above did not pass std::endl to the Inform object,
  // so no flush() was ever done.  Later on, the output context
  // for the Inform object is set to all contexts.  When std::endl
  // is finally passed to the Inform object, all the contexts print
  // out everything in the buffer since the last flush().  The second
  // output statement below had a std::endl in it already, causing a
  // a flush() of the buffer and output only by context 0.  I went 
  // ahead and added a std::endl to the first output statement above,
  // just to avoid any confusion or future problems.  The barrier seems
  // to work just fine, with or without these output statements.

  tester.out() << "Running on " << Pooma::contexts() 
  	       << " contexts." << std::endl;

  Grid<1> g;

  if (Pooma::context() == 0)
    {
      Range<1> r(0, 16, 2);

      g = Grid<1>(r);
    }

  RemoteProxy<Grid<1> > broadcast(g);
  Grid<1> ans = broadcast.value();

  BARRIER;

  tester.out().setOutputContext(-1);
  tester.out() << ans << std::endl;

  BARRIER;

  tester.check(ans == Grid<1>(Range<1>(0,16,2)));

  int ret = tester.results("GridBroadcast Test");
  Pooma::finalize();

  return ret;
}

