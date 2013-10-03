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

  tester.out() << "Testing vector broadcasting . . ." << std::endl;

  tester.out() << "Running on " << Pooma::contexts() 
  	       << " contexts." << std::endl;

  std::vector<int> v(10);

  if (Pooma::context() == 0)
    {
      for (int i = 0; i < 10; i++)
	v[i] = i;
    }

  RemoteProxy<std::vector<int> > broadcast(v);
  std::vector<int> ans = broadcast.value();

  BARRIER;

  tester.out().setOutputContext(-1);
  tester.out() << ans[3] << std::endl;
  tester.out() << ans.size() << std::endl;
  BARRIER;

  tester.check(ans[3] == 3);

  int ret = tester.results("Vector Broadcast Test");
  Pooma::finalize();

  return ret;
}

