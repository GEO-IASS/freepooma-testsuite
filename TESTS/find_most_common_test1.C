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
// ----------------------------------------------------------------------
// Delete test #1.
//-----------------------------------------------------------------------------

#include "Utilities/algorithms.h"
#include "Utilities/Tester.h"
#include "Pooma/Pooma.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

int randN(int n);

void test(Pooma::Tester &t, int numElements);

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  Pooma::Tester tester(argc, argv);

  test(tester, 0);
  test(tester, 1);
  test(tester, 10);
  test(tester, 15);

  int res = tester.results("find_most_common_test1 " );
  Pooma::finalize();

  return res;
}

// Generate a random number in the range [1,...,N]:

int randN(int n)
{
  return int(1 + n*(double(std::rand())/RAND_MAX));
}

// Test function.

void test(Pooma::Tester &t, int numElements)
{  
  std::vector<int> v;
  std::map<int, int> m;
  
  int i;
  for (i = 0; i < numElements; i++)
    {
      int e = randN(numElements);
      v.push_back(e);
      std::map<int, int>::iterator ei = m.find(e);
      if (ei != m.end())
        (*ei).second += 1;
      else
        m[e] = 1;
    }
  
  std::sort(v.begin(), v.end());
  for (i = 0; i < numElements; i++)
    t.out() << v[i] << " ";
  t.out() << std::endl;
      
  std::vector<int>::iterator mc = 
    Pooma::Algorithms::find_most_common(v.begin(), v.end());

  if (numElements == 0)
    t.check ("zero length", mc == v.end());
  else
    {  
      int cmc, nmc = -1;
      std::map<int, int>::const_iterator j = m.begin();
      while (j != m.end())
        {
          if ((*j).second > nmc)
            {
              nmc = (*j).second;
              cmc = (*j).first;
            }
          t.out() << (*j).first << ":" << (*j).second << " ";
          ++j;
        }
  
      t.out() << std::endl;  
      t.check("most common", *mc, cmc);
    }
}
