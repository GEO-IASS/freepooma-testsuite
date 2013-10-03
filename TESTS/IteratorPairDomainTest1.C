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

#include "Domain/IteratorPairDomain.h"
#include "Utilities/Tester.h"
#include "Pooma/Pooma.h"

#include <iostream>
#include <vector>
#include <list>

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  using std::vector;
  using std::list;

  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  tester.out() << "Starting IteratorPairDomain test.\n" << endl;
  tester.out() << "First testing with std::vector..." << endl;

  {
    vector<int> klist(7); // klist(0) ... klist(6)

    klist.assign(klist.size(),1);

    for (int i = 1; i < 7; ++i)
      klist[i] = klist[i-1] + i;
    klist[2] = 3;
    klist[5] = 12;
    klist[6] = 20;
  
    typedef vector<int>::iterator Iter_t;
    typedef Pooma::IteratorPairDomain<Iter_t> IPDomain_t;
  
    IPDomain_t dom(klist.begin(),klist.end());
  
    tester.out() << "dom   = " << dom << endl;
    tester.out() << "max   = " << dom.max() << endl;
    tester.out() << "min   = " << dom.min() << endl;
    tester.out() << "first = " << dom.first() << endl;
    tester.out() << "last  = " << dom.last() << endl;
    tester.out() << "size  = " << dom.size() << endl;

    IPDomain_t::iterator pos = dom.begin();
    for (int i = 0; i < 7; ++i)
      {
        tester.check(klist[i] == *pos++);
      }

    IPDomain_t dom2 = dom;

    pos = dom2.begin();
    int i = 0;
    while (pos < dom2.end())
      {
        tester.check(klist[i++] == *pos++);
      }

    IPDomain_t dom3, tmp;

    tester.check(!dom3.initialized());
    tester.check(!tmp.initialized());
    tester.check(dom3.size() == 0);
    tester.check(tmp.size() == 0);
    tester.check(dom3.length() == 0);
    tester.check(tmp.length() == 0);

    dom3 = (tmp = dom);
    for (int i = 0; i < dom3.size(); ++i)
      {
        tester.check(klist[i] == dom3(i));
      }

    dom3(3) = 100;
    pos = dom3.begin() + 5;
    *pos = -201;

    tester.out() << "dom   = " << dom << endl;
  }

  tester.out() << "\nRepeating same test with a std::list..." << endl;

  {
    vector<int> vlist(7); // vlist(0) ... vlist(6)

    vlist.assign(vlist.size(),1);

    for (int i = 1; i < 7; ++i)
      vlist[i] = vlist[i-1] + i;
    vlist[2] = 3;
    vlist[5] = 12;
    vlist[6] = 20;

    list<int> klist; 

    for (int i = 0; i < 7; ++i)
      klist.push_back(vlist[i]);

    typedef list<int>::iterator Iter_t;
    typedef Pooma::IteratorPairDomain<Iter_t> IPDomain_t;
  
    IPDomain_t dom(klist.begin(),klist.end());
  
    tester.out() << "dom   = " << dom << endl;
    tester.out() << "max   = " << dom.max() << endl;
    tester.out() << "min   = " << dom.min() << endl;
    tester.out() << "first = " << dom.first() << endl;
    tester.out() << "last  = " << dom.last() << endl;
    tester.out() << "size  = " << dom.size() << endl;

    IPDomain_t::iterator pos = dom.begin();
    list<int>::iterator lpos = klist.begin();
    for (int i = 0; i < 7; ++i)
      {
        tester.check(*lpos++ == *pos++);
      }

    IPDomain_t dom2 = dom;

    pos = dom2.begin();
    lpos = klist.begin();
    while (pos != dom2.end())
      {
        tester.check(*lpos++ == *pos++);
      }

    IPDomain_t dom3, tmp;

    tester.check(!dom3.initialized());
    tester.check(!tmp.initialized());
    tester.check(dom3.size() == 0);
    tester.check(tmp.size() == 0);
    tester.check(dom3.length() == 0);
    tester.check(tmp.length() == 0);

    dom3 = (tmp = dom);
    lpos = klist.begin();
    for (int i = 0; i < dom3.size(); ++i)
      {
        tester.check(*lpos++ == dom3(i));
      }

    dom3(3) = 100;
    pos = dom3.begin();
    std::advance(pos,5);
    *pos = -201;

    tester.out() << "dom   = " << dom << endl;
  }

  tester.out() << "Finished IteratorPairDomain test 1.\n" << endl;

  int res = tester.results("IteratorPairDomainTest1");
  Pooma::finalize();
  return res;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IteratorPairDomainTest1.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
