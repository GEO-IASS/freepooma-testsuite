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
// MetaTokenIteratorTest1: Test the MetaTokenIterator class.
//-----------------------------------------------------------------------------

#include "IO/MetaTokenIterator.h"
#include "Utilities/Tester.h"
using namespace Pooma;

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>  // copy
#include <iterator>   // back_inserter
#include <stdlib.h>   // atoi

using std::string;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  Pooma::Tester tester(argc, argv);

  // Test lines

  string s0 = "Type = double = 1 2 3 4 #This should be ignored";
  string s1 = "Type= int=";
  string s2 = "Type =float";
  string s3 = "Type long";

  MetaTokenIterator pend; // end iterator

  // First test that the line is parsed correctly. We only do this for
  // the first line.

  MetaTokenIterator pt(s0);
  tester.check(*pt++ == "Type");
  tester.check(pt != pend);
  tester.check(*pt++ == "double");
  tester.check(pt != pend);
  tester.check(*pt++ == "=");
  tester.check(pt != pend);
  tester.check(*pt++ == "1");
  tester.check(pt != pend);
  tester.check(*pt++ == "2");
  tester.check(pt != pend);
  tester.check(*pt++ == "3");
  tester.check(pt != pend);
  tester.check(*pt++ == "4");
  tester.check(pt == pend);
  
  // Test the iterator functionality

  MetaTokenIterator pv(s0);
  std::vector<std::string> v;
  std::copy(pv, pend, std::back_inserter(v));
  tester.check(v.size() == 7);
  tester.check(v[0] == "Type");
  tester.check(v[1] == "double");
  tester.check(v[2] == "=");
  tester.check(v[3] == "1");
  tester.check(v[4] == "2");
  tester.check(v[5] == "3");
  tester.check(v[6] == "4");

  // Another test of iterator functionality
  // (I'm not sure this will work in parallel - what does
  // Inform::stream() return on a context that isn't doing I/O???)

  MetaTokenIterator ps(s0);
  std::copy(pv, pend, 
            std::ostream_iterator<std::string>(tester.out().stream(), " "));
  tester.out() << std::endl;

  // For the other lines, just print out the details.

  MetaTokenIterator p0(s0);
  while (p0 != pend) tester.out() << *p0++ << " ";
  tester.out() << endl;
  MetaTokenIterator p1(s1);
  while (p1 != pend) tester.out() << *p1++ << " ";
  tester.out() << endl;
  MetaTokenIterator p2(s2);
  while (p2 != pend) tester.out() << *p2++ << " ";
  tester.out() << endl;
  MetaTokenIterator p3(s3);
  while (p3 != pend) tester.out() << *p3++ << " ";
  tester.out() << endl;

  MetaTokenIterator pn(s0);
  tester.out() << "Testing operator->" << endl;
  tester.out() << "Skipping words: ";
  while (*pn != "1") 
    {
      tester.out() << *pn++ << " ";
    }
  tester.out() << endl;
  tester.out() << "Found first number." << endl;
  tester.out() << "Numbers are: ";
  while (pn != pend)
    {
      int i = atoi(pn++->c_str());
      tester.out() << i << " ";
    }
  tester.out() << endl;

  int ret = tester.results("MetaTokenIteratorTest1");
  
  return ret;
}

