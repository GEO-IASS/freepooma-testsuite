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

#include <iostream>
using std::cout;
using std::endl;

#include <stdlib.h>
#include <vector>

int main(int argc, char *argv[])
{
  cout << "sizeof(short)     = " << sizeof(short) << endl;
  cout << "sizeof(int)       = " << sizeof(int) << endl;
  cout << "sizeof(long)      = " << sizeof(long) << endl;
  cout << "sizeof(long long) = " << sizeof(long long) << endl;
  cout << "sizeof(float)     = " << sizeof(float) << endl;
  cout << "sizeof(double)    = " << sizeof(double) << endl;
  cout << "sizeof(size_t)    = " << sizeof(size_t) << endl;

  std::vector<int> a(10);
  cout << "sizeof(a.size())  = " << sizeof(a.size()) << endl;
}
