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
//   DataBlockPtr test code.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/DataBlockPtr.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"

#include <iostream>

typedef DataBlockPtr<double,true> RCBlock_t;
typedef DataBlockPtr<double,false> RCFBlock_t;

void err_report(const char *what, int n,Pooma::Tester &);
void foo(RCBlock_t, Pooma::Tester &);
void bar(RCBlock_t&, Pooma::Tester &);
void recurse(RCBlock_t a, Pooma::Tester &);

template <bool Check>
void print(const DataBlockPtr<double,Check> &b, Pooma::Tester &);

int main(int argc, char* argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int return_status = 0;

  int test_number = 0;

#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTesting DataBlockPtr." << std::endl;

    RCBlock_t p(10);

    test_number++;

    PAssert(!p.isShared());

    test_number++;

    for (int i = 0; i < 10; i++)
      p[i] = (i-5)*(i-5);

    test_number++;

    print(p,tester);

#if POOMA_EXCEPTIONS
    test_number++;

    try {

      for (int i = 0; i < 11; i++)
	p[i] = -p[i];

      throw "Bounds checking failed!";
    }
    catch(const Pooma::Assertion &) 
      { 
	tester.out() << "Bounds check worked." << std::endl; 
      }
#endif

    test_number++;

    for (int i = 0; i < 10; i++)
      PInsist( p[i] == *(p+i), "p[i] != *(p+i)" );

    for (int i = 0; i < 10; i++)
      PInsist( p[i] == *(p+i), "p[i] != *(p+i)" );

    test_number++;

    PAssert(!p.isShared());

    test_number++;

    foo(p,tester);
    
    test_number++;

    PAssert(!p.isShared());

    test_number++;

    bar(p,tester);

    PAssert(!p.isShared());

    test_number++;

    print(p,tester);

    RCBlock_t a;
    a++;

#if POOMA_EXCEPTIONS
    test_number++;

    try {
      tester.out() << a[4];
      throw "Bounds checking failed!";
    }
    catch(const Pooma::Assertion &) 
      { 
	tester.out() << "Bounds check worked." << std::endl; 
      }
#endif

    RCBlock_t q1 = p;

#if POOMA_EXCEPTIONS
    test_number++;

    try {
      RCBlock_t q2; q2 = p;

      PAssert(q1 == p);
      PAssert(q2 == p);
      PAssert(q1 == q2);

      PAssert(p.isShared());
      PAssert(q1.isShared());
      PAssert(q2.isShared());

      for (int i = 0; i < 10; i++)
	PAssert(q1[i] == q2[i]);
    }
    catch (...) { tester.out() << "Something is very wrong!" << std::endl; }
#endif

    PAssert(p.isShared());
    PAssert(q1.isShared());

    p[1] = -999;
    PAssert(q1[1] == -999);

    test_number++;

    p.invalidate();

    PAssert(!p.isValid());

#if POOMA_EXCEPTIONS
    try {
      tester.out() << p[3];
      throw "Bounds checking failed!";
    }
    catch(const Pooma::Assertion &) 
      { 
	tester.out() << "Bounds check worked." << std::endl; 
      }
#endif

    PAssert(!q1.isShared());

    test_number++;

    recurse(q1,tester);

    PAssert(!q1.isShared());
    tester.out() << "q1.isShared = " << q1.isShared() << std::endl;

    print(q1,tester);

    test_number++;

    {
      const RCBlock_t r = q1;

      PAssert(r.isShared());

      print(r,tester);

      for (int i = 0; i < 10; i++)
	tester.out() << *(r+i) << " ";

      tester.out() << std::endl;

      p = r;

      PAssert(p.isShared());
    }

    PAssert(p.isShared());

    test_number++;

    q1.invalidate();

    PAssert(!p.isShared());

    test_number++;

    tester.out() << "\nTesting conversions to non-boundschecked" << std::endl;
    RCFBlock_t s = p;

    PAssert(s.isShared());
    PAssert(p.isShared());
    PAssert(s == p);

    print(s,tester);

    // Hmmm. CW is not performing the implicit conversion here.

#ifdef __MWERKS__
    recurse(RCBlock_t(s),tester);
#else
    recurse(s,tester);
#endif

    PAssert(s.isShared());

    test_number++;

    s.makeOwnCopy();
    PAssert(!s.isShared());
    PAssert(!p.isShared());
    PAssert(s != p);

    for (int i = 0; i < 10; i++)
      s[i] = i*i;

    tester.out() << "These should not be the same." << std::endl;

    for (int i = 0; i < 10; i++)
      tester.out() << p[i] << " ";

    tester.out() << std::endl;
    for (int i = 0; i < 10; i++)
      tester.out() << s[i] << " ";

    tester.out() << std::endl;

    tester.out() << "printed ok that time." << std::endl;

    print(s,tester);
    print(p,tester);

    s.invalidate();
    
    PAssert(!p.isShared());

    p.invalidate();

#if POOMA_EXCEPTIONS
  }
  catch(const char *err) 
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
  catch(const Pooma::Assertion &err)
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
#endif    

  tester.out() << "All Done!" << std::endl;
  int res = tester.results(" dbptr_test1 " );
  Pooma::finalize();  
  return res;
}

void foo(RCBlock_t b,Pooma::Tester &tester)
{
  PAssert(b.isShared());

  for (int i = 0; i < 10; i++)
    tester.out() << *(b+i) << " ";

  tester.out() << std::endl;

  for (int i = 0; i < 10; i++)
    tester.out() << *b++ << " ";

  tester.out() << std::endl;

  for (int i = -10; i < 0; i++)
    tester.out() << b[i] << " ";

  tester.out() << std::endl;

#if POOMA_EXCEPTIONS
  try {
    double a = b[0];
    throw "Bounds check failed!";
  }
  catch(const Pooma::Assertion &a)
    { tester.out() << "Bounds check worked." << std::endl; }
#endif
}

void bar(RCBlock_t &b,Pooma::Tester &tester)
{
  PAssert(!b.isShared());

  for (int i = 0; i < 10; i++)
    tester.out() << *(b+i) << " ";

  tester.out() << std::endl;

  for (int i = 0; i < 10; i++)
    tester.out() << *b++ << " ";

  tester.out() << std::endl;

  for (int i = -10; i < 0; i++)
    tester.out() << b[i] << " ";

  tester.out() << std::endl;

#if POOMA_EXCEPTIONS
  try {
    double a = b[0];
    throw "Bounds check failed!";
  }
  catch(const Pooma::Assertion &a)
    { tester.out() << "Bounds check worked." << std::endl; }
#endif
  b -= 10;

}

void recurse(RCBlock_t a,Pooma::Tester &tester)
{
  static int depth = 0;
  if (depth++ < 10) 
    {
      tester.out() << "Depth = " << depth 
		   << ";\t a.isShared = " << a.isShared() 
		   << ";\t *a++ = " << *a++ << std::endl;
      recurse(a,tester);
    }
  depth--;
}

template<bool Checked>
void print(const DataBlockPtr<double,Checked> &b,Pooma::Tester &tester)
{
  DataBlockPtr<double,Checked> pi = b.begin();

  DataBlockPtr<double,Checked> end = b.end();

  while ( pi != end )
    tester.out() << *pi++ << " ";

  tester.out() << std::endl;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dbptr_test1.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
