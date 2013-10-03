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
// DataBlockPtr test code.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/DataBlockPtr.h"
#include "Utilities/ElementProperties.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"

#include <iostream>

// This specialized version of ElementProperties makes a deep copy
// of a DataBlockPtr.

template <class T, bool b>
struct ElementProperties< DataBlockPtr<T,b> >
  : MakeOwnCopyProperties< DataBlockPtr<T,b> >
{ };

typedef DataBlockPtr<double,true> RCBlock_t;
typedef DataBlockPtr<double,false> RCFBlock_t;
typedef DataBlockPtr<RCBlock_t,true> RCBlock2D_t;

void print(const RCBlock_t &b,Pooma::Tester &);


int main(int argc, char* argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int return_status = 0;
  int test_number = 0;

#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTest that Block<Block<T> > works" 
		 << std::endl;

    RCBlock_t p(10);
    for (int i = 0; i < 10; i++)
      p[i] = (i-5)*(i-5);

    print(p,tester);

    tester.out() << std::endl;

    PAssert(!p.isShared());

    test_number = 1;

    RCBlock2D_t a(5,p);

    test_number = 2;

    PAssert(!p.isShared());

    test_number = 3;

    for (int i = 0; i < 5; i++)
      { PAssert(!a[i].isShared()); }

    test_number = 4;

    for (int i = 0; i < 5; i++)
      print(a[i],tester);

    test_number = 5;

    for (int i = 0; i < 5; i++)
      a[i][3] = -1;

    test_number = 6;

    for (int i = 0; i < 5; i++)
      print(a[i],tester);

    {
      RCBlock2D_t b = a;

      PAssert(b.isShared());
      PAssert(a.isShared());
      
      // Hmmm. The subarrays themselves will NOT show up
      // as being shared. This is because they are NOT. Only
      // the outer block controller was copied. Each copy of
      // the outer controller has a pointer to a single
      // inner block controller for each row. This is
      // a bit confusing, but can it actually cause problems?

      for (int i = 0; i < 5; i++)
        { 
	  PAssert(!b[i].isShared());
	}

      for (int i = 0; i < 5; i++)
	{ 
	  PAssert(!a[i].isShared()); 
	}
      
      b.makeOwnCopy();
      
      PAssert(!b.isShared());
      PAssert(!a.isShared());
      
      b[0][0] = 0;
      b[0][1] = 0;
      b[1][1] = 0;
      b[1][2] = 0;
      b[2][2] = 0;
      b[2][3] = 0;
      b[3][3] = 0;
      b[3][4] = 0;
      b[4][4] = 0;
      b[4][5] = 0;
      
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 10; j++)
          if (j > i) b[i][j] = 0;

#if POOMA_EXCEPTIONS      
      try {
        b[5][5] = 0;
        throw "Bounds checking failed!";
      } 
      catch(const Pooma::Assertion &) 
	{ 
	  tester.out() << "Bounds check worked." << std::endl; 
	}
#endif

      for (int i = 0; i < 5; i++)
        print(a[i],tester);
      
      for (int i = 0; i < 5; i++)
        print(b[i],tester);
        
      RCBlock2D_t c = a;
      
      PAssert(a.isShared());
      PAssert(c.isShared());   
      
    }

    PAssert(!a.isShared());
     
    for (int i = 0; i < 5; i++)
      { 
	PAssert(!a[i].isShared()); 
      }

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
  int res = tester.results(" dbptr_test3 " );
  Pooma::finalize();
  return res;
}


void print(const RCBlock_t &b,Pooma::Tester & tester)
{
  for (int i = 0; i < 10; i++)
    tester.out() << b[i] << " ";

  tester.out() << std::endl;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: dbptr_test3.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
