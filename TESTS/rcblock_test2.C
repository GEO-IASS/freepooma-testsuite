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
// RefCountedBlockPtr test code.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/RefCountedBlockPtr.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"
#include <iostream>

class SharedInt : public RefCounted
{
public:
  SharedInt(int i) : d_m(i) {};
  SharedInt(const SharedInt &model) : d_m(model.d_m) { }

  SharedInt & operator=(const SharedInt &model)
    { 
      if (&model == this) return *this;
      d_m = model.d_m;
      return *this;
    }

  SharedInt & operator=(int i) { d_m = i; return *this; }

  bool operator==(const SharedInt &rhs) const
    { return d_m == rhs.d_m; }

  bool operator!=(const SharedInt &rhs) const
    { return d_m != rhs.d_m; }

  int val() const {return d_m;}

private:

  int d_m;
};


typedef RefCountedBlockPtr<SharedInt,true> SBlock_t;

void err_report(const char *what, int n, Pooma::Tester &);

int main(int argc, char* argv[])
{
  // initialize Pooma
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTest that blocks work if T has no T()" 
		 << std::endl;

    SBlock_t foo(10,SharedInt(3));

    foo[2] = 2;
    foo[6] = 8;

#ifdef TEST_3B
    {
      SBlock_t bar(2); // illegal - SharedInt has no default constructor.
    
      foo[8] = bar[0];
      foo[9] = bar[1];
    }
#endif
    
    for(int i = 0; i < 10; i++)
      tester.out()<<  "Value = " <<foo[i].val() << std::endl;

    SBlock_t bar = foo;
    
    PAssert(foo.isShared());
    PAssert(bar.isShared());
    
    for(int i = 0; i < 10; i++)
      tester.out()<<  "Value = "<< bar[i].val() << std::endl;

	bar.makeOwnCopy();
	
	PAssert(!foo.isShared());
	PAssert(!bar.isShared());
	
	bar[0] = -111;
	bar[1] = -222;

    for(int i = 0; i < 10; i++)
     tester.out()<<  "Value = "<< bar[i].val() << std::endl;
    
    for(int i = 0; i < 10; i++)
     tester.out()<<  "Value = "<< foo[i].val() << std::endl;
    
    bar.invalidate();
    foo.invalidate();
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
  int ret = tester.results("rcblock_test2");
  Pooma::finalize();  
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: rcblock_test2.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:17:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
