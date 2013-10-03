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
// RefCountedPtr test code.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/RefCounted.h"
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

  //  void print() const { tester.out() << "Value = " << d_m << std::endl; }

  int val() const {return d_m;}

private:

  int d_m;
};

typedef RefCountedPtr<SharedInt> RCIntPtr_t;

int main(int argc, char* argv[])
{
  // initialize Pooma
  Pooma::initialize(argc,argv);

  Pooma::Tester tester(argc, argv);

  int return_status = 0;

  int test_number = 0;
#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\nTesting RefCountedPtr." << std::endl;

    RCIntPtr_t pn;
    pn = new SharedInt(2);

    tester.out() << "pn->val() = " << pn->val() << std::endl;

    *pn = 5;

    tester.out() << "pn->val() = " << pn->val() << std::endl;

    RCIntPtr_t p1(new SharedInt(1));
    RCIntPtr_t p2(new SharedInt(2));
    RCIntPtr_t p3(new SharedInt(3));

    tester.out() << p1->val() << " " 
		 << p2->val() << " "
		 << p3->val() << " "
		 << std::endl;

    *p1 = *p2 = *p3 = -777;

    tester.out() << p1->val() << " " 
		 << p2->val() << " "
		 << p3->val() << " "
		 << std::endl;

    PAssert(*p1 == *p2 && *p1 == *p3 && *p2 == *p3);
    PAssert(p1 != p2 && p1 != p3 && p2 != p3);

    PAssert(!p1.isShared());
    PAssert(!p2.isShared());
    PAssert(!p3.isShared());
    PAssert(!pn.isShared());

    RCIntPtr_t pc = pn;

    PAssert(pn.isShared());
    PAssert(pc.isShared());

    PAssert(pn == pc);
    PAssert(*pn == *pc);

    pn.invalidate();

    PAssert(!pc.isShared());

    tester.out()<< "Value = " << pc->val() << std::endl;

    {
      RCIntPtr_t pn = pc;
      PAssert(pn == pc);
      PAssert(*pn == *pc);
      PAssert(pn.isShared());
      PAssert(pc.isShared());

      pn.makeOwnCopy();
      
      PAssert(pn != pc);
      PAssert(*pn == *pc);
      PAssert(!pn.isShared());
      PAssert(!pc.isShared());
      
      tester.out() << "Making copy and modifying. "
		   << "Next two shouldn't be the same" << std::endl;
      RCIntPtr_t pv = pn;
      pv.makeOwnCopy();
      
      *pv = pv->val() - 999;
      tester.out() << "Value = " << pv->val() << std::endl;
      tester.out() << "Value = " << pn->val() << std::endl;
    
      PAssert(*pv != *pn);
      tester.out() << std::endl;

      RCIntPtr_t p1 = pc;
      PAssert(p1.isShared());
      PAssert(pc.isShared());
      PAssert(*p1 == *pn);
    }

    PAssert(!pc.isShared());

    pc.invalidate();

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
  tester.results(" rcptr_test1 " );
  Pooma::finalize();  
  return return_status;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: rcptr_test1.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:17:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
