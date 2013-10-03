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

typedef RefCountedPtr<Shared<int> > RCIntPtr_t;

int main(int argc, char* argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\n\nTesting RefCountedPtr with Shared<int>.\n" 
		 << std::endl;

    RCIntPtr_t pn;
    pn = new Shared<int>(2);

    tester.out() << "pn->data() = " << pn->data() << std::endl;

    *pn = 5;

    tester.out() << "pn->data() = " << pn->data() << std::endl;

    RCIntPtr_t p1(new Shared<int>(1));
    RCIntPtr_t p2(new Shared<int>(2));
    RCIntPtr_t p3(new Shared<int>(3));

    tester.out() << p1->data() << " " 
         << p2->data() << " "
         << p3->data() << " "
         << std::endl;

    *p1 = *p2 = *p3 = -777;

    tester.out() << p1->data() << " " 
         << p2->data() << " "
         << p3->data() << " "
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

    tester.out() << "Value = " << pc->data() << std::endl;

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
      
      *pv = pv->data() - 999;
      tester.out() << "Value = " << pv->data() << std::endl;
      tester.out() << "Value = " << pn->data() << std::endl;
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
  int retval = tester.results("rcptr_test2");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: rcptr_test2.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
