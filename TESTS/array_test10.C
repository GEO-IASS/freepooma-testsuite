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
// Array test 10: compressible brick data objects.
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Engine/BrickEngine.h"
#include "Engine/CompressibleBrick.h"
#include "Array/Array.h"

static bool OK = true;

inline void check(bool ans, bool correct, Pooma::Tester &tester)
{
  OK = (OK && ans == correct);
  tester.check(OK);
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  Pooma::blockingExpressions(true);

  // Checking compressible brick arrays.

  typedef Array<1, double, CompressibleBrick> Array_t;
  typedef Array<1, double, BrickView> View_t;

  // We're going to be tracking the compressed status and the data
  // object for the CBC contained by several views of a compressible
  // brick array.

  Pooma::DataObject_t *obj, *tstobj, *bobj, *cobj;
  bool status;

  // Allocate a compressible array:

  Array_t a(100);

  // Assign some data:

  a = 3.0;

  // See if it is compressed.

  status = a.engine().compressed();

  check(status, true, tester);

  // Get the compressed data

  double & a_ref = a.engine().compressedReadWrite();
  double * a_ptr = &a_ref;

  check(a_ref == 3.0, true, tester);

  // Get a's data object pointer and save it for later comparisons.

  obj = a.engine().dataObject();

  check(obj != 0, true, tester);

  // Now make a copy of a:

  Array_t b(a);

  // See if it is compressed.

  status = b.engine().compressed();

  check(status, true, tester);

  // Get the compressed data

  double & b_ref = b.engine().compressedReadWrite();
  double * b_ptr = &b_ref;

  check(b_ref == 3, true, tester);
  check(b_ptr == a_ptr, true, tester);

  // Get b's data object pointer and compare to a's.

  tstobj = b.engine().dataObject();

  check(tstobj != 0, true, tester);
  check(tstobj == obj, true, tester);

  {
    // Cause a to uncompress.

    a(10) = 5;

    status = a.engine().compressed();

    check(status, false, tester);

    status = b.engine().compressed();

    check(status, false, tester);

    tstobj = a.engine().dataObject();

    check(tstobj == obj, true, tester);

    tstobj = b.engine().dataObject();

    check(tstobj == obj, true, tester);

    // Now make a brick-view and test it.

    View_t av(a);

    status = a.engine().compressed();

    check(status, false, tester);

    tstobj = a.engine().dataObject();

    check(tstobj == obj, true, tester);

    tstobj = av.engine().dataObject();

    check(tstobj == obj, true, tester);
  }

  // Should still be compressed.

  status = a.engine().compressed();

  check(status, false, tester);

  tstobj = a.engine().dataObject();

  check(tstobj == obj, true, tester);

  {
    // Make it compressible again.

    a(10) = 3.0;

    status = a.engine().compressed();
    check(status, false, tester);

    tstobj = a.engine().dataObject();
    check(tstobj == obj, true, tester);

    // Take another view. When it goes out of scope, the array should
    // recompress.

    View_t bv(b);

    tstobj = b.engine().dataObject();
    check(tstobj == obj, true, tester);

    tstobj = bv.engine().dataObject();
    check(tstobj == obj, true, tester);
  }

  status = a.engine().compressed();
  check(status, true, tester);

  tstobj = a.engine().dataObject();
  check(tstobj == obj, true, tester);

  tstobj = b.engine().dataObject();
  check(tstobj == obj, true, tester);

  {
    View_t bv(b);

    status = b.engine().compressed();
    check(status, false, tester);

    tstobj = b.engine().dataObject();
    check(tstobj == obj, true, tester);

    tstobj = bv.engine().dataObject();
    check(tstobj == obj, true, tester);
  }

  status = b.engine().compressed();
  check(status, true, tester);

  tstobj = b.engine().dataObject();
  check(tstobj == obj, true, tester);

  tstobj = a.engine().dataObject();
  check(tstobj == obj, true, tester);

  // Now for the makeOwnCopy test....
  // First test it with the existing compressed data.

  b.makeOwnCopy();

  bobj = b.engine().dataObject();
  check(bobj != obj, true, tester);

  status = b.engine().compressed();
  check(status, true, tester);

  // Take a BrickView and make sure that the resulting view has the
  // same DataObject.

  {
    View_t bv(b);

    status = b.engine().compressed();
    check(status, false, tester);

    tstobj = b.engine().dataObject();
    check(tstobj == bobj, true, tester);

    tstobj = bv.engine().dataObject();
    check(tstobj == bobj, true, tester);
  }

  // Now try it with an uncompressed array.

  a(7) = 45;

  status = a.engine().compressed();
  check(status, false, tester);

  tstobj = a.engine().dataObject();
  check(tstobj == obj, true, tester);

  Array_t c(a);

  status = c.engine().compressed();
  check(status, false, tester);

  tstobj = c.engine().dataObject();
  check(tstobj == obj, true, tester);

  tstobj = a.engine().dataObject();
  check(tstobj == obj, true, tester);

  // Now make our own copy:

  c.makeOwnCopy();

  cobj = c.engine().dataObject();
  check(cobj != obj, true, tester);

  status = c.engine().compressed();
  check(status, false, tester);

  // And do the view test again...
  
  {
    View_t cv(c);

    status = c.engine().compressed();
    check(status, false, tester);

    tstobj = c.engine().dataObject();
    check(tstobj == cobj, true, tester);

    tstobj = cv.engine().dataObject();
    check(tstobj == cobj, true, tester);
  }

  status = c.engine().compressed();
  check(status, false, tester);

  tstobj = c.engine().dataObject();
  check(tstobj == cobj, true, tester);

  // Force c to recompress.
  {
    View_t cv(c);

    cv(7) = 3.0;
  }

  status = c.engine().compressed();
  check(status, true, tester);

  tstobj = c.engine().dataObject();
  check(tstobj == cobj, true, tester);

  // They should be disconnected, so a should
  // be left uncompressed.

  check(a.read(7) == 45.0, true, tester);

  status = a.engine().compressed();
  check(status, false, tester);

  tstobj = a.engine().dataObject();
  check(tstobj == obj, true, tester);

  bool etst = Array_t::Engine_t::hasDataObject;

  check(etst, true, tester);

  int ret = tester.results( "array_test10" );
  Pooma::finalize();
  return ret; 

}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test10.cpp,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
