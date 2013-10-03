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
// Class:
// AssertEquals
//-----------------------------------------------------------------------------

/** @file
 * @ingroup PETE
 * @brief
 * AssertEquals is a handy class for asserting conformance of an integer value
 * in expressions.
 */

#ifndef POOMA_PETE_ASSERTEQUALS_H
#define POOMA_PETE_ASSERTEQUALS_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * AssertEquals is a handy class for asserting conformance of an integer value
 * in expressions.  An example is the numPatches function that arrays have.
 * For the patch function to make sense on an array containing an expression,
 * all the arrays in the expression must have the same number of patches.
 * Scalars have zero patches, but we wish to ignore them, so AssertEquals
 * allows you to set an integer value that you wish to ignore.
 */

struct AssertEquals
{
  AssertEquals(int ignore = 0) : ignore_m(ignore) { }
  int ignore_m;
};

template<class Op>
struct Combine2<int, int, Op, AssertEquals>
{
  typedef int Type_t;
  inline static
  Type_t combine(const int &a, const int &b, const AssertEquals &ae)
  {
    int ret = a;
    if ((a != ae.ignore_m) && (b != ae.ignore_m))
    {
      PAssert(a == b);
    }
    else
    {
      if (b != ae.ignore_m) return ret = b;
    }
    return ret;
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PETE_ASSERTEQUALS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: AssertEquals.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:05 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
