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
// Function:
//   bool findLeftCommonEndpoint(a0,a1,s,b0,b1,t,endpoint)
//   bool findIntersectionEndpoints(a0,a1,s,b0,b1,t,i0,i1,is)
//-----------------------------------------------------------------------------

// include files
#include "Domain/DomainCalculus.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
//
// findLeftCommonEndpoint is used by the domain calculus routines to
// find the leftmost common endpoint for two domains defined by
// [a0:a1:s] and [b0:b1:t].  If an endpoint is possible, this returns
// true and the endpoint in the final argument.  If one is not possible,
// due to incompatible striding, this returns false and leaves the final
// argument unchanged.
//
//-----------------------------------------------------------------------------

bool findLeftCommonEndpoint(int a0, int a1, int s, int b0, int b1, int t,
			    int &endpoint)
{
  // We should have both domains moving in positive directions, although
  // we might have to reverse the sign on the strides.
  PAssert(a0 <= a1 && b0 <= b1);
  if (s < 0)  s = -s;
  if (t < 0)  t = -t;

  // Consider the values here as defining two domains
  //    a = [a0:a1:s]
  //    b = [b0:b1:t]
  // We must find out if there is at least one point in a that is also in b.
  // To do this, we do the following:
  //   1) Adjust the endpoints so that they always have a0 < a1, b0 < b1,
  //      with positive stride, and so that a0 < b0.
  //   2) Starting as close to (but not greater than) the b0 endpoint as
  //      possible, start two test values growing toward larger values,
  //      increasing by the strides, until the test values are the same
  //      or we determine they can never be the same.
  //   3) The domains touch if the test values are the same at a point
  //      which is less than a1 and b1.

  // Find the proper endpoints and stride to get a0 < a1, b0 < b1, a0 <= b0,
  // positive strides
  if (a0 > b0) {
    int tmp;
    tmp = a0;  a0 = b0;  b0 = tmp;
    tmp = a1;  a1 = b1;  b1 = tmp;
    tmp =  s;   s =  t;   t = tmp;
  }

  // start out the two test probes, one just below b0, the other at b0.
  int i1 = b0 - ((b0 - a0) % s);
  int i2 = b0;

  // grow both by s or t until they match, get too big, or we determine we
  // will never get a match
  int maxdiff = 0;
  int minright = a1;
  if (b1 < a1)
    minright = b1;
  while (i1 <= minright && i2 <= minright) {
    // move the first probe value forward, until it passes the second
    while (i1 < i2)
      i1 += s;

    // examine the difference between the probe values
    int newdiff = i1 - i2;
    if (i1 == i2 || newdiff == maxdiff) {
      // if the probes match, break out of this loop.  Or, if the maximum
      // difference between i1 and i2 when i1 gets larger than i2 is seen
      // more than once, then we know that we'll never be able to find
      // values of i1 and i2 that are equal, so we can quit now.
      break;
    } else if (newdiff > maxdiff) {
      maxdiff = newdiff;
    }

    // move the second probe value forward
    i2 += t;
  }

  // at this point, either the probes match, or they don't.  If they do,
  // the domains touch.
  if (i1 == i2 && i1 <= minright) {
    endpoint = i1;
    return true;
  } else {
    return false;
  }
}


//-----------------------------------------------------------------------------
//
// findLCM calculates the least-common-multiple of the two arguments.
//
//-----------------------------------------------------------------------------

int findLCM(int s, int t)
{
  // Both of the values must be positive
  PAssert(s > 0 && t > 0);

  // This code works for s < t.  If the opposite is true, swap 'em.
  int i1 = s, i2 = t;
  if (s > t) {
    s = t;
    t = i1;
    i1 = s;
    i2 = t;
  }

  // Start advancing the two probes i1 and i2 in steps of s and t, until they
  // both get to the same value.
  while (i1 != i2) {
    while (i1 < i2)
      i1 += s;
    if (i1 > i2)
      i2 += t;
  }

  // At this point, we should have i1 == i2 == LCM.
  return i1;
}


//-----------------------------------------------------------------------------
//
// findIntersectionEndpoints is used by the domain calculus routines to
// find the endpoints and stride of an intersection domain given
// two other strided domains defined by [a0:a1:s] and [b0:b1:t].  If
// this is possible, this return true and the min, max, stride of the
// resulting domain [i0,i1,is] in the final arguments.  If no intersection
// is possible, due to incompatible striding, this returns false and
// leaves the final arguments unchanged.  Note that if a domain is returned,
// it will always be true that i0 <= i1, is > 0.
//
//-----------------------------------------------------------------------------

bool findIntersectionEndpoints(int a0, int a1, int s, int b0, int b1, int t,
                               int &left, int &right, int &stride)
{
  // We should have both domains moving in positive directions, although
  // we might have to reverse the sign on the strides.
  PAssert(a0 <= a1 && b0 <= b1);
  if (s < 0)  s = -s;
  if (t < 0)  t = -t;

  // find the left endpoint first.  If we cannot, we're done.
  if (!findLeftCommonEndpoint(a0, a1, s, b0, b1, t, left))
    return false;

  // find the least-common-multiple of the strides s and t.  This will
  // be the stride of the resulting domain, if there is one possible.  But
  // first make sure the strides are positive, since we're assuming in this
  // routine that it was called with a0 <= a1, b0 <= b1.
  stride = findLCM(s, t);

  // find the minimum of the right endpoints, and then how far we must
  // move in from this endpoint to get a point in both domains.
  int m = a1;
  if (b1 < a1)
    m = b1;
  right = m - ((m - left) % stride);

  // we were able to find an intersection
  return true;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainCalculus.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
