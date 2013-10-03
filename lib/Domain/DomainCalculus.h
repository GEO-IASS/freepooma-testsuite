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

#ifndef POOMA_DOMAIN_DOMAIN_CALCULUS_H
#define POOMA_DOMAIN_DOMAIN_CALCULUS_H

//-----------------------------------------------------------------------------
// Function:
//   bool findLeftCommonEndpoint(a0,a1,s,b0,b1,t,endpoint)
//   bool findIntersectionEndpoints(a0,a1,s,b0,b1,t,i0,i1,is)
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file 
 * @ingroup Domain
 * @brief
 * This file defines the prototypes for routines used in the domain
 * calculus computations.  These routines are not part of the main user
 * API, they are mainly for the domain calculus implementation.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/// findLeftCommonEndpoint is used by the domain calculus routines to
/// find the leftmost common endpoint for two domains defined by
/// [a0:a1:s] and [b0:b1:t].  If an endpoint is possible, this returns
/// true and the endpoint in the final argument.  If one is not possible,
/// due to incompatible striding, this returns false and leaves the final
/// argument unchanged.

extern
bool findLeftCommonEndpoint(int a0, int a1, int s, int b0, int b1, int t,
                            int &endpoint);


/// findIntersectionEndpoints is used by the domain calculus routines to
/// find the endpoints and stride of an intersection domain given
/// two other strided domains defined by [a0:a1:s] and [b0:b1:t].  If
/// this is possible, this return true and the min, max, stride of the
/// resulting domain [i0,i1,is] in the final arguments.  If no intersection
/// is possible, due to incompatible striding, this returns false and
/// leaves the final arguments unchanged.  Note that if a domain is returned,
/// it will always be true that i0 <= i1, is > 0.

extern
bool findIntersectionEndpoints(int a0, int a1, int s, int b0, int b1, int t,
                               int &i0, int &i1, int &is);


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_CALCULUS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainCalculus.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
