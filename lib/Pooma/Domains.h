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

#ifndef POOMA_POOMA_DOMAINS_H
#define POOMA_POOMA_DOMAINS_H

/** @file
 * @ingroup Pooma
 * @brief
 * A one-stop-shopping header file that sets up everything one needs to use
 * all Pooma II domain objects and functions.
 */

// Include files for the basic user domain types

#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Grid.h"
#include "Domain/Region.h"

// Include files for the wildcard domain types

#include "Domain/AllDomain.h"
#include "Domain/LeftDomain.h"
#include "Domain/RightDomain.h"

// Include files for the domain calculus routines

#include "Domain/Split.h"
#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/Touches.h"
#include "Domain/EquivSubset.h"

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Domains.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
