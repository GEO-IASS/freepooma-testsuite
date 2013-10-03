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

/** @file
 * @ingroup Field
 * @brief
 * functor for getting the nth patch out of a field.
 */

#ifndef POOMA_FIELD_FIELDENGINEPATCH_H
#define POOMA_FIELD_FIELDENGINEPATCH_H

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

#include "Engine/EnginePatch.h"

/**
 * Fields know about two domains, a physical and total domain.  To get a patch
 * view of a field, we need to pass the physical domain along, since the engine
 * will typically extend beyond the domain when there are guard cells.
 */

template<int Dim>
struct FieldEnginePatch
{
  FieldEnginePatch(int patch, Interval<Dim> domain)
    : patch_m(patch), domain_m(domain)
  { }

  int patch_m;
  Interval<Dim> domain_m;
};

/////////////////////////////////////////////////////////////////////

#endif     // POOMA_FIELD_FIELDENGINEPATCH_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldEnginePatch.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:45 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
