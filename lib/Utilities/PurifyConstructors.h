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
// macro POOMA_PURIFY_CONSTRUCTORS()
//-----------------------------------------------------------------------------

#ifndef POOMA_UTILITIES_PURIFYCONSTRUCTORS_H
#define POOMA_UTILITIES_PURIFYCONSTRUCTORS_H

#include "Pooma/Configuration.h"

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// POOMA_PURIFY_CONSTRUCTORS:
//
// POOMA_PURIFY_CONSTRUCTORS is a macro that can be used to generate an empty
// constructor, copy constructor and destructor for tag classes.  Because
// empty classes need to be implemented with some storage, copying empty
// classes can lead to accessing uninitialized memory, which leads purify
// to report errors.
//
// This macro generates no code if POOMA_PURIFY is not enabled, since some
// compilers (egcs) can generate more optimized code if the compiler generates
// the default constructors.
//
//-----------------------------------------------------------------------------

#if POOMA_PURIFY

#define POOMA_PURIFY_CONSTRUCTORS(CLASS)  \
  CLASS() { }   \
  ~CLASS() { }   \
  CLASS(const CLASS &) { }

#else

#define POOMA_PURIFY_CONSTRUCTORS(CLASS)

#endif

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_PURIFYCONSTRUCTORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PurifyConstructors.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
