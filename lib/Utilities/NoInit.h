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
//   Pooma::NoInit
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * Pooma::NoInit tag class.
 */

#ifndef POOMA_UTILITIES_NOINIT_H
#define POOMA_UTILITIES_NOINIT_H

#include "Pooma/Configuration.h"
#include "Utilities/PurifyConstructors.h"

//////////////////////////////////////////////////////////////////////

namespace Pooma {

/**
 * Tag class used to signal constructors that certain initializations are
 * to be skipped. 
 */

class NoInit
{
public:
  POOMA_PURIFY_CONSTRUCTORS(NoInit)
};

} // namespace Pooma

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_NOINIT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NoInit.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
