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

//-----------------------------------------------------------------------------
// Unique non-template definitions.
//-----------------------------------------------------------------------------

// include files

#include "Threads/PoomaMutex.h"
#include "Utilities/Unique.h"

//-----------------------------------------------------------------------------
// The static counter used in Unique
//-----------------------------------------------------------------------------

Unique::Value_t Unique::next_s = 0;


//-----------------------------------------------------------------------------
// The static mutexs used in Unique
//-----------------------------------------------------------------------------

Pooma::Mutex_t Unique::mutex_s;


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Unique.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
