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

// POOMA include files

#include "Connect/Paws/PawsAppPointer.h"
#include "Utilities/PAssert.h"
#include "Pooma/Pooma.h"
#include "Pooma/Configuration.h"

// PAWS C++ API include files

#if POOMA_PAWS
#include "Paws/PawsApplication.h"
#endif


//-----------------------------------------------------------------------------
// static data members of PawsAppPointer
//-----------------------------------------------------------------------------

// The Paws application object.  There is only one of these at a time.  If
// you create a PawsAppPointer and this pointer is null, the application
// will initialize the Paws connection and contact the Paws controller.

PawsApplication *PawsAppPointer::paws_s  = 0;
int              PawsAppPointer::users_s = 0;


//----------------------------------------------------------------------
//
// PawsAppPointer constructor
// The constructor takes a string name; the type is "paws".
//
//----------------------------------------------------------------------

PawsAppPointer::PawsAppPointer(const char *conname, int argc, char *argv[])
  : connected_m(false)
{
  // if we do not have a Paws connection yet, create it now

  if (paws_s == 0)
    {
      PAssert(conname != 0);
#if POOMA_PAWS
      paws_s = new PawsApplication(conname, argc, argv,
				   Pooma::context(), Pooma::contexts());
#endif
    }

  // indicate we're one more object using Paws

  if (connected_m = (paws_s != 0))
    users_s++;
}


//----------------------------------------------------------------------
//
// PawsAppPointer destructor
//
//----------------------------------------------------------------------

PawsAppPointer::~PawsAppPointer()
{
  close();
}


//----------------------------------------------------------------------
//
// Return the Paws connection object.
//
//----------------------------------------------------------------------

PawsApplication &PawsAppPointer::paws() const
{
  PAssert(paws_s != 0);
  PAssert(connected());
  return *paws_s;
}


//----------------------------------------------------------------------
//
// Perform a Paws poll
//
//----------------------------------------------------------------------

void PawsAppPointer::poll()
{
#if POOMA_PAWS
  paws().poll();
#endif
}


//----------------------------------------------------------------------
//
// Wait for a ready signal from the Paws controller
//
//----------------------------------------------------------------------

void PawsAppPointer::ready()
{
#if POOMA_PAWS
  paws().ready();
#endif
}


//----------------------------------------------------------------------
//
// Shut down the connection, and put ourselves in an "unconnected" state.
//
//----------------------------------------------------------------------

void PawsAppPointer::close()
{
  // close paws connection if we're the last one

  if (connected() && --users_s == 0)
    {
#if POOMA_PAWS
      delete paws_s;
#endif
      paws_s = 0;
    }

  connected_m = false;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsAppPointer.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo



