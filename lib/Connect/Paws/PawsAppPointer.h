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

#ifndef POOMA_CONNECT_PAWS_PAWS_APP_POINTER_H
#define POOMA_CONNECT_PAWS_PAWS_APP_POINTER_H

//-----------------------------------------------------------------------------
// Classes:
// PawsAppPointer
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// PawsAppPointer holds a pointer to a PawsApplication object, and makes
// sure it is a singleton object.  It maintains a count on objects that
// have requested to use it, and will delete the PawsApplication object
// when the count goes back to zero.  To use it, create an instance in
// your own class.  Use the "paws()" method to get
// a reference to the PawsApplication object that it stores.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

class PawsApplication;


//-----------------------------------------------------------------------------
//
// Full Description:
//
// PawsAppPointer is a wrapper around an instance of a PawsApplication
// class.  PawsApplication must be a "singleton" type instance, and this
// class uses a static pointer to a PawsApplication instance and a counter
// to tell when the instance must be created or deleted.  If more than
// one PawsAppPointer is created, each one will refer to the same
// PawsApplication object, until the last PawsAppPointer goes away.
//
// PawsAppPointer is created with the arguments needed to create
// PawsApplication, namely a string with a name for the application, and
// the command-line arguments.  It will instantiate a new PawsApplication
// if none has yet been created.  The destructor will call the "close()"
// method, that will delete the PawsApplication if it is the last object
// around of this type.
//
// PawsAppPointer also provides methods "poll()", "ready()", and others
// that are just deferred to the PawsApplication.  The implementation for
// these are put in the .cmpl.cpp file so that we can encapsulate the
// details of the Paws code in a separate file.
//
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

class PawsAppPointer
{
public:
  //============================================================
  // Constructor
  //============================================================

  // This class is constructed with a string name, and the argc,argv
  // values used to initialize PAWS.

  PawsAppPointer(const char *conname, int argc, char *argv[]);


  //============================================================
  // Destructor
  //============================================================

  // When destroyed, disconnect from Paws.

  virtual ~PawsAppPointer();


  //============================================================
  // Accessors
  //============================================================

  // Return whether we are connected.

  bool connected() const
  {
    return connected_m;
  }

  // Return the Paws connection object.

  PawsApplication &paws() const;


  //============================================================
  // Operations
  //============================================================

  // Perform a Paws poll.

  void poll();

  // Wait for a ready signal from the Paws controller.

  void ready();

  // Shut down the connection, and put ourselves in an "unconnected" state.

  void close();

private:
  // The Paws connection object.

  static PawsApplication *paws_s;

  // The number of instances using Paws

  static int users_s;

  // Are we connected, and using the app pointer?

  bool connected_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  PawsAppPointer();
  PawsAppPointer(const PawsAppPointer &);
  PawsAppPointer &operator=(const PawsAppPointer &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_PAWS_PAWS_APP_POINTER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsAppPointer.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
