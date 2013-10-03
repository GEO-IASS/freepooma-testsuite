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

#ifndef POOMA_CONNECT_PAWS_PAWS_CONNECTION_H
#define POOMA_CONNECT_PAWS_PAWS_CONNECTION_H

//-----------------------------------------------------------------------------
// Classes:
// Connection<Paws>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connection<Paws> is a ConnectionBase subclass that manages a connection
// to a another program, using the Paws library for inter-app communication.
// When a Connection<Paws> object is created, it will initialize the Paws
// library and connect to the Paws controller.  This is a specialization
// of the Connection class, using the simple tag class "Paws" to specialize.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connection.h"
#include "Connect/Paws/PawsAppPointer.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

class Paws;


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connection<Paws> is a specialization of Connection<T> that maintains
// a connection for this application to a Paws controller and other 
// Paws apps.  It stores (via PawsAppPointer) a PawsApplication instance
// used to establish and access this connection.
// 
// When created, Connection<Paws> will create a PawsApplication instance
// (or use an existing one if available), and register this app with the
// provided string name with the Paws controller.  Command-line parameters
// must be provided that include special Paws flags to indicate the
// location of the Paws controller and other settings.
//
// Once registered with the controller, connectors can be created that will
// register individual data objects for sharing with other Paws apps.
//
// Connection<Paws> includes some extra methods not in the base class:
//   ready() indicates this application will wait for a "ready" signal
//           from the controller.  ready() will not return until that signal
//           has been given (or if the application is not connected at all).
//   poll() gives Paws a chance to check for Paws-specific events and
//          process them.  This is generally not something the user ever
//          has to call.
//
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

template<>
class Connection<Paws> : public ConnectionBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef Paws                         ConnectionTag_t;
  typedef Connection<Paws>             Connection_t;


  //============================================================
  // Connection<Paws> Constructor
  //============================================================

  // The constructor takes a string name; the type is "paws".  It also
  // takes the argc, argv of the program, for use in initializing PAWS.

  Connection(const char *conname, int argc, char *argv[])
    : ConnectionBase(conname, "paws"), paws_m(conname, argc, argv)
  {
  }


  //============================================================
  // Connection<Paws> Destructor
  //============================================================

  // Connection<Paws> destructor.  Remove all connectors and close Paws link.

  virtual ~Connection()
  {
    close();
  }


  //============================================================
  // Connection<Paws> accessors
  //============================================================

  // Return whether we are connected properly right now.

  virtual bool connected() const
  {
    return paws_m.connected();
  }

  // Return the Paws connection object.

  inline PawsApplication &paws() const
  {
    return paws_m.paws();
  }


  //============================================================
  // Connection<Paws> operations
  //============================================================

  // NOTE: The basic operations are all inherited from the base class.
  // The extra routines here are specific to a Connection<Paws>.

  // Perform a Paws poll.

  inline void poll()
  {
    paws_m.poll();
  }

  // Wait for a ready signal from the Paws controller.

  inline void ready()
  {
    paws_m.ready();
  }

  // Completely close the connection and remove all connectors.  This is
  // the equivalent of running the destructor, using this routine you can
  // control when the connection is closed.  This removes all connectors
  // and shuts down link to Paws controller.

  virtual void close()
  {
    disconnectConnectors();
    paws_m.close();
  }

  // Connect in an object of type T.  This will create a new Connector
  // instance and add it in to our list of connected items.

  template<class T>
  ConnectorBase *connect(const char *cname, const T &obj, int mode)
  {
    ConnectorBase *cb = new Connector<T,Paws>(cname, obj, *this, mode);
    return ConnectionBase::connect(cb);
  }

  // Connect in an object of type T, by providing a non-const ref.  The
  // connection will use this reference.  This should be used for
  // types that cannot be specified with a const-ref, like regular
  // scalars.

  template<class T>
  ConnectorBase *connectScalar(const char *cname, T &obj, int mode)
  {
    ConnectorBase *cb = new Connector<T,Paws>(cname, obj, *this, mode);
    return ConnectionBase::connect(cb);
  }

private:
  // The Paws connection object, managed by PawsAppPointer.

  PawsAppPointer paws_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Connection();
  Connection(const Connection_t &);
  Connection_t &operator=(const Connection_t &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_PAWS_PAWS_CONNECTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsConnection.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
