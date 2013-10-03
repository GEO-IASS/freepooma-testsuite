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

#ifndef POOMA_CONNECT_LUX_LUX_CONNECTION_H
#define POOMA_CONNECT_LUX_LUX_CONNECTION_H

//-----------------------------------------------------------------------------
// Classes:
// Connection<Lux>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connection<Lux> is a ConnectionBase subclass that manages a connection
// to a another program, using the Lux package for run-time visualization.
// When a Connection<Lux> object is created, it will initialize the Lux
// library and display environment.  This is a specialization
// of the Connection class, using the simple tag class "Lux" to specialize.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connection.h"
#include "Connect/ConnectPair.h"
#include "Connect/Connector.h"
#include "Connect/Lux/LuxAppPointer.h"
#include "Utilities/PAssert.h"


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

// A simple tag class for Lux; we just need a forward declaration, since
// we never instantiate anything of this type.

class Lux;


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connection<Lux> is a specialization of Connection<T> that maintains
// a connection for this application to the Lux runtime interactive
// visualization package.  This tool is used to display field and
// particle data in a visualization window.  Lux is a toolkit from the
// visualization team in the Advanced Computing Laboratory.
// 
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

template<>
class Connection<Lux> : public ConnectionBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef Lux                          ConnectionTag_t;
  typedef Connection<Lux>              Connection_t;


  //============================================================
  // Connection<Lux> Constructor
  //============================================================

  // The constructor takes a string name; the type is "lux".  The name
  // may be used someday for a display window title.

  Connection(const char *conname)
    : ConnectionBase(conname, "lux"), lux_m(conname)
  {
  }


  //============================================================
  // Connection<Lux> Destructor
  //============================================================

  // Connection<Lux> destructor.  Remove all connectors and close Lux window.

  virtual ~Connection()
  {
    close();
  }


  //============================================================
  // Connection<Lux> accessors
  //============================================================

  // Return the main Lux tool object

  inline LuxAppPointer &lux()
  {
    return lux_m;
  }

  // Return whether we are connected properly right now.

  virtual bool connected() const
  {
    return lux_m.connected();
  }


  //============================================================
  // Connection<Lux> operations
  //============================================================

  // NOTE: The basic operations are all inherited from the base class.
  // The extra routines here are specific to a Connection<Lux>.

  // Make sure everything is ready after initialization and setup.  Here,
  // this means to do an update and then an interact.

  inline void ready()
  {
    update();
    interact();
  }

  // Perform a Lux poll for user interaction.

  inline void poll()
  {
    lux().poll();
  }

  // Completely close the connection and remove all connectors.  This is
  // the equivalent of running the destructor, using this routine you can
  // control when the connection is closed.

  virtual void close()
  {
    disconnectConnectors();
    lux().close();
  }

  // Connect in an object of type T.  This will create a new Connector
  // instance and add it in to our list of connected items.

  template<class T>
  ConnectorBase *connect(const char *cname, const T &obj,
			 int mode = ConnectionBase::out)
  {
    ConnectorBase *cb = new Connector<T,Lux>(cname, obj, *this, mode);
    return ConnectionBase::connect(cb);
  }

  // Connect a pair of objects of types T1 and T2.

  template<class T1, class T2>
  ConnectorBase *connect(const char *cname, const T1 &obj1, const T2 &obj2,
			 int mode)
  {
    typedef ConnectPair<T1,T2> Pair_t;
    ConnectorBase *cb = new Connector<Pair_t,Lux>(cname, Pair_t(obj1, obj2),
						  *this, mode);
    return ConnectionBase::connect(cb);
  }

  // Allow for user interaction with the display.  This may result in
  // handing control to the Lux routines until the user is finished with
  // interaction.  This is the same as "poll".

  virtual void interact(const char * = 0)
  {
    poll();
  }

private:
  // The visualization tool object, managed by LuxAppPointer.

  LuxAppPointer lux_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Connection();
  Connection(const Connection_t &);
  Connection_t &operator=(const Connection_t &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_LUX_LUX_CONNECTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LuxConnection.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
