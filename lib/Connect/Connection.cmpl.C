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
// Classes:
// ConnectionBase
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connection.h"
#include "Connect/Connector.h"
#include "Utilities/PAssert.h"
#include <string>
#include <vector>


//----------------------------------------------------------------------
//
// ConnectionBase Constructor
// The constructor takes the strings for the name and type.
//
//----------------------------------------------------------------------

ConnectionBase::ConnectionBase(const char *conname, const char *contype)
  : name_m(conname), type_m(contype)
{
}


//----------------------------------------------------------------------
//
// ConnectionBase Destructor
//
//----------------------------------------------------------------------

ConnectionBase::~ConnectionBase()
{
  // At this point, the subclasses should have removed all
  // connections

  PAssert(size() == 0);
}


//----------------------------------------------------------------------
//
// Return whether we are connected properly right now.  By default,
// we are.
//
//----------------------------------------------------------------------

bool ConnectionBase::connected() const
{
  return true;
}


//----------------------------------------------------------------------
//
// Connect a new Connector.  The Connector maps from a single data object
// (such as a Field) to this ConnectionBase.  The base class pointer is used
// to store the data, since Connector is a templated class.  The user must
// 'new' thie Connector object and give the allocated instance to this
// class, which will delete the instance when the connection closes.
// Return the ConnectorBase pointer just connected.
//
//----------------------------------------------------------------------

ConnectorBase *ConnectionBase::connect(ConnectorBase *cb)
{
  PAssert(connected());

  // make sure we do not already have this connector; if we do, just return

  for (iterator a = begin(); a != end(); ++a)
    if (cb == *a)
      return cb;

  // this is a new one; add it to our list

  PAssert(cb != 0);
  connectors_m.push_back(cb);
  return cb;
}


//----------------------------------------------------------------------
//
// Disconnect an existing Connector, but do not delete it.  This removes
// the connector from our list, and informs it it should not use the
// ConnectionBase anymore.
// ConnectorBase must implement a 'notify(int)' virtual function.
// Return the ConnectorBase pointer just disconnected.
//
//----------------------------------------------------------------------

ConnectorBase *ConnectionBase::disconnect(ConnectorBase *cb)
{
  PAssert(connected());

  // find the connection; if we do not have it, it is an error

  for (iterator a = begin(); a != end(); ++a)
    {
      if (cb == *a)
	{
	  // we have the item; notify it and remove it from our list

	  (*a)->notify(*this, ConnectionBase::disconnectEvent);
	  connectors_m.erase(a);

	  return cb;
	}
    }

  // if we're here, its an error

  PAssert(false);
  return 0;
}


//----------------------------------------------------------------------
//
// Update the connections - tell the connections to update, or just an
// individual one.
//
//----------------------------------------------------------------------

void ConnectionBase::update(ConnectorBase *cb)
{
  PAssert(connected());

  // loop through all the connectors; update them all, or just the
  // one mentioned.

  for (iterator a = begin(); a != end(); ++a)
    if (cb == 0 || cb == *a)
      (*a)->update();
}


//----------------------------------------------------------------------
//
// Allow for interaction with the connection.  An optional string
// can be provided to tell how to do the interaction.
//
//----------------------------------------------------------------------

void ConnectionBase::interact(const char *)
{
  // by default, this does nothing
}


//----------------------------------------------------------------------
//
// Disconnect any remaining connectors, and delete them.
//
//----------------------------------------------------------------------

void ConnectionBase::disconnectConnectors()
{
  while (size() > 0)
    {
      ConnectorBase *cb = *(begin());
      disconnect(cb);
      delete cb;
    }
}


//----------------------------------------------------------------------
//
// Completely close the connection and remove all connectors.  This is
// the equivalent of running the destructor, using this routine you can
// control when the connection is closed.  By default it just removes
// all connectors.
//
//----------------------------------------------------------------------

void ConnectionBase::close()
{
  disconnectConnectors();
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Connection.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:16 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
