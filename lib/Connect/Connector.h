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

#ifndef POOMA_CONNECT_CONNECTOR_BASE_H
#define POOMA_CONNECT_CONNECTOR_BASE_H

//-----------------------------------------------------------------------------
// Classes:
// ConnectorBase
// Connector<Data,Con> default type
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Connect
 * @brief
 * Base class for all objects that are used to take
 * data from a single item and connect it to some connection.
 *
 * ConnectorBase is a base class for all objects that are used to take
 * data from a single item such as an Array and connect it to some connection
 * such as a file or another program.  Subclasses of ConnectorBase should
 * be specializations of Connector<Data,Con> where 'Data' is the type of
 * data being connected to the connection of type 'Con'.  ConnectorBase
 * defines the abstract interface to all Connectors, including the 'update'
 * and 'notify' methods.  A single Connector represents the connection between
 * just one ConnectionBase object and one data object.
 *
 * Each ConnectorBase has a name and a ConnectionBase base class pointer.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connection.h"
#include "Utilities/PAssert.h"
#include <string>


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/** declaration of the Connector class, that users will need to specialize.
 * Specialized versions of Connector should inherit from ConnectorBase, but
 * this general declaration does not have to.
 */

template<class D, class C>
class Connector { };


/**
 * ConnectionBase and Connection<Tag> are part of the POOMA external data
 * connection API.  The other part of this API are the classes ConnectorBase
 * and Connector<T,Tag>.
 *
 * Connection<Tag> represents the information about how a "connection" is
 * made to some external agency.  Examples: files, Paws-based programs,
 * visualization facilities.  The parameter "Tag" is a simple tag class used
 * to partially specialize Connection.  This same tag is used in Connector.
 *
 * Connector<T,Tag> represents the information about how to map data from
 * a single instance of a type T to some Connection<Tag> object.  For example,
 * if data must be serialized into a buffer and an API routine for the "Tag"
 * connection type called, Connector<T,Tag> will do that, for a given instance
 * of T.
 *
 * There can be several Connector<T,Tag> objects associated with a given
 * Connection<Tag>.  Each object maintains the "channel" for a data object to
 * the connection.  The connection object stores a list of connectors and
 * provides an API for "updating" and "interacting" with the external agency.
 *
 * ConnectorBase is a non-templated base class for all Connector<T,Tag>
 * objects, used to store common information and provide a non-templated
 * base class for heterogenous lists of connectors.  It provides:
 *   -# A pointer to the ConnectionBase class for the connection this
 *      connector is associated with.
 *   -# Provides an API to:
 *      - update the connection channel
 *      - disconnect
 *      .
 *   -# Provides storage for a user-supplied name string for the connector.
 *   -# Provides storage for a transfer mode setting.
 *
 * New connector's are added via a "connect" method in ConnectionBase.
 * You can create your own Connector<T,Tag> instance, and add that in via
 * connect(), or you can call a specialized (and generally templated)
 * connect() method in the Connection<Tag> class.  Most Connection classes
 * provide the specialized connect in order to make it nicer to add in
 * new connectors without having to instantiate a Connector<T,Tag> directly
 * in user code.
 *
 * When you update the connections, you ask each connector to do an update.
 * Updates are done in the order they were connected.  Each connector has
 * a transfer mode, one of the following:
 *  - ConnectionBase::in ..... data is imported from the external agency
 *  - ConnectionBase::out .... data is exported to the external agency
 *  - ConnectionBase::inout .. data is imported and exported to external agency
 *
 * The transfer mode is established when the connector is added.  An update
 * causes new values to be imported or exported based on the current state
 * of the data and the external agency.  For a visualization connection, an
 * update would result in new data being sent to the renderer and hopefully
 * a new image drawn in the visualization program's display windows(s).
 * You can update all connectors attached to a connection, or just one (by
 * providing a pointer to the ConnectorBase of the connection you want to
 * update).
 *
 * When a connector is disconnected, it is told via a call to its virtual
 * "notify" method that a disconnect has occurred.  This lets the
 * connector note that it will no longer be able to update, even though it
 * may not yet have been deleted.
 */

class ConnectorBase
{
public:
  //============================================================
  // ConnectorBase Constructor
  //============================================================

  /// The constructor takes a string for the name, the connection to use,
  /// and the transfer mode.

  ConnectorBase(const char *conname, ConnectionBase &c, int mode)
    : connection_m(&c), name_m(conname), mode_m(mode)
    {
      PAssert(mode == ConnectionBase::in ||
	      mode == ConnectionBase::out ||
	      mode == ConnectionBase::inout);
    }


  //============================================================
  // ConnectorBase Destructor
  //============================================================

  /// When destructed, this disconnects from its ConnectionBase (if any).
  /// By this point, we may no longer have a connection, so make sure
  /// by checking the pointer.

  virtual ~ConnectorBase()
    {
      if (connected())
	connection().disconnect(this);
    }


  //============================================================
  // ConnectorBase operations
  //============================================================

  /// Do special activities to disconnect ourselves from the ConnectionBase.

  virtual void disconnect() = 0;

  /// Update our connection, for example, transfer data or read/write a
  /// file.  This must be provided by derived classes.

  virtual void update() = 0;

  /// Allow for interaction with the connection.  An optional string
  /// can be provided to tell how to do the interaction.

  virtual void interact(const char * = 0) = 0;

  /// Get notified of something by our ConnectionBase object.  If the event
  /// code is 0, it means the ConnectionBase is going away so we should no
  /// longer use it.  Subclasses may need to be able to handle other
  /// event codes, so don't complain if we see one we don't know about.

  virtual void notify(ConnectionBase &c, int event)
    {
      if (event == ConnectionBase::disconnectEvent)
	{
	  // make sure this is the right connection, then unset our pointer

	  PAssert(&c == connection_m);
	  disconnect();
	  connection_m = 0;
	}
    }


  //============================================================
  // ConnectorBase accessors
  //============================================================

  /// Return the name of this connection.

  const std::string &name() const
    {
      return name_m;
    }

  /// Are we connected to something?

  bool connected() const
    {
      return (connection_m != 0);
    }

  /// Return the connection we are using

  ConnectionBase &connection() const
    {
      PAssert(connected());
      return *connection_m;
    }

  /// Return the data transfer mode

  int transferMode() const
    {
      return mode_m;
    }

private:
  /// Our connection

  ConnectionBase *connection_m;

  /// The name of this connector

  std::string name_m;

  /// The transfer mode (in, out or both)

  int mode_m;

  /// The default and copy constructors are made private and undefined
  /// since they should not be used

  ConnectorBase();
  ConnectorBase(const ConnectorBase &);
  ConnectorBase &operator=(const ConnectorBase &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_CONNECTOR_BASE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Connector.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:16 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
