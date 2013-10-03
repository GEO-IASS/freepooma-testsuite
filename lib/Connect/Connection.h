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

#ifndef POOMA_CONNECT_CONNECTION_H
#define POOMA_CONNECT_CONNECTION_H

//-----------------------------------------------------------------------------
// Classes:
// ConnectionBase
// Connection
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Connect
 * @brief base class for specialized classes that manage a
 * connection or "channel" to some external agency.
 *
 * ConnectionBase is a base class for specialized classes that manage a
 * connection or "channel" to some external agency, such as a visualization
 * window, a file, another program, etc.  It allows other "Connector"
 * objects to connect it to some data object such as an Array.
 * For each ConnectionBase, there can be several Connector's (stored
 * via a ConnectorBase pointer) registered with it.  When a ConnectionBase is
 * deleted or closed, it informs all registered observers that it is being
 * deleted.
 *
 * Each ConnectionBase also has a string name and string connection type.
 *
 * ConnectionBase provides a set of virtual functions that can be used to
 * provide the standard behavior, with the option to override them to do
 * specialized operations.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <string>
#include <vector>


//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

class ConnectorBase;


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * declaration of the Connection class, that users will need to specialize.
 * Specialized versions of Connection should inherit from ConnectionBase, but
 * this general declaration does not have to.  The type T should be a simple
 * tag used to specify the general connection type.
 */

template<class T>
class Connection { };


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
 * ConnectionBase is a non-templated base class that provides some general
 * facilities for Connection<Tag> objects.  It does the following:
 *   -# Stores a (heterogeneous) list of connectors, by storing a list of
 *      pointers to the base class for all connectors, ConnectorBase.
 *   -# Provides an API to:
 *      - iterate over connectors
 *      - update all connectors
 *      - interact with the connection
 *      - disconnect connectors
 *      - close the connection
 *      .
 *   -# Provides storage for a user-supplied name string for the connection.
 *   -# Provides storage for a user-supplied type string for the connection.
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
 * You can disconnect all or specific connectors via the disconnect() method.
 * This will remove that connector's data from the external agency, and
 * future update() calls will not use that connector.
 *
 * To completely shut down the connection, either delete the connector
 * instance or call the close() method.  It will disconnect all connectors,
 * and shut down the connection (for a visualizer, this would close the
 * windows; for a Paws connection, it will close the connection to the Paws
 * controller).
 *
 * Some external agencies may need to occasionally take control, to do some
 * event processing or something.  The interact() method allows this to
 * happen.  It may not be necessary for some connections.
 */

class ConnectionBase
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  typedef std::vector<ConnectorBase *> List_t;
  typedef List_t::size_type            Size_t;
  typedef List_t::iterator             iterator;
  typedef List_t::const_iterator       const_iterator;

  /// Enumeration of event types used to notify connectors.

  enum { disconnectEvent = 0, connectEvent = 1 };

  /// Enumeration of data transfer directions.

  enum { in = 0, out = 1, inout = 2 };


  //============================================================
  // ConnectionBase Constructor
  //============================================================

  /// The constructor takes the strings for the name and type.

  ConnectionBase(const char *conname, const char *contype);


  //============================================================
  // ConnectionBase Destructor
  //============================================================

  /// When destructed, a ConnectionBase disconnects all Connectors and then
  /// deletes them.

  virtual ~ConnectionBase();


  //============================================================
  // ConnectionBase operations
  //============================================================

  /// Establish a connection from the first argument to the second.  The
  /// first should be
  /// Connect a new Connector.  The Connector maps from a single data object
  /// (such as a Field) to this ConnectionBase.  The base class pointer is used
  /// to store the data, since Connector is a templated class.  The user must
  /// 'new' thie Connector object and give the allocated instance to this
  /// class, which will delete the instance when the connection closes.
  /// Return the ConnectorBase pointer just connected.

  virtual ConnectorBase *connect(ConnectorBase *);

  /// Disconnect an existing Connector, but do not delete it.  This removes
  /// the connector from our list, and informs it it should not use the
  /// ConnectionBase anymore.
  /// ConnectorBase must implement a 'notify(int)' virtual function.
  /// Return the ConnectorBase pointer just disconnected.

  virtual ConnectorBase *disconnect(ConnectorBase *);

  /// Update the connections - tell the connections to update, or just an
  /// individual one.

  virtual void update(ConnectorBase * = 0);

  /// Allow for interaction with the connection.  An optional string
  /// can be provided to tell how to do the interaction.

  virtual void interact(const char * = 0);

  /// Completely close the connection and remove all connectors.  This is
  /// the equivalent of running the destructor, using this routine you can
  /// control when the connection is closed.  By default it just removes
  /// all connectors.

  virtual void close();


  //============================================================
  // ConnectionBase accessors
  //============================================================

  /// Return whether we are connected properly right now.

  virtual bool connected() const;

  /// Return the name of this connection.

  const std::string &name() const
    {
      return name_m;
    }

  /// Return the type string of this connection.

  const std::string &type() const
    {
      return type_m;
    }

  /// @name Begin/end iterators for the list of connectors
  //@{

  iterator begin()
    {
      return connectors_m.begin();
    }

  iterator end()
    {
      return connectors_m.end();
    }
  //@}

  /// @name Const begin/end iterators for the list of connectors
  //@{

  const_iterator begin() const
    {
      return connectors_m.begin();
    }

  const_iterator end() const
    {
      return connectors_m.end();
    }
  //@}

  /// Return the number of registered connectors

  Size_t size() const
    {
      return connectors_m.size();
    }

protected:
  /// Disconnect any remaining connectors.

  void disconnectConnectors();

private:
  /// The list of connectors

  List_t connectors_m;

  /// The name and type of this connection

  std::string name_m;
  std::string type_m;

  /// The default and copy constructors are made private and undefined
  /// since they should not be used

  ConnectionBase();
  ConnectionBase(const ConnectionBase &);
  ConnectionBase &operator=(const ConnectionBase &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_CONNECTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Connection.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:16 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
