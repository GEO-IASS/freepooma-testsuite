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

#ifndef POOMA_CONNECT_LUX_LUX_CONNECTOR_ARRAY_H
#define POOMA_CONNECT_LUX_LUX_CONNECTOR_ARRAY_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<Array<>, Lux>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<Array<>, Lux> is a specialization
// of the general Connector class used to connect a POOMA II Array to 
// the Lux run-time visualization package.  It can only be used for output
// of data.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connector.h"
#include "Connect/Lux/LuxConnection.h"
#include "Connect/Lux/LuxAppPointer.h"
#include "Pooma/Arrays.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<Array<>, Lux> is a specialization of Connector<T,Tag> that can
// map data from a POOMA II Array of type T to the Lux run-time visualization
// package.  The Array can be of any time and structure, but will be
// limited in the maximum dimension (3).
//
// Connector<Array, Lux> should be created with a name for the data object and
// a non-const reference to the Array to share with another program.  This
// will store a shallow copy of that Array, and use its value to update the
// connection when update() is called.  Since this is a copy, the provided
// Array could be a view of some other Array so that you can share portions
// of a larger Array with another application.
//
// This can only work with a Connection<Lux> connection object.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag>
class Connector<Array<Dim, T, EngineTag>, Lux> : public ConnectorBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef Array<Dim, T, EngineTag> Array_t;
  typedef T                        Scalar_t;
  typedef ReadFieldTool            LuxData_t;
  typedef Connection<Lux>          Connection_t;
  typedef Connector<Array_t,Lux>   Connector_t;
  typedef Lux                      ConnectionTag_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout).

  Connector(const char *conname, const Array_t &a, Connection_t &c,
	    int mode = ConnectionBase::out)
    : ConnectorBase(conname, c, mode), array_m(a), data_m(0)
  {
    // This will only work for output connections

    PAssert(mode == ConnectionBase::out);

    // Compute the sizes of the domain and geometry.

    findsize();

    // Add the new array to the visualization tool, and save a reference
    // to this connection object.

    data_m = luxConnection().lux().createArray(name());
  }


  //============================================================
  // Connector Destructor
  //============================================================

  // When destroyed, disconnect from Lux.  This will eventually call
  // disconnect, where the data object will get deleted.

  virtual ~Connector()
  {
    if (connected())
      connection().disconnect(this);
  }


  //============================================================
  // Connector accessors
  //============================================================

  // Return the connection, cast as a Connection<Lux>

  Connection_t &luxConnection() const
  {
    PAssert(connected());
    return dynamic_cast<Connection_t &>(connection());
  }

  // Return a reference to the array itself

  Array_t &array() const
  {
    return array_m;
  }

  //============================================================
  // Connector operations
  //============================================================

  // Retarget this connector to a new data object.  For some items, the
  // data may be of different size, for others it will be the same size
  // always.

  void resize(const Array_t &newarray)
  {
    // For an array, re-initialize the internal array.

    array_m.initialize(newarray);

    // Compute the sizes of the domain and geometry.

    findsize();
  }


  //============================================================
  // ConnectorBase operations
  //============================================================

  // Do special activities to disconnect ourselves from the Connection<Lux>.

  virtual void disconnect()
  {
    // Disconnect by having lux delete the item.

    if (data_m != 0)
      luxConnection().lux().destroyArray(data_m, name());
    data_m = 0;
  }

  // Update our connection.  For Lux, this results in a data transfer
  // operation, either send or receive, based on the connection method.

  virtual void update()
  {
    if (connected() && data_m != 0)
      {
	// Get the data type

	int datatype = LuxDataType<T>::datatype;

	// Create storage for this data for each element

	float data[LuxDataType<T>::dimensions];

	// Prepare for data

	luxConnection().lux().beginArray(data_m, datatype, size_m,
					 spacing_m, origin_m);

	// For all the elements in the array, get them and put them
	// in the storage.

	typedef typename Array_t::Domain_t Domain_t;
	typedef typename Domain_t::iterator iterator_t;
	iterator_t domiter = array_m.domain().begin();
	for (int indx=0; indx < total_m; ++indx, ++domiter)
	  {
	    // Get the next value of the array, and put it in 'data'

	    LuxDataType<T>::copy(array_m(*domiter), data);

	    // Put this data in the Lux storage

	    luxConnection().lux().insertArray(data_m, datatype, indx, data);
	  }

	// Finish up with data and do update

	luxConnection().lux().endArray(data_m, name());
      }
  }

  // Allow for interaction with the connection.  An optional string
  // can be provided to tell how to do the interaction.  Here, does nothing
  // except poll.

  virtual void interact(const char *s = 0)
  {
    if (connected())
      luxConnection().interact();
  }

protected:
  //============================================================
  // Protected data storing mesh information
  //============================================================

  int total_m;
  int size_m[3];
  float origin_m[3];
  float spacing_m[3];

private:
  //============================================================
  // Private data
  //============================================================

  // The array we're connecting

  Array_t array_m;

  // The Lux connected-data object we're using

  LuxData_t *data_m;

  //============================================================
  // Private methods
  //============================================================

  // Find the size, origin, and spacing of this array and save them.

  void findsize()
  {
    total_m = 0;
    for (int d=0; d < 3; ++d)
      {
	if (d < Dim)
	  {
	    size_m[d] = array_m.domain()[d].length();
	    origin_m[d] = array_m.domain()[d].first();
	    spacing_m[d] = 1.0;
	  }
	else
	  {
	    size_m[d] = 1;
	    origin_m[d] = 0.0;
	    spacing_m[d] = 1.0;
	  }

	if (d > 0)
	  total_m *= size_m[d];
	else
	  total_m = size_m[d];
      }
  }

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Connector();
  Connector(const Connector_t &);
  const Connector_t &operator=(const Connector_t &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_LUX_LUX_CONNECTOR_ARRAY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LuxConnector.Array.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
