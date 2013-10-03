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

#ifndef POOMA_CONNECT_LUX_LUX_CONNECTOR_PARTICLES_H
#define POOMA_CONNECT_LUX_LUX_CONNECTOR_PARTICLES_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<ConnectPair<DynamicArray, DynamicArray>, Lux>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<ConnectPair<DynamicArray, DynamicArray>, Lux> is a specialization
// of the general Connector class used to connect a POOMA II DynamicArray
// pair that both use the same layout with the Lux run-time visualization
// package.  The first item should be a position attribute, and the second
// should be the associated data attribute.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connector.h"
#include "Connect/ConnectPair.h"
#include "Connect/Lux/LuxConnection.h"
#include "Connect/Lux/LuxAppPointer.h"
#include "Pooma/DynamicArrays.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<..., Lux> is a specialization of Connector<T,Tag> that can
// map data from two POOMA II DynamicArray's of type T1 and T2 to the Lux
// run-time visualization package.  The two arrays are a position and data
// attribute from a Particles object, and are plotted as particle items
// in Lux.  Both attributes must use the same layout, since they are
// used together.
//
// Connector<..., Lux> should be created with a name for the data object and
// a non-const reference to the DA's to share with another program.  This
// will store a shallow copy of those arrays, and use the values to update the
// connection when update() is called.  Since this is a copy, the provided
// Array could be a view of some other Array so that you can share portions
// of a larger Array with another application.
//
// This can only work with a Connection<Lux> connection object.
//
//-----------------------------------------------------------------------------

template<class T1, class T2, class E>
class Connector<ConnectPair<DynamicArray<T1,E>, DynamicArray<T2,E> >, Lux>
  : public ConnectorBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef DynamicArray<T1, E>                   PosAttrib_t;
  typedef DynamicArray<T2, E>                   ValAttrib_t;
  typedef ConnectPair<PosAttrib_t, ValAttrib_t> Pair_t;
  typedef ReadParticleTool                      LuxData_t;
  typedef Connection<Lux>                       Connection_t;
  typedef Connector<Pair_t, Lux>                Connector_t;
  typedef Lux                                   ConnectionTag_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout).

  Connector(const char *conname, const Pair_t &a, Connection_t &c,
	    int mode = ConnectionBase::out)
    : ConnectorBase(conname, c, mode),
      pos_m(a.first()), val_m(a.second()), data_m(0)
  {
    // This will only work for output connections

    PAssert(mode == ConnectionBase::out);

    // Add the new array to the visualization tool, and save a reference
    // to this connection object.

    data_m = luxConnection().lux().createParticles(name());
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

  // Return a reference to the position array itself

  PosAttrib_t &posAttrib() const
  {
    return pos_m;
  }

  // Return a reference to the value array itself

  ValAttrib_t &valAttrib() const
  {
    return val_m;
  }


  //============================================================
  // Connector operations
  //============================================================

  // Retarget this connector to a new data object.  For some items, the
  // data may be of different size, for others it will be the same size
  // always.

  void resize(const Pair_t &newpair)
  {
    // For an array, re-initialize the internal array.  The copy array
    // will be reset later.

    pos_m.initialize(newpair.first());
    val_m.initialize(newpair.second());
  }


  //============================================================
  // ConnectorBase operations
  //============================================================

  // Do special activities to disconnect ourselves from the Connection<Lux>.

  virtual void disconnect()
  {
    // Disconnect by having lux delete the item.

    if (data_m != 0)
      luxConnection().lux().destroyParticles(data_m, name());
    data_m = 0;
  }

  // Update our connection.  For Lux, this results in a data transfer
  // operation, either send or receive, based on the connection method.

  virtual void update()
  {
    if (connected() && data_m != 0)
      {
	// Get the data type for the value attribute.  Position must
	// be a vector.

	int datatype = LuxDataType<T2>::datatype;
	PAssert(LuxDataType<T1>::datatype == LuxAppPointer::vector);

	// Create storage for this data for each element

	float pos[3];
	float data[LuxDataType<T2>::dimensions];

	// Get the size of the attributes, make sure they're equal

	int totsize = pos_m.domain().size();
	PAssert(totsize == val_m.domain().size());

	// Prepare for data

	luxConnection().lux().beginParticles(data_m, datatype, totsize);

	// For all the elements in the array, get them and put them
	// in the storage.

	for (int indx=0; indx < totsize; ++indx)
	  {
	    // Get the next value of the array, and put it in 'data'

	    LuxDataType<T1>::copy(pos_m(indx), pos);
	    LuxDataType<T2>::copy(val_m(indx), data);

	    // Put this data in the Lux storage

	    luxConnection().lux().insertParticles(data_m, datatype, indx,
						  pos, data, indx);
	  }

	// Finish up with data and do update

	luxConnection().lux().endParticles(data_m, name());
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

private:
  //============================================================
  // Private data
  //============================================================

  // The position attribute

  PosAttrib_t pos_m;

  // The value attribute

  ValAttrib_t val_m;

  // The Lux connected-data object we're using

  LuxData_t *data_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Connector();
  Connector(const Connector_t &);
  const Connector_t &operator=(const Connector_t &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_LUX_LUX_CONNECTOR_PARTICLES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LuxConnector.Particles.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
