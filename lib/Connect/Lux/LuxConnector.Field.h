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

#ifndef POOMA_CONNECT_LUX_LUX_CONNECTOR_FIELD_H
#define POOMA_CONNECT_LUX_LUX_CONNECTOR_FIELD_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<Field<>, Lux>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<Field<>, Lux> is a specialization
// of the general Connector class used to connect a POOMA II Field to 
// the Lux run-time visualization package.  It can only be used for output
// of data.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connector.h"
#include "Connect/Lux/LuxConnector.Array.h"
#include "Pooma/Fields.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<Field<>, Lux> is a specialization of Connector<T,Tag> that can
// map data from a POOMA II Field of type T to the Lux run-time visualization
// package.  The Field can be of any time and structure, but will be
// limited in the maximum dimension (3).
//
// Since POOMA II works with one context right now, the easiest way to
// implement this is to make a copy of the data into a single Brick engine
// Field, and use that to give the data to Lux.
//
// Connector<Field, Lux> should be created with a name for the data object and
// a non-const reference to the Field to share with another program.  This
// will store a shallow copy of that Field, and use its value to update the
// connection when update() is called.  Since this is a copy, the provided
// Field could be a view of some other Field so that you can share portions
// of a larger Field with another application.
//
// This can only work with a Connection<Lux> connection object.
//
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag>
class Connector<Field<Mesh, T, EngineTag>, Lux>
  : public Connector<Array<FieldEngine<Mesh, T, EngineTag>::dimensions,
                           T, EngineTag>,
                     Lux>
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef Field<Mesh, T, EngineTag>                      Field_t;
  typedef typename Field_t::Domain_t                     Domain_t;
  typedef Array<FieldEngine<Mesh, T, EngineTag>::dimensions,
                T, EngineTag>                            Array_t;
  typedef Connector<Array_t, Lux>                        Base_t;
  typedef Connection<Lux>                                Connection_t;
  typedef Connector<Field_t, Lux>                        Connector_t;
  typedef Lux                                            ConnectionTag_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout).

  Connector(const char *conname, const Field_t &a, Connection_t &c,
	    int mode = ConnectionBase::out)
    : Base_t(conname, a.array(), c, mode)
  {
    // Initialize the base classes mesh info

    setupMeshInfo(a);
  }


  //============================================================
  // Connector Destructor
  //============================================================

  // When destroyed, disconnect from Lux.  This will eventually call
  // disconnect, where the data object will get deleted.

  virtual ~Connector()
  {
  }


  //============================================================
  // Connector operations
  //============================================================

  // Retarget this connector to a new data object.  For some items, the
  // data may be of different size, for others it will be the same size
  // always.

  void resize(const Field_t &newfield)
  {
    // For an array, re-initialize the internal array, and re-initialize
    // the copy array.  Use the base classes resize, and then reset the
    // mesh info.
    
    Base_t::resize(newfield.array());
    setupMeshInfo(newfield);
  }

private:
  //============================================================
  // Private methods
  //============================================================

  // Reset the mesh data in the base class to the data for the given
  // Field's mesh.

  void setupMeshInfo(const Field_t &f)
  {
    const int Dim = Field_t::dimensions;
    for (int d=0; d < 3; ++d)
      if (d < Dim)
	{
	  this->origin_m[d] = f.geometry().mesh().origin()(d);
	  this->spacing_m[d] =
	    f.geometry().boundingBox(f.geometry().totalDomain())[d].length();
	  this->spacing_m[d] /= this->size_m[d];
	}
      else
	{
	  this->origin_m[d] = 0.0;
	  this->spacing_m[d] = 0.0;
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

#endif     // POOMA_CONNECT_LUX_LUX_CONNECTOR_FIELD_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LuxConnector.Field.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
