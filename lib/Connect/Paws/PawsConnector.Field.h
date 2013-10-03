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

#ifndef POOMA_CONNECT_PAWS_PAWS_CONNECTOR_FIELD_H
#define POOMA_CONNECT_PAWS_PAWS_CONNECTOR_FIELD_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<Field<>, Paws>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<Field<>, Paws> is a specialization
// of the general Connector class used to connect a POOMA II Field to another
// application via Paws.  This class will use the PAWS API to register this
// Field for sharing, describe its storage and layout, and send/receive it.
// Only the storage data for Field will be communicated; the mesh and bc
// information is not transferred.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Paws/PawsConnection.h"
#include "Connect/Paws/PawsConnector.Array.h"
#include "Field/Field.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<Field<>, Paws> is a specialization of Connector<T,Tag> that can
// map data from a POOMA II Field of type T to another Paws application.  The
// Field can be any type and any dimension.
//
// Since POOMA II works with one context right now, the easiest way to
// implement this is to make a copy of the data into a single Brick engine
// Array, and use that to send/receive the actual data.  This is the fallback
// situation, for all types of Fields that are not Brick's or MultiPatches
// of Brick's.  For the special case of Brick-based engines, we could create
// specializations that work with the blocks directly instead of with a copy.
//
// Connector<Field, Paws> should be created with a name for the data object and
// a non-const reference to the Field to share with another program.  This
// will store a shallow copy of that Field, and use its value to update the
// connection when update() is called.  Since this is a copy, the provided
// Field could be a view of some other Field so that you can share portions
// of a larger Field with another application.
//
// This can only work with a Connection<Paws> connection object.
//
//-----------------------------------------------------------------------------

template<class Geom, class T, class EngineTag>
class Connector<Field<Geom, T, EngineTag>, Paws>
  : public Connector<typename ArrayView<Field<Geom,T,EngineTag>,
                     typename Field<Geom,T,EngineTag>::Domain_t>::Type_t,
                     Paws>
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef Field<Geom, T, EngineTag>                      Field_t;
  typedef typename Field_t::Domain_t                     Domain_t;
  typedef typename ArrayView<Field_t, Domain_t>::Type_t  Array_t;
  typedef Connector<Array_t, Paws>                       Base_t;
  typedef Connection<Paws>                               Connection_t;
  typedef Connector<Field_t,Paws>                        Connector_t;
  typedef Paws                                           ConnectionTag_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout)

  Connector(const char *conname, const Field_t &a, Connection_t &c, int mode)
    : Base_t(conname, a.array(), c, mode)
  {
  }


  //============================================================
  // Connector Destructor
  //============================================================

  // When destroyed, disconnect from Paws.  This will eventually call
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
    // the copy array.  Use the base classes resize.
    
    Base_t::resize(newfield.array());
  }

private:
  //============================================================
  // Private methods
  //============================================================

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Connector();
  Connector(const Connector_t &);
  const Connector_t &operator=(const Connector_t &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_PAWS_PAWS_CONNECTOR_FIELD_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsConnector.Field.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
