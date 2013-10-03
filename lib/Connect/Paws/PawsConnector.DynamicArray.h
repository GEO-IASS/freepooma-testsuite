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

#ifndef POOMA_CONNECT_PAWS_PAWS_CONNECTOR_DYNAMICARRAY_H
#define POOMA_CONNECT_PAWS_PAWS_CONNECTOR_DYNAMICARRAY_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<DynamicArray<>, Paws>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<DynamicArray<>, Paws> is a specialization
// of the general Connector class used to connect a POOMA II DArray to another
// application via Paws.  This class will use the PAWS API to register this
// array for sharing, describe its storage and layout, and send/receive it.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Paws/PawsConnection.h"
#include "Connect/Paws/PawsConnector.Array.h"
#include "DynamicArray/DynamicArray.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<DynamicArray<>, Paws> is a specialization of Connector<T,Tag>
// that can map data from a POOMA II DynamicArray of type T to another Paws
// application.  The DynamicArray can be any type.
//
// This can only work with a Connection<Paws> connection object.
//
// The resulting connection will be fully dynamic, and both sides will
// need to do the extra steps required when data distributions change.
//
//-----------------------------------------------------------------------------

template<class T, class EngineTag>
class Connector<DynamicArray<T, EngineTag>, Paws>
  : public Connector<typename DynamicArray<T, EngineTag>::Base_t, Paws>
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef DynamicArray<T, EngineTag>                     DynamicArray_t;
  typedef typename DynamicArray_t::Domain_t              Domain_t;
  typedef typename DynamicArray_t::Base_t                Array_t;
  typedef Connector<Array_t, Paws>                       Base_t;
  typedef Connection<Paws>                               Connection_t;
  typedef Connector<DynamicArray_t,Paws>                 Connector_t;
  typedef Paws                                           ConnectionTag_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout)
  // The 'true' in the final arg to the base class indicates this should
  // always be a dynamic connection.

  Connector(const char *conname, const DynamicArray_t &a, Connection_t &c,
	    int mode)
    : Base_t(conname, a.array(), c, mode, true)
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

  void resize(const DynamicArray_t &newdata)
  {
    // For an array, re-initialize the internal array, and re-initialize
    // the copy array.  Use the base classes resize.
    
    Base_t::resize(newdata.array());
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

#endif     // POOMA_CONNECT_PAWS_PAWS_CONNECTOR_DYNAMICARRAY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsConnector.DynamicArray.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
