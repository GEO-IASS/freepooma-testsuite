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

#ifndef POOMA_CONNECT_PAWS_PAWS_CONNECTOR_H
#define POOMA_CONNECT_PAWS_PAWS_CONNECTOR_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<T, Paws>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<T, Paws> is a specialization
// of the general Connector class used to connect a scalar object to an
// application via Paws.  This is the default case when you try to connect
// an object via Paws - it is treated as a scalar.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connector.h"
#include "Connect/Paws/PawsConnection.h"
#include "Utilities/PAssert.h"
#include "Pooma/Configuration.h"

#if POOMA_PAWS
#include "Paws/Paws.h"
#include "Paws/PawsScalarData.h"
#endif


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

class Paws;
template<class T> class PawsScalarData;


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<T, Paws> is a specialization of Connector<T,Tag> that can
// map data from an object of type T to another Paws application.  The
// general version of this is for the case where T is a scalar; for other
// data types, like Array, other specializations must be created.
//
// Connector<T, Paws> should be created with a name for the data object and
// a non-const reference to the scalar to share with another program.  This
// will store a reference to that scalar, and use its value to update the
// connection when update() is called.
//
// This can only work with a Connection<Paws> connection object.
//
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

template<class T>
class Connector<T, Paws> : public ConnectorBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef T                  Scalar_t;
  typedef Paws               ConnectionTag_t;
  typedef PawsScalarData<T>  PawsData_t;
  typedef Connection<Paws>   Connection_t;
  typedef Connector<T,Paws>  Connector_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout)

  Connector(const char *conname, Scalar_t &a, Connection_t &c, int mode)
    : ConnectorBase(conname, c, mode), userScalar_m(&a), data_m(0)
  {
#if POOMA_PAWS
    // Determine the transfer mode flag for PAWS
    
    int pawsmode = PAWS_IN;
    if (transferMode() == ConnectionBase::out)
      pawsmode = PAWS_OUT;
    else if (transferMode() == ConnectionBase::inout)
      pawsmode = PAWS_INOUT;

    // Create a new PawsScalarData object to manage sending/receiving
    // this scalar.

    data_m = new PawsData_t(name().c_str(),
			    &scalar_m,
			    pawsmode,
			    PAWS_SYNC,
			    pawsConnection().paws());
#endif
  }


  //============================================================
  // Connector Destructor
  //============================================================

  // When destroyed, disconnect from Paws.

  virtual ~Connector()
  {
    if (connected())
      connection().disconnect(this);
  }


  //============================================================
  // Connector accessors
  //============================================================

  // Return the connection, cast as a Connection<Paws>

  Connection_t &pawsConnection() const
  {
    PAssert(connected());
    return dynamic_cast<Connection_t &>(connection());
  }

  // Return the scalar paws data object

  PawsData_t &pawsData() const
  {
    PAssert(data_m != 0);
    return *data_m;
  }

  // Return the scalar itself

  Scalar_t &scalar() const
  {
    return *userScalar_m;
  }

  //============================================================
  // Connector operations
  //============================================================

  // Retarget this connector to a new data object.  For some items, the
  // data may be of different size, for others it will be the same size
  // always.

  void resize(Scalar_t &newscalar)
  {
    // For a scalar, this is simple, just use the given pointer.

    userScalar_m = &newscalar;
  }


  //============================================================
  // ConnectorBase operations
  //============================================================

  // Do special activities to disconnect ourselves from the Connection<Paws>.

  virtual void disconnect()
  {
    // Disconnect by calling the finalize method, and then delete the
    // data object

#if POOMA_PAWS
    PAssert(data_m != 0);
    data_m->finalize();
    delete data_m;
    data_m = 0;
#endif
  }

  // Update our connection.  For Paws, this results in a data transfer
  // operation, either send or receive, based on the connection method.

  virtual void update()
  {
#if POOMA_PAWS
    if (connected())
      {
	// Either send or receive, based on the transfer mode
	
	if (transferMode() == ConnectionBase::in ||
	    transferMode() == ConnectionBase::inout)
	  {
	    pawsData().receive();
	    *userScalar_m = scalar_m;
	  }

	if (transferMode() == ConnectionBase::out ||
	    transferMode() == ConnectionBase::inout)
	  {
	    scalar_m = *userScalar_m;
	    pawsData().send();
	  }
      }
#endif
  }

  // Allow for interaction with the connection.  An optional string
  // can be provided to tell how to do the interaction.  Here, does nothing
  // except poll.

  virtual void interact(const char * = 0)
  {
    if (connected())
      pawsConnection().poll();
  }

private:
  // A pointer to the scalar we're connecting

  Scalar_t *userScalar_m;

  // A scalar actually used in the send/receive calls

  Scalar_t scalar_m;

  // The Paws scalar data object we're using

  PawsData_t *data_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  Connector();
  Connector(const Connector_t &);
  const Connector_t &operator=(const Connector_t &);
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_PAWS_PAWS_CONNECTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsConnector.h,v $   $Author: richard $
// $Revision: 1.8 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
