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

#ifndef POOMA_CONNECT_PAWS_PAWS_CONNECTOR_STRING_H
#define POOMA_CONNECT_PAWS_PAWS_CONNECTOR_STRING_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<std::string, Paws>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<std::string, Paws> is a specialization
// of the general Connector class used to connect a possibly-changing string to
// an application via Paws.  The string will be resized to accomodate the
// string if it is being received.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connector.h"
#include "Connect/Paws/PawsConnection.h"
#include "Utilities/PAssert.h"
#include "Pooma/Configuration.h"
#include <string>

#if POOMA_PAWS
#include "Paws/Paws.h"
#include "Paws/PawsScalarData.h"
#include "Paws/PawsStringData.h"
#endif


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

class Paws;


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<std::string, Paws> is a specialization of Connector<T,Tag> that
// can map data from a string to another Paws application.  The string
// length can vary from transfer to transfer; the receiver will automatically
// resize its string storage to accomodate a string that changes size.
//
// Connector<std::string, Paws> should be created with a name for the string
// and a non-const reference to the string to share with another program.  This
// will store a reference to that string, and use its value to update the
// connection when update() is called.  The string will be resized if needed.
//
// This can only work with a Connection<Paws> connection object.
//
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

template<>
class Connector<std::string, Paws> : public ConnectorBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef std::string                  Scalar_t;
  typedef Paws                         ConnectionTag_t;
  typedef PawsStringData<char>         PawsData_t;
  typedef Connection<Paws>             Connection_t;
  typedef Connector<std::string,Paws>  Connector_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout)

  Connector(const char *conname, Scalar_t &a, Connection_t &c, int mode)
    : ConnectorBase(conname, c, mode), scalar_m(&a), data_m(0),
      buffer_m(0), buflen_m(0)
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
    return *scalar_m;
  }

  // Another version of the same routine to return the scalar (here a string)

  Scalar_t &string() const
  {
    return *scalar_m;
  }

  //============================================================
  // Connector operations
  //============================================================

  // Retarget this connector to a new data object.  For some items, the
  // data may be of different size, for others it will be the same size
  // always.

  void resize(Scalar_t &newscalar)
  {
    // For a string, this is simple.

    scalar_m = &newscalar;
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
	    // When receiving, first find out the size of the sent string,
	    // then make room for it, then do the receive.

	    // First wait for updated size info

	    if (pawsData().resizeWait() != PAWS_ERROR)
	      {
		// Get the new string length

		int newsize;
		pawsData().size(&newsize);

		// Resize our string to this new size, and our internal buffer.

		scalar_m->resize(newsize);
		resizeBuffer(newsize);

		// Update our internal info and receive

		pawsData().update(buffer_m, scalar_m->length());
		pawsData().receive();

		// Copy back the buffer to the string

		*scalar_m = buffer_m;
	      }
	  }

	if (transferMode() == ConnectionBase::out ||
	    transferMode() == ConnectionBase::inout)
	  {
	    // Copy the string into our buffer

	    resizeBuffer(scalar_m->length());
	    strcpy(buffer_m, scalar_m->c_str());

	    // When sending, first let the other side know you're sending
	    // a newly changed string size, then send the string.

	    // First do the resize

	    pawsData().resize(buffer_m, scalar_m->length());

	    // And then the send, which will not actually proceed
	    // until the other side has finished its 'update' action.

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
  // A pointer to the string we're connecting

  Scalar_t *scalar_m;

  // The Paws string data object we're using

  PawsData_t *data_m;

  // A char buffer, and its length

  char *buffer_m;
  int buflen_m;

  // Resize our internal buffer to have space for a string of the given
  // size.  Does not preserve contents if the new size is > old size.

  void resizeBuffer(int newsize)
  {
    if (newsize > buflen_m || buffer_m == 0)
      {
	if (buffer_m != 0)
	  delete [] buffer_m;

	buflen_m = newsize;
	buffer_m = new char[newsize + 2];
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

#endif     // POOMA_CONNECT_PAWS_PAWS_CONNECTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsConnector.String.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
