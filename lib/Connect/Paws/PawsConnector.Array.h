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

#ifndef POOMA_CONNECT_PAWS_PAWS_CONNECTOR_ARRAY_H
#define POOMA_CONNECT_PAWS_PAWS_CONNECTOR_ARRAY_H

//-----------------------------------------------------------------------------
// Classes:
// Connector<Array<>, Paws>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Connector<Array<>, Paws> is a specialization
// of the general Connector class used to connect a POOMA II Array to another
// application via Paws.  This class will use the PAWS API to register this
// Array for sharing, describe its storage and layout, and send/receive it.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Connect/Connector.h"
#include "Connect/Paws/PawsConnection.h"
#include "Connect/Paws/Resize.h"
#include "Pooma/Configuration.h"
#include "Pooma/Arrays.h"
#include "Utilities/PAssert.h"

#if POOMA_PAWS
#include "Paws/Paws.h"
#include "Paws/PawsArrayData.h"
#endif


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

class Paws;
template<class T> class PawsArrayData;


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
// Connector<Array<>, Paws> is a specialization of Connector<T,Tag> that can
// map data from a POOMA II Array of type T to another Paws application.  The
// Array can be any type and any dimension.
//
// Since POOMA II works with one context right now, the easiest way to
// implement this is to make a copy of the data into a single Brick engine
// Array, and use that to send/receive the actual data.  This is the fallback
// situation, for all types of Arrays that are not Brick's or MultiPatches
// of Brick's.  For the special case of Brick-based engines, we could create
// specializations that work with the blocks directly instead of with a copy.
//
// Connector<Array, Paws> should be created with a name for the data object and
// a non-const reference to the Array to share with another program.  This
// will store a shallow copy of that Array, and use its value to update the
// connection when update() is called.  Since this is a copy, the provided
// Array could be a view of some other Array so that you can share portions
// of a larger Array with another application.
//
// This can only work with a Connection<Paws> connection object.
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag>
class Connector<Array<Dim, T, EngineTag>, Paws> : public ConnectorBase
{
public:
  //============================================================
  // Public typedefs and enums
  //============================================================

  typedef Array<Dim, T, EngineTag> Array_t;
  typedef Array<Dim, T, Brick>     CopyArray_t;
  typedef T                        Scalar_t;
  typedef PawsArrayData<T>         PawsData_t;
  typedef Connection<Paws>         Connection_t;
  typedef Connector<Array_t,Paws>  Connector_t;
  typedef Paws                     ConnectionTag_t;


  //============================================================
  // Connector Constructor
  //============================================================

  // The constructor takes a string name, the data to connect, and the
  // data transfer mode (ConnectionBase::in, out, or inout)

  Connector(const char *conname, const Array_t &a, Connection_t &c, int mode,
	    bool dynamic = false)
    : ConnectorBase(conname, c, mode), array_m(a), data_m(0), ptr_m(0),
      dynamic_m(dynamic), resize_m(dynamic)
  {
#if POOMA_PAWS
    // Determine the transfer mode flag for PAWS

    int pawsmode = PAWS_IN;
    if (transferMode() == ConnectionBase::out)
      pawsmode = PAWS_OUT;
    else if (transferMode() == ConnectionBase::inout)
      pawsmode = PAWS_INOUT;

    // Set up an internal Array that uses a Brick engine, of the same
    // total domain.

    copy_m.initialize(array_m.domain());
    copy_m = array_m;		// make sure to initialize to get mem locality
    Pooma::blockAndEvaluate();

    // We just register this copied Array's brick for sharing,
    // and describe it has having one block.  First set up info to use
    // in PAWS API calls.
    
    int d, first[Dim], last[Dim], stride[Dim], localblocks[3*Dim];
    for (d=0; d < Dim; ++d)
      {
	first[d]  = localblocks[3*d]     = copy_m.first(d);
	last[d]   = localblocks[3*d + 1] = copy_m.last(d);
	stride[d] = localblocks[3*d + 2] = 1;
      }

    // Get a pointer to the beginning of the brick of data

    ptr_m = &(copy_m(copy_m.firsts()));

    // Create a new PawsArrayData object to manage sending/receiving
    // this array.

    data_m = new PawsData_t(name().c_str(),
			    &ptr_m,
			    Dim,
			    first,
			    last,
			    stride,
			    1,
			    localblocks,
			    pawsmode,
			    PAWS_SYNC,
			    PAWS_COLUMN,
			    pawsConnection().paws());
#endif
  }


  //============================================================
  // Connector Destructor
  //============================================================

  // When destroyed, disconnect from Paws.  This will eventually call
  // disconnect, where the data object will get deleted.

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
    // For an array, re-initialize the internal array, and indicate
    // the next send/rec must do a resize.

    array_m.initialize(newarray);
    resize_m = true;
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
	    bool doreceive = true;

	    // If we are receiving, and need to resize, do so now

	    if (resize_m || dynamic_m)
	      {
		// Wait for updated size info ... if nothing is
		// connected to us, though, we can just quit, and this
		// is signaled by having resizeWait() return an error.

		if (doreceive = (pawsData().resizeWait() != PAWS_ERROR))
		  {

		    // Get new total domain info

		    int d, first[Dim], last[Dim], stride[Dim];
		    int localblocks[3*Dim];
		    pawsData().domain(first, last, stride);
		    Interval<Dim> newdomain;
		    for (d=0; d < Dim; ++d)
		      {
			localblocks[3*d]     = first[d];
			localblocks[3*d + 1] = last[d];
			localblocks[3*d + 2] = stride[d];
			newdomain[d] = Interval<1>(first[d], last[d]);
		      }

		    // Resize the existing array, if it has the wrong size.

		    if (copy_m.domain().size() != newdomain.size())
		      {
			copy_m.initialize(newdomain);
			ptr_m = &(copy_m(copy_m.firsts()));
		      }

		    // Update the paws object with the new storage.

		    pawsData().update(&ptr_m,
				      Dim,
				      first,
				      last,
				      stride,
				      1,
				      localblocks);
		  }
	      }

	    // Receive the data now, if we still need to

	    if (doreceive)
	      {
		pawsData().receive();

		// Check the total size - are they equal?  If not, we must
		// change the size of the destination array.
		// Resize changes patches in array_m to get the total sizes
		// equal.  Most reasonable solution - spread total
		// size over all array_m patches as equally as possible.
		// If the engine type for array_m does not allow this,
		// we will throw a run-time error.

		if (array_m.domain().size() != copy_m.domain().size())
		  Resize<Array_t>::resize(array_m, copy_m.domain());

		// Finally, copy data over to the destination array
		array_m = copy_m;
		Pooma::blockAndEvaluate();
	      }
	  }

	if (transferMode() == ConnectionBase::out ||
	    transferMode() == ConnectionBase::inout)
	  {
	    // If we are sending, and need to resize, do so now

	    if (resize_m || dynamic_m)
	      {
		// Resize the temporary.

		copy_m.initialize(array_m.domain());

		// Tell PAWS about the new size and location of data

		int d, first[Dim], last[Dim], stride[Dim], localblocks[3*Dim];
		for (d=0; d < Dim; ++d)
		  {
		    first[d]  = localblocks[3*d]     = copy_m.first(d);
		    last[d]   = localblocks[3*d + 1] = copy_m.last(d);
		    stride[d] = localblocks[3*d + 2] = 1;
		  }

		// Get a pointer to the beginning of the brick of data

		ptr_m = &(copy_m(copy_m.firsts()));

		// Call the resize function now after setting up args.

		pawsData().resize(&ptr_m,
				  Dim,
				  first,
				  last,
				  stride,
				  1,
				  localblocks);
	      }

	    // Send the data now, by copying it to the temporary and
	    // then sending data from the temporary.

	    copy_m = array_m;
	    Pooma::blockAndEvaluate();
	    pawsData().send();
	  }

	// We don't need a resize any more, unless we're dynamic.

	resize_m = false;
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
  //============================================================
  // Private data
  //============================================================

  // The array we're connecting

  Array_t array_m;

  // The internal copy of the data we use to collect data into a contiguous
  // brick

  CopyArray_t copy_m;

  // The Paws connected-data object we're using

  PawsData_t *data_m;

  // A pointer to the beginning of the data, used by the data_m object.

  Scalar_t *ptr_m;

  // A flag indicating whether the send/receives should assume
  // dynamic data objects

  bool dynamic_m;

  // A flag indicating if the 'resize' method has been called recently.

  bool resize_m;

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

#endif     // POOMA_CONNECT_PAWS_PAWS_CONNECTOR_ARRAY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PawsConnector.Array.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
