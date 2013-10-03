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

#ifndef POOMA_CONNECT_LUX_LUX_APP_POINTER_H
#define POOMA_CONNECT_LUX_LUX_APP_POINTER_H

//-----------------------------------------------------------------------------
// Classes:
// LuxAppPointer
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// LuxAppPointer holds a pointer to a VizTool object, and makes
// sure it is a singleton object.  It maintains a count on objects that
// have requested to use it, and will delete the VizTool object
// when the count goes back to zero.  To use it, create an instance in
// your own class.  Use the "lux()" method to get
// a reference to the VizTool object that it stores.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include <string>


//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

class VizTool;
class ReadFieldTool;
class ReadParticleTool;

template<int D, class T, class E> class Vector;
template<int D, class T, class E> class Tensor;


//-----------------------------------------------------------------------------
//
// Full Description:
//
// LuxAppPointer is a wrapper around an instance of a VizTool
// class.  VizTool must be a "singleton" type instance, and this
// class uses a static pointer to a VizTool instance and a counter
// to tell when the instance must be created or deleted.  If more than
// one LuxAppPointer is created, each one will refer to the same
// VizTool object, until the last LuxAppPointer goes away.
//
// LuxAppPointer is created with the arguments needed to create a
// VizTool, namely a string with a name for the display and script file.
// The destructor will call the "close()"
// method, that will delete the VizTool if it is the last object
// around of this type.
//
// LuxAppPointer also provides methods "poll()" and others
// that are just deferred to the VizTool.  The implementation for
// these are put in the .cmpl.cpp file so that we can encapsulate the
// details of the Lux code in a separate file.
//
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

class LuxAppPointer
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // Enumeration of data types

  enum { scalar = 1, vector = 2, tensor = 3 };


  //============================================================
  // Constructor
  //============================================================

  // This class is constructed with a string name, used to
  // initialize LUX.

  LuxAppPointer(const char *conname);


  //============================================================
  // Destructor
  //============================================================

  // When destroyed, disconnect from Lux.

  virtual ~LuxAppPointer();


  //============================================================
  // Accessors
  //============================================================

  // Return whether we are connected.

  bool connected() const
    {
      return connected_m;
    }

  // Return the Lux connection object.

  VizTool &lux() const;


  //============================================================
  // Operations
  //============================================================

  // Perform a Lux poll.

  void poll();

  // Shut down the connection, and put ourselves in an "unconnected" state.

  void close();


  //============================================================
  // Create/destroy/manipulate tools for Arrays
  //============================================================

  // Create a new array tool, and add it in.  Return a ReadFieldTool handle.

  ReadFieldTool *createArray(const std::string &nm);

  // Start working on a new set of values for an array.

  void beginArray(ReadFieldTool *tool, int datatype,
		  int *size, float *spacing, float *origin);

  // Add in an element of the array

  void insertArray(ReadFieldTool *tool, int datatype,
		   int indx, float *val);

  // Finish up work on a set of values for an array, and update Lux.

  void endArray(ReadFieldTool *tool, const std::string &nm);

  // Remove and delete an existing array tool

  void destroyArray(ReadFieldTool *tool, const std::string &nm);


  //============================================================
  // Create/destroy/manipulate tools for Particles
  //============================================================

  // Create a new particles tool, and add it in.  Return ReadParticleTool.

  ReadParticleTool *createParticles(const std::string &nm);

  // Start working on a new set of values for a set of particles.

  void beginParticles(ReadParticleTool *tool, int datatype, int totsize);

  // Add in an element of the particles

  void insertParticles(ReadParticleTool *tool, int datatype,
		       int indx, float *pos, float *val, int id);

  // Finish up work on a set of values for particles, and update Lux.

  void endParticles(ReadParticleTool *tool, const std::string &nm);

  // Remove and delete an existing particles tool

  void destroyParticles(ReadParticleTool *tool, const std::string &nm);

private:
  // The Lux connection object.

  static VizTool *lux_s;

  // The number of instances using Lux

  static int users_s;

  // Are we connected, and using the app pointer?

  bool connected_m;

  // The default and copy constructors are made private and undefined
  // since they should not be used

  LuxAppPointer();
  LuxAppPointer(const LuxAppPointer &);
  LuxAppPointer &operator=(const LuxAppPointer &);
};


//-----------------------------------------------------------------------------
//
// A simple partially-specialized class to tell if something is
// a scalar, vector, or tensor
//
//-----------------------------------------------------------------------------

template<class T>
struct LuxDataType
{
  enum { datatype   = LuxAppPointer::scalar };
  enum { dimensions = 1 };
  static void copy(const T &val, float data[])
  {
    data[0] = static_cast<float>(val);
  }
};

template<int Dim, class T, class E>
struct LuxDataType< Vector<Dim, T, E> >
{
  enum { datatype = LuxAppPointer::vector };
  enum { dimensions = 3 };
  static void copy(const Vector<Dim,T,E> &val, float data[])
  {
    for (int d=0; d < 3; ++d)
      data[d] = (d < Dim ? static_cast<float>(val(d)) : 0.0);
  }
};

template<int Dim, class T, class E>
struct LuxDataType< Tensor<Dim, T, E> >
{
  enum { datatype = LuxAppPointer::tensor };
  enum { dimensions = 9 };
  static void copy(const Tensor<Dim,T,E> &val, float data[])
  {
    int sz = 0;
    for (int d1=0; d1 < 3; ++d1)
      for (int d2=0; d2 < 3; ++d2, ++sz)
	data[sz] = (d1<Dim && d2<Dim ? static_cast<float>(val(d1,d2)):0.0);
  }
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_LUX_LUX_APP_POINTER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LuxAppPointer.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
