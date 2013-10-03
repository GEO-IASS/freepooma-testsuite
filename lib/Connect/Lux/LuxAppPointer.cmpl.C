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

// POOMA include files

#include "Connect/Lux/LuxAppPointer.h"
#include "Pooma/Configuration.h"
#include "Utilities/PAssert.h"
#include "Utilities/Inform.h"

// LUX include files

#if POOMA_LUX
# include "interface/script.h"
#endif


//-----------------------------------------------------------------------------
// static data members of LuxAppPointer
//-----------------------------------------------------------------------------

// The VizTool object.  There is only one of these at a time.  If
// you create a LuxAppPointer and this pointer is null, the application
// will initialize the Lux display.

VizTool *LuxAppPointer::lux_s   = 0;
int      LuxAppPointer::users_s = 0;


//----------------------------------------------------------------------
//
// LuxAppPointer constructor
// The constructor takes a string name; the type is "lux".
//
//----------------------------------------------------------------------

LuxAppPointer::LuxAppPointer(const char *conname)
  : connected_m(false)
{
  // if we do not have a Lux connection yet, create it now

  if (lux_s == 0)
    {
#if POOMA_LUX
      PAssert(conname != 0);
      lux_s = lux_init(const_cast<char *>(conname));
#endif
    }

  // indicate we're one more object using Lux

  if (connected_m = (lux_s != 0))
    users_s++;
}


//----------------------------------------------------------------------
//
// LuxAppPointer destructor
//
//----------------------------------------------------------------------

LuxAppPointer::~LuxAppPointer()
{
  close();
}


//----------------------------------------------------------------------
//
// Return the Lux connection object.
//
//----------------------------------------------------------------------

VizTool &LuxAppPointer::lux() const
{
  PAssert(lux_s != 0);
  PAssert(connected());
  return *lux_s;
}


//----------------------------------------------------------------------
//
// Perform a Lux interaction
//
//----------------------------------------------------------------------

void LuxAppPointer::poll()
{
#if POOMA_LUX
  PAssert(connected());
  lux_interact();
#endif
}


//----------------------------------------------------------------------
//
// Shut down the connection, and put ourselves in an "unconnected" state.
//
//----------------------------------------------------------------------

void LuxAppPointer::close()
{
  // close lux connection if we're the last one

  if (connected() && --users_s == 0)
    {
#if POOMA_LUX
      delete lux_s;
#endif
      lux_s = 0;
    }

  connected_m = false;
}


//----------------------------------------------------------------------
//
// Create a new array tool, and add it in.  Return a ReadFieldTool handle.
//
//----------------------------------------------------------------------

ReadFieldTool *LuxAppPointer::createArray(const std::string &nm)
{
  ReadFieldTool *tool = 0;

#if POOMA_LUX
  PAssert(connected());

  // Create the new tool

  tool = new ReadFieldTool;

  // Then give it to Lux

  lux().connect(static_cast<void *>(tool), const_cast<char *>(nm.c_str()),
		tool, FieldDataType);
#endif

  // Return the tool, or zero if nothing happened

  return tool;
}


//----------------------------------------------------------------------
//
// Initialize array tool to start getting new values.
//
//----------------------------------------------------------------------

void LuxAppPointer::beginArray(ReadFieldTool *tool, int datatype,
			       int *size, float *spacing, float *origin)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);
  PAssert(size != 0);
  PAssert(spacing != 0);
  PAssert(origin != 0);

  // Determine the data type

  int luxtype = vizStructuredFieldDataType::ACLVIS_SCALAR;
  if (datatype == LuxAppPointer::vector)
    luxtype = vizStructuredFieldDataType::ACLVIS_VECTOR;
  else if (datatype == LuxAppPointer::tensor)
    luxtype = vizStructuredFieldDataType::ACLVIS_TENSOR;

  // Compute the total size

  int totsize = 1;
  for (int d=0; d < 3; ++d)
    totsize *= size[d];

  // Tell Lux what type and how much data to expect, and of what type

  tool->GetVizData()->InitData(totsize, luxtype);

  // Set up the mesh information

  tool->SetDimensions(size);
  tool->SetAspectRatio(spacing);
  tool->SetOrigin(origin);
#endif
}


//----------------------------------------------------------------------
//
// Add in the value of the Nth array element.
// 'val' has one element for datatype == scalar, 3 elements for
// datatype == vector, etc.
//
//----------------------------------------------------------------------

void LuxAppPointer::insertArray(ReadFieldTool *tool, int datatype,
				int indx, float *val)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);
  PAssert(val != 0);
  PAssert(indx >= 0);

  // Set value based on the data type

  if (datatype == LuxAppPointer::scalar)
    tool->GetVizData()->AddScalar(indx, *val);
  else if (datatype == LuxAppPointer::vector)
    tool->GetVizData()->AddVector(indx, val);
#endif
}


//----------------------------------------------------------------------
//
// Finish and update existing array tool with new data
//
//----------------------------------------------------------------------

void LuxAppPointer::endArray(ReadFieldTool *tool,
			     const std::string &nm)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);

  // Finish up and give data to Lux for display

  tool->PrepareFinishedData();
  lux().update(const_cast<char *>(nm.c_str()));
#endif
}


//----------------------------------------------------------------------
//
// Remove and delete an existing array tool
//
//----------------------------------------------------------------------

void LuxAppPointer::destroyArray(ReadFieldTool *tool,
				 const std::string &nm)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);

  // Disconnect the named item, and delete the tool

  lux().disconnect(const_cast<char *>(nm.c_str()));
  delete tool;
#endif
}


//----------------------------------------------------------------------
//
// Create a new particles tool, and add it in.  Return ReadParticleTool.
//
//----------------------------------------------------------------------

ReadParticleTool *LuxAppPointer::createParticles(const std::string &nm)
{
  ReadParticleTool *tool = 0;

#if POOMA_LUX
  PAssert(connected());

  // Create the new tool

  tool = new ReadParticleTool;

  // Then give it to Lux

  lux().connect(static_cast<void *>(tool), const_cast<char *>(nm.c_str()),
		tool, ParticleDataType);
#endif

  // Return the tool, or zero if nothing happened

  return tool;
}


//----------------------------------------------------------------------
//
// Initialize existing particles tool for new data
//
//----------------------------------------------------------------------

void LuxAppPointer::beginParticles(ReadParticleTool *tool, int datatype,
				   int totsize)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);
  PAssert(totsize >= 0);

  // Determine the data type

  int luxtype = vizStructuredFieldDataType::ACLVIS_SCALAR;
  if (datatype == LuxAppPointer::vector)
    luxtype = vizStructuredFieldDataType::ACLVIS_VECTOR;
  else if (datatype == LuxAppPointer::tensor)
    luxtype = vizStructuredFieldDataType::ACLVIS_TENSOR;

  // Initialize the data tool

  tool->GetVizData()->InitData(totsize, luxtype, 1);
#endif
}


//----------------------------------------------------------------------
//
// Add in the value of the Nth particle. 'pos' must be a 3-vector.
// 'val' has one element for datatype == scalar, 3 elements for
// datatype == vector, etc.
//
//----------------------------------------------------------------------

void LuxAppPointer::insertParticles(ReadParticleTool *tool, int datatype,
				    int indx, float *pos, float *val, int id)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);
  PAssert(pos != 0);
  PAssert(val != 0);
  PAssert(indx >= 0);

  // Set position for this particle

  tool->GetVizData()->AddPoint(indx, pos);

  // Set value for this particle

  if (datatype == LuxAppPointer::scalar)
    tool->GetVizData()->AddScalar(indx, *val);
  else if (datatype == LuxAppPointer::vector)
    tool->GetVizData()->AddVector(indx, val);

  // Set ID for the particle

  tool->GetVizData()->AddIdInfoVal(indx, id);
#endif
}


//----------------------------------------------------------------------
//
// Finish and update existing particles tool with new data
//
//----------------------------------------------------------------------

void LuxAppPointer::endParticles(ReadParticleTool *tool,
				 const std::string &nm)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);

  // Finish up and give data to Lux for display

  tool->PrepareFinishedData();
  lux().update(const_cast<char *>(nm.c_str()));
#endif
}


//----------------------------------------------------------------------
//
// Remove and delete an existing particles tool
//
//----------------------------------------------------------------------

void LuxAppPointer::destroyParticles(ReadParticleTool *tool,
				     const std::string &nm)
{
#if POOMA_LUX
  PAssert(connected());
  PAssert(tool != 0);

  // Disconnect the named item, and delete the tool

  lux().disconnect(const_cast<char *>(nm.c_str()));
  delete tool;
#endif
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LuxAppPointer.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
