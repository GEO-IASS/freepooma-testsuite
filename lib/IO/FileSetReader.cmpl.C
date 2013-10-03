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

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

// Pooma

#include "IO/FileSetReader.h"   // Class definition
#include "Pooma/Pooma.h"        // Pooma::context()
#include "Tulip/RemoteProxy.h"  // RemoteProxy template
#include "Utilities/PAssert.h"  // PInsist

// System

#include <string>

#include <iostream>

#if !POOMA_NO_IOS_HEADER
# include <ios>
#endif

//-----------------------------------------------------------------------------
// Implementation of FileSetReader non-template member functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor
//
// Currently this does a PInsist if it encounters problems opening or
// reading files. If exceptions are enabled, this is OK.
// Unfortunately, exceptions are usually disabled and failure to open
// a file will result in an abort. We should probably change to a
// separate "open" function that returns success or failure in order
// to make this software more flexible. 
//-----------------------------------------------------------------------------

template <int Dim>
FileSetReader<Dim>::
FileSetReader(const char *pfileset)
  : basename_m(pfileset),
    diskLayout_m(pfileset), 
    pmetafile_m(0),
    currentRecord_m(0),
    currentField_m(0),
    mycontext_m(Pooma::context()), 
    iocontext_m(0), // For now!!!
    bytesReversed_m(false)
{ 
  if (mycontext_m == iocontext_m)
    {
      pmetafile_m = new DiskMeta(pfileset);
    }
}

//-----------------------------------------------------------------------------
// open(bool abortOnError = false)
//
// Open the fileset files. By default, this returns false if any
// problems are encountered, like files not existing. If abortOnError
// is true, then PInsist is used instead.
//-----------------------------------------------------------------------------

#define IOInsist(c,str) \
if (!abortOnError && !(c)) { errorMsg_m = str; return false; } \
else PInsist(c,str)

template <int Dim>
bool FileSetReader<Dim>::readMeta()
{
  PAssert(mycontext_m == iocontext_m);
  PAssert(pmetafile_m);

  // Use IOInsist, but never abort since it won't be sync'd up

  const bool abortOnError = false;

  // Open and parse the metafile

  bool success;
  success = pmetafile_m->open(false);
  IOInsist(success, "Couldn't open .meta file");

  success = pmetafile_m->read(false);
  IOInsist(success, "Couldn't read .meta file");

  numRecords_m = pmetafile_m->numRecords();
  fieldsPerRecord_m = pmetafile_m->fieldsPerRecord();

  // Check consistency with the meta file

  IOInsist(pmetafile_m->dimension() == Dim,
           "File set has wrong dimensionality");

  IOInsist(pmetafile_m->numFileSets() == 1,
           "Multiple filesets not supported (YET)!");

  return true;
}

template <int Dim>
bool FileSetReader<Dim>::open(bool abortOnError)
{
  // Read and check the metafile

  bool success = true;

  if (mycontext_m == iocontext_m)
    {
      // Read the .meta file on the IO context

      success = readMeta();
    }

  // Broadcast success/failure to all contexts; return false if 
  // readMeta failed.

  bool iosuccess = RemoteProxy<bool>(success);
  if (!iosuccess) return false;

  // Make local copies of certain fields to braodcast. 

  int nrecs, fprec;
  Domain_t dom;

  if (mycontext_m == iocontext_m)
    {
      // local copies

      nrecs = numRecords_m;
      fprec = fieldsPerRecord_m;

      for (int d = 0; d < Dim; ++d)
        {
          dom[d] = pmetafile_m->domain(d);
        }
    }

  // Broadcast the number of records and the number of fields per record.

  numRecords_m = RemoteProxy<int>(nrecs);
  fieldsPerRecord_m = RemoteProxy<int>(fprec);

  // Broadcast the global domain

  domain_m = RemoteProxy<Domain_t>(dom).value();

  // Open the disk layout file (DL broadcasts success internally).

  iosuccess = diskLayout_m.open();

  if (!iosuccess) return false;

  // Assume that all files have the same byte ordering

  bytesReversed_m = diskLayout_m.bytesReversed();

  // Open the .offset and .data files

  const char *errmsg = 0;

  if (mycontext_m == iocontext_m)
    {
      std::string offsetname = basename_m + ".offset";
      foffset_m.open(offsetname.c_str(), std::ios::binary);
      if (!foffset_m) 
        {
          errmsg = "Couldn't open .offset file";
          success = false;
        }
      else
        {
          std::string dataname = basename_m + ".data";
          fdata_m.open(dataname.c_str(), std::ios::binary);
          if (!fdata_m) 
            {
              errmsg = "Couldn't open .data file";
              success = false;
            }
        }
    }

  // Broadcast error status

  iosuccess = RemoteProxy<bool>(success);
  IOInsist(iosuccess, errmsg);

  return true;
}

//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

// Trivial since fstreams close their files on destruction.

template <int Dim>
FileSetReader<Dim>::~FileSetReader()
{ 
  if (pmetafile_m) delete pmetafile_m;
}

//-----------------------------------------------------------------------------
// readFieldID
// 
// This one read operation is independent of T or any other template
// parameters, so I separated it out to where it can be
// preinstantiated.
//-----------------------------------------------------------------------------

template <int Dim>
int FileSetReader<Dim>::readFieldID()
{
  PAssert(mycontext_m == iocontext_m);

  // Read the field ID from the .offset file. 
  // This is only done on the IO context. 

  int fieldID;

  foffset_m.read((char*)&fieldID, sizeof(int));
  if (foffset_m.gcount() != sizeof(int)) return -1;
  if (bytesReversed_m) ::reverseBytes(fieldID);

  return fieldID;
}

//-----------------------------------------------------------------------------
// Template preinstantiations
//-----------------------------------------------------------------------------

template class FileSetReader<1>;
template class FileSetReader<2>;
template class FileSetReader<3>;

// Uncomment this if you need to read DiskFields with Dim > 3.

#if 0
template class FileSetReader<4>;
template class FileSetReader<5>;
template class FileSetReader<6>;
template class FileSetReader<7>;
#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetReader.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
