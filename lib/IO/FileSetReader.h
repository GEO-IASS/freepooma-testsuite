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

/** @file
 * @ingroup IO
 * @brief
 * FileSetReader<Dim> manages the reading of Arrays and Fields from
 * "DiscField" format files.
 */

#ifndef POOMA_IO_FILESETREADER_H
#define POOMA_IO_FILESETREADER_H

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

// Pooma

#include "IO/DiskLayout.h"
#include "IO/DiskMeta.h"
#include "Domain/Interval.h"

// System

#include <fstream> // ifstream

//-----------------------------------------------------------------------------
// Forward Declarations
//-----------------------------------------------------------------------------

template <int Dim, class T, class ETag> class Array;
template <class Mesh, class T, class ETag> class Field;
template <int Dim, class T, class ETag> class Engine;
struct CompressibleBrick;

//-----------------------------------------------------------------------------
// FileSetReader
//-----------------------------------------------------------------------------

/**
 * FileSetReader<Dim> manages the reading of Arrays and Fields from
 * "DiscField" format files (as designed for Pooma r1 by Bill
 * Humphrey). It uses a DiskMeta object to read the fileset's .meta
 * file and a DiskLayout object to read the .layout file. It then
 * reads the .offset and .data files to extract the data for each
 * field.
 *
 * NOTE: We do not currently provide random access to the fields
 * within a record, as r1 did. Rather, reading an Array just advances
 * the field counter until fieldsPerRecord is hit, and then it
 * advances the record counter. Reading a Field assumes that there is
 * enough room remaining in the record to get all the "sub-Fields" -
 * i.e. a Pooma Field can't span multiple records.
 */

template <int Dim>
class FileSetReader
{
public:

  //---------------------------------------------------------------------------
  // Typedefs
  //---------------------------------------------------------------------------

  typedef Interval<Dim> Domain_t;

  //---------------------------------------------------------------------------
  // Constructors and destructor
  //---------------------------------------------------------------------------

  // The constructor takes a fileset base-name, "base", and
  // initializes the various member data. It does not open any files -
  // see open below. 
  
  // NOTE: We probably need to add a constructor that takes a
  // DiskConfig object, or some such thing, or switch over to using
  // DiskConfig internally.

  FileSetReader(const char *filesetname);

  ~FileSetReader();

  //---------------------------------------------------------------------------
  // FileSetReader interface
  //---------------------------------------------------------------------------

  // Open the file - return false if there is a problem
  // If abortOnError is true, errors cause a PInsist instead, which
  // will either throw an exception or abort depending on the
  // configuration (whether exceptions are enabled).

  bool open(bool abortOnError = false);

  // Read the next Field or Array in the file.

  template <class T, class ETag>
  bool read(Array<Dim,T,ETag> &a);

  template <class Mesh, class T, class ETag>
  bool read(Field<Mesh,T,ETag> &f);

  //---------------------------------------------------------------------------
  // Accessors to various state
  //---------------------------------------------------------------------------

  // Global domain. This is read from the .meta file and broadcast
  // to all nodes. 

  const Domain_t &domain() const { return domain_m; }

  // Access to the DiskLayout object

  const DiskLayout<Dim> &diskLayout() const { return diskLayout_m; }

  // Next record and field to be read

  int nextRecord() const { return currentRecord_m; }
  int nextField() const { return currentField_m; }

  // Access to error message, if open returned false

  std::string errorMessage() const {  return errorMsg_m; }

  // Do we have to reverse the bytes? 

  bool bytesReversed() const { return bytesReversed_m; }

  // Access to the DiskMeta object (IO context only - all others return 0).

  const DiskMeta *diskMeta() const { return pmetafile_m; }

private:
  
  //---------------------------------------------------------------------------
  // Data
  //---------------------------------------------------------------------------

  // Fileset base name

  std::string basename_m;

  // Error message, if error happens.

  std::string errorMsg_m;

  // File streams for reading and writing.

  std::ifstream foffset_m;
  std::ifstream fdata_m;

  // DiskLayout reader

  DiskLayout<Dim> diskLayout_m;

  // DiskMeta reader (only exists on IO context)

  DiskMeta *pmetafile_m;

  // Global domain

  Domain_t domain_m;

  // Info about where we are in the file

  int currentRecord_m;
  int currentField_m;

  int numRecords_m;
  int fieldsPerRecord_m;

  // Local context and IO context

  int mycontext_m;
  int iocontext_m;

  // Do we need to reverse the bytes? 

  bool bytesReversed_m;

  //---------------------------------------------------------------------------
  // Offset data structure
  //---------------------------------------------------------------------------

  // Local structure for reading from the OffsetData file.

  // If we were build with long long defined, use long long for the
  // type of the offset information in the .offset file.

  // Using long long is wacky since all of the normal IO functions
  // (read, fread, fstream::read) are all going to convert it to
  // size_t, which is probably not long long, the latter being
  // non-standard. Unless IRIX 6.5 defined long to be 64 bits and you
  // have to use long long to read the files on 32 bit machines???

  // This is really a bad way to go about things since it requires
  // that we know at compile time what the format of our input files
  // are. It would be nice to be able to determine this dynamically!!!

  // I believe we can actually figure out the length of Offset_t from
  // the length of the file and implement a version that will do the
  // right thing for any file. TODO

public:

  typedef POOMA_INT64 Offset_t;

private:

  // Struct used to read the .offset records.

  template <class T>
  struct OffsetData
  {
    void reverseBytes();

    int nodedata[6*Dim];  // domain data (same format as .layout)
    union {
      bool isCompressed;    // Is the data compressed
      char pad[8];
    } u;
    Offset_t offset;      // offset in sizeof(T) units
    T compressedValue;    // Data value, if compressed
  };

  //---------------------------------------------------------------------------
  // Utility functions
  //---------------------------------------------------------------------------

  // Try to cut down on code bloat by performing tasks that don't
  // really care about T, ETag, Mesh, or whatever, in separate
  // functions. Some of this functionality could even be moved into a
  // non-template base-class.

  // Read and check the metafile info (iocontext only)

  bool readMeta();

  // Read the field ID from the .offset file.

  int readFieldID();

  // Unfortunately, this function has to know about T in order to do
  // the byte reversal. 

  template <class T>
  bool readOffsetData(OffsetData<T> &odata, Interval<Dim> &dom);

  // Read the field data into a compressible brick engine of the
  // appropriate type. 

  template <class T>
  bool readFieldData(const DiskNode<Dim> &dnode,
                     const Engine<Dim, T, CompressibleBrick> &e, 
                     OffsetData<T> &odata);

  // Read a subfield - simplifies the logic of read(Field_t&)

  template <class Mesh, class T, class ETag>
  bool readSubField(Field<Mesh,T,ETag> &sf);
};

// Include .cpp file to get out-of-line functions.

#include "IO/FileSetReader.cpp"

#endif // POOMA_IO_FILESETREADER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetReader.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
