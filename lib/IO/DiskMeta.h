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
 * This .meta format is based on the Pooma r1 "DiscField" .meta
 * format, designed by Bill Humphrey and others.
 */

#ifndef POOMA_IO_DISKMETA_H
#define POOMA_IO_DISKMETA_H

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"  // field domain information

#include <string>
#include <vector>
#include <fstream>

/**
 * DiskMeta - reads in information in an r1 "DiscField" meta file,
 * which is a file in the format
 * <PRE>
 *   # comment line
 *   keyword [=] value
 *   keyword [=] value
 *   ...
 * </PRE>
 *
 * An optional = may separate the single-word keyword and the value.  This
 * class parses the file, storing the state, and allowing access to
 * that state via the various accessors. 
 * 
 * NOTE: The DiskMeta class is only useful on an "IO Context". The
 * mutators check and are protected so that they'll only run on an IO
 * context, but the accessors are not. 
 */

class DiskMeta
{
public:

  //---------------------------------------------------------------------------
  // Constructor and destructor
  //---------------------------------------------------------------------------

  DiskMeta(const char *basename);

  ~DiskMeta();

  //---------------------------------------------------------------------------
  // Mutators
  //---------------------------------------------------------------------------

  // Open the .meta file - return false if the file can't be open.

  bool open(bool abortOnError = false);

  // Read and parse the .meta file. 
  //   If abortOnError is true, read will PInsist on all errors.
  //   If abortOnError is false, read will return false on error. 

  bool read(bool abortOnError = false);

  //---------------------------------------------------------------------------
  // Accessor functions
  //---------------------------------------------------------------------------

  // Return the .meta filename
  
  const std::string &filename() const { return filename_m; }

  // Return the type string

  const std::string &type() const { return type_m; }

  // Return the dimension of the stored field

  int dimension() const { return dim_m; }

  // Return the domain of the stored fields. 
  // For complex centerings, this is the vert domain.

  const Interval<1> &domain(int d) const;

  // Return the number of fields in each record.

  int fieldsPerRecord() const { return fieldsPerRecord_m; }

  // Return the number of records in the fileset.

  int numRecords() const { return numRecords_m; }

  // Return the number of filesets used to store the field. 

  int numFileSets() const { return numFileSets_m; }

  // Return a list of the number of patches in each record of this
  // fileset. 

  const std::vector<int> &patchesPerRecord() const
  { return patchesPerRecord_m; }

  // Return a list of the number of patches written in previous
  // records. 

  const std::vector<int> &patchTally() const
  { return patchTally_m; }

  // Return the error message, if read or open fail.

  std::string errorMessage() const { return errorMsg_m; }

private:

  // My context

  int mycontext_m;

  // IO context

  int iocontext_m;

  // The name of the meta file

  std::string filename_m;

  // Error string

  std::string errorMsg_m;

  // istream for reading file

  std::ifstream fin_m;
  
  // The type field from the file

  std::string type_m;

  // The dimension of the field

  int dim_m;

  // The domain of the field

  Interval<1> domain_m[7];

  // The number of fields per record

  int fieldsPerRecord_m;

  // The number of records

  int numRecords_m;

  // The number of filesets (aka SMPs)

  int numFileSets_m;

  // A list of the number of patches in each record.

  std::vector<int> patchesPerRecord_m;

  // A list of the number of patches written in previous records

  std::vector<int> patchTally_m;

};

#endif // POOMA_IO_DISKMETA_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DiskMeta.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
