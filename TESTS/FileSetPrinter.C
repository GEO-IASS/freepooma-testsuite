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
// FileSetPrinter: ASCII dump of 3D double file sets. 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

// Class to test

#include "IO/FileSetReader.h"

// Basic infrastructure

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"

// Pooma stuff

#include "Domain/Interval.h"
#include "Array/Array.h"
#include "Engine/BrickEngine.h"
#include "Engine/RemoteEngine.h"

// Standard includes

#include <iostream>

//-----------------------------------------------------------------------------
// Code that does the actual work
//-----------------------------------------------------------------------------

void printFileSet(const char *basename)
{
  typedef Remote<Brick>                   PatchTag_t;
  typedef Array<3, double, PatchTag_t>    Array_t;

  // Output

  Inform pout;

  // Open the files for reading

  FileSetReader<3> reader(basename);

  // Error checking...

  bool success = reader.open();
  if (!success) 
    {
      pout << "Couldn't open file set" << std::endl;
      exit(1);
    }

  // Unfortunately, there is no way to easily check that the data type
  // in the file is double. Sigh.

  
  // Print the file set...

  pout << "Reading fileset " << basename << "\n";
  if (reader.bytesReversed()) 
    pout << "File set has bytes reversed\n";

  const Interval<3> &dom = reader.domain();
  
  Array_t a(dom);

  int rec = reader.nextRecord();
  int field = reader.nextField();

  while (reader.read(a))
    {
      if (rec == 0 && field == 0) 
        pout << "Global domain = " << dom << std::endl;

      pout << 
        "=============================================================="
        "========\n\n";
      pout << "Record " << rec << "; Field " << field << "\n";

      pout << a << std::endl;

      rec = reader.nextRecord();
      field = reader.nextField();

      a = 0;
    }

  pout << 
    "=============================================================="
    "========" << std::endl;
}

//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  if (argc == 2 && *argv[1] != '-')
    {
      printFileSet(argv[1]);
    }
  else
    {
      Pooma::Tester tester(argc, argv);
      tester.out() << "Usage: FileSetPrinter basename" << std::endl;
      int ret = tester.results("FileSetPrinter");
      return ret;
    }

  Pooma::finalize();

  return 0;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetPrinter.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
