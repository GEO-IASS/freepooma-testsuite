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
// FileSetReaderTest1: Test the FileSetReader functionality for a single
//   fileset, including the ability to dynamically detect the need to
//   fix byte ordering. 
// NOTE: You must run FileSetReaderTest0 before running this test.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

// Class to test

#include "IO/FileSetReader.h"

// Basic infrastructure

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"

// Pooma stuff needed for test

#include "Domain/Interval.h"

#include "Array/Array.h"
#include "Engine/BrickEngine.h"
#include "Engine/RemoteEngine.h"

#include "Field/Field.h"
#include "Field/FieldCentering.h"
#include "Field/Mesh/UniformRectilinearMesh.h"

//-----------------------------------------------------------------------------
// Main program
//  -- run the test
//-----------------------------------------------------------------------------

typedef Remote<Brick>                     PatchTag_t;
typedef Array<3, double, PatchTag_t>      Array_t;
typedef DomainLayout<3>                   Layout_t;
typedef UniformRectilinearMesh<3>         Mesh_t;
typedef Field<Mesh_t, double, PatchTag_t> Field_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Open the files for reading

  FileSetReader<3> reader("TestData");

  bool success = reader.open();
  tester.check(success);

  if (!success)
    {
      using Pooma::perr;
      perr << "You must create the TestData file set first.\n"
           << "This is done by running FileSetReaderTest0 and moving the\n" 
           << "file set to the location where this test will be run. " 
           << std::endl;
      return 1;
    }

  tester.out() << "Bytes are reversed? " 
               << (reader.bytesReversed() ? "yes" : "no") << std::endl;

  // Check the domain

  typedef Interval<1> Dom1_t;
  Interval<3> dom(Dom1_t(4),Dom1_t(5),Dom1_t(6));

  tester.check(reader.domain() == dom);

  // Check metafile information

  if (Pooma::context() == 0)
    {
      tester.check(reader.diskMeta()->numRecords() == 1);
      tester.check(reader.diskMeta()->fieldsPerRecord() == 2);
      tester.check(reader.diskMeta()->dimension() == 3);
    }

  // Create an array to read the first field

  Array_t a(dom);

  // Read the array.

  reader.read(a);

  tester.check(success);

  tester.out() << "a = \n" << a << std::endl;
  
  // Now construct a simple Field and read the second record with it.

  // First we need a layout . . .

  Layout_t layout(dom);

  // . . . and some centerings . . .
  
  Centering<3> vert = canonicalCentering<3>(VertexType, Continuous);

  // . . . and finally a Field

  Field_t f(vert, layout);

  // No read the next record . . .

  success = reader.read(f);

  tester.check(success);

  // . . . and print it out the field

  tester.out() << "f = \n" << f << std::endl;
  
  // Try it again - this should fail.

  success = reader.read(f);
  tester.check(!success);

  int ret = tester.results("FileSetReaderTest1");
  
  Pooma::finalize();

  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetReaderTest1.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
