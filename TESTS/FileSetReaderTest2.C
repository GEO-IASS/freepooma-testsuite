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
// FileSetReaderTest2: Test the FileSetReader functionality for a single
//   fileset, including the ability to dynamically detect the need to
//   fix byte ordering. Same as FileSetReaderTest1, except that it
//   reads the data into a multipatch array. 
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
#include "Engine/MultiPatchEngine.h"
#include "Layout/GridLayout.h"

#include "Field/Field.h"
#include "Field/FieldCentering.h"
#include "Field/Mesh/UniformRectilinearMesh.h"

// Why???

#include "Partition/BisectionMapper.h"
#include "Partition/ContiguousMapper.h"

//-----------------------------------------------------------------------------
// Main program
//  -- run the test
//-----------------------------------------------------------------------------

typedef Remote<Brick>                     PatchTag_t;
typedef MultiPatch<GridTag, PatchTag_t>   MP_t;
typedef Array<3, double, MP_t>            Array_t;
typedef DomainLayout<3>                   Layout_t;
typedef UniformRectilinearMesh<3>         Mesh_t;
typedef Field<Mesh_t, double, MP_t>       Field_t;

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

  // Create two layouts, each with two patches, but different
  // partitionings. (a is for array, f is for field).

  Loc<3> ablocks;
  int d;
  for (d = 0; d < 3; d++) 
    {
      ablocks[d] = Loc<1>(d == 2 ? 2 : 1);
    }

  GridLayout<3> alayout(dom, ablocks, DistributedTag()); 

  // Create an array to read the first field

  Array_t a(alayout);

  // Read the array.

  reader.read(a);

  tester.check(success);

  tester.out() << "a = \n" << a << std::endl;
  
  // Now construct a simple Field and read the second record with it.

  // First create a layout with some guards . . .

  Loc<3> fblocks;
  for (d = 0; d < 3; d++) 
    {
      fblocks[d] = Loc<1>(d == 1 ? 2 : 1);
    }

  GridLayout<3> flayout(dom, fblocks, GuardLayers<3>(2), DistributedTag());

  // . . . and some centerings . . .
  
  Centering<3> vert = canonicalCentering<3>(VertexType, Continuous);

  // . . . and finally a Field

  Field_t f(vert, flayout);

  // No read the next record . . .

  success = reader.read(f);

  tester.check(success);

  // . . . and print it out the field

  tester.out() << "f = \n" << f << std::endl;
  
  // Try it again - this should fail.

  success = reader.read(f);
  tester.check(!success);

  int ret = tester.results("FileSetReaderTest2");
  
  Pooma::finalize();

  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetReaderTest2.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
