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
// FileSetWriter operations test.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

#include "IO/FileSetWriter.h"

//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

const int dim = 3;

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
#endif

typedef UniformRectilinearMesh<dim> Mesh_t;
typedef MultiPatch<GridTag, BrickTag_t> MP_t;
typedef Field<Mesh_t, double, MP_t> Field_t;
typedef Array<dim, double, MP_t> Array_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int d;
  Interval<dim> physicalVertexDomain;
  
  for (d = 0; d < dim; d++) 
    {
      physicalVertexDomain[d] = Interval<1>(d + 4);
    }

  // Set up some arrays.
  
  Vector<dim, double> origin;
  Vector<dim, double> spacings;
  Loc<dim> blocks;
  for (d = 0; d < dim; d++) 
    {
      origin(d) = d;
      spacings(d) = d + 1;
      blocks[d] = Loc<1>(d == 2 ? 1 : 2);
    }
      
  // Make the layout.

  GridLayout<dim> layout(physicalVertexDomain, blocks, 
    GuardLayers<dim>(2), LayoutTag_t());

  // ..and some centerings..
  
  Centering<dim> vert = canonicalCentering<dim>(VertexType, Continuous);

  // ..and, finally, a field and an array.
  
  Field_t f(vert, layout, origin, spacings);
  Array_t a(layout);
        
  Pooma::blockAndEvaluate();
  
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 5; j++)
      for (int k = 0; k < 6; k++)
        {
          f(i, j, k) = i + j + k;
          a(i, j, k) = i + j + k;
        }
        
  Pooma::blockAndEvaluate();
  
  FileSetWriter<dim> w("fset", 2);
  
  w.write(f);
  w.write(a);
  
  int ret = tester.results("FileSetWriter");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetWriterTest1.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
