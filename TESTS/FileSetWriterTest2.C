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

#include "IO/FileSetWriter.h"
#include "IO/FileSetReader.h"

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

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
typedef Field<Mesh_t, Vector<dim,double>, MP_t> Field_t;

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
  
  typedef Vector<dim, double> Vector_t;

  Vector_t origin;
  Vector_t spacings;
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
  
  Field_t x(vert, layout, origin, spacings);
        
  Pooma::blockAndEvaluate();
  
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 5; j++)
      for (int k = 0; k < 6; k++)
        {
          x(i, j, k) = origin + 
            Vector_t(i*spacings(0), j*spacings(1), k*spacings(2));
        }
        
  {
    FileSetWriter<dim> w("xset", 1);
  
    w.write(x);
  }

  {
    typedef Array<dim, Vector_t, Brick> Array_t;
    Array_t a(physicalVertexDomain);

    FileSetReader<dim> w("xset");
    bool success = w.open();
    tester.check(success);

    if (success)
      {
        w.read(a);

        Pooma::blockAndEvaluate();

        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 5; j++)
            for (int k = 0; k < 6; k++)
              {
                double 
                  dx = i * spacings(0), 
                  dy = j * spacings(1),
                  dz = k * spacings(2);

                tester.check(a(i, j, k) == origin + Vector_t(dx,dy,dz));
              }

        tester.out() << "a = \n" << a << std::endl;
      }
  }    

  int ret = tester.results("FileSetWriterTest2");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetWriterTest2.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
