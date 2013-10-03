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
// A tour of the new Field class.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
#endif

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  // To declare a field, you first need to set up a layout. This requires
  // knowing the physical vertex-domain and the number of external guard
  // cell layers. Vertex domains contain enough points to hold all of the
  // rectilinear centerings that POOMA is likely to support for quite
  // awhile. Also, it means that the same layout can be used for all
  // fields, regardless of centering.
  
  Interval<2> physicalVertexDomain(4, 4);  // 0..3 x 0..3
  
  Loc<2> blocks(1, 2);  // x-direction has one block, y-dir has two blocks
  UniformGridPartition<2> partition(blocks, GuardLayers<2>(1));   // add one layer of guard cells in each direction
  UniformGridLayout<2> layout(physicalVertexDomain, partition, LayoutTag_t());
  
  std::cout << layout << std::endl;
  std::cout << layout.domain() << std::endl;
  
  // Now, we can declare a field.

  Centering<2> allFace = canonicalCentering<2>(FaceType, Continuous);

  typedef 
    Field<UniformRectilinearMesh<2>, double,
    MultiPatch<UniformTag, BrickTag_t> > Field_t;
  Field_t f(allFace, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));

  std::cout << f.centering() << std::endl;
  std::cout << f[0].centering() << std::endl;
  std::cout << f[1].centering() << std::endl;

  // Ask for the field's physical cell domain.
  
  std::cout << f.physicalCellDomain() << std::endl;

  // If we ask for the physical domain, we should get the physical cell
  // domain back because of the all-face centering. We can get the
  // face-domains by specifying the sub-fields.

  std::cout << f.physicalDomain() << std::endl;  // cell orientation
  std::cout << f.physicalDomain(0) << std::endl; // x face orientation
  std::cout << f.physicalDomain(1) << std::endl; // y face orientation

  // Total domains work similarly.

  std::cout << f.totalDomain() << std::endl;
  std::cout << f.totalDomain(0) << std::endl;
  std::cout << f.totalDomain(1) << std::endl;

  // We can do a similar sort of thing by taking sub-field views.

  std::cout << f[0].physicalDomain() << std::endl; // x faces
  std::cout << f[1].physicalDomain() << std::endl; // y faces

  // Total domains work similarly. Note: taking a sub-field view doesn't
  // remove the guard layers.

  std::cout << f[0].totalDomain() << std::endl;
  std::cout << f[1].totalDomain() << std::endl;
    
  // We can actually index fields after taking a sub-field view. The
  // indices refer to the actual domain.
  
  f[0](1, 2) = 3.0;
  f[1](1, 2) = f[0](1, 2) + 1.2;
  
  std::cout << f[0](1, 2) << std::endl;
  std::cout << f[1](1, 2) << std::endl;

  // Same thing after taking domain & sub-field views.
  
  Interval<1> I(1,2);
  std::cout << f[0](I, I)(0, 1) << std::endl;
  std::cout << f(I, I)[1](0, 1) << std::endl;
  
  // The guard layers are removed when you take a domain view.

  std::cout << f(I, I).physicalDomain() << std::endl;
  std::cout << f(I, I).totalDomain() << std::endl;
  std::cout << f(I, I).physicalDomain(0) << std::endl;
  std::cout << f(I, I).totalDomain(0) << std::endl;
  
  // Check assignment of a scalar.
  
  f = -1.0;  // assign physical domain
  f(I, I) = -2.0;
  std::cout << f << std::endl;

  // Declare another field. Note how we can reuse the layout for a field
  // with a different centering.

  Centering<2> face1 = canonicalCentering<2>(FaceType, Continuous, YDim);
  
  Field_t g(face1, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));

  g = -3.0;
  g(I, I) = -4.0;
  f[1] = g;

  std::cout << f.all() << std::endl;
  std::cout << g.all() << std::endl;

  /*
  typedef 
    Field<Lagrangian<2>, double,
    MultiPatch<UniformTag, BrickTag_t> > 
     LagrField_t;
   LagrField_t h(allFace, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
  
  std::cout << h.fieldEngine().vertexPositions() << std::endl;
  
  // Try assigning to a field with a Lagrangian mesh.
  
  h.all() = 3.0; // .all means also set guards as well as physical
  h = -6.0;
  std::cout << h.all() << std::endl;
  */
   
  // Check out the patch function.
  
  f.all() = 1.0;
  f = 2.0;

  f[0](1, 1) = 3.0;
  f[1](1, 1) = 3.0;

  int nLocal = f[0].numPatchesLocal();
  std::cout << "context " << Pooma::context() << " has "
	    << nLocal << " patches" << std::endl;
  if (nLocal > 0)
    {
      std::cout << "context " << Pooma::context()
		<< " local patch 0: " << f[0].patchLocal(0) << std::endl;
    }
  if (nLocal > 1)
    {
      std::cout << "context " << Pooma::context()
		<< " local patch 1: " << f[0].patchLocal(1) << std::endl;
    }

  // Play with relations.

  Pooma::addAllPosReflectFaceBC(f);

  std::cout << f.all() << std::endl;
  
  // Try to create a vector field.

  typedef 
    Field<UniformRectilinearMesh<2>, Vector<2>,
    MultiPatch<UniformTag, BrickTag_t> > 
     VectorField_t;
  VectorField_t l(allFace, layout, Vector<2>(0.0), Vector<2>(1.0, 2.0));
   
  l.all() = Vector<2>(-1.0, 2.0);
  l = Vector<2>(4.0, 6.0);

  std::cout << l.all().comp(0) << std::endl;

  Pooma::finalize();
  return 0; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldTour1.cpp,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
