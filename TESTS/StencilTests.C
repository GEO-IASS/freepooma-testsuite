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
// Test the use of some field stencils.
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Pooma/Fields.h"

#include "Field/DiffOps/Div.h"
#include "Field/DiffOps/Div.UR.h"

#include <iostream>
#include <cmath>

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
#endif

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<2> physicalVertexDomain(10, 10);
  
  Loc<2> blocks(2, 2);
  UniformGridPartition<2> partition(blocks, GuardLayers<2>(1));   
  UniformGridLayout<2> layout(physicalVertexDomain, partition,
			      LayoutTag_t());
  
  // Now, we can declare fields.

  Centering<2> cell = canonicalCentering<2>(CellType, Continuous);
  Centering<2> vertex = canonicalCentering<2>(VertexType, Continuous);
  Centering<2> allFace = canonicalCentering<2>(FaceType, Continuous);

  typedef UniformRectilinearMesh<2> Geometry_t;

  typedef 
    Field<Geometry_t, double, MultiPatch<UniformTag, BrickTag_t> >
    Field_t;

  typedef 
    Field<Geometry_t, Vector<2>, MultiPatch<UniformTag, BrickTag_t> >
    VField_t;

  Vector<2> origin(0.0, 0.0);
  Vector<2> spacings(1.0, 1.0);

  VField_t vfield(vertex, layout, origin, spacings);
  Field_t cfield(cell, layout, origin, spacings);
  Field_t facefield(allFace, layout, origin, spacings);

  vfield = positions(vfield);

  tester.out() << vfield << std::endl;

  cfield = divVertToCell(vfield);

  tester.out() << cfield << std::endl;

  tester.check("divergence is 2", sum(cfield -2) == 0);

  int ret = tester.results("StencilTests");
  Pooma::finalize();
  return ret;
}




// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: StencilTests.cpp,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
