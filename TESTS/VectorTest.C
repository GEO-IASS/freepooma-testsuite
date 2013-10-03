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
// Test some vector concepts with new field.
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
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Engine/UserFunction.h"
#include "Engine/Stencil.h"
#include "Tiny/Vector.h"
#include "Pooma/FunctorResult.h"
#include "Pooma/Fields.h"

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
  DomainLayout<2> layoutDom(physicalVertexDomain, GuardLayers<2>(1));
  
  Loc<2> blocks(2, 2);
  UniformGridPartition<2> partition(blocks, GuardLayers<2>(1));   
  UniformGridLayout<2> layout(physicalVertexDomain, partition,
  			      LayoutTag_t());
  
  tester.out() << "layout domain: " << layoutDom.domain() << std::endl;
  tester.out() << "layout domain: " << layout.domain() << std::endl;
  
  // Now, we can declare a field.

  Centering<2> allFace = canonicalCentering<2>(FaceType, Continuous);

  typedef UniformRectilinearMesh<2> Geometry_t;
  typedef MultiPatch<UniformTag, BrickTag_t> EngineTag_t;
  //  typedef BrickTag_t EngineTag_t;

  typedef Field<Geometry_t, double, EngineTag_t > Field_t;
  typedef Field<Geometry_t, Vector<2>, EngineTag_t > VField_t;

  Vector<2> origin(0.0, 0.0);
  Vector<2> spacings(1.0, 1.0);

  Field_t a(allFace, layout, origin, spacings);
  VField_t b(allFace, layout, origin, spacings);
  VField_t c(allFace, layout, origin, spacings);
  //  Field_t a(allFace, layoutDom, origin, spacings);
  //  VField_t b(allFace, layoutDom, origin, spacings);
  //  VField_t c(allFace, layoutDom, origin, spacings);

  b = positions(b);
  c = Vector<2>(1.0, 2.0);

  a = dot(b, c);

  tester.out() << a << std::endl;

  int ret = tester.results("VectorTest");
  Pooma::finalize();
  return ret;
}




// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: VectorTest.cpp,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
