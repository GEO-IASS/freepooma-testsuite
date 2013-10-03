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
// Some simple tests of evaluation using message passing with the new field.
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

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<2> physicalVertexDomain(10, 10);

  // Note, the layout uses the DistributedTag, since we are using
  // Remote engines in the fields.
  
  Loc<2> blocks(2, 2);
  UniformGridPartition<2> partition(blocks, GuardLayers<2>(1));   
  UniformGridLayout<2> layout(physicalVertexDomain, partition,
			      DistributedTag());
  
  tester.out() << "layout domain: " << layout.domain() << std::endl;
  
  // Now, we can declare a field.

  Centering<2> allFace = canonicalCentering<2>(FaceType, Continuous);

  typedef UniformRectilinearMesh<2> Mesh_t;
  typedef MultiPatch<UniformTag, Remote<Brick> > EngineTag_t;

  typedef Field<Mesh_t, double, EngineTag_t > Field_t;

  typedef Field<Mesh_t, Vector<2>, EngineTag_t > VField_t;

  Vector<2> origin(0.0, 0.0);
  Vector<2> spacings(1.0, 1.0);

  Field_t a(allFace, layout, origin, spacings);
  Field_t b(allFace, layout, origin, spacings);
  Field_t c(allFace, layout, origin, spacings);

  // Should really figure out how to repackage these three lines:

  DomainLayout<2> layoutDom(physicalVertexDomain, GuardLayers<2>(1));
  XField<Mesh_t>::Type_t x(allFace, layoutDom, origin, spacings);
  setXField(x);

  b = 0.0;
  c = 0.0;

  Vector<2> line(1.0, 1.0);

  a = where(dot(x, line) > 8.0, x.comp(0), x.comp(1));

  tester.out() << a << std::endl;

  tester.check("sum a[0]", sum(a[0]), 423.0);
  tester.check("sum a[0]*x[0](0)", sum(a[0] * x[0].comp(0)), 2397.0);
  tester.check("sum a[0]*x[0](1)", sum(a[0] * x[0].comp(1)), 2083.5);
  tester.check("sum a[1]", sum(a[1]), 387.0);
  tester.check("sum a[1]*x[1](0)", sum(a[1] * x[1].comp(0)), 2161.5);
  tester.check("sum a[1]*x[1](1)", sum(a[1] * x[1].comp(1)), 1990.5);

  int ret = tester.results("CrossBox");
  Pooma::finalize();
  return ret;
}




// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CrossBox.cpp,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
