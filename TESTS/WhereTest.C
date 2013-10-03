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
// Test the uses of where() with the new field.
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
  
  Loc<2> blocks(2, 2);
  UniformGridPartition<2> partition(blocks, GuardLayers<2>(1));   
  UniformGridLayout<2> layout(physicalVertexDomain, partition,
			      LayoutTag_t());
  
  tester.out() << "layout domain: " << layout.domain() << std::endl;
  
  // Now, we can declare a field.

  Centering<2> allFace = canonicalCentering<2>(FaceType, Continuous);
  Centering<2> allCell = canonicalCentering<2>(CellType, Continuous);

  typedef UniformRectilinearMesh<2> Geometry_t;

  typedef 
    Field<Geometry_t, double, MultiPatch<UniformTag, BrickTag_t> >
    Field_t;

  typedef 
    Field<Geometry_t, Vector<2>, MultiPatch<UniformTag, BrickTag_t> >
    VField_t;

  Vector<2> origin(0.0, 0.0);
  Vector<2> spacings(1.0, 1.0);

  Field_t a(allFace, layout, origin, spacings);
  Field_t b(allFace, layout, origin, spacings);
  Field_t c(allFace, layout, origin, spacings);
  Field_t d(allCell, layout, origin, spacings);
  Field_t e(allCell, layout, origin, spacings);
  Field_t f(allCell, layout, origin, spacings);

  PositionsTraits<Geometry_t>::Type_t x = positions(a);

  b = 0.0;
  c = 0.0;

  Vector<2> line(1.0, 1.0);


  // 3-arg where

  a = where(dot(x, line) > 8.0, x.comp(0), x.comp(1));

  // equivalent to:
  //  a = where(x.comp(0) + x.comp(1) > 8.0, x.comp(0), x.comp(1));

  tester.out() << "where(dot(x, line) > 8.0, x.comp(0), x.comp(1))\n" << a << std::endl;

  // Should verify these results by hand.
  // These are basically regression tests.

  tester.check("sum a[0]", sum(a[0]), 423.0);
  tester.check("sum a[0]*x[0](0)", sum(a[0] * x[0].comp(0)), 2397.0);
  tester.check("sum a[0]*x[0](1)", sum(a[0] * x[0].comp(1)), 2083.5);
  tester.check("sum a[1]", sum(a[1]), 387.0);
  tester.check("sum a[1]*x[1](0)", sum(a[1] * x[1].comp(0)), 2161.5);
  tester.check("sum a[1]*x[1](1)", sum(a[1] * x[1].comp(1)), 1990.5);


  // 2-arg where

  b = where(dot(x, line) > 8.0, x.comp(0));
  c = where(dot(x, line) <= 8.0, x.comp(1));

  tester.out() << "where(dot(x, line) > 8.0, x.comp(0))" << b << std::endl;
  tester.out() << "where(dot(x, line) <= 8.0, x.comp(1))" << c << std::endl;

  // verify using 3-arg where verified above

  tester.check("twoarg where result 0.0 part, centering zero",
	       all(where(dot(x.subField(0, 0), line) > 8.0,
                   c.subField(0, 0), b.subField(0, 0)) == 0.0));
  tester.check("twoarg where result 0.0 part, centering one",
	       all(where(dot(x.subField(0, 1), line) > 8.0,
                   c.subField(0, 1), b.subField(0, 1)) == 0.0));
  tester.check("twoarg where result dirtied part, centering zero",
               all(where(dot(x.subField(0, 0), line) > 8.0,
                   b.subField(0, 0), c.subField(0, 0)) == a.subField(0, 0)));
  tester.check("twoarg where result dirtied part, centering one",
               all(where(dot(x.subField(0, 1), line) > 8.0,
                   b.subField(0, 1), c.subField(0, 1)) == a.subField(0, 1)));

  // 2-arg where reduction

  d = 1.0;
  e = positions(e).read(e.physicalDomain()).comp(0);
  tester.check("reduction over twoarg where",
	       sum(where(e(e.physicalDomain()) < 4.0, d)) == 4.0*9.0);

  // 3-arg where reduction

  d = 1.0;
  f = 0.0;
  e = positions(e).read(e.physicalDomain()).comp(0);
  tester.check("reduction over twoarg where",
	       sum(where(e(e.physicalDomain()) < 4.0, d, f)) == 4.0*9.0);

  // 2-arg where with scalar expression and reduction variant thereof

  d = where(e(e.physicalDomain()) >= 4.0, 0.0);
  tester.check("counting reduction",
	       sum(where(d(d.physicalDomain()) != 0.0, 1)) == 4*9);

  // 2-arg where with scalar test and reduction variant thereof

  d = where(true, f);
  tester.check("simple where", all(d(d.physicalDomain()) == 0.0));
  tester.check("simple where reduction", prod(where(true, d)) == 0.0);
 
  // note that where with both expression and test being scalar does not
  // work because of
  // src/Pooma/PETE/ExpressionTraits.h:121:
  // error: no type named `Type_t' in
  //        `struct CombineExpressionTraits<ExpressionIsScalar, ExpressionIsScalar>'
  // and that is probably not the only reason.
  //
  // d = where(true, 1.0);
  // tester.check("even more simple where", all(d(d.physicalDomain()) == 1.0));

  int ret = tester.results("WhereTest");
  Pooma::finalize();
  return ret;
}




// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: WhereTest.cpp,v $   $Author: richi $
// $Revision: 1.6 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
