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
// Test of the patchLocal function and demonstration of how to use it to
// perform local SPMD computations on fields.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
  typedef Remote<CompressibleBrick> CompressibleBrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
  typedef CompressibleBrick CompressibleBrickTag_t;
#endif

// perform some nonsense on some memory containing doubles.

void nonsense(double *data, int size)
{
  for (int i = 0; i < size; ++i)
  {
    data[i] += i;
  }
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // To declare a field, you first need to set up a layout. This requires
  // knowing the physical vertex-domain and the number of external guard
  // cell layers. Vertex domains contain enough points to hold all of the
  // rectilinear centerings that POOMA is likely to support for quite
  // awhile. Also, it means that the same layout can be used for all
  // fields, regardless of centering.
  
  Interval<2> physicalVertexDomain(14, 14);
  Loc<2> blocks(3, 3);
  GridLayout<2> layout1(physicalVertexDomain, blocks, GuardLayers<2>(1),
                        LayoutTag_t());
  GridLayout<2> layout0(physicalVertexDomain, blocks, GuardLayers<2>(0),
                        LayoutTag_t());

  Centering<2> cell = canonicalCentering<2>(CellType, Continuous, AllDim);
  Centering<2> vert = canonicalCentering<2>(VertexType, Continuous, AllDim);
  Centering<2> yedge = canonicalCentering<2>(EdgeType, Continuous, YDim);

  Vector<2> origin(0.0);
  Vector<2> spacings(1.0, 2.0);

  // First basic test verifies that we're assigning to the correct areas
  // on a brick.

  typedef
    Field<UniformRectilinearMesh<2>, double,
    MultiPatch<GridTag, BrickTag_t> > Field_t;
  Field_t b0(cell, layout1, origin, spacings);
  Field_t b1(vert, layout1, origin, spacings);
  Field_t b2(yedge, layout1, origin, spacings);
  Field_t b3(yedge, layout1, origin, spacings);
  Field_t bb0(cell, layout0, origin, spacings);
  Field_t bb1(vert, layout0, origin, spacings);
  Field_t bb2(yedge, layout0, origin, spacings);

  b0.all() = 0.0;
  b1.all() = 0.0;
  b2.all() = 0.0;

  b0 = 1.0;
  b1 = 1.0;
  b2 = 1.0;

  bb0.all() = 0.0;
  bb1.all() = 0.0;
  bb2.all() = 0.0;

  bb0 = 1.0;
  bb1 = 1.0;
  bb2 = 1.0;

  // SPMD code follows.
  // Note, SPMD code will work with the evaluator if you are careful
  // to perform assignment on all the relevant contexts.  The patchLocal
  // function creates a brick on the local context, so you can just perform
  // the assignment on that context.

  int i;

  for (i = 0; i < b0.numPatchesLocal(); ++i)
  {
    Patch<Field_t>::Type_t patch = b0.patchLocal(i);
    //    tester.out() << "context " << Pooma::context() << ":  assigning to patch " << i
    //              << " with domain " << patch.domain() << std::endl;
    patch += 1.5;
  }

  // This is safe to do since b1 and b2 are built with the same layout.
  for (i = 0; i < b1.numPatchesLocal(); ++i)
  {
    b1.patchLocal(i) += 1.5;
    b2.patchLocal(i) += 1.5;
  }

  for (i = 0; i < bb0.numPatchesLocal(); ++i)
  {
    Patch<Field_t>::Type_t patch = bb0.patchLocal(i);
    //    tester.out() << "context " << Pooma::context() << ":  assigning to patch on bb0 " << i
    //              << " with domain " << patch.domain() << std::endl;
    patch += 1.5;
  }

  // This is safe to do since bb1 and bb2 are built with the same layout.
  for (i = 0; i < bb1.numPatchesLocal(); ++i)
  {
    bb1.patchLocal(i) += 1.5;
    bb2.patchLocal(i) += 1.5;
  }

  tester.check("cell centered field is 2.5", all(b0 == 2.5));
  tester.check("vert centered field is 2.5", all(b1 == 2.5));
  tester.check("edge centered field is 2.5", all(b2 == 2.5));

  tester.out() << "b0.all():" << std::endl << b0.all() << std::endl;
  tester.out() << "b1.all():" << std::endl << b1.all() << std::endl;
  tester.out() << "b2.all():" << std::endl << b2.all() << std::endl;

  tester.check("didn't write into b0 boundary",
               sum(b0.all()) == 2.5 * b0.physicalDomain().size());
  tester.check("didn't write into b1 boundary",
               sum(b1.all()) == 2.5 * b1.physicalDomain().size());
  tester.check("didn't write into b2 boundary",
               sum(b2.all()) == 2.5 * b2.physicalDomain().size());

  tester.check("cell centered field is 2.5", all(bb0 == 2.5));
  tester.check("vert centered field is 2.5", all(bb1 == 2.5));
  tester.check("edge centered field is 2.5", all(bb2 == 2.5));

  tester.out() << "bb0:" << std::endl << bb0 << std::endl;
  tester.out() << "bb1:" << std::endl << bb1 << std::endl;
  tester.out() << "bb2:" << std::endl << bb2 << std::endl;

  typedef
    Field<UniformRectilinearMesh<2>, double,
    MultiPatch<GridTag, CompressibleBrickTag_t> > CField_t;
  CField_t c0(cell, layout1, origin, spacings);
  CField_t c1(vert, layout1, origin, spacings);
  CField_t c2(yedge, layout1, origin, spacings);
  CField_t cb0(cell, layout0, origin, spacings);
  CField_t cb1(vert, layout0, origin, spacings);
  CField_t cb2(yedge, layout0, origin, spacings);

  c0.all() = 0.0;
  c1.all() = 0.0;
  c2.all() = 0.0;

  c0 = 1.0;
  c1 = 1.0;
  c2 = 1.0;

  cb0.all() = 0.0;
  cb1.all() = 0.0;
  cb2.all() = 0.0;

  cb0 = 1.0;
  cb1 = 1.0;
  cb2 = 1.0;

  // SPMD code follows.
  // Note, SPMD code will work with the evaluator if you are careful
  // to perform assignment on all the relevant contexts.  The patchLocal
  // function creates a brick on the local context, so you can just perform
  // the assignment on that context.

  for (i = 0; i < c0.numPatchesLocal(); ++i)
  {
    Patch<CField_t>::Type_t patch = c0.patchLocal(i);
    tester.out() << "context " << Pooma::context() << ":  assigning to patch " << i
                 << " with domain " << patch.domain() << std::endl;
    patch += 1.5;
  }

  // This is safe to do since c1 and c2 are built with the same layout.
  for (i = 0; i < c1.numPatchesLocal(); ++i)
  {
    c1.patchLocal(i) += 1.5;
    c2.patchLocal(i) += 1.5;
  }

  for (i = 0; i < cb0.numPatchesLocal(); ++i)
  {
    Patch<CField_t>::Type_t patch = cb0.patchLocal(i);
    tester.out() << "context " << Pooma::context() << ":  assigning to patch on cb0 " << i
                 << " with domain " << patch.domain() << std::endl;
    patch += 1.5;
  }

  // This is safe to do since cb1 and cb2 are cuilt with the same layout.
  for (i = 0; i < cb1.numPatchesLocal(); ++i)
  {
    cb1.patchLocal(i) += 1.5;
    cb2.patchLocal(i) += 1.5;
  }

  tester.check("cell centered field is 2.5", all(c0 == 2.5));
  tester.check("vert centered field is 2.5", all(c1 == 2.5));
  tester.check("edge centered field is 2.5", all(c2 == 2.5));

  tester.out() << "c0.all():" << std::endl << c0.all() << std::endl;
  tester.out() << "c1.all():" << std::endl << c1.all() << std::endl;
  tester.out() << "c2.all():" << std::endl << c2.all() << std::endl;

  tester.check("didn't write into c0 boundary",
               sum(c0.all()) == 2.5 * c0.physicalDomain().size());
  tester.check("didn't write into c1 boundary",
               sum(c1.all()) == 2.5 * c1.physicalDomain().size());
  tester.check("didn't write into c2 boundary",
               sum(c2.all()) == 2.5 * c2.physicalDomain().size());

  tester.check("cell centered field is 2.5", all(cb0 == 2.5));
  tester.check("vert centered field is 2.5", all(cb1 == 2.5));
  tester.check("edge centered field is 2.5", all(cb2 == 2.5));

  tester.out() << "cb0:" << std::endl << cb0 << std::endl;
  tester.out() << "cb1:" << std::endl << cb1 << std::endl;
  tester.out() << "cb2:" << std::endl << cb2 << std::endl;

  //------------------------------------------------------------------
  // Scalar code example:
  //

  c0 = iota(c0.domain()).comp(0);
  c1 = iota(c1.domain()).comp(1);

  // Make sure all the data-parallel are done:

  Pooma::blockAndEvaluate();

  for (i = 0; i < c0.numPatchesLocal(); ++i)
  {
    Patch<CField_t>::Type_t local0 = c0.patchLocal(i);
    Patch<CField_t>::Type_t local1 = c1.patchLocal(i);
    Patch<CField_t>::Type_t local2 = c2.patchLocal(i);

    Interval<2> domain = local2.domain();  // physical domain of local y-edges

    // --------------------------------------------------------------
    // I believe the following is probably the most efficient approach
    // for sparse computations.  For data-parallel computations, the
    // evaluator will uncompress the patches and take brick views, which
    // provide the most efficient access.  If you are only performing
    // the computation on a small portion of cells, then the gains would
    // be outweighed by the act of copying the compressed value to all the
    // cells.
    //
    // The read function is used on the right hand side, because
    // operator() is forced to uncompress the patch just in case you want
    // to write to it.

    for(Interval<2>::iterator pos = domain.begin(); pos != domain.end(); ++pos)
    {
      Loc<2> edge = *pos;
      Loc<2> rightCell = edge;  // cell to right is same cell
      Loc<2> leftCell = edge - Loc<2>(1,0);
      Loc<2> topVert = edge + Loc<2>(0, 1);
      Loc<2> bottomVert = edge;

      local2(edge) =
        local0.read(rightCell) + local0.read(leftCell) +
        local1.read(topVert) + local1.read(bottomVert);
    }

    // This statement is optional, it tries to compress the patch after
    // we're done computing on it.  Since I used .read() for the local0 and 1
    // they remained in their original state. compress() can be expensive, so
    // it may not be worth trying unless space is really important.

    compress(local2);
  }

  tester.out() << "c0" << std::endl << c0 << std::endl;
  tester.out() << "c1" << std::endl << c1 << std::endl;
  tester.out() << "c2" << std::endl << c2 << std::endl;

  //------------------------------------------------------------------
  // Interfacing with a c-function:
  //
  // This example handles the corner cases, where the patches from a
  // cell centered field with no guard layers actually contain some
  // extra data.

  Pooma::blockAndEvaluate();

  for (i = 0; i < cb0.numPatchesLocal(); ++i)
  {
    Patch<CField_t>::Type_t local0 = cb0.patchLocal(i);
    Interval<2> physicalDomain = local0.physicalDomain();
    double *data;
    int size = physicalDomain.size();

    if (physicalDomain == local0.totalDomain())
    {
      uncompress(local0);
      data = &local0(physicalDomain.firsts());
      nonsense(data, size);
    }
    else
    {
      // In this case, the engine has extra storage even though the
      // field has the right domain. We copy it to a brick engine,
      // call the function and copy it back.  No uncompress is required,
      // since the assignment will copy the compressed value into the
      // brick.

      // arrayView is a work-around.  Array = Field doesn't work at
      // the moment.

      Array<2, double, Brick> brick(physicalDomain);
      Array<2, double, CompressibleBrick> arrayView(local0.engine());
      brick = arrayView(physicalDomain);
      Pooma::blockAndEvaluate();
      data = &brick(Loc<2>(0));
      nonsense(data, size);
      arrayView(physicalDomain) = brick;

      // Note that we don't need a blockAndEvaluate here, since an iterate has
      // been spawned to perform the copy.
    }

    // If you want to try compress(local0) here, you should do blockAndEvaluate
    // first in case the local0 = brick hasn't been executed yet.
  }
      
  tester.out() << "cb0.all()" << std::endl << cb0 << std::endl;

  b2 = positions(b2).comp(0);

  RefCountedBlockPtr<double> block = pack(b2);

  // The following functions give you access to the raw data from pack.
  // Note that the lifetime of the data is managed by the RefCountedBlockPtr,
  // so when "block" goes out of scope, the data goes away.  (i.e. Don't write
  // a function where you return block.beginPointer().)

  double *start = block.beginPointer();  // start of the data
  double *end = block.endPointer();      // one past the end
  int size = block.size();               // size of the data

  tester.out() << Pooma::context() << ":" << block.size() << std::endl;

  unpack(b3, block);

  tester.out() << "b2" << std::endl << b2 << std::endl;
  tester.out() << "b3" << std::endl << b3 << std::endl;

  tester.check("pack, unpack", all(b2 == b3));

  int ret = tester.results("LocalPatch");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: LocalPatch.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
