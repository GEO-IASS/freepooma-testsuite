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
// Illustration of the ScalarCode evaluator.  This example computes the
// average value of the faces of a cell in a face centered field and stores
// the result in a cell centered field.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"
#include "Engine/RemoteEngine.h"
#include <cmath>
#include <vector>

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
#endif

//-----------------------------------------------------------------------------
// Example provided by John Hall
//-----------------------------------------------------------------------------

template<int Dim>
struct EdgeQPressureInfo
{
  void scalarCodeInfo(ScalarCodeInfo &info) const
  {
    info.arguments(5);

    info.write(0, true );
    info.write(1, false);
    info.write(2, false);
    info.write(3, false);
    info.write(4, false);

    // Does this operation index neighboring
    // cells?  (Do we need to update the internal
    // guard layers.)

    info.useGuards(0, false);
    info.useGuards(1, false);
    info.useGuards(2, false);
    info.useGuards(3, false);
    info.useGuards(4, false);

    info.dimensions(Dim);

    for (int i = 0; i < Dim; ++i)
    {
      info.lowerExtent(i) = 0;
      info.upperExtent(i) = 0;
    }
  }
};

typedef double Real;

template<int Dim>
struct ScalarEdgeQPressure
  : public EdgeQPressureInfo<Dim>
{
  ScalarEdgeQPressure(const Real& inLinearQ )
      : EdgeQPressureInfo<Dim>(),
    linearQ(inLinearQ)
  {
  }

  template<class F1, class F2, class F3, class F4, class F5>
  void operator()(const F1& EdgeQPressure,     const F2& EdgeGammaConstant,
		  const F3& EdgeSoundSpeed,    const F4& EdgeVelocity,
		  const F5& EdgePsiLimiter,    const 
		  Loc<Dim> &loc) const
       {
	 if( EdgePsiLimiter(loc) < 0.2 ) {  // epsilon
	   EdgeQPressure(loc) = 0.0;
	   return;
	 }

	 Real edgeVelocityMagnitude = 
	   sqrt(dot(EdgeVelocity(loc),EdgeVelocity(loc)));
	 
	 EdgeQPressure(loc) = edgeVelocityMagnitude * 
	   EdgePsiLimiter(loc) *
	   (EdgeGammaConstant(loc) * 
	    edgeVelocityMagnitude +
	    sqrt( linearQ * linearQ * 
		  EdgeSoundSpeed(loc) * EdgeSoundSpeed(loc) +
		  EdgeGammaConstant(loc) * 
		  EdgeGammaConstant(loc) *
		  edgeVelocityMagnitude * 
		  edgeVelocityMagnitude));
       }

private:
  Real linearQ;
};



//-----------------------------------------------------------------------------
// Test function
//-----------------------------------------------------------------------------

// This example averages the values from an all-face field and puts the
// result in a cell-centered field.
// The example is interesting for 2 reasons:
// -useGuards for the input field is true, because we will need to
// index into the guard layers---you need the same face on two different
// processors, so one of the values comes from the guard layer.
// -on the other hand, we've written the code to operate on views of
// the all-face field, so the extents are actually 0.  The view of the
// all-face field take the cell-based domain and include all the relevant
// faces.

struct AllFaceToCellInfo
{
  AllFaceToCellInfo(int dimensions)
    : dimensions_m(dimensions)
  {
  }

  void scalarCodeInfo(ScalarCodeInfo &info) const
  {
    info.arguments(2);

    info.write(0, true );
    info.write(1, false);
    info.useGuards(0, false);
    info.useGuards(1, true);

    info.dimensions(dimensions_m);

    int i;
    for (i = 0; i < dimensions_m; ++i)
    {
      info.lowerExtent(i) = 0;
      info.upperExtent(i) = 0;
    }
  }

  int dimensions_m;
};

template<int Dim>
struct AllFaceToCellAverage
  : public AllFaceToCellInfo
{
  AllFaceToCellAverage()
    : AllFaceToCellInfo(Dim)
  {
    int i;
    for (i = 0; i < Dim; ++i)
    {
      off_m[i] = Loc<Dim>(0);
      off_m[i][i] = 1;
    }
    factor_m = 1.0 / (2.0 * Dim);
  }

  template<class F1, class F2>
  inline
  void operator()(F1 &f1, const F2 &f2, const Loc<Dim> &loc) const
  {
    int i;
    store_m = 0.0;
    for (i = 0; i < Dim; ++i)
    {
      store_m += f2[i](loc) + f2[i](loc + off_m[i]);
    }
    f1(loc) = factor_m * store_m;
  }

  mutable double store_m;
  Loc<Dim> off_m[Dim];
  double factor_m;
};

//-----------------------------------------------------------------------------
// Gradient example.  (This example exists to test operations that use
// guard layers and have a left extent to make sure we compute on the
// correct region.)
//-----------------------------------------------------------------------------

struct EdgeFromCenterDerivativeInfo
{
  EdgeFromCenterDerivativeInfo(int dimensions)
    : dimensions_m(dimensions)
  {
  }

  void scalarCodeInfo(ScalarCodeInfo &info) const
  {
    info.arguments(2);

    info.write(0, true );
    info.write(1, false);
    info.useGuards(0, false);
    info.useGuards(1, true);

    info.dimensions(dimensions_m);

    int i;
    for (i = 0; i < dimensions_m; ++i)
    {
      info.lowerExtent(i) = 1;
      info.upperExtent(i) = 0;
    }
  }

  int dimensions_m;
};

template<int Dim>
struct EdgeFromCenterDerivative
  : public  EdgeFromCenterDerivativeInfo
{
  EdgeFromCenterDerivative()
    : EdgeFromCenterDerivativeInfo(Dim),
      off_m(0)
  {
    off_m[0] = 1;
  }

  template<class F1, class F2>
  inline
  void operator()(F1 &f1, const F2 &f2, const Loc<Dim> &loc) const
  {
    if (f2(loc) > 4.0)
    {
      f1(loc) = f2(loc) - 2.0 * f2(loc - off_m);
    }
    else
    {
      f1(loc) = f2(loc) - 1.1 * f2(loc - off_m);
    }
  }

  Loc<Dim> off_m;
};

//-----------------------------------------------------------------------------
// Main program
//-----------------------------------------------------------------------------

// set the problem dimension here:  must be >= 2

const int dim = 2;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  int NX = 5;
  int d;

  Interval<1> I(NX);
  Interval<dim> physicalVertexDomain;
  
  for (d = 0; d < dim; d++) 
  {
    physicalVertexDomain[d] = I;
  }

  // Create the mesh.
  
  Vector<dim, double> origin;
  Vector<dim, double> spacings;
  for (d = 0; d < dim; d++) 
  {
    origin(d) = d;
    spacings(d) = d + 1;
  }
      
  // Make a Brick-Engine-based field.

  DomainLayout<dim> layout1(physicalVertexDomain,
			    GuardLayers<dim>(0));
  Loc<dim> blocks(2);
  GridLayout<dim> layout2(physicalVertexDomain, blocks, 
			  GuardLayers<dim>(1), GuardLayers<dim>(0),
			  LayoutTag_t());

  typedef UniformRectilinearMesh<dim> Mesh_t;
  typedef Field<Mesh_t, double, Brick> FieldBrick_t;
  typedef MultiPatch<GridTag, BrickTag_t> MP2_t;
  typedef Field<Mesh_t, double, MP2_t> Field_t;
  typedef Field<Mesh_t, Vector<dim>, MP2_t> FieldV_t;

  Centering<dim> cell = canonicalCentering<dim>(CellType, Continuous);
  Centering<dim> allFace = canonicalCentering<dim>(FaceType, Continuous);

  FieldBrick_t f(allFace, layout1, origin, spacings);
  Pooma::addAllConstantFaceBC(f, 4.0, true);

  PositionsTraits<Mesh_t>::Type_t x = positions(f);

  for(d = 0; d < dim; ++d)
  {
    f[d] = x[d].comp(d);
  }

  tester.out() << "input field:" << std::endl
	       << f << std::endl;

  // Make a Nonuniform Multipatch-Engine-based field.

  Field_t fg(cell, layout2, origin, spacings);

  ScalarCode<AllFaceToCellAverage<dim> > faceToCell;

  faceToCell(fg, f);
    
  tester.out() << "result:" << std::endl
	       << fg << std::endl;  

  Field_t fgCheck(cell, layout2, origin, spacings);

  Interval<dim> cellDomain = f.physicalCellDomain();

  fgCheck = 0.0;
  for(d = 0; d < dim; ++d)
  {
    Loc<dim> off(0);
    off[d] = 1;
    fgCheck(cellDomain) += f[d](cellDomain) + f[d](cellDomain + off);
  }
  fgCheck /= (2.0 * dim);

  tester.out() << "input field:" << std::endl
	       << f << std::endl;

  tester.out() << "check:" << std::endl
	       << fgCheck << std::endl;  

  tester.check("scalar code differs from explicit computation",
	       sum(fgCheck - fg) < 0.001);

  // Now try a problem relevant to Blanca:

  GridLayout<dim> layout3(physicalVertexDomain, blocks, 
			  GuardLayers<dim>(1), GuardLayers<dim>(1),
			  LayoutTag_t());

  Centering<dim> edge = canonicalCentering<dim>(EdgeType, Continuous, YDim);
  Field_t EdgeQPressure(edge, layout3, origin, spacings);
  Field_t EdgeGammaConstant(edge, layout3, origin, spacings);
  Field_t EdgeSoundSpeed(edge, layout3, origin, spacings);
  FieldV_t EdgeVelocity(edge, layout3, origin, spacings);
  Field_t EdgePsiLimiter(edge, layout3, origin, spacings);

  DomainLayout<dim> layout4(physicalVertexDomain, GuardLayers<dim>(1));
  FieldBrick_t fEdge(edge, layout4, origin, spacings);

  EdgeGammaConstant = 1.4;
  EdgeSoundSpeed = 42.0 + positions(fEdge).comp(1);
  EdgeVelocity = 3.0 * positions(fEdge);
  EdgePsiLimiter = 1.0;

  typedef ScalarEdgeQPressure<dim> SEQP_t;
  SEQP_t sEQP(3.4);
  ScalarCode<SEQP_t> edgeQcompute(sEQP);

  edgeQcompute(EdgeQPressure,
	       EdgeGammaConstant, EdgeSoundSpeed,
	       EdgeVelocity, EdgePsiLimiter);

  tester.out() << "EdgeQPressure" << EdgeQPressure << std::endl;
  tester.out() << "EdgeGamma" << EdgeGammaConstant << std::endl;
  tester.out() << "EdgeSound" << EdgeSoundSpeed << std::endl;
  tester.out() << "EdgeV" << EdgeVelocity << std::endl;
  tester.out() << "EdgePsi" << EdgePsiLimiter << std::endl;

  // 2) scalar code with extents (lower extent in particular)

  Centering<dim> edgeY = canonicalCentering<dim>(EdgeType, Continuous, YDim);

  Field_t edgeValues(edgeY, layout3, origin, spacings);
  Field_t cellValues(cell, layout3, origin, spacings);
  FieldBrick_t fEdgeY(edgeY, layout4, origin, spacings);

  edgeValues.all() = 42.0;
  cellValues.all() = 5.0;
  cellValues = positions(fEdgeY).comp(0);

  tester.out() << "starting cell values: " << std::endl
	       << cellValues << std::endl
	       << "all" << std::endl
	       << cellValues.all() << std::endl;
  tester.out() << "starting edge values: " << std::endl
	       << edgeValues << std::endl
	       << "all" << std::endl
	       << edgeValues.all() << std::endl;

  ScalarCode<EdgeFromCenterDerivative<dim> > edgeFromCenter;

  edgeFromCenter(edgeValues, cellValues);

  tester.out() << "final edge values: " << std::endl
	       << edgeValues << std::endl
	       << "all" << std::endl
	       << edgeValues.all() << std::endl;

  double check2 = sum(edgeValues * edgeValues);
  tester.out() << "check value: " << check2 << std::endl;

  tester.check("value from derivative computation", std::abs(check2 - 134.8) < 0.2);

  // final cases to consider:
  // 1) replicated fields
  // 2) Lagrangian fields

  int ret = tester.results("ScalarCode");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ScalarCode.cpp,v $   $Author: richi $
// $Revision: 1.5 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo


