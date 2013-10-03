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

// include files
#include "Field/FieldCentering.h"


//----------------------------------------------------------------------
// Construct canonical centerings so the programmer need not do this.
//----------------------------------------------------------------------

template <int Dim>
CanonicalCentering<Dim>::CanonicalCentering()
{
  Centering<Dim> centering;
  typename Centering<Dim>::Orientation orientation;
  typename Centering<Dim>::Position position;
  typename Centering<Dim>::Orientations orientations[Dim][2];
  typename Centering<Dim>::Positions positions[Dim][2];
  enum { x = 0, y, z };

  // Create the tables if necessary.
  if (class_count_m == 0) {
    centering_table_m = new Centering<Dim>**[CellType+1];
    for (int i = 0; i <= CellType; ++i) {
      centering_table_m[i] = new Centering<Dim>*[2];
      for (int j = 0; j < 2; ++j)
	centering_table_m[i][j] = new Centering<Dim>[1<<Dim];
    }
  }
  ++class_count_m;

  // Add the cell centerings.
  centering = Centering<Dim>(CellType, Continuous);
  orientation = 1;
  position = 0.5;
  centering.addValue(orientation, position);
  centering_table_m[CellType][Continuous][AllDim%(1<<Dim)] = centering;

  // Add the edge centerings.
  // X Edge
  orientation = 0; orientation[0] = 1;
  position = 0.0; position(0) = 0.5;
  addValue(orientations[x][Continuous],
	   positions[x][Continuous],
	   orientation, position);
  orientations[x][Discontinuous] =
    orientations[x][Continuous];
  positions[x][Discontinuous] =
    positions[x][Continuous];
  if (Dim > 1) {
    position(1) = 1.0;
    addValue(orientations[x][Discontinuous],
	     positions[x][Discontinuous],
	     orientation, position);
    if (Dim > 2) {
      position(2) = 1.0;
      addValue(orientations[x][Discontinuous],
	       positions[x][Discontinuous],
	       orientation, position);
      position(1) = 0.0;
      addValue(orientations[x][Discontinuous],
	       positions[x][Discontinuous],
	       orientation, position);
    }
  }

  // Y Edge
  if (Dim > 1) {
    orientation = 0; orientation[1] = 1;
    position = 0.0; position(1) = 0.5;
    addValue(orientations[y][Continuous],
	     positions[y][Continuous],
	     orientation, position);
    orientations[y][Discontinuous] =
      orientations[y][Continuous];
    positions[y][Discontinuous] =
      positions[y][Continuous];
    position(0) = 1.0;
    addValue(orientations[y][Discontinuous],
	     positions[y][Discontinuous],
	     orientation, position);
    if (Dim > 2) {
      position(2) = 1.0;
      addValue(orientations[y][Discontinuous],
	       positions[y][Discontinuous],
	       orientation, position);
      position(0) = 0.0;
      addValue(orientations[y][Discontinuous],
	       positions[y][Discontinuous],
	       orientation, position);
    }
  }

  // Z Edge
  if (Dim > 2) {
    orientation = 0; orientation[2] = 1;
    position = 0.0; position(2) = 0.5;
    addValue(orientations[z][Continuous],
	     positions[z][Continuous],
	     orientation, position);
    orientations[z][Discontinuous] =
      orientations[z][Continuous];
    positions[z][Discontinuous] =
      positions[z][Continuous];
    position(0) = 1.0;
    addValue(orientations[z][Discontinuous],
	     positions[z][Discontinuous],
	     orientation, position);
    position(1) = 1.0;
    addValue(orientations[z][Discontinuous],
	     positions[z][Discontinuous],
	     orientation, position);
    position(0) = 0.0;
    addValue(orientations[z][Discontinuous],
	     positions[z][Discontinuous],
	     orientation, position);
  }

  for (int cont = 0; cont < 2; ++cont) {
    centering_table_m[EdgeType][cont][XDim] =
      Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		     orientations[x][cont], positions[x][cont]);
    if (Dim > 1) {
      centering_table_m[EdgeType][cont][YDim] =
	Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		       orientations[y][cont], positions[y][cont]);
      centering_table_m[EdgeType][cont][XDim|YDim] =
	Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[x][cont],orientations[y][cont]),
		       combine(positions[x][cont], positions[y][cont]));
    }
    if (Dim > 2) {
      centering_table_m[EdgeType][cont][ZDim] =
	Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		       orientations[z][cont], positions[z][cont]);
      centering_table_m[EdgeType][cont][XDim|ZDim] =
	Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[x][cont],orientations[z][cont]),
		       combine(positions[x][cont], positions[z][cont]));
      centering_table_m[EdgeType][cont][YDim|ZDim] =
	Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[y][cont],orientations[z][cont]),
		       combine(positions[y][cont], positions[z][cont]));
      centering_table_m[EdgeType][cont][XDim|YDim|ZDim] =
	Centering<Dim>(EdgeType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[x][cont],
			       combine(orientations[y][cont],orientations[z][cont])),
		       combine(positions[x][cont],
			       combine(positions[y][cont], positions[z][cont])));
    }
  }  

  for (int dim = 0; dim < Dim; ++dim)
    for (int cont = 0; cont < 2; ++cont)
      {
	orientations[dim][cont].clear();
	positions[dim][cont].clear();
      }

  // Add the face centerings.
  // X Face
  orientation = 1; orientation[0] = 0;
  position = 0.5; position(0) = 0.0;
  addValue(orientations[x][Continuous],
	   positions[x][Continuous],
	   orientation, position);
  orientations[x][Discontinuous] =
    orientations[x][Continuous];
  positions[x][Discontinuous] =
    positions[x][Continuous];
  position(0) = 1.0;
  addValue(orientations[x][Discontinuous],
	   positions[x][Discontinuous],
	   orientation, position);

  // Y Face
  if (Dim > 1) {
    orientation = 1; orientation[1] = 0;
    position = 0.5; position(1) = 0.0;
    addValue(orientations[y][Continuous],
	     positions[y][Continuous],
	     orientation, position);
    orientations[y][Discontinuous] =
      orientations[y][Continuous];
    positions[y][Discontinuous] =
      positions[y][Continuous];
    position(1) = 1.0;
    addValue(orientations[y][Discontinuous],
	     positions[y][Discontinuous],
	     orientation, position);
  }

  // Z Face
  if (Dim > 2) {
    orientation = 1; orientation[2] = 0;
    position = 0.5; position(2) = 0.0;
    addValue(orientations[z][Continuous],
	     positions[z][Continuous],
	     orientation, position);
    orientations[z][Discontinuous] =
      orientations[z][Continuous];
    positions[z][Discontinuous] =
      positions[z][Continuous];
    position(2) = 1.0;
    addValue(orientations[z][Discontinuous],
	     positions[z][Discontinuous],
	     orientation, position);
  }

  for (int cont = 0; cont < 2; ++cont) {
    centering_table_m[FaceType][cont][XDim] =
      Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		     orientations[x][cont], positions[x][cont]);
    if (Dim > 1) {
      centering_table_m[FaceType][cont][XDim|YDim] =
	Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[x][cont],orientations[y][cont]),
		       combine(positions[x][cont], positions[y][cont]));
      centering_table_m[FaceType][cont][YDim] =
	Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		       orientations[y][cont], positions[y][cont]);
    }
    if (Dim > 2) {
      centering_table_m[FaceType][cont][ZDim] =
	Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		       orientations[z][cont], positions[z][cont]);
      centering_table_m[FaceType][cont][XDim|ZDim] =
	Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[x][cont],orientations[z][cont]),
		       combine(positions[x][cont], positions[z][cont]));
      centering_table_m[FaceType][cont][YDim|ZDim] =
	Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[y][cont],orientations[z][cont]),
		       combine(positions[y][cont], positions[z][cont]));
      centering_table_m[FaceType][cont][XDim|YDim|ZDim] =
	Centering<Dim>(FaceType, static_cast<enum ContinuityType>(cont),
		       combine(orientations[x][cont],
			       combine(orientations[y][cont],orientations[z][cont])),
		       combine(positions[x][cont],
			       combine(positions[y][cont], positions[z][cont])));
    }
  }  
  for (int dim = 0; dim < Dim; ++dim)
    for (int cont = 0; cont < 2; ++cont)
      {
	orientations[dim][cont].clear();
	positions[dim][cont].clear();
      }

  // Add the vertex centerings.
  centering = Centering<Dim>(VertexType, Continuous);
  orientation = 0;
  position = 0.0;
  centering.addValue(orientation, position);
  centering_table_m[VertexType][Continuous][AllDim%(1<<Dim)] =
    centering;

  centering = Centering<Dim>(VertexType, Discontinuous);
  orientation = 0;
  position = 0.0;
  centering.addValue(orientation, position);
  position(0) = 1.0; centering.addValue(orientation, position);
  if (Dim > 1) {
    position(1) = 1.0; centering.addValue(orientation, position);
    position(0) = 0.0; centering.addValue(orientation, position);
    if (Dim > 2) {
      position(2) = 1.0; centering.addValue(orientation, position);
      position(0) = 1.0; centering.addValue(orientation, position);
      position(1) = 0.0; centering.addValue(orientation, position);
      position(0) = 0.0; centering.addValue(orientation, position);
    }
  }
  centering_table_m[VertexType][Discontinuous][AllDim%(1<<Dim)] =
    centering;

  return;
}


// The canonical centering objects
const CanonicalCentering<1> canonicalCenteringOne_g;
const CanonicalCentering<2> canonicalCenteringTwo_g;
const CanonicalCentering<3> canonicalCenteringThree_g;

template <int Dim>
const Centering<Dim> canonicalCentering
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension); // Unsupported dimension

template <>
const Centering<1> canonicalCentering<1>
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension)
{
  return canonicalCenteringOne_g(type, discontinuous, dimension);
}

// Instantiate the constructors we care about:

template <>
const Centering<2> canonicalCentering<2>
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension)
{
  return canonicalCenteringTwo_g(type, discontinuous, dimension);
}


template <>
const Centering<3> canonicalCentering<3>
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension)
{
  return canonicalCenteringThree_g(type, discontinuous, dimension);
}


// Pre-instantiate canonical centerings

template class CanonicalCentering<1>;
template class CanonicalCentering<2>;
template class CanonicalCentering<3>;


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldCentering.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:42 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
