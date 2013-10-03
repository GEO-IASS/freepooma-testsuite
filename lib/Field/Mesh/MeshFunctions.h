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

/** @file
 * @ingroup Mesh
 * @brief
 * Mesh functions for querying mesh properties.
 *
 * Functions: 
 *   - positions() returns the centering point locations for a Field. 
 *   - outwardNormals() returns outward-facing normals for a Field.
 *   - coordinateNormals() returns coordinate normals for a Field.
 *   - cellVolumes() returns coordinate normals for a Field.
 *   - faceAreas() returns coordinate normals for a Field.
 *   - edgeLengths() returns coordinate normals for a Field.
 */

#ifndef POOMA_FIELD_MESH_MESHFUNCTIONS_H
#define POOMA_FIELD_MESH_MESHFUNCTIONS_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Field/Mesh/NoMesh.h"                  // Returned class

//-----------------------------------------------------------------------------
// Forward declarartions:
//-----------------------------------------------------------------------------

template<class Mesh, class T, class EngineTag> class Field;


template<class Mesh>
struct PositionsTraits {
  typedef Field<NoMesh<Mesh::dimensions>, 
    typename Mesh::PointType_t, 
    typename Mesh::PositionsEngineTag_t> Type_t;
};

/// This function returns the centering point locations for a Field f.
/// These are returned in a field with no mesh, but the same centering
/// as the original Field.
    
template<class Mesh, class T, class EngineTag>
typename PositionsTraits<Mesh>::Type_t
positions(const Field<Mesh, T, EngineTag> &f)
{
  NoMesh<Mesh::dimensions> m(f.layout());
  
  typename PositionsTraits<Mesh>::Type_t 
    of(f.numMaterials(), f.centering(), f.layout(), m);
  
  for (int i = 0; i < of.numMaterials(); i++)
    for (int j = 0; j < of.centeringSize(); j++)
      f.mesh().initializePositions(of.subField(i, j).engine(), of.centering(j));
         
  return of;
}

template<class Mesh>
struct NormalsTraits {
  typedef Field<NoMesh<Mesh::dimensions>, 
    typename Mesh::VectorType_t, 
    typename Mesh::NormalsEngineTag_t> Type_t;
};

/// This function returns outward-facing normals for a Field f.
/// These are returned in a discontinuous face-centered field with 
/// no mesh.

template<class Mesh, class T, class EngineTag>
typename NormalsTraits<Mesh>::Type_t
outwardNormals(const Field<Mesh, T, EngineTag> &f)
{
  NoMesh<Mesh::dimensions> m(f.layout());
  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(FaceType, Discontinuous);
    
  typename NormalsTraits<Mesh>::Type_t 
    of(f.numMaterials(), c, f.layout(), m);
  
  for (int i = 0; i < of.numMaterials(); i++)
    for (int j = 0; j < of.centeringSize(); j++)
      f.mesh().initializeNormals(of.subField(i, j).engine(), of.centering(j), true);
         
  return of;
}

/// This function returns coordinate normals for a Field f.
/// These are returned in a continuous face-centered field with 
/// no mesh.

template<class Mesh, class T, class EngineTag>
typename NormalsTraits<Mesh>::Type_t
coordinateNormals(const Field<Mesh, T, EngineTag> &f)
{
  NoMesh<Mesh::dimensions> m(f.layout());
  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(FaceType, Continuous);
    
  typename NormalsTraits<Mesh>::Type_t 
    of(f.numMaterials(), c, f.layout(), m);
  
  for (int i = 0; i < of.numMaterials(); i++)
    for (int j = 0; j < of.centeringSize(); j++)
      f.mesh().initializeNormals(of.subField(i, j).engine(), of.centering(j), false);
         
  return of;
}


template<class Mesh>
struct CellVolumesTraits {
  typedef Field<NoMesh<Mesh::dimensions>, 
    typename Mesh::T_t, 
    typename Mesh::CellVolumesEngineTag_t> Type_t;
};

/// This function returns the cell volumes for a Field f.
/// These are returned in a cell-centered field with no mesh.

template<class Mesh, class T, class EngineTag>
typename CellVolumesTraits<Mesh>::Type_t
cellVolumes(const Field<Mesh, T, EngineTag> &f)
{
  NoMesh<Mesh::dimensions> m(f.layout());
  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(CellType, Continuous);
    
  typename CellVolumesTraits<Mesh>::Type_t
    of(f.numMaterials(), c, f.layout(), m);
  
  for (int i = 0; i < of.numMaterials(); i++)
    f.mesh().initializeCellVolumes(of.subField(i, 0).engine(), of.centering(0));
         
  return of;
}


template<class Mesh>
struct FaceAreasTraits {
  typedef Field<NoMesh<Mesh::dimensions>, 
    typename Mesh::T_t, 
    typename Mesh::FaceAreasEngineTag_t> Type_t;
};

/// This function returns the face areas for a Field f.
/// These are returned in a continuous face-centered field 
/// with no mesh.

template<class Mesh, class T, class EngineTag>
typename FaceAreasTraits<Mesh>::Type_t 
faceAreas(const Field<Mesh, T, EngineTag> &f)
{
  NoMesh<Mesh::dimensions> m(f.layout());
  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(FaceType, Continuous);
    
  typename FaceAreasTraits<Mesh>::Type_t 
    of(f.numMaterials(), c, f.layout(), m);
  
  for (int i = 0; i < of.numMaterials(); i++)
    for (int j = 0; j < of.centeringSize(); j++)
      f.mesh().initializeFaceAreas(of.subField(i, j).engine(), of.centering(j));
         
  return of;
}


template<class Mesh>
struct EdgeLengthsTraits {
  typedef Field<NoMesh<Mesh::dimensions>, 
    typename Mesh::T_t, 
    typename Mesh::EdgeLengthsEngineTag_t> Type_t;
};

/// This function returns the edge lengths for a Field f.
/// These are returned in a continuous edge-centered field 
/// with no mesh.

template<class Mesh, class T, class EngineTag>
typename EdgeLengthsTraits<Mesh>::Type_t 
edgeLengths(const Field<Mesh, T, EngineTag> &f)
{
  NoMesh<Mesh::dimensions> m(f.layout());
  Centering<Mesh::dimensions> c = 
    canonicalCentering<Mesh::dimensions>(EdgeType, Continuous);
    
  typename EdgeLengthsTraits<Mesh>::Type_t
    of(f.numMaterials(), c, f.layout(), m);
  
  for (int i = 0; i < of.numMaterials(); i++)
    for (int j = 0; j < of.centeringSize(); j++)
      f.mesh().initializeEdgeLengths(of.subField(i, j).engine(), of.centering(j));
         
  return of;
}

#endif // POOMA_FIELD_MESH_MESHFUNCTIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MeshFunctions.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:46 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
