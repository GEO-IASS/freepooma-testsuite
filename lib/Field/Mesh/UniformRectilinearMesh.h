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
 * A rectilinear mesh with uniform spacing between vertices.
 */

#ifndef POOMA_FIELD_MESH_UNIFORMRECTILINEARMESH_H
#define POOMA_FIELD_MESH_UNIFORMRECTILINEARMESH_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/ConstantFunctionEngine.h"         // Used in functors
#include "Engine/IndexFunctionEngine.h"            // Used in functors
#include "Layout/INode.h"                          // Used in ctors
#include "Field/FieldEngine/FieldEnginePatch.h" // Used in ctors
#include "Field/Mesh/NoMesh.h"                  // Base class
#include "Field/FieldCentering.h"               // Centering<Dim> inline
#include "Tiny/Vector.h"                        // Class member

/**
 * Holds the data for a uniform rectilinear mesh. That class has a ref-counted
 * instance of this class.
 */

template <int Dim, class T>
class UniformRectilinearMeshData : public NoMeshData<Dim>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  /// We provide a default constructor that creates the object with empty
  /// domains. To be useful, this object must be replaced by another 
  /// version via assignment.
  
  UniformRectilinearMeshData()
    { 
      // This space intentionally left blank.
    }

  /// This constructor fully constructs the object. It uses the layout to
  /// compute domains and also initializes the origin and the spacings in each
  /// coordinate direction. The indices in the layout refer to VERTEX
  /// positions.

  template<class Layout>  
  UniformRectilinearMeshData(
    const Layout &layout,
    const Vector<Dim, T> &origin,
    const Vector<Dim, T> &spacings)
  : NoMeshData<Dim>(layout), 
    origin_m(origin),
    spacings_m(spacings)
    {
    }
    
  /// Copy constructor.

  UniformRectilinearMeshData(const UniformRectilinearMeshData<Dim, T> &model)
  : NoMeshData<Dim>(model), 
    origin_m(model.origin_m),
    spacings_m(model.spacings_m)
    {
      // This space intentionally left blank.
    } 
    
  /// @name View constructors.
  //@{
  
  /// Interval view. This means that we simply need to adjust the
  /// origin by the amount the view is offset from the model's physical
  /// cell domain. We rely on the base class to do the heavy lifting
  /// with respect to figuring out the domains correctly.
  ///
  /// The Interval supplied must refer to VERTEX positions.
  
  UniformRectilinearMeshData(const UniformRectilinearMeshData<Dim, T> &model, 
                             const Interval<Dim> &d)
  : NoMeshData<Dim>(d),
    origin_m(model.origin_m),
    spacings_m(model.spacings_m)
    {
      for (int i = 0; i < Dim; i++)
        origin_m(i) += 
          spacings_m(i) * 
            (d[i].first() - model.physicalCellDomain()[i].first());
    }
  
  /// FieldEnginePatch view. We don't fiddle with the origin because we are not
  /// making the domain zero-based.
  ///
  /// The domain supplied by the FieldEnginePatch must refer to VERTEX
  /// positions.
  
  UniformRectilinearMeshData(const UniformRectilinearMeshData<Dim, T> &model, 
                             const FieldEnginePatch<Dim> &p)
  : NoMeshData<Dim>(model, p),
    origin_m(model.origin_m),
    spacings_m(model.spacings_m)
    {
    }
  
  //@}

  //---------------------------------------------------------------------------
  /// Copy assignment operator.
  
  UniformRectilinearMeshData<Dim, T> &
  operator=(const UniformRectilinearMeshData<Dim, T> &rhs)
    {
      if (this != &rhs)
        {
          NoMeshData<Dim>::operator=(rhs);
          origin_m = rhs.origin_m;
          spacings_m = rhs.spacings_m;
        }
        
      return *this;
    }

  //---------------------------------------------------------------------------
  /// Empty destructor is fine. Note, however, that NoMeshData does not have
  /// a virtual destructor. We must be careful to delete these puppies as
  /// UniformRectilinearMeshData.

  ~UniformRectilinearMeshData() { }

  //---------------------------------------------------------------------------
  /// @name General accessors.
  //@{

  /// The mesh spacing.
  
  inline const Vector<Dim, T> &spacings() const 
    { 
      return spacings_m; 
    }

  /// The mesh origin.

  inline const Vector<Dim, T> &origin() const 
    { 
      return origin_m; 
    }

  //@}

private:

  /// Origin of mesh and spacings between vertices.

  Vector<Dim, T> origin_m, spacings_m;

};


/**
 * UniformRectilinearMesh is the simplest form of rectilinear mesh in that it
 * has uniform spacing between vertices. This spacing can be different in
 * each coordinate direction.
 */

template<int Dim, class T = POOMA_DEFAULT_POSITION_TYPE>
class UniformRectilinearMesh
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.

  /// The type of mesh points.
    
  typedef Vector<Dim, T> PointType_t;

  /// The type of vectors used to represent, for example, normals.
  
  typedef Vector<Dim, T> VectorType_t;
  
  /// The type T, used to represent, for example, volumes & areas, etc.
  
  typedef T T_t;
    
  /// The number of indices required to select a point in this mesh.
  
  enum { dimensions = Dim };

  //---------------------------------------------------------------------------
  // Constructors.
  
  /// We supply a default constructor, but it doesn't generate a useful mesh.
  /// This is accomplished through assignment.
  
  UniformRectilinearMesh() 
  : data_m(new UniformRectilinearMeshData<Dim, T>)
    { 
      // This space intentionally left blank.
    }

  /// This constructor fully constructs the object. It uses the layout to
  /// compute domains and also initializes the origin and the spacings in each
  /// coordinate direction.
  ///
  /// The Layout supplied must refer to VERTEX positions.
  
  template<class Layout>
  inline UniformRectilinearMesh(const Layout &layout, 
                                const PointType_t &origin,
                                const PointType_t &spacings)
  : data_m(new UniformRectilinearMeshData<Dim, T>(layout, origin, spacings))
    { 
    }
    
  template<class Layout>
  inline UniformRectilinearMesh(const Layout &layout)
  : data_m(new UniformRectilinearMeshData<Dim, T>(layout, 
                                                  PointType_t(0), 
                                                  PointType_t(1)))
    { 
    }
    
  /// Copy constructor. 
  
  inline UniformRectilinearMesh(const UniformRectilinearMesh<Dim, T> &model)
  : data_m(model.data_m)
    {
    }
    
  /// @name View constructors
  /// These are the only possible views of this
  /// mesh. Other views will make a NoMesh.
  //@{ 
  
  /// Interval view.
  ///
  /// The Interval supplied must refer to VERTEX positions.
  
  inline UniformRectilinearMesh(const UniformRectilinearMesh<Dim, T> &model, 
                                const Interval<Dim> &d)
  : data_m(new UniformRectilinearMeshData<Dim, T>(*model.data_m, d))
    {
    }
  
  /// INode view.
  ///
  /// The INode supplied must refer to VERTEX positions.
  
  inline UniformRectilinearMesh(const UniformRectilinearMesh<Dim, T> &model, 
                                const INode<Dim> &i)
  : data_m(new UniformRectilinearMeshData<Dim, T>(*model.data_m, i.domain()))
    {
    }
    
  /// FieldEnginePatch view.
  ///
  /// The FieldEnginePatch supplied must refer to VERTEX positions.
  
  inline UniformRectilinearMesh(const UniformRectilinearMesh<Dim, T> &model, 
                                const FieldEnginePatch<Dim> &p)
  : data_m(new UniformRectilinearMeshData<Dim, T>(*model.data_m, p))
    {
    }

  //@}

  //---------------------------------------------------------------------------
  /// Copy assignment operator.
  
  inline UniformRectilinearMesh<Dim, T> &
  operator=(const UniformRectilinearMesh<Dim, T> &rhs)
    {
      if (&rhs != this)
        {
          data_m = rhs.data_m;
        }
      
      return *this;
    }

  //---------------------------------------------------------------------------
  /// Empty destructor is fine. The pointer to the data is ref-counted so its
  /// lifetime is correctly managed.
  
  ~UniformRectilinearMesh()
    {
    }
  
  //---------------------------------------------------------------------------
  /// @name Domain functions.
  //@{
  
  /// The vertex domain, as the mesh was constructed with.

  inline const Interval<Dim> &physicalVertexDomain() const
    {
      return data_m->physicalVertexDomain(); 
    }

  /// Function that returns a domain adjusted to give the indices of the cells.

  inline const Interval<Dim> &physicalCellDomain() const
    {
      return data_m->physicalCellDomain(); 
    }

  /// The total vertex domain, including mesh guard vertices.

  inline const Interval<Dim> &totalVertexDomain() const
    {
      return data_m->totalVertexDomain(); 
    }

  /// The total cell domain, including mesh guard cells.

  inline const Interval<Dim> &totalCellDomain() const
    {
      return data_m->totalCellDomain(); 
    }

  //@}

  //---------------------------------------------------------------------------
  /// @name General accessors.
  //@{
  
  /// The mesh spacing.
  
  inline const Vector<Dim, T> &spacings() const 
    { 
      return data_m->spacings(); 
    }

  /// The mesh origin.

  inline const Vector<Dim, T> &origin() const 
    { 
      return data_m->origin(); 
    }

  /// The cell containing a particular point.
  
  inline Loc<Dim> cellContaining(const Vector<Dim, T> &point) const
    {
      Loc<Dim> loc((0, Pooma::NoInit()));	// Avoid a g++ parse error.
      
      for (int i = 0; i < Dim; i++)
        loc[i] = 
          Loc<1>(static_cast<int>((point(i) - origin()(i)) / spacings()(i)));
      
      return loc;
    }

  /// The lower-left vertex associated with a given cell location.
    
  inline Vector<Dim, T> vertexPosition(const Loc<Dim> &loc) const
    {
      Vector<Dim, T> point;
      
      for (int i = 0; i < Dim; i++)
        point(i) = origin()(i) + spacings()(i) * 
          (loc[i].first() - physicalCellDomain()[i].first());

      return point;
    }

  //@}

  //---------------------------------------------------------------------------
  /// Support for the positions() function. We need to provide a functor for
  /// use with IndexFunction-engine. We also need to export the
  /// PositionsEngineTag_t typedef and the positionsFunctor() member function,
  /// which computes the positions using the centering point positions.
  /// The indices passed in refer to cells.
  
  class PositionsFunctor {
  public:
  
    /// Need to be able to default construct since we fill in the details
    /// after the fact.
    
    PositionsFunctor() { }
    
    PositionsFunctor(const UniformRectilinearMesh<Dim, T> &m, 
                     const Centering<Dim> &c)
      : origin_m(m.origin()), spacings_m(m.spacings())
      {
        for (int i = 0; i < Dim; i++)
          origin_m(i) += spacings_m(i) * 
            (c.position(0)(i) - m.physicalCellDomain()[i].first());
      }
      
    inline PointType_t operator()(int i0) const
      {
        return origin_m + PointType_t(i0) * spacings_m;
      }
      
    inline PointType_t operator()(int i0, int i1) const
      {
        return origin_m + PointType_t(i0, i1) * spacings_m;
      }

    inline PointType_t operator()(int i0, int i1, int i2) const
      {
        return origin_m + PointType_t(i0, i1, i2) * spacings_m;
      }

  private:
  
    PointType_t origin_m, spacings_m;

  };
  
  typedef IndexFunction<PositionsFunctor> PositionsEngineTag_t;
  
  void initializePositions(
    Engine<Dim, PointType_t, PositionsEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      e.setFunctor(PositionsFunctor(*this, c));
    }
  
  //---------------------------------------------------------------------------
  /// Support for the outwardNormals() and coordinateNormals() functions. 
  /// We also need to export the NormalsEngineTag_t typedef and the 
  /// initializeNormals() member function, which sets the appropriate constant 
  /// value (since the normals exactly align with the coordinate axes).
  /// The boolean value passed is true if we are asking for outward normals,
  /// as opposed to coordinate normals. The indices passed in refer to cells.

  typedef ConstantFunction NormalsEngineTag_t;
  
  void initializeNormals(
    Engine<Dim, VectorType_t, NormalsEngineTag_t> &e, 
    const Centering<Dim> &c,
    bool outward = true) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be face-centered.
      
      PAssert(c.size() == 1);
      PAssert(c.centeringType() == FaceType);
      
      // Generate the normals. The coordinate normals are computed from
      // 1 - orientation. Then, if we are on the near face, indicated by
      // position == 0.0, we need to multiply by -1.0 if we are doing
      // outward normals.
      
      VectorType_t normal;
      for (int i = 0; i < Dim; i++)
        {
          normal(i) = static_cast<T_t>(1 - c.orientation(0)[i].first());
          if (outward && c.position(0)(i) == 0.0)
            normal(i) *= static_cast<T_t>(-1);
        }
        
      e.setConstant(normal);
    }
  
  //---------------------------------------------------------------------------
  /// Support for the cellVolumes() function. We also need to export the 
  /// CellVolumesEngineTag_t typedef and the initializeCellVolumes() member 
  /// function, which sets the appropriate constant value for the volume. 
  /// The indices passed in refer to cells.

  typedef ConstantFunction CellVolumesEngineTag_t;
  
  void initializeCellVolumes(
    Engine<Dim, T_t, CellVolumesEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be cell-centered.
      
      PAssert(c.size() == 1);
      PAssert(c.centeringType() == CellType);
      
      // Use the general function to do the job.
      
      initializeGeneralVolume(e, c);
    }
  
  //---------------------------------------------------------------------------
  /// Support for the faceAreas() function. We also need to export the 
  /// FaceAreasEngineTag_t typedef and the initializeFaceAreas() member 
  /// function, which sets the appropriate constant face area value.
  /// The indices passed in refer to cells.

  typedef ConstantFunction FaceAreasEngineTag_t;
  
  void initializeFaceAreas(
    Engine<Dim, T_t, FaceAreasEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be face-centered.
      
      PAssert(c.size() == 1);
      PAssert(c.centeringType() == FaceType);
      
      // Use the general function to do the job.
      
      initializeGeneralVolume(e, c);
    }
  
  //---------------------------------------------------------------------------
  /// Support for the edgeLengths() function. We also need to export the 
  /// EdgeLengthsEngineTag_t typedef and the initializeEdgeLengths() member 
  /// function, which sets the appropriate constant edge length value.
  /// The indices passed in refer to cells.

  typedef ConstantFunction EdgeLengthsEngineTag_t;
  
  void initializeEdgeLengths(
    Engine<Dim, T_t, EdgeLengthsEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be edge-centered.
      
      PAssert(c.size() == 1);
      PAssert(c.centeringType() == EdgeType);
      
      // Use the general function to do the job.
      
      initializeGeneralVolume(e, c);
    }
  
private:

  /// General "volume" computation: works for edges, faces, and cells.
  
  void initializeGeneralVolume(
    Engine<Dim, T_t, ConstantFunction> &e, 
    const Centering<Dim> &c) const
    {
      // Generate the volume by multiplying the spacings in the
      // directions where the orientation is non-zero.
      
      T_t volume = static_cast<T_t>(1);
      for (int i = 1; i < Dim; i++)
        {
          if (c.orientation(0)[i].first() != 0)
            volume *= spacings()(i);
        }
        
      e.setConstant(volume);
    }

  /// Our data, stored as a ref-counted pointer to simplify memory management.
  
  RefCountedPtr<UniformRectilinearMeshData<Dim, T> > data_m;
};

#endif // POOMA_FIELD_MESH_UNIFORMRECTILINEARMESH_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformRectilinearMesh.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:46 $
// ----------------------------------------------------------------------
// ACL:rcsinfo

