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
 * A rectilinear mesh without uniform spacing between vertices.
 */

#ifndef POOMA_FIELD_MESH_RECTILINEARMESH_H
#define POOMA_FIELD_MESH_RECTILINEARMESH_H

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


//-----------------------------------------------------------------------------
/// Holds the data for a rectilinear mesh. That class has a ref-counted
/// instance of this class.
//-----------------------------------------------------------------------------

template <int Dim, class T>
class RectilinearMeshData : public NoMeshData<Dim>
{
public:

  typedef Array<1, T> SpacingsType_t[Dim];

  //---------------------------------------------------------------------------
  // Constructors.

  /// We provide a default constructor that creates the object with empty
  /// domains. To be useful, this object must be replaced by another 
  /// version via assignment.
  
  RectilinearMeshData()
    { 
      // This space intentionally left blank.
    }

  /// This constructor fully constructs the object. It uses the layout to
  /// compute domains and also initializes the origin and the spacings in each
  /// coordinate direction. The indices in the layout refer to VERTEX
  /// positions.

  template<class Layout, class EngineTag>
  RectilinearMeshData(
    const Layout &layout,
    const Vector<Dim, T> &origin,
    const SpacingsType_t &spacings)
  : NoMeshData<Dim>(layout), 
    origin_m(origin)
    //spacings_m(spacings)
    {
      for (int i=0; i<Dim; i++) {
	spacings_m[i].engine() = spacings[i].engine(); // init
	spacings_m[i].engine().makeOwnCopy(); // FIXME? Do we want this?
	Interval<1> I(layout.domain()[i]);
	positions_m[i].engine() = Engine<1, T, Brick>(I);
	positions_m[i](0) = origin_m(i);
	// initialize from origin downward the ghost cells
	for (int j=-1; j>=I.min(); j--)
	  positions_m[i](j) = positions_m[i].read(j+1) - spacings_m[i].read(j);
	// initialize from origin upward
	for (int j=1; j<=I.max(); j++)
	  positions_m[i](j) = positions_m[i].read(j-1) + spacings_m[i].read(j-1);
      }
    }

  /// Constructor for constructing evenly spaced rectilinear meshes just
  /// like UniformRectilinearMesh does.

  template<class Layout>
  RectilinearMeshData(
    const Layout &layout,
    const Vector<Dim, T> &origin,
    const Vector<Dim, T> &spacings)
  : NoMeshData<Dim>(layout), 
    origin_m(origin)
    {
      // for each dimension we allocate engines for spacings & positions
      // and initialize them according to origin/spacings
      for (int i=0; i<Dim; i++) {
	Interval<1> I(layout.domain()[i]);
	// allocate and assign spacings
	spacings_m[i].engine() = Engine<1, T, Brick>(I);
	spacings_m[i](I) = spacings(i); // no Array.all()
	Pooma::blockAndEvaluate();
	// allocate positions, assign origin
	positions_m[i].engine() = Engine<1, T, Brick>(I);
	positions_m[i](0) = origin_m(i);
	// initialize from origin downward the ghost cells
	for (int j=-1; j>=I.min(); j--)
	  positions_m[i](j) = positions_m[i].read(j+1) - spacings_m[i].read(j);
	// initialize from origin upward
	for (int j=1; j<=I.max(); j++)
	  positions_m[i](j) = positions_m[i].read(j-1) + spacings_m[i].read(j-1);
      }
    }
    
  /// Copy constructor.

  RectilinearMeshData(const RectilinearMeshData<Dim, T> &model)
  : NoMeshData<Dim>(model), 
    origin_m(model.origin_m)
    //spacings_m(model.spacings_m),
    //positions_m(model.positions_m)
    {
      for (int i=0; i<Dim; i++) {
	spacings_m[i].engine() = model.spacings_m[i].engine();
	positions_m[i].engine() = model.positions_m[i].engine();
      }
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
  
  RectilinearMeshData(const RectilinearMeshData<Dim, T> &model, 
		      const Interval<Dim> &d)
  : NoMeshData<Dim>(d)
    {
      for (int i = 0; i < Dim; i++) {
	// FIXME: Wheeee ;) (we cant store a BrickView...
	// and still dont want to copy)
	spacings_m[i].engine() = Engine<1, T, Brick>(&model.spacings_m[i](d[i])(0), d[i]);
	positions_m[i].engine() = Engine<1, T, Brick>(&model.positions_m[i](d[i])(0), d[i]);
	origin_m(i) = positions_m[i](d[i].min());
      }
    }
  
  /// FieldEnginePatch view. We don't fiddle with the origin because we are not
  /// making the domain zero-based.
  ///
  /// The domain supplied by the FieldEnginePatch must refer to VERTEX
  /// positions.
  
  RectilinearMeshData(const RectilinearMeshData<Dim, T> &model, 
		      const FieldEnginePatch<Dim> &p)
  : NoMeshData<Dim>(model, p),
    origin_m(model.origin_m),
    spacings_m(model.spacings_m),
    positions_m(model.spacings_m)
    {
      // FIXME: what does FieldEnginePatch do???
      for (int i=0; i<Dim; i++) {
	spacings_m[i].engine() = model.spacings_m[i].engine();
	positions_m[i].engine() = model.positions_m[i].engine();
      }
    }
  
  //@}

  //---------------------------------------------------------------------------
  /// Copy assignment operator.
  
  RectilinearMeshData<Dim, T> &
  operator=(const RectilinearMeshData<Dim, T> &rhs)
    {
      if (this != &rhs)
        {
          NoMeshData<Dim>::operator=(rhs);
          origin_m = rhs.origin_m;
	  for (int i=0; i<Dim; i++) {
	    spacings_m[i].engine() = rhs.spacings_m[i].engine();
	    positions_m[i].engine() = rhs.positions_m[i].engine();
	  }
        }
        
      return *this;
    }

  //---------------------------------------------------------------------------
  /// Empty destructor is fine. Note, however, that NoMeshData does not have
  /// a virtual destructor. We must be careful to delete these puppies as
  /// RectilinearMeshData.

  ~RectilinearMeshData() { }

  //---------------------------------------------------------------------------
  /// @name General accessors.
  //@{

  /// The mesh spacing.
  
  inline const SpacingsType_t& spacings() const 
    { 
      return spacings_m; 
    }

  /// The mesh vertex positions.
  
  inline const SpacingsType_t& positions() const 
    { 
      return positions_m; 
    }

  /// The mesh origin.

  inline const Vector<Dim, T> &origin() const 
    { 
      return origin_m; 
    }

  //@}

private:

  /// Origin of mesh (coordinate vector of first vertex).

  Vector<Dim, T> origin_m;

  /// Spacings between vertices.

  SpacingsType_t spacings_m;

  /// Vertex positions.

  SpacingsType_t positions_m;

};


///
/// RectilinearMesh is a rectilinear mesh sometimes called a 
/// "cartesian product" or "tensor product" mesh. Each dimension has a
/// spacing value between every pair of vertices along that dimension;
/// these spacings can all be different.
///
template<int Dim, class T = POOMA_DEFAULT_POSITION_TYPE>
class RectilinearMesh
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.

  /// The type of mesh points.
    
  typedef Vector<Dim, T> PointType_t;

  /// The type of vectors used to represent, for example, normals.
  
  typedef Vector<Dim, T> VectorType_t;

  /// The type used to store spacings.

  typedef Array<1, T> SpacingsType_t[Dim];
  
  /// The type T, used to represent, for example, volumes & areas, etc.
  
  typedef T T_t;

  /// The number of indices required to select a point in this mesh.
  
  enum { dimensions = Dim };

  //---------------------------------------------------------------------------
  // Constructors.
  
  /// We supply a default constructor, but it doesn't generate a useful mesh.
  /// This is accomplished through assignment.
  
  RectilinearMesh() 
  : data_m(new RectilinearMeshData<Dim, T>)
    { 
      // This space intentionally left blank.
    }

  /// This constructor fully constructs the object. It uses the layout to
  /// compute domains and also initializes the origin and the spacings in each
  /// coordinate direction.
  ///
  /// The Layout supplied must refer to VERTEX positions.
  
  template<class Layout, class EngineTag>
  inline RectilinearMesh(const Layout &layout, 
			 const PointType_t &origin,
			 const SpacingsType_t &spacings)
  : data_m(new RectilinearMeshData<Dim, T>(layout, origin, spacings))
    { 
    }

  /// Constructor compatible to UniformRectilinearMesh.

  template<class Layout>
  inline RectilinearMesh(const Layout &layout,
			 const PointType_t &origin,
			 const PointType_t &spacings)
  : data_m(new RectilinearMeshData<Dim, T>(layout, origin, spacings))
    { 
    }

  template<class Layout>
  inline RectilinearMesh(const Layout &layout)
  : data_m(new RectilinearMeshData<Dim, T>(layout, 
					   PointType_t(0), 
					   PointType_t(1)))
    { 
    }

  /// Copy constructor. 
  
  inline RectilinearMesh(const RectilinearMesh<Dim, T> &model)
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
  
  inline RectilinearMesh(const RectilinearMesh<Dim, T> &model, 
			 const Interval<Dim> &d)
  : data_m(new RectilinearMeshData<Dim, T>(*model.data_m, d))
    {
    }
  
  /// INode view.
  ///
  /// The INode supplied must refer to VERTEX positions.
  
  inline RectilinearMesh(const RectilinearMesh<Dim, T> &model, 
			 const INode<Dim> &i)
  : data_m(new RectilinearMeshData<Dim, T>(*model.data_m, i.domain()))
    {
    }
    
  /// FieldEnginePatch view.
  ///
  /// The FieldEnginePatch supplied must refer to VERTEX positions.
  
  inline RectilinearMesh(const RectilinearMesh<Dim, T> &model, 
			 const FieldEnginePatch<Dim> &p)
  : data_m(new RectilinearMeshData<Dim, T>(*model.data_m, p))
    {
    }

  //@}

  //---------------------------------------------------------------------------
  /// Copy assignment operator.
  
  inline RectilinearMesh<Dim, T> &
  operator=(const RectilinearMesh<Dim, T> &rhs)
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
  
  ~RectilinearMesh()
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

  inline const SpacingsType_t &spacings() const 
    { 
      return data_m->spacings();
    }

  /// The mesh positions.

  inline const SpacingsType_t &positions() const 
    { 
      return data_m->positions();
    }

  /// The mesh origin.

  inline const Vector<Dim, T> &origin() const 
    { 
      return data_m->origin();
    }

  /// The cell containing a particular point.

  inline Loc<Dim> cellContaining(const Vector<Dim, T> &point) const
    {
      /// FIXME
      Loc<Dim> loc((0, Pooma::NoInit()));	// Avoid a g++ parse error.
      for (int i = 0; i < Dim; i++)
	{
	  const T *start = &positions()[i](0);
	  const T *finish = start + positions()[i].physicalDomain()[i].length();
	  const T *p = std::lower_bound(start, finish, point(i));
#if POOMA_BOUNDS_CHECK
	  PInsist(p != finish,
		  "Rectilinear::cellContaining(): point is outside mesh.");
#endif
	  // The lower_bound function returns the first element that is not
	  // less than the point we're searching for.
	  int j = static_cast<int>(std::distance(start, p));
	  if (*p == point(i))
	    loc[i] = j;
	  else
	    loc[i] = j-1;
	}

      return loc;
    }

  /// The lower-left vertex associated with a given cell location.
    
  inline Vector<Dim, T> vertexPosition(const Loc<Dim> &loc) const
    {
      Vector<Dim, T> point;
      for (int i = 0; i < Dim; i++)
        point(i) = positions()[i](loc[i]); 
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

    // WARNING! For Arrays to be initialized (copy constructed, assigned,
    //          etc.) correctly, even in the case of uninitialized targets
    //          we need to copy the engines explicitly rather than rely
    //          on the compiler generating correct copy constructors and
    //          assignment operators.
    // FIXME! Technically we either can dump the copy constructor or the
    //        assignment operator.

    PositionsFunctor() { }
    
    PositionsFunctor(const RectilinearMesh<Dim, T> &m, 
                     const Centering<Dim> &c)
      : centering_m(c.position(0))
      {
	for (int i=0; i<Dim; i++) {
	  positions_m[i].engine() = m.positions()[i].engine();
	  spacings_m[i].engine() = m.spacings()[i].engine();
	}
      }

    PositionsFunctor(const PositionsFunctor &m)
      :	centering_m(m.centering_m)
    {
      for (int i=0; i<Dim; i++) {
	positions_m[i].engine() = m.positions_m[i].engine();
	spacings_m[i].engine() = m.spacings_m[i].engine();
      }
    }

    PositionsFunctor& operator=(const PositionsFunctor &m)
    {
      centering_m = m.centering_m;
      for (int i=0; i<Dim; i++) {
	positions_m[i].engine() = m.positions_m[i].engine();
	spacings_m[i].engine() = m.spacings_m[i].engine();
      }

      return *this;
    }

    inline PointType_t operator()(int i0) const
      {
        return PointType_t(positions_m[0].read(i0) + spacings_m[0].read(i0)*centering_m(0));
      }
      
    inline PointType_t operator()(int i0, int i1) const
      {
        return PointType_t(positions_m[0].read(i0) + spacings_m[0].read(i0)*centering_m(0),
			   positions_m[1].read(i1) + spacings_m[1].read(i1)*centering_m(1));
      }

    inline PointType_t operator()(int i0, int i1, int i2) const
      {
        return PointType_t(positions_m[0].read(i0) + spacings_m[0].read(i0)*centering_m(0),
			   positions_m[1].read(i1) + spacings_m[1].read(i1)*centering_m(1),
			   positions_m[2].read(i2) + spacings_m[2].read(i2)*centering_m(2));
      }

  private:

    SpacingsType_t positions_m;
    SpacingsType_t spacings_m;
    typename Centering<Dim>::Position centering_m;

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


  /// General "volume" functor: works for edges, faces and cells.

  class GeneralVolumesFunctor {
  public:
  
    /// Need to be able to default construct since we fill in the details
    /// after the fact.
    
    GeneralVolumesFunctor() { }
    
    GeneralVolumesFunctor(const RectilinearMesh<Dim, T> &m, 
			  const Centering<Dim> &c)
      : orientation_m(c.orientation(0))
      {
	for (int i=0; i<Dim; i++)
	  spacings_m[i].engine() = m.spacings()[i].engine();
      }

    GeneralVolumesFunctor(const GeneralVolumesFunctor &m)
      : orientation_m(m.orientation_m)
      {
	for (int i=0; i<Dim; i++)
	  spacings_m[i].engine() = m.spacings_m[i].engine();	
      }

    GeneralVolumesFunctor& operator=(const GeneralVolumesFunctor &m)
    {
      orientation_m = m.orientation_m;
      for (int i=0; i<Dim; i++)
	spacings_m[i].engine() = m.spacings_m[i].engine();
      return *this;
    }

    inline T operator()(int i0) const
      {
	// It does not make sense to have a zero orientation for 1D
	return spacings_m[0].read(i0);
      }
      
    inline T operator()(int i0, int i1) const
      {
	// It does not make sense to have all zero orientations for 2D
	if (orientation_m[0].first() == 0)
	  return spacings_m[1].read(i1);
	else if (orientation_m[1].first() == 0)
	  return spacings_m[0].read(i0);
	else
	  return spacings_m[0].read(i0) * spacings_m[1].read(i1);
      }

    inline T operator()(int i0, int i1, int i2) const
      {
	// Could optimize as above
	T volume = static_cast<T>(1);
	if (orientation_m[0].first() != 0)
	  volume *= spacings_m[0].read(i0);
	if (orientation_m[1].first() != 0)
	  volume *= spacings_m[1].read(i1);
	if (orientation_m[2].first() != 0)
	  volume *= spacings_m[2].read(i2);
	return volume;
      }

  private:

    SpacingsType_t spacings_m;
    typename Centering<Dim>::Orientation orientation_m;

  };

  
  //---------------------------------------------------------------------------
  /// Support for the cellVolumes() function. We also need to export the 
  /// CellVolumesEngineTag_t typedef and the initializeCellVolumes() member 
  /// function, which sets the appropriate constant value for the volume. 
  /// The indices passed in refer to cells.
  
  typedef IndexFunction<GeneralVolumesFunctor> CellVolumesEngineTag_t;
  
  void initializeCellVolumes(
    Engine<Dim, T, CellVolumesEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be cell-centered.

      PAssert(c.size() == 1);
      PAssert(c.centeringType() == CellType);

      // Use the general functor to do the job.

      e.setFunctor(GeneralVolumesFunctor(*this, c));
    }

  //---------------------------------------------------------------------------
  /// Support for the faceAreas() function. We also need to export the 
  /// FaceAreasEngineTag_t typedef and the initializeFaceAreas() member 
  /// function, which sets the appropriate constant face area value.
  /// The indices passed in refer to cells.

  typedef IndexFunction<GeneralVolumesFunctor> FaceAreasEngineTag_t;
  
  void initializeFaceAreas(
    Engine<Dim, T, FaceAreasEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be face-centered.

      PAssert(c.size() == 1);
      PAssert(c.centeringType() == FaceType);

      // Use the general functor to do the job.
      
      e.setFunctor(GeneralVolumesFunctor(*this, c));
    }
  
  //---------------------------------------------------------------------------
  /// Support for the edgeLengths() function. We also need to export the 
  /// EdgeLengthsEngineTag_t typedef and the initializeEdgeLengths() member 
  /// function, which sets the appropriate constant edge length value.
  /// The indices passed in refer to cells.

  typedef IndexFunction<GeneralVolumesFunctor> EdgeLengthsEngineTag_t;
  
  void initializeEdgeLengths(
    Engine<Dim, T, EdgeLengthsEngineTag_t> &e, 
    const Centering<Dim> &c) const
    {
      // Check some pre-conditions. We need there to be a single centering
      // point and it must be edge-centered.
      
      PAssert(c.size() == 1);
      PAssert(c.centeringType() == EdgeType);
      
      // Use the general functor to do the job.
      
      e.setFunctor(GeneralVolumesFunctor(*this, c));
    }
  
private:

  /// Our data, stored as a ref-counted pointer to simplify memory management.
  
  RefCountedPtr<RectilinearMeshData<Dim, T> > data_m;
};

#endif // POOMA_FIELD_MESH_RECTILINEARMESH_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RectilinearMesh.h,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/23 16:20:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo

