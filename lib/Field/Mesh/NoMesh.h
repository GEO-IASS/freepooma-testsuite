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
// Classes: 
//   NoMesh
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Mesh
 * @brief
 * NoMesh is an extremely lightweight class that indicates a Field cannot
 * answer mesh-type questions
 */

#ifndef POOMA_FIELD_MESH_NOMESH_H
#define POOMA_FIELD_MESH_NOMESH_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"          // NoMeshData<> class member
#include "Domain/Shrink.h"            // shrinkRight() in NoMeshData<> ctor
#include "Layout/INode.h"             // Used in ctors
#include "Field/FieldEngine/FieldEnginePatch.h" 
                                      // Used in ctors
#include "Utilities/RefCountedPtr.h"  // Reference counting used
#include "Utilities/RefCounted.h"     // Reference counting used

/**
 * Holds the data for a NoMesh mesh. That class has a ref-counted
 * instance of this class
 */

template <int Dim>
class NoMeshData : public RefCounted
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  /// We provide a default constructor that creates the object with empty
  /// domains. To be useful, this object must be replaced by another version 
  /// via assignment.
  
  NoMeshData() 
    { 
    }

  /// This constructor fully constructs the object. It uses the layout to set
  /// up the appropriate domains.
  ///
  /// The Layout supplied must refer to VERTEX positions.

  template<class Layout>  
  explicit NoMeshData(const Layout &layout)
  : physicalVertexDomain_m(layout.innerDomain()),
    physicalCellDomain_m(shrinkRight(physicalVertexDomain_m, 1)),
    totalVertexDomain_m(layout.domain()),
    totalCellDomain_m(shrinkRight(totalVertexDomain_m, 1))
    {
    }
    
  /// Copy constructor.

  NoMeshData(const NoMeshData<Dim> &model)
  : physicalVertexDomain_m(model.physicalVertexDomain_m),
    physicalCellDomain_m(model.physicalCellDomain_m),
    totalVertexDomain_m(model.totalVertexDomain_m),
    totalCellDomain_m(model.totalCellDomain_m)
    {
    }
    
  /// @name View constructors.
  //@{
  
  /// Interval view. For now, we simply make the zero-based 
  /// total domain == physical domain.
  ///
  /// The Interval supplied must refer to VERTEX positions.
  
  NoMeshData(const Interval<Dim> &d)
  : physicalVertexDomain_m(d - d.firsts()),
    physicalCellDomain_m(shrinkRight(physicalVertexDomain_m, 1)),
    totalVertexDomain_m(d - d.firsts()),
    totalCellDomain_m(shrinkRight(totalVertexDomain_m, 1))
    {
    }
    
  /// FieldEnginePatch constructor.
  ///
  /// The FieldEnginePatch supplied must refer to VERTEX positions.

  NoMeshData(const NoMeshData<Dim> &model, const FieldEnginePatch<Dim> &p)
  : physicalVertexDomain_m(p.domain_m),
    physicalCellDomain_m(shrinkRight(physicalVertexDomain_m, 1)),
    totalVertexDomain_m(p.domain_m),
    totalCellDomain_m(shrinkRight(totalVertexDomain_m, 1))
    {
    }

  //@}

  //---------------------------------------------------------------------------
  /// Copy assignment operator.
  
  NoMeshData<Dim> &operator=(const NoMeshData<Dim> &rhs)
    {
      if (this != &rhs)
        {
          physicalVertexDomain_m = rhs.physicalVertexDomain_m;
          physicalCellDomain_m = rhs.physicalCellDomain_m;
          totalVertexDomain_m = rhs.totalVertexDomain_m;
          totalCellDomain_m = rhs.totalCellDomain_m;
        }
      return *this;
    }

  //---------------------------------------------------------------------------
  /// Empty destructor is fine. However, note that it is not virtual. So,
  /// even though we are inheriting implementation from this class, we must
  /// take care not to delete through a pointer to this base class.

  ~NoMeshData() 
    { 
      // This space intentionally left blank.
    }

  //---------------------------------------------------------------------------
  // General accessors.
  
  /// @name Domains.
  //@{
  
  inline const Interval<Dim> &physicalVertexDomain() const 
    { 
      return physicalVertexDomain_m; 
    }
  
  inline const Interval<Dim> &physicalCellDomain() const 
    { 
      return physicalCellDomain_m; 
    }
  
  inline const Interval<Dim> &totalVertexDomain() const 
    { 
      return totalVertexDomain_m; 
    }
  
  inline const Interval<Dim> &totalCellDomain() const 
    { 
      return totalCellDomain_m; 
    }

  //@}

private:

  /// Domains.

  Interval<Dim> physicalVertexDomain_m, physicalCellDomain_m;
  Interval<Dim> totalVertexDomain_m, totalCellDomain_m;

};

/**
 * NoMesh is an extremely lightweight class that indicates a Field cannot
 * answer mesh-type questions. When a Field has a NoMesh, it has the flavor
 * of a "multi-array"; that is, an array with multiple engines.
 */

template<int Dim>
class NoMesh
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.

  enum { dimensions = Dim };

  //---------------------------------------------------------------------------
  // Constructors.

  /// We provide a default constructor that creates the object with empty
  /// domains. To be useful, this object must be replaced by another version 
  /// via assignment.
  
  inline NoMesh() 
  : data_m(new NoMeshData<Dim>)
    { 
    }

  /// This constructor fully constructs the object using the layout to
  /// compute domains.
  ///
  /// The Layout supplied must refer to VERTEX positions.
  
  template<class Layout>
  inline explicit NoMesh(const Layout &layout)
  : data_m(new NoMeshData<Dim>(layout))
    { 
    }
    
  /// Copy constructor. 
  
  inline NoMesh(const NoMesh<Dim> &model)
  : data_m(model.data_m)
    {
    }
    
  /** @name View constructors.
   */
  //@{
  
  /// Interval view.
  ///
  /// The Interval supplied must refer to VERTEX positions.
  
  inline NoMesh(const NoMesh<Dim> &model, const Interval<Dim> &d)
  : data_m(new NoMeshData<Dim>(d))
    {
    }
  
  /// INode view.
  ///
  /// The INode supplied must refer to VERTEX positions.
  
  NoMesh(const NoMesh<Dim> &model, const INode<Dim> &i)
  : data_m(new NoMeshData<Dim>(i.domain()))
    {
    }
  
  /// FieldEnginePatch view.
  ///
  /// The FieldEnginePatch supplied must refer to VERTEX positions.
  
  inline NoMesh(const NoMesh<Dim> &model, const FieldEnginePatch<Dim> &p)
  : data_m(new NoMeshData<Dim>(*model.data_m, p))
    {
    }
    
  /// General view. Made, for instance, by taking a Range-view of some
  /// other mesh.
  ///
  /// The Domain supplied must refer to VERTEX positions.
  
  template<class Mesh, class Domain>
  inline NoMesh(const Mesh &, const Domain &d)
    {
      Interval<Dim> dom;
      for (int i = 0; i < Dim; i++)
        dom[i] = d[i].size();
      data_m = new NoMeshData<Dim>(dom);
    }
  
  //@}

  //---------------------------------------------------------------------------
  /// Empty destructor is fine. The pointer to the data is ref-counted so its
  /// lifetime is correctly managed.
  
  ~NoMesh() { }

  //---------------------------------------------------------------------------
  /// Copy assignment operator.
  
  NoMesh<Dim> &
  operator=(const NoMesh<Dim> &rhs)
    {
      if (&rhs != this)
        {
          data_m = rhs.data_m;
        }
      
      return *this;
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

private:

  /// Our data, stored as a ref-counted pointer to simplify memory management.
  
  RefCountedPtr<NoMeshData<Dim> > data_m;
   
};

#endif // POOMA_FIELD_MESH_NOMESH_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NoMesh.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:46 $
// ----------------------------------------------------------------------
// ACL:rcsinfo


