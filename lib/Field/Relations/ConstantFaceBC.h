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
// Class:
//   ConstantFaceBC
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_RELATIONS_CONSTANTFACEBC_H
#define POOMA_FIELD_RELATIONS_CONSTANTFACEBC_H

/** @file
 * @ingroup Relations
 * @brief
 * Relation functor class setting all guard layers beyond a
 * specified (logically) rectilinear mesh face to a  
 * constant value.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

#include "Field/Relations/Relations.h"

/** ConstantFaceBC is an Relation functor class.
 * 
 * It represents a Dirichlet boundary condition on a logically rectilinear
 * domain where the value on that face is a constant. The setConstant()
 * function provides a means to have a time-dependent BC.  A constructor switch
 * allows the BC to enforce that the mesh-boundary value is zero; this affects
 * only vertex-centered Field values/components because the boundary is defined
 * to be the last vertex.
 * 
 * The T template parameter is the type of the constant value.
 */

template<int Dim, class T>
class ConstantFaceBC
{
public:

  //---------------------------------------------------------------------------
  // Constructors. 
  
  ConstantFaceBC(int face, const T &constant, 
    bool enforceConstantBoundary = false) 
  : domain_m(Pooma::NoInit()),
    face_m(face),
    constant_m(constant),
    enforceConstantBoundary_m(enforceConstantBoundary)  
    { }
    
  ConstantFaceBC(const ConstantFaceBC<Dim, T> &model) 
  : domain_m(model.domain_m),
    face_m(model.face_m), 
    constant_m(model.constant_m),
    enforceConstantBoundary_m(model.enforceConstantBoundary_m)  
    { }

  template<class Target>    
  ConstantFaceBC(const ConstantFaceBC<Dim, T> &init, const Target &t) 
  : domain_m(t.totalDomain()),
    face_m(init.face_m), 
    constant_m(init.constant_m),
    enforceConstantBoundary_m(init.enforceConstantBoundary_m)  
  {
    // This only makes sense if the target has no sub-fields.
	  
    PAssert(t.numSubFields() == 0);
	 
    // Get the direction.
	 
    int d = face_m / 2;
	 
    // The other directions span the subject's total domain. 
    // Therefore, we just chop out the guard layers, taking care to 
    // handle the case where we are enforcing a constant boundary 
    // (appropriate only for vert-centering).

    // FIXME: what about discontinuous fields?
	 
    int adjust;
    if (enforceConstantBoundary_m &&
        t.centering().orientation(0)[d].min() == 0)
      adjust = 0;
    else
      adjust = 1;
        
    // Select the high or low face.

    if (face_m & 1) 
      {
        // High face.
        // Get the number of guard layers in the upper direction.
        
        int nGuards = t.fieldEngine().guardLayers().upper(d);
          
        // Adjust the domain.
                   
        domain_m[d] = Interval<1>(domain_m[d].max() - nGuards + adjust, 
			          domain_m[d].max());
      } 
    else 
      {
        // Low face.   
	// Get the number of guard layers in the lower direction.

        int nGuards = t.fieldEngine().guardLayers().lower(d);
          
        // Adjust the domain.
                   	
        domain_m[d] = Interval<1>(domain_m[d].min(), 
				  domain_m[d].min() + nGuards - adjust);
      }
  }    

  //---------------------------------------------------------------------------
  // Assignment operator. Does deep assignment.
  
  ConstantFaceBC<Dim, T> &operator=(const ConstantFaceBC<Dim, T> &rhs)
  {
    domain_m = rhs.domain_m;
    face_m = rhs.face_m;
    constant_m = rhs.constant_m;
    enforceConstantBoundary_m = rhs.enforceConstantBoundary_m; 

    return *this;
  }

  //---------------------------------------------------------------------------
  /// Constant we set the boundary to.

  T constant() const { return constant_m; }

  /// User may want to change the constant's value, e.g., for time-dependence.

  void setConstant(T newConstant) { constant_m = newConstant; }

  /// Face we operate on.

  int face() const { return face_m; }

  //---------------------------------------------------------------------------
  // Update function.

  template<class Target>
  void operator()(const Target &t) const
  {  
    t(domain_m) = constant_m;
  }
    
private:

  Interval<Dim> domain_m;
  int face_m;
  T constant_m;
  bool enforceConstantBoundary_m;
};


//-----------------------------------------------------------------------------
// Override the default priority so that boundary conditions get executed
// last.
//-----------------------------------------------------------------------------

template<int Dim, class T>
struct RelationFunctorTraits<ConstantFaceBC<Dim, T> > {
  
  enum { defaultPriority = 100 }; 

};



namespace Pooma {

  /// addConstantFaceBC installs ConstantFace boundary conditions on the
  /// specified face of every subfield of the Target.

  template<class Target, class T>
  void addConstantFaceBC(const Target &f, int face, const T &constant,
    bool enforceConstantBoundary = false)
  {
    newRelation(ConstantFaceBC<Target::dimensions, T>
      (face, constant, enforceConstantBoundary), f);
  }

  /// addAllConstantFaceBC installs ConstantFace boundary conditions on all
  /// of the faces of every subfield of the Target.

  template<class Target, class T>
  void addAllConstantFaceBC(const Target &f, const T &constant, 
    bool enforceConstantBoundary = false)
  {
    for (int i = 0; i < 2 * Target::dimensions; i++)
      {
        addConstantFaceBC(f, i, constant, enforceConstantBoundary);
      }
  }
  
} // namespace Pooma

#endif // POOMA_FIELD_RELATIONS_CONSTANTFACEBC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ConstantFaceBC.h,v $   $Author: richi $
// $Revision: 1.5 $   $Date: 2004/11/10 22:00:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
