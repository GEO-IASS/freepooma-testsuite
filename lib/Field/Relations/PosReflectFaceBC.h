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
//   PosReflectFaceBC
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_RELATIONS_POSREFLECTFACEBC_H
#define POOMA_FIELD_RELATIONS_POSREFLECTFACEBC_H

/** @file
 * @ingroup Relations
 * @brief
 * Relation functor class setting all guard layers beyond a
 * specified (logically) rectilinear mesh face to a
 * positively reflected value.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

#include "Field/Relations/Relations.h"

/** PosReflectFaceBC is an Relation functor class.
 * 
 * It represents a Dirichlet boundary condition on a logically rectilinear
 * domain where the value on that face is obtained by positive reflection.
 * A constructor switch allows the BC to enforce that the mesh-boundary value
 * is zero; this affects only vertex-centered Field values/components because
 * the boundary is defined to be the last vertex.
 */ 

template<int Dim>
class PosReflectFaceBC
{
public:

  //---------------------------------------------------------------------------
  // Constructors. 
  
  PosReflectFaceBC(int face, bool enforceZeroBoundary = false) 
  : domain_m(Pooma::NoInit()),
    vertFaceDomain_m(Pooma::NoInit()),
    srcRange_m(Pooma::NoInit()),
    face_m(face), 
    enforceZeroBoundary_m(enforceZeroBoundary)  
    { }
    
  PosReflectFaceBC(const PosReflectFaceBC<Dim> &model) 
  : domain_m(model.domain_m),
    vertFaceDomain_m(model.vertFaceDomain_m),
    srcRange_m(model.srcRange_m),
    face_m(model.face_m), 
    enforceZeroBoundary_m(model.enforceZeroBoundary_m)
    { }

  template<class Target>    
  PosReflectFaceBC(const PosReflectFaceBC<Dim> &init, const Target &t) 
  : domain_m(t.totalDomain()),
    vertFaceDomain_m(t.totalDomain()),
    srcRange_m(Pooma::NoInit()),
    face_m(init.face_m), 
    enforceZeroBoundary_m(init.enforceZeroBoundary_m)
  {
    // This only makes sense if the target has no sub-fields.
	  
    PAssert(t.numSubFields() == 0);
	 
    // Note: a convertor from Interval<Dim> Range<Dim> would be handy here.
    // srcRange is used to get data to copy from.

    for (int dd = 0; dd < Dim; ++dd) 
      {
        srcRange_m[dd] = 
          Range<1>(domain_m[dd].min(), domain_m[dd].max(), 1);
      }
	 
    // Get the direction.
	 
    int d = face_m / 2;
	 
    // The other directions span the subject's total domain. 
    // Therefore, we just chop out the guard layers, taking care to 
    // handle the case where we are enforcing a zero boundary 
    // (appropriate only for vert-centering).

    // FIXME: what about discontinuous fields?
	 
    int adjust = 1 - t.centering().orientation(0)[d].min();

    // Select the high or low face.

    if (face_m & 1) 
      {
        // High face.
        // Get the number of guard layers in the upper direction.
        
        int nGuards = t.fieldEngine().guardLayers().upper(d);
	if (enforceZeroBoundary_m && adjust == 1) 
	  {
            vertFaceDomain_m[d] =
	      Interval<1>(t.physicalDomain()[d].max(), 
			  t.physicalDomain()[d].max());
          }

        // Adjust the domains.
	              
        srcRange_m[d] = 
          Range<1>(t.physicalDomain()[d].max() - adjust, 
		   t.physicalDomain()[d].max() - adjust - (nGuards - 1),
		   -1);

        domain_m[d] = Interval<1>(domain_m[d].max() - (nGuards - 1), 
				  domain_m[d].max());
      } 
    else 
      {
        // Low face.   
	// Get the number of guard layers in the lower direction.
	
        int nGuards = t.fieldEngine().guardLayers().lower(d);
	if (enforceZeroBoundary_m && adjust == 0) 
	  {
            vertFaceDomain_m[d] =
              Interval<1>(t.physicalDomain()[d].min(), 
			  t.physicalDomain()[d].min());
	   }
	     
        // Adjust the domains.
        
        srcRange_m[d] = 
          Range<1>(t.physicalDomain()[d].min() + adjust +
		   (nGuards - 1), 
		   t.physicalDomain()[d].min() + adjust, -1);
	domain_m[d] = Interval<1>(domain_m[d].min(), 
				  domain_m[d].min() + (nGuards - 1));
      }
  }    

  //---------------------------------------------------------------------------
  // Assignment operator. Does deep assignment.
  
  PosReflectFaceBC<Dim> &operator=(const PosReflectFaceBC<Dim> &rhs)
  {
    domain_m = rhs.domain_m;
    vertFaceDomain_m = rhs.vertFaceDomain_m;
    srcRange_m = rhs.srcRange_m;
    face_m = rhs.face_m;
    enforceZeroBoundary_m = rhs.enforceZeroBoundary_m;

    return *this;
  }

  /// Face this operates on.

  int face() const { return face_m; }

  /// Whether we enforce zero boundary by setting it so.

  bool enforceZeroBoundary() const { return enforceZeroBoundary_m; }

  //---------------------------------------------------------------------------
  // Update function.

  template<class Target>
  void operator()(const Target &t) const
  {  
    t(domain_m) = t(srcRange_m);

    if (enforceZeroBoundary_m &&
        t.centering().orientation(0)[face_m / 2].min() == 0)
      {
        typedef typename Target::Element_t T;
        t(vertFaceDomain_m) = T(0.0);
      }
  }
    
private:

  Interval<Dim> domain_m, vertFaceDomain_m;
  Range<Dim> srcRange_m;
  int face_m;
  bool enforceZeroBoundary_m;
};


//-----------------------------------------------------------------------------
// Override the default priority so that boundary conditions get executed
// last.
//-----------------------------------------------------------------------------

template<int Dim>
struct RelationFunctorTraits<PosReflectFaceBC<Dim> > {
  
  enum { defaultPriority = 100 }; 

};


namespace Pooma {

  /// addPosReflectFaceBC installs PosReflectFace boundary conditions on the
  /// specified face of every subfield of the Target.

  template<class Target>
  void addPosReflectFaceBC(const Target &f, int face, 
    bool enforceZeroBoundary = false)
  {
      newRelation(PosReflectFaceBC<Target::dimensions>
        (face, enforceZeroBoundary), f);
  }

  /// addAllPosReflectFaceBC installs PosReflectFace boundary conditions on all
  /// of the faces of every subfield of the Target.

  template<class Target>
  void addAllPosReflectFaceBC(const Target &f, bool enforceZeroBoundary = false)
  {
    for (int i = 0; i < 2 * Target::dimensions; i++)
      {
        addPosReflectFaceBC(f, i, enforceZeroBoundary);
      }
  }

} // namespace Pooma

#endif // POOMA_FIELD_RELATIONS_POSREFLECTFACEBC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PosReflectFaceBC.h,v $   $Author: richi $
// $Revision: 1.5 $   $Date: 2004/11/10 22:00:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
