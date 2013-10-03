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
//   PeriodicFaceBC
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_RELATIONS_PERIODICFACEBC_H
#define POOMA_FIELD_RELATIONS_PERIODICFACEBC_H

/** @file
 * @ingroup Relations
 * @brief
 * Updater setting all guard layers beyond a
 * specified (logically) rectilinear mesh face to the value
 * from the non-guard element symmetrically across the face
 * (the face is defined at the last vertex).
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"

#include "Field/Relations/Relations.h"

/** PeriodicFaceBC is an Relation functor class.
 * 
 * It represents a periodic boundary condition in one direction (perpendicular
 * to the specified face) of a logically rectilinear domain.
 */

template<int Dim>
class PeriodicFaceBC
{
public:

  //---------------------------------------------------------------------------
  // Constructors. 
  
  PeriodicFaceBC(int face) 
  : domain_m(Pooma::NoInit()),
    srcDomain_m(Pooma::NoInit()),
    face_m(face)
    { }
    
  PeriodicFaceBC(const PeriodicFaceBC<Dim> &model) 
  : domain_m(model.domain_m),
    srcDomain_m(model.srcDomain_m),
    face_m(model.face_m)
    { }

  template<class Target>    
  PeriodicFaceBC(const PeriodicFaceBC<Dim> &init, const Target &t) 
  : domain_m(t.totalDomain()),
    srcDomain_m(t.totalDomain()),
    face_m(init.face_m)
  {
    // This only makes sense if the target has no sub-fields.
	  
    PAssert(t.numSubFields() == 0);
	 
    // Get the direction.
	 
    int d = face_m / 2;

    // Check if we're on a vertex in the current direction.
	 	 
    int adjust = 1 - t.centering().orientation(0)[d].min();
        
    // Select the high or low face.

    if (face_m & 1) 
      {
        // High face.
        // Get the number of guard layers in the upper direction.
        
        int nGuards = t.fieldEngine().guardLayers().upper(d);
          
        // Adjust the domain.
                   
	domain_m[d] = 
          Interval<1>(domain_m[d].max() - (nGuards - 1 + adjust), 
	              domain_m[d].max());

	// The source domain is just the destination domain offset by the
	// periodicity length (number of *cells*):
	
	srcDomain_m[d] = 
	  Interval<1>(domain_m[d].min() - 
		      (t.physicalDomain()[d].length() - adjust),
		       domain_m[d].max() - 
		       (t.physicalDomain()[d].length() - adjust));

      } 
    else 
      {
        // Low face.   
	// Get the number of guard layers in the lower direction.

        int nGuards = t.fieldEngine().guardLayers().lower(d);
          
        // Adjust the domain.
                   	
	domain_m[d] = Interval<1>(domain_m[d].min(), 
				  domain_m[d].min() + (nGuards - 1));

	// The source domain is just the destination domain offset by the
	// periodicity length (number of *cells*):
	
	srcDomain_m[d] = 
	  Interval<1>(domain_m[d].min() + 
		      (t.physicalDomain()[d].length() - adjust),
		       domain_m[d].max() + 
		       (t.physicalDomain()[d].length() - adjust));
      }
  }    

  //---------------------------------------------------------------------------
  // Assignment operator. Does deep assignment.
  
  PeriodicFaceBC<Dim> &operator=(const PeriodicFaceBC<Dim> &rhs)
  {
    if (&rhs != this)
      {
        domain_m = rhs.domain_m;
        srcDomain_m = rhs.srcDomain_m;
        face_m = rhs.face_m;
      }

    return *this;
  }

  /// Face this operates on.

  int face() const { return face_m; }

  //---------------------------------------------------------------------------
  // Update function.

  template<class Target>
  void operator()(const Target &t) const
  {  
    t(domain_m) = t(srcDomain_m);
  }
    
private:

  Interval<Dim> domain_m, srcDomain_m;
  int face_m;
};


//-----------------------------------------------------------------------------
// Override the default priority so that boundary conditions get executed
// last.
//-----------------------------------------------------------------------------

template<int Dim>
struct RelationFunctorTraits<PeriodicFaceBC<Dim> > {
  
  enum { defaultPriority = 100 }; 

};


namespace Pooma {

  /// addPeriodicFaceBC installs PeriodicFace boundary conditions on the
  /// specified face of every subfield of the Target.

  template<class Target>
  void addPeriodicFaceBC(const Target &f, int face)
  {
    newRelation(PeriodicFaceBC<Target::dimensions>(face), f);
  }

  /// addAllPeriodicFaceBC installs PeriodicFace boundary conditions on all
  /// of the faces of every subfield of the Target.

  template<class Target>
  void addAllPeriodicFaceBC(const Target &f)
  {
    for (int i = 0; i < 2 * Target::dimensions; i++)
      {
        addPeriodicFaceBC(f, i);
      }
  }

} // namespace Pooma

#endif // POOMA_FIELD_RELATIONS_PERIODICFACEBC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PeriodicFaceBC.h,v $   $Author: richi $
// $Revision: 1.8 $   $Date: 2004/11/10 22:00:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
