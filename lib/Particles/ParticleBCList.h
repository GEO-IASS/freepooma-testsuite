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
//   ParticleBCList
//-----------------------------------------------------------------------------

#ifndef POOMA_PARTICLES_PARTICLEBCLIST_H
#define POOMA_PARTICLES_PARTICLEBCLIST_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * ParticleBCList holds a list of ParticleBCItems for the Particles class,
 * storing them in an STL vector of ParticleBCItem pointers.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Particles/ParticleBCItem.h"
#include "Utilities/PAssert.h"

#include <vector>
#include <iosfwd>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

class ParticleBCList
{
public:

  // Typedefs.
  typedef std::vector<ParticleBCItem*>  BCContainer_t;
  typedef BCContainer_t::size_type      Size_t;

  // Constructors.
  ParticleBCList();
  ParticleBCList(const ParticleBCList&);

  // Destructor.
  ~ParticleBCList();

  // Methods.

  // Retrun number of boundary conditions.
  Size_t size() const { return bcList_m.size(); }

  // Access the ith boundary condition.
  const ParticleBCItem* operator()(Size_t i) const
  {
    PAssert(i < size());
    return bcList_m[i];
  }
  ParticleBCItem* operator()(Size_t i)
  {
    PAssert(i < size());
    return bcList_m[i];
  }

  // Add a boundary condition to the end of the list.
  // Returns the index number of the newly added ParticleBC.
  template <class Subject, class Object, class BCType>
  Size_t addBC(const Subject& s, const Object& o, const BCType& bc)
  {
    // BCType creates ParticleBC, which is added to list.
    bcList_m.push_back( bc.create(s, o) );
    return (size() - 1);
  }
  // Same, except one argument acts as both subject and object.
  template <class Subject, class BCType>
  Size_t addBC(const Subject& s, const BCType& bc)
  {
    bcList_m.push_back( bc.create(s) );
    return (size() - 1);
  }

  // Remove boundary condition with given index number from the list.
  void removeBC(Size_t);

  // Print each ParticleBC in the list to an ostream.
  void print(std::ostream&) const;

private:

  // The actual vector of ParticleBCItem pointers.
  BCContainer_t bcList_m;
};


// declare operator<< for a std::ostream and a ParticleBCList

std::ostream&
operator<<(std::ostream&, const ParticleBCList&);


#endif     // POOMA_PARTICLES_PARTICLEBCLIST_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ParticleBCList.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
