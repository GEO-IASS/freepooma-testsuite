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
//   ParticleBCItem
//-----------------------------------------------------------------------------

#ifndef POOMA_PARTICLES_PARTICLEBCITEM_H
#define POOMA_PARTICLES_PARTICLEBCITEM_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * ParticleBCItem is an abstract base class for any boundary condition 
 * applied to a Particles object.
 *
 * The ParticleBCItems are stored in a
 * ParticleBCList.  ParticleBCItem provides the virtual interface for
 * applying a particle boundary condition.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include <iosfwd>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

class ParticleBCItem
{
public:

  // Constructors.
  ParticleBCItem() {}
  ParticleBCItem(const ParticleBCItem&) {}

  // Destructor.
  virtual ~ParticleBCItem() {}

  // Apply boundary condition, either to all the patches or just
  // to a specified one (if pid < 0, apply to all).
  virtual void applyBoundaryCondition(int pid) = 0;

  inline void applyBoundaryCondition()
    {
      applyBoundaryCondition(-1);
    }

  // Print out to a stream
  virtual void print(std::ostream&) const = 0;
};


// Specialization of operator<< for ParticleBCItem, which invokes print method.

inline
std::ostream&
operator<<(std::ostream& o, const ParticleBCItem& bcitem)
{
  bcitem.print(o);
  return o;
}


#endif     // POOMA_PARTICLES_PARTICLEBCITEM_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ParticleBCItem.h,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
