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
//   KillBC
//-----------------------------------------------------------------------------

#ifndef POOMA_PARTICLES_KILLBC_H
#define POOMA_PARTICLES_KILLBC_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * A "kill" boundary condition for Particles.
 *
 * When Attribute value goes outside boundary by some amount,
 * KillBC destroys this particle.
 * In this case, the subject is an Attribute, while the object
 * is the containing Particles object.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Particles/ParticleBC.h"
#include <iostream>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Here is an example of adding a new particle boundary condition.
 * First, add a BCType that inherits from ParticleBCType, describes 
 * a boundary condition, and stores any data that is needed.
 * Then, we specialize ParticleBC for this BCType.
 */

template <class T>
class KillBC : public ParticleBCType< KillBC<T> >
{
public:

  // Constructors.
  KillBC(T min, T max)
    : MinVal_m(min), MaxVal_m(max) {}
  KillBC(const KillBC<T>& model)
    : MinVal_m(model.min()), MaxVal_m(model.max()) {}

  // Destructor.
  ~KillBC() {}

  // Methods.

  // Assignment operator does deep copy.
  KillBC<T>& operator=(const KillBC<T>& rhs)
  {
    MinVal_m = rhs.min();
    MaxVal_m = rhs.max();
    return *this;
  }

  // Accessor functions
  T min() const { return MinVal_m; }
  T max() const { return MaxVal_m; }
  T& min() { return MinVal_m; }
  T& max() { return MaxVal_m; }

private:

  // extents of KillBC region
  T MinVal_m, MaxVal_m;
};



template <class Subject, class Object, class T>
class ParticleBC< Subject, Object, KillBC<T> >
  : public ParticleBCItem
{
public:

  // Typedefs.
  typedef Subject Subject_t;
  typedef Object Object_t;
  typedef KillBC<T> BCType_t;
  typedef ParticleBC<Subject_t,Object_t,BCType_t> This_t;

  // Constructors.
  ParticleBC(const Subject_t& s, const Object_t& o,
             const BCType_t& bc)
    : subject_m(s),
      object_m(const_cast<Object_t&>(o)),
      bc_m(bc) { }
  ParticleBC(const This_t& model)
    : subject_m(model.subject()),
      object_m(const_cast<Object_t&>(model.object())),
      bc_m(model.bc()) { }

  // Destructor.
  ~ParticleBC() { }

  // Methods.

  // Accessors for subject and object.

  const Subject_t& subject() const { return subject_m; }
  const Object_t& object() const { return object_m; }

  // Access boundary condition.
  const BCType_t& bc() const { return bc_m; }

  // Apply boundary condition.
  void applyBoundaryCondition(int pid);

  // Print out to an ostream.
  void print(std::ostream& o) const
    {
      o << "BC Type: Kill, Range: (" << bc_m.min()
        << "," << bc_m.max() << ")" << std::endl;
    }

private:

  // Can't call default constructor.
  ParticleBC();

  // Subject of the boundary condition. 
  Subject_t subject_m;

  // Object of the boundary condition.
  // Note that we store the object by reference here.
  Object_t& object_m;

  // Our boundary condition.
  BCType_t bc_m;

};

#include "Particles/KillBC.cpp"

#endif     // POOMA_PARTICLES_KILLBC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: KillBC.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
