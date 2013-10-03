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
//   ParticleBCType, ParticleBC, ParticleCompBC
//-----------------------------------------------------------------------------

#ifndef POOMA_PARTICLES_PARTICLEBC_H
#define POOMA_PARTICLES_PARTICLEBC_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * ParticleBCType is a tag-like base class for a category of boundary
 * condition applied to an Attribute of a Particles object.
 *
 * It provides
 * a factory method for creating the desired boundary condition object.
 * Each subclass of ParticleBCType gives itself as the template parameter
 * to ParticleBCType, so that this class knows what to create.
 *
 * ParticleBC is the class representing a generalized particle boundary
 * condition.  It has a subject type, object type, and a ParticleBCType.
 * The subject can be a DynamicArray or a view or expression involving
 * DynamicArrays.  The subject is examined to determine if its value is 
 * outside the bounds of the ParticleBC.  If the subject is out of bounds,
 * the object of the ParticleBC is modified accordingly.  Thus, the object
 * must be modifiable in the appropriate way.  The object of a ParticleBC
 * can never be an expression, and for a KillBC, the object must be a 
 * DynamicArray or a Particles object.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Particles/ParticleBCItem.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// General template for ParticleBC.  We will specialize on BCType.
//-----------------------------------------------------------------------------

template <class Subject, class Object, class BCType>
class ParticleBC : public ParticleBCItem
{
};


//-----------------------------------------------------------------------------
// Class ParticleBCType:
//-----------------------------------------------------------------------------

template <class BCType>
class ParticleBCType
{
public:

  // Constructors.
  ParticleBCType() {}
  ParticleBCType(const ParticleBCType<BCType>&) {}

  // Destructor.
  virtual ~ParticleBCType() {}

  // Methods.

  // Assignment operator does nothing.
  ParticleBCType<BCType>& operator=(const ParticleBCType<BCType>&)
  {
    return *this;
  }

  // Factory method to create the boundary condition.
  template <class Subject, class Object>
  ParticleBCItem* create(const Subject& s, const Object& o) const
  {
    typedef ParticleBC<Subject,Object,BCType> PBC_t;
    return new PBC_t(s, o, static_cast<const BCType&>(*this)); 
  }

  // Factory method to create the boundary condition with just a subject.
  // We assume that the argument acts as both subject and object.
  template <class Subject>
  ParticleBCItem* create(const Subject& s) const
  {
    typedef ParticleBC<Subject,Subject,BCType> PBC_t;
    return new PBC_t(s, s, static_cast<const BCType&>(*this)); 
  }
};

//-----------------------------------------------------------------------------
// ParticleCompBC is used to construct boundary conditions that work on
// particular components of a multi-component Particle attribute. The idea
// is to, for example, build a ReflectBC object and then use this to create
// a ParticleCompBC<1, ReflectBC> object, which builds a boundary condition
// by taking a one-dimensional component view of the specified subject and
// object. This means that particle  boundary condition objects need only
// work on scalar fields.
//-----------------------------------------------------------------------------

// General template

template <int Dim, class BCType>
class ParticleCompBC { };

// Partial specialization for one-dimensional elements

template <class BCType>
class ParticleCompBC<1,BCType>
{
public:

  // Constructors. Need to store a copy of the bc type and which component
  // we're dealing with.
  
  ParticleCompBC(const BCType& bc, int c1) 
    : bc_m(bc), c1_m(c1) { }
  ParticleCompBC(const ParticleCompBC<1,BCType>& model)
    : bc_m(model.bc()), c1_m(model.comp1()) { }
  
  // Assignment operator. Does deep assignment.
  
  ParticleCompBC<1,BCType>& operator=(const ParticleCompBC<1,BCType>& rhs)
  {
    bc_m = rhs.bc();
    c1_m = rhs.comp1();
    return *this;
  }
  
  // Accessors for data members.
  
  const BCType& bc() const { return bc_m; }
  int comp1() const { return c1_m; }
  
  // Factory method. Actually makes the boundary condition.
  // Takes a component view of the subject and object. 
  
  template <class Subject, class Object>
  ParticleBCItem* create(const Subject& s, const Object& o) const
  {
    return create2(s.comp(comp1()),o.comp(comp1()));
  }
  
  // Factory method. Actually makes the boundary condition.
  // Takes a component view of the subject and object. 
  // This one takes only a subject.
  
  template <class Subject>
  ParticleBCItem* create(const Subject& s) const
  {
    return create2(s.comp(comp1()),s.comp(comp1()));
  }
  
private:

  // Super-secret factory method. Actually makes the boundary condition. Called
  // by create(). We use this two-step approach so we don't have to worry about
  // explicitly specifying the subject type and object type coming from the
  // component-view operations.
  
  template <class ComponentSubject, class ComponentObject>
  ParticleBCItem* create2(const ComponentSubject& s,
                          const ComponentObject& o) const
  {
    typedef ParticleBC<ComponentSubject,ComponentObject,BCType> PBC_t;
    return new PBC_t(s, o, bc());
  }  

  BCType bc_m;
  int c1_m;
};


// Partial specialization for two-dimensional elements

template <class BCType>
class ParticleCompBC<2,BCType>
{
public:

  // Constructors. Need to store a copy of the bc type and which component
  // we're dealing with.
  
  ParticleCompBC(const BCType& bc, int c1, int c2) 
    : bc_m(bc), c1_m(c1), c2_m(c2) { }
  ParticleCompBC(const ParticleCompBC<2,BCType>& model)
    : bc_m(model.bc()), c1_m(model.comp1()), c2_m(model.comp2()) { }
  
  // Assignment operator. Does deep assignment.
  
  ParticleCompBC<2,BCType>& operator=(const ParticleCompBC<2,BCType>& rhs)
  {
    bc_m = rhs.bc();
    c1_m = rhs.comp1();
    c2_m = rhs.comp2();
    return *this;
  }
  
  // Accessors for data members.
  
  const BCType& bc() const { return bc_m; }
  int comp1() const { return c1_m; }
  int comp2() const { return c2_m; }
  
  // Factory method. Actually makes the boundary condition.
  // Takes a component view of the subject and object. 
  
  template <class Subject, class Object>
  ParticleBCItem* create(const Subject& s, const Object& o) const
  {
    return create2(s.comp(comp1(),comp2()),o.comp(comp1(),comp2()));
  }
  
  // Factory method. Actually makes the boundary condition.
  // Takes a component view of the subject and object. 
  // This one takes only a subject.
  
  template <class Subject>
  ParticleBCItem* create(const Subject& s) const
  {
    return create2(s.comp(comp1(),comp2()),s.comp(comp1(),comp2()));
  }
  
private:

  // Super-secret factory method. Actually makes the boundary condition. Called
  // by create(). We use this two-step approach so we don't have to worry about
  // explicitly specifying the subject type and object type coming from the
  // component-view operations.
  
  template <class ComponentSubject, class ComponentObject>
  ParticleBCItem* create2(const ComponentSubject& s,
                          const ComponentObject& o) const
  {
    typedef ParticleBC<ComponentSubject,ComponentObject,BCType> PBC_t;
    return new PBC_t(s, o, bc());
  }  

  BCType bc_m;
  int c1_m, c2_m;
};


#endif     // POOMA_PARTICLES_PARTICLEBC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ParticleBC.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
