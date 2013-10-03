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

#ifndef POOMA_PARTICLES_PARTICLES_H
#define POOMA_PARTICLES_PARTICLES_H

//-----------------------------------------------------------------------------
// Classes:
//   Particles<Traits>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Particles
 * @brief
 * The Particles class is the base class to be used for development of
 * user-defined particles classes that describe application-specific 
 * particle populations.
 */

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

#include "Particles/AttributeList.h"
#include "Particles/ParticleBCList.h"
#include "DynamicArray/DynamicArray.h"
#include "Engine/DynamicEngine.h"
#include "Layout/DynamicEvents.h"
#include <iosfwd>

//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * The Particles class is the base class to be used for development of
 * user-defined particles classes that describe application-specific 
 * particle populations.
 *
 * Particles is templated on a ParticleTraits
 * class that provides typedefs for:
 *   - the particle layout strategy type (ParticleLayout_t)
 *   - the engine tag specifying an attrib engine type (AttributeEngineTag_t)
 *
 * From these types, Particles can determine the geometry of the
 * system, and the type of layout that should be used by the attributes.
 * Particles provides the basic interface for collecting a set of
 * heterogeneous Attributes, each of which describes a particular
 * characteristic of all the particles.  With this interface, the user 
 * can create or destroy particles, apply boundary conditions and 
 * update the particle data distribution.
 *
 * Here is an example of what the ParticleTraits class might look like:
 *
 * @code
 * template <class FieldLayoutType, class GeometryType>
 * struct StdParticleTraits {
 *   typedef SpatialLayout<FieldLayoutType, GeometryType>  ParticleLayout_t;
 *   typedef MultiPatch AttributeEngineTag_t;
 * };
 * @endcode
 */

template <class ParticleTraits>
class Particles
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  /// Typedefs for this object, and its trait class type.

  typedef Particles<ParticleTraits>                     This_t;
  typedef ParticleTraits                                Traits_t;

  /// The data type used to refer to the number of particles

  typedef long                                          Size_t;

  /// The tag type for the engine used to store attributes.
  /// All attributes should use the same engine, since they will use the
  /// same layout object.

  typedef typename ParticleTraits::AttributeEngineTag_t AttributeEngineTag_t;

  /// The full particle layout type (not just a tag).  The particle layout
  /// determines what particles go where based on some algorithm.  The
  /// layout should have the type of geometry, dimensions, etc., for the
  /// simulation domain of the particles.

  typedef typename ParticleTraits::ParticleLayout_t     ParticleLayout_t;

  /// The attribute layout type.  To get it, we need the type of engine
  /// in a typical attribute, so we use a temporary engine type (TempEngine_t
  /// has no other purpose but to be used in the next line to get
  /// its layout type).

  typedef Engine<1, double, AttributeEngineTag_t>       TempEngine_t;
  typedef typename TempEngine_t::Layout_t               AttributeLayout_t;

  /// Type of data used to specify a patch number in dynamic operations

  typedef typename AttributeLayout_t::PatchID_t         PatchID_t;

  /// The type of domain used in the attribute layout

  typedef typename AttributeLayout_t::Domain_t          AttributeDomain_t;


  //============================================================
  // Constructors and initializers
  //============================================================

  /// Construct a Particles object, using the given particle layout object.
  /// Let the particle layout object initialize the attribute data layout.

  Particles(const ParticleLayout_t& pl);

  /// Default constructor: this should be used if
  ///   - The particle layout has a useful default constructor, and
  ///   - The attribute layout has a useful default constructor, and
  ///   - The initialize() method will be used later.

  Particles();

  /// Copy constructor: only copy the layout data and destroy flag, do
  /// NOT copy the attributes.  The new particle object starts out with
  /// just layouts and no attributes (except for the internal destroy list).

  Particles(const This_t &p);

  /// Initialize this object with the given particle layout.

  void initialize(const ParticleLayout_t& pl);


  //============================================================
  // Destructor
  //============================================================

  ~Particles();


  //============================================================
  /** @name Attribute accessors */
  //============================================================
  //@{

  /// Has this object been initialized?

  bool initialized() const
  {
    return particleLayout_m.initialized();
  }

  /// Return the total number of particles.

  Size_t size() const
  {
    return attributeLayout_m.domain().size();
  }

  /// Return the domain of the attributes of the Particles

  const AttributeDomain_t& domain() const
  {
    return attributeLayout().domain();
  }

  /// Return a reference to the particle layout object

  ParticleLayout_t& particleLayout()
  {
    return particleLayout_m;
  }

  const ParticleLayout_t& particleLayout() const
  {
    return particleLayout_m;
  }

  /// Return a reference to the attribute layout object

  AttributeLayout_t& attributeLayout()
  {
    return attributeLayout_m;
  }

  const AttributeLayout_t& attributeLayout() const
  {
    return attributeLayout_m;
  }

  /// Return the number of attributes

  AttributeList::Size_t attributes() const { return attributes_m.size(); }

  /// Return the Nth attribute

  Attribute* attribute(AttributeList::Size_t n)
  {
    return attributes_m.attribute(n);
  }

  const Attribute* attribute(AttributeList::Size_t n) const
  {
    return attributes_m.attribute(n);
  }

  //@}

  /// Return our current destroy method code

  int destroyMethod() const { return destroyMethod_m; }

  /// Change the current destroy method.  This changes a flag, which
  /// is examined during the next sync() or performDestroy() call to
  /// determine how to perform the cached destroy requests.  So changing
  /// this flag will affect all cached destroys, not just the ones issued
  /// after this call.

  template <class DestroyMethod>
  void setDestroyMethod(const DestroyMethod &)
  {
    destroyMethod_m = DestroyMethod::code;
  }

  /// Return the number of particles we have queued to destroy in the
  /// next performDestroy call, for all the patches (if pid < 0), or just
  /// for the specified patch.

  Size_t deferredDestroyAmount(PatchID_t pid = (-1)) const;

  /// Return a reference to the destroy list for local patch pid

  DynamicArray<int,Dynamic>& destroyList(int pid)
  {
    PAssert(pid >=0 && pid < attributeLayout_m.sizeLocal());
    return destroyList_m[pid];
  }

  const DynamicArray<int,Dynamic>& destroyList(int pid) const
  {
    PAssert(pid >=0 && pid < attributeLayout_m.sizeLocal());
    return destroyList_m[pid];
  }

  //============================================================
  /** @name Attribute modifier methods */
  //============================================================
  //@{

  /// Add a new attribute to this object ... it will be initialized
  /// with the proper attribute layout, and put in our list.  It must
  /// have the proper engine type.  The attribute can be given an
  /// optional name for use in later look-ups of the attribute.
  /// Return the index of the newly added attribute.

  template <class T>
  AttributeList::Size_t addAttribute(T &attrib)
  {
    // Initialize the attribute, with our layout.  It will be an
    // error if the layout type that the attrib expects does not
    // match the attribute layout type that we own and give to it.

    attrib.initialize(attributeLayout_m);

    // Add this attribute to our list.  AttributeList will create
    // an AttributeWrapper object to wrap around the special type T.

    return attributes_m.add(attrib);
  }

  /// Remove an Attribute from our list by index number.
  /// Return success.

  bool removeAttribute(AttributeList::Size_t index)
  {
    return attributes_m.remove(index);
  }
  //@}

  //============================================================
  /** @name Dynamic particle interface (sync, swap, create, destroy, etc.) */
  //============================================================
  //@{

  /// Sync up this object, by letting all contexts know what the
  /// other contexts have created and destroyed, plus move particles
  /// around to their proper place.  The no-argument version of sync()
  /// should be used for particle layouts that do not need to know the
  /// position or other information to do their job.
  /// Sync will also apply all cached delete operations and apply boundary
  /// conditions.  The order of operations is:
  ///  - apply boundary conditions
  ///  - perform cached destroys
  ///  - swap particles to proper patches

  void sync() { particleLayout_m.sync(*this); }

  /// Sync up this object, but also specify an object that can be used
  /// as an attribute for the particle layout's distribution algorithm.
  /// For example, this could be a position attribute when using a
  /// spatial layout.  The object should have the same interface and
  /// layout as the other particle attributes.  If more than one 
  /// attribute is needed (e.g., separate components of the position),
  /// they must be packaged up as a single attribute using, for example,
  /// the global function "vector".
  /// Sync will also apply all cached delete operations and apply boundary
  /// conditions.

  template <class T>
  void sync(const T& attrib) { particleLayout_m.sync(*this, attrib); }

  /// Just perform the swap operation, which invokes the corresponding
  /// routine in the particle layout object to move particles around
  /// to other patches.  No other operation is done (destroys, BC's, etc).
  /// To do ALL the update operations, call sync() (which in turn will call
  /// swap(), among other things).  We have versions of swap() that takes
  /// the same types of arguments as sync().

  void swap() { particleLayout_m.swap(*this); }

  template <class T>
  void swap(const T& attrib) { particleLayout_m.swap(*this, attrib); }

  /// Recompute the global domain.  No particle swapping or application
  /// of boundary conditions is done here.

  void renumber()
  {
    // Just tell the Attribute layout to resync ... this does not
    // adjust where particles are located, but it will make the
    // attributes the proper size so that storage can be initialized.

    attributeLayout_m.sync();
  }

  /// Create the specified number of particles on the local context,
  /// and add them to the last local patch, or the specified patch.
  /// The domain is renumbered after the create operation if the last
  /// argument is true.  Note that if you are renumbering, you really
  /// must call this function in SPMD style on every context.

  void create(Size_t np, PatchID_t patch = (-1), bool renum = true)
  {
    if (np > 0)
    {
      if (patch >= 0)
        attributeLayout_m.create(np, patch);
      else
        attributeLayout_m.create(np);
    }

    if (renum)
      renumber();
  }

  /// Create the specified number of particles, total, spread across
  /// all patches.  An equal (or close to) number of particles are added
  /// to each patch on each context.  This must be called in SPMD-style.
  /// The domain is renumbered after the create operation if the second
  /// argument is true.

  void globalCreate(Size_t np, bool renum = true);

  /// Destroy set of particles specified by DomainType, which may be a
  /// 1D Range of particle index numbers or a list of index numbers.
  /// This version does the destroy immediately, instead of deferring
  /// it until later.  If the patchID is < 0, the domain should be
  /// a global domain.  If the patchID is >= 0, the domain should be
  /// a "local" one, where index = 0 refers to the first particle in
  /// that patch.  If the final argument is true (the default), the
  /// domain is renumbered after the destroy is complete.

  template <class DomainType>
  void destroy(const DomainType& dom, PatchID_t pid = (-1), bool renum = true);

  /// Same as above, but use begin/end iterator pair to specify domain.

  template <class Iter>
  void destroy(Iter begin, Iter end, PatchID_t pid = (-1), bool renum = true);

  /// Destroy the particles in the given global domain by storing them
  /// a list for use in a later performDestroy call.  This lets you
  /// perform several different destroy operations without going through
  /// all the shifting-around work until destroys are done.  This does
  /// NOT renumber the domain.  If the patchID is < 0, the domain should
  /// be a global domain. If the patchID is >= 0, the domain should
  /// be a "local" one, where index = 0 refers to the first particle in
  /// that patch.

  template <class DomainType>
  void deferredDestroy(const DomainType& domain, PatchID_t pid = (-1));

  /// Same as above, but use begin/end iterator pair to specify domain.

  template <class Iter>
  void deferredDestroy(Iter begin, Iter end, PatchID_t pid = (-1));

  /// Carry out all the previously requested destroy commands, and
  /// renumber the domain if requested.  If the requested patch ID is
  /// < 0, do all deferred destroys, otherwise do it for just the requested
  /// patch.

  void performDestroy(PatchID_t pid = (-1), bool renum = true);

  //@}

  //============================================================
  /** @name Boundary Condition methods */
  //============================================================
  //@{

  /// Add a new boundary condition of given type, acting on info from the
  /// Subject and affecting values in the Object.

  template <class Subject, class Object, class BCType>
  void addBoundaryCondition(const Subject& s, const Object& o,
			    const BCType& bc);

  /// Add a new boundary condition of given type, acting on info from the
  /// Subject and affecting values in the same Subject.

  template <class Subject, class BCType>
  void addBoundaryCondition(const Subject& s, const BCType& bc);

  /// Remove the ith boundary condition in this Particles object.

  void removeBoundaryCondition(ParticleBCList::Size_t i);

  /// Remove all the boundary conditions

  void removeBoundaryConditions();

  /// Apply the boundary conditions to the existing particles for all
  /// patches, or just one (if all, the value for the patch ID is < 0).

  void applyBoundaryConditions(PatchID_t pid = (-1));

  /// Return number of boundary conditions

  ParticleBCList::Size_t boundaryConditions() const { return bcList_m.size(); }

  /// Return the ith boundary condition in the list

  ParticleBCItem* boundaryCondition(ParticleBCList::Size_t i)
  {
    return bcList_m(i);
  }

  const ParticleBCItem* boundaryCondition(ParticleBCList::Size_t i) const
  {
    return bcList_m(i);
  }

  //@}

  //============================================================
  // I/O
  //============================================================

  template <class Out>
  void print(Out& o) const;

private:

  /// A list of pointers to our Attributes, so we can loop over them.

  AttributeList attributes_m;

  /// Our particle layout object, the item that embodies the algorithm
  /// used to determine what particle goes where.

  ParticleLayout_t particleLayout_m;

  /// Our attribute layout object, the item that manages the information about
  /// what particle is located where.

  AttributeLayout_t attributeLayout_m;

  /// Our current destroy method code.

  int destroyMethod_m;

  /// A list of deferred particle destroy requests.  The list is processed
  /// whenever the sync() function is called.  Stored as a set of Dynamic
  /// Arrays (one per local patch) containing local particle indices.

  DynamicArray<int,Dynamic> *destroyList_m;

  /// A list of particle boundary conditions to be applied.  The list is
  /// processed whenever the sync() function is called.

  ParticleBCList bcList_m;

};


//-----------------------------------------------------------------------------
//
/// A specialization of the Inform traits used to say that Particles has
/// a print method.
//
//-----------------------------------------------------------------------------

template <class PT>
std::ostream &operator<<(std::ostream &o, const Particles<PT> &particles)
{
  particles.print(o);
  return o;
}


// Include out-of-line definitions

#include "Particles/Particles.cpp"


// } // namespace POOMA


//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_PARTICLES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Particles.h,v $   $Author: richard $
// $Revision: 1.36 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
