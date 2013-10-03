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

#ifndef POOMA_PARTICLES_INTERPOLATORNGP_H
#define POOMA_PARTICLES_INTERPOLATORNGP_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * Specialization of Interpolator class template for nearest-grid-point
 * (NGP) interpolation between Particle Attributes and Fields.
 *
 * Specialization of Interpolator class template for 
 * nearest-grid-point interpolation.  Interpolation is performed
 * using a PatchFunction that spawns threads to work on each patch
 * and loop over the particles on that patch.  This functor will
 * store a copy of the Field to be gathered from or scattered to
 * and take the appropriate view of the Field for each patch.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Particles/Interpolation.h"
#include "Domain/Loc.h"
#include "Utilities/PAssert.h"
#include "Utilities/ElementProperties.h"

#include <iostream>


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

// NGP Functor classes for gather and scatter

template <class FC, int Dim> struct NGPGather;
template <class FC, int Dim> struct NGPScatter;
template <class FC, int Dim, class ValueT> struct NGPScatterValue;

template <class FC, int Dim> struct NGPGatherFillCache;
template <class FC, int Dim> struct NGPScatterFillCache;
template <class FC, int Dim, class ValueT> struct NGPScatterValueFillCache;

template <class FC, int Dim> struct NGPGatherUseCache;
template <class FC, int Dim> struct NGPScatterUseCache;
template <class FC, int Dim, class ValueT> struct NGPScatterValueUseCache;


//-----------------------------------------------------------------------------
// Interpolator tag class for NGP interpolation.
//-----------------------------------------------------------------------------

struct NGP
{
};


//-----------------------------------------------------------------------------
// Struct to store cached interpolation data
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct NGPCacheData
{
  // static variables and typedefs
  static const int dimensions = Dim;
  typedef T AxisType_t;

  // data members
  Loc<Dim> ngp;

  // send to stream
  template <class Out>
  void print(Out& o) const
  {
    o << ngp;
  }
};


//-----------------------------------------------------------------------------
// Global function for sending Interpolator cache data to an ostream
//-----------------------------------------------------------------------------

template <int Dim, class T>
std::ostream&
operator<<(std::ostream& o, const NGPCacheData<Dim,T>& cache)
{
  cache.print(o);
  return o;
}

//-----------------------------------------------------------------------------
// Specialization of ElementProperties struct for NGPCacheData.
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct ElementProperties< NGPCacheData<Dim,T> > 
  : public TrivialElementProperties< NGPCacheData<Dim,T> >
{ };

//-----------------------------------------------------------------------------
// Interpolator class template specialized for NGP interpolation.
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct Interpolator<Dim,T,NGP>
{
  // typedef for the CacheData struct

  typedef NGPCacheData<Dim,T> Cache_t;

  // gather/scatter using particle position attribute

  template <class PA, class FC, class PPos>
  static
  void gather(const PA& attrib, const FC& field, const PPos& pos)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == attrib.layout().sizeGlobal(),
        "Field and Particle Attribute must have same number of patches!");
      PInsist(layout.sizeGlobal() == pos.layout().sizeGlobal(),
        "Field and Particle Position must have same number of patches!");

      // Create functor and give it the field to store
      NGPGather<FC,Dim> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< NGPGather<FC,Dim>,
                     PatchParticle2<true,false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute
      patchfun.block(attrib,pos);
    }

  template <class PA, class FC, class PPos>
  static
  void scatter(const PA& attrib, const FC& field, const PPos& pos)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == attrib.layout().sizeGlobal(),
        "Field and Particle Attribute must have same number of patches!");
      PInsist(layout.sizeGlobal() == pos.layout().sizeGlobal(),
        "Field and Particle Position must have same number of patches!");

      // Create functor and give it the field to store
      NGPScatter<FC,Dim> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< NGPScatter<FC,Dim>,
                     PatchParticle2<false,false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute
      patchfun.block(attrib,pos);
    }

  template <class ValueT, class FC, class PPos>
  static
  void scatterValue(const ValueT& value, const FC& field,
                    const PPos& pos)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == pos.layout().sizeGlobal(),
        "Field and Particle Position must have same number of patches!");

      // Create functor and give it the field and value to store
      NGPScatterValue<FC,Dim,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< NGPScatterValue<FC,Dim,ValueT>,
                     PatchParticle1<false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute
      patchfun.block(pos);
    }

  // gather/scatter using particle position attribute and cache
  // interpolation data

  template <class PA, class FC, class PPos, class EngineTag>
  static
  void gatherCache(const PA& attrib, const FC& field, const PPos& pos,
                   const DynamicArray<Cache_t,EngineTag>& cache)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == attrib.layout().sizeGlobal(),
        "Field and Particle Attribute must have same number of patches!");
      PInsist(layout.sizeGlobal() == pos.layout().sizeGlobal(),
        "Field and Particle Position must have same number of patches!");
      PInsist(layout.sizeGlobal() == cache.layout().sizeGlobal(),
        "Field and Particle CacheData must have same number of patches!");

      // Create functor and give it the field to store
      NGPGatherFillCache<FC,Dim> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< NGPGatherFillCache<FC,Dim>,
                     PatchParticle3<true,false,true> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute and caching the interpolation data
      patchfun.block(attrib,pos,cache);
    }

  template <class PA, class FC, class PPos, class EngineTag>
  static
  void scatterCache(const PA& attrib, const FC& field, const PPos& pos,
                    const DynamicArray<Cache_t,EngineTag>& cache)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == attrib.layout().sizeGlobal(),
        "Field and Particle Attribute must have same number of patches!");
      PInsist(layout.sizeGlobal() == pos.layout().sizeGlobal(),
        "Field and Particle Position must have same number of patches!");
      PInsist(layout.sizeGlobal() == cache.layout().sizeGlobal(),
        "Field and Particle CacheData must have same number of patches!");

      // Create functor and give it the field to store
      NGPScatterFillCache<FC,Dim> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< NGPScatterFillCache<FC,Dim>,
                     PatchParticle3<false,false,true> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute and caching the interpolation data
      patchfun.block(attrib,pos,cache);
    }

  template <class ValueT, class FC, class PPos, class EngineTag>
  static
  void scatterValueCache(const ValueT& value, const FC& field, const PPos& pos,
                         const DynamicArray<Cache_t,EngineTag>& cache)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == pos.layout().sizeGlobal(),
        "Field and Particle Position must have same number of patches!");
      PInsist(layout.sizeGlobal() == cache.layout().sizeGlobal(),
        "Field and Particle CacheData must have same number of patches!");

      // Create functor and give it the field and value to store
      NGPScatterValueFillCache<FC,Dim,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< NGPScatterValueFillCache<FC,Dim,ValueT>,
                     PatchParticle2<false,true> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute and caching the interpolation data
      patchfun.block(pos,cache);
    }

  // gather/scatter using cached interpolation data

  template <class PA, class FC, class EngineTag>
  static
  void gatherCache(const PA& attrib, const FC& field,
                   const DynamicArray<Cache_t,EngineTag>& cache)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == attrib.layout().sizeGlobal(),
        "Field and Particle Attribute must have same number of patches!");
      PInsist(layout.sizeGlobal() == cache.layout().sizeGlobal(),
        "Field and Particle CacheData must have same number of patches!");

      // Create functor and give it the field to store
      NGPGatherUseCache<FC,Dim> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< NGPGatherUseCache<FC,Dim>,
                     PatchParticle2<true,false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // cached interpolation data
      patchfun.block(attrib,cache);
    }

  template <class PA, class FC, class EngineTag>
  static
  void scatterCache(const PA& attrib, const FC& field,
                    const DynamicArray<Cache_t,EngineTag>& cache)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == attrib.layout().sizeGlobal(),
        "Field and Particle Attribute must have same number of patches!");
      PInsist(layout.sizeGlobal() == cache.layout().sizeGlobal(),
        "Field and Particle CacheData must have same number of patches!");

      // Create functor and give it the field to store
      NGPScatterUseCache<FC,Dim> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< NGPScatterUseCache<FC,Dim>,
                     PatchParticle2<false,false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // cached interpolation data
      patchfun.block(attrib,cache);
    }

  template <class ValueT, class FC, class EngineTag>
  static
  void scatterValueCache(const ValueT& value, const FC& field,
                         const DynamicArray<Cache_t,EngineTag>& cache)
    {
      // Check that arguments have the same number of patches
      typedef typename FC::Engine_t Engine_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field.engine().layout();
      PInsist(layout.sizeGlobal() == cache.layout().sizeGlobal(),
        "Field and Particle CacheData must have same number of patches!");

      // Create functor and give it the field and value to store
      NGPScatterValueUseCache<FC,Dim,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< NGPScatterValueUseCache<FC,Dim,ValueT>,
                     PatchParticle1<false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // cached interpolation data
      patchfun.block(cache);
    }
};


//-----------------------------------------------------------------------------
// Definitions for NGPGather, NGPScatter and NGPScatterValue functors
//-----------------------------------------------------------------------------

template <class FC, int Dim>
struct NGPGather
{
  // Typedefs

  typedef NGPGather<FC,Dim> This_t;
  typedef FC                Field_t;
  typedef int               PatchID_t;
  typedef int               Size_t;

  // Constructor

  NGPGather(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  NGPGather(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~NGPGather()
  {
  }

  // apply method implements NGP gather

  template <class Patch1, class Patch2>
  void apply(const Patch1& attrib, const Patch2& pos, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = attrib.domain().size();
      if (n == 0)
        return;

      // Get the global patch ID for this patch

      typedef typename Field_t::Engine_t Engine_t;
      typedef typename Field_t::Mesh_t Mesh_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field_m.engine().layout();
      PatchID_t gid = layout.nodeListLocal()[pid]->globalID();

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and gather
      // field values into the particle attribute.

      Size_t i;
      Loc<Dim> indx;
      const Mesh_t& mesh = field_m.mesh();
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now get the value at this position in the field
          // and put it into the particle attribute

          attrib(i) = fpatch(indx);
        }
    }

  // Copy of Field to be gathered from.

  Field_t field_m;
};

template <class FC, int Dim>
struct NGPScatter
{
  // Typedefs

  typedef NGPScatter<FC,Dim> This_t;
  typedef FC                 Field_t;
  typedef int                PatchID_t;
  typedef int                Size_t;

  // Constructor

  NGPScatter(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  NGPScatter(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~NGPScatter()
  {
  }

  // apply method implements NGP scatter

  template <class Patch1, class Patch2>
  void apply(const Patch1& attrib, const Patch2& pos, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = attrib.domain().size();
      if (n == 0)
        return;

      // Get the global patch ID for this patch

      typedef typename Field_t::Engine_t Engine_t;
      typedef typename Field_t::Mesh_t Mesh_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field_m.engine().layout();
      PatchID_t gid = layout.nodeListLocal()[pid]->globalID();

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and scatter
      // particle attribute into the field elements.

      Size_t i;
      Loc<Dim> indx;
      const Mesh_t& mesh = field_m.mesh();
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now add the particle attribute value into the field
          
          fpatch(indx) += attrib(i);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class ValueT>
struct NGPScatterValue
{
  // Typedefs

  typedef NGPScatterValue<FC,Dim,ValueT> This_t;
  typedef FC                             Field_t;
  typedef int                            PatchID_t;
  typedef int                            Size_t;
  typedef ValueT                         Value_t;

  // Constructor

  NGPScatterValue(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  NGPScatterValue(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~NGPScatterValue()
  {
  }

  // apply method implements NGP scatter

  template <class Patch1>
  void apply(const Patch1& pos, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = pos.domain().size();
      if (n == 0)
        return;

      // Get the global patch ID for this patch

      typedef typename Field_t::Engine_t Engine_t;
      typedef typename Field_t::Mesh_t Mesh_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field_m.engine().layout();
      PatchID_t gid = layout.nodeListLocal()[pid]->globalID();

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and scatter
      // value into the field elements.

      Size_t i;
      Loc<Dim> indx;
      const Mesh_t& mesh = field_m.mesh();
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now add the value into the field
          
          fpatch(indx) += value_m;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};

//-----------------------------------------------------------------------------
// Definitions for NGPGatherFillCache, NGPScatterFillCache and
// NGPScatterValueFillCache functors
//-----------------------------------------------------------------------------

template <class FC, int Dim>
struct NGPGatherFillCache
{
  // Typedefs

  typedef NGPGatherFillCache<FC,Dim> This_t;
  typedef FC                         Field_t;
  typedef int                        PatchID_t;
  typedef int                        Size_t;

  // Constructor

  NGPGatherFillCache(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  NGPGatherFillCache(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~NGPGatherFillCache()
  {
  }

  // apply method implements NGP gather and fills cache

  template <class Patch1, class Patch2, class Patch3>
  void apply(const Patch1& attrib, const Patch2& pos,
             const Patch3& cache, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = attrib.domain().size();
      if (n == 0)
        return;

      // Get the global patch ID for this patch.

      typedef typename Field_t::Engine_t Engine_t;
      typedef typename Field_t::Mesh_t Mesh_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field_m.engine().layout();
      PatchID_t gid = layout.nodeListLocal()[pid]->globalID();

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and gather
      // field values into the particle attribute.

      Size_t i;
      Loc<Dim> indx;
      const Mesh_t& mesh = field_m.mesh();
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));

          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now get the value at this position in the field
          // and put it into the particle attribute

          attrib(i) = fpatch(indx);
          
          // Cache this interpolation data

          cache(i).ngp = indx;
        }
    }

  // Copy of Field to be gathered from

  Field_t field_m;
};

template <class FC, int Dim>
struct NGPScatterFillCache
{
  // Typedefs

  typedef NGPScatterFillCache<FC,Dim> This_t;
  typedef FC                          Field_t;
  typedef int                         PatchID_t;
  typedef int                         Size_t;

  // Constructor

  NGPScatterFillCache(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  NGPScatterFillCache(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~NGPScatterFillCache()
  {
  }

  // apply method implements NGP scatter and fills cache

  template <class Patch1, class Patch2, class Patch3>
  void apply(const Patch1& attrib, const Patch2& pos,
             const Patch3& cache, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = attrib.domain().size();
      if (n == 0)
        return;

      // Get the global patch ID for this patch.

      typedef typename Field_t::Engine_t Engine_t;
      typedef typename Field_t::Mesh_t Mesh_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field_m.engine().layout();
      PatchID_t gid = layout.nodeListLocal()[pid]->globalID();

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and scatter
      // particle attribute into the field elements.

      Size_t i;
      Loc<Dim> indx;
      const Mesh_t& mesh = field_m.mesh();
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now add the particle attribute value into the field
          
          fpatch(indx) += attrib(i);

          // Cache this interpolation data.

          cache(i).ngp = indx;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class ValueT>
struct NGPScatterValueFillCache
{
  // Typedefs

  typedef NGPScatterValueFillCache<FC,Dim,ValueT> This_t;
  typedef FC                                      Field_t;
  typedef int                                     PatchID_t;
  typedef int                                     Size_t;
  typedef ValueT                                  Value_t;

  // Constructor

  NGPScatterValueFillCache(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  NGPScatterValueFillCache(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~NGPScatterValueFillCache()
  {
  }

  // apply method implements NGP scatter and fills cache

  template <class Patch1, class Patch2>
  void apply(const Patch1& pos, const Patch2& cache, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = cache.domain().size();
      if (n == 0)
        return;

      // Get the global patch ID for this patch.

      typedef typename Field_t::Engine_t Engine_t;
      typedef typename Field_t::Mesh_t Mesh_t;
      typedef typename Engine_t::Layout_t Layout_t;
      const Layout_t& layout = field_m.engine().layout();
      PatchID_t gid = layout.nodeListLocal()[pid]->globalID();

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and scatter
      // value into the field elements.

      Size_t i;
      Loc<Dim> indx;
      const Mesh_t& mesh = field_m.mesh();
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now add the value into the field
          
          fpatch(indx) += value_m;

          // Cache this interpolation data.

          cache(i).ngp = indx;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};


//-----------------------------------------------------------------------------
// Definitions for NGPGatherUseCache, NGPScatterUseCache and
// NGPScatterValueUseCache functors
//-----------------------------------------------------------------------------

template <class FC, int Dim>
struct NGPGatherUseCache
{
  // Typedefs

  typedef NGPGatherUseCache<FC,Dim> This_t;
  typedef FC                        Field_t;
  typedef int                       PatchID_t;
  typedef int                       Size_t;

  // Constructor

  NGPGatherUseCache(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  NGPGatherUseCache(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~NGPGatherUseCache()
  {
  }

  // apply method implements NGP gather using cached data

  template <class Patch1, class Patch2>
  void apply(const Patch1& attrib, const Patch2& cache, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = attrib.domain().size();
      if (n == 0)
        return;

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and gather
      // field values into the particle attribute.

      Size_t i;
      for (i=0; i<n; ++i)
        {
          // get the value at this position in the field
          // using cached data and put it into the particle attribute

          attrib(i) = fpatch(cache(i).ngp);
        }
    }

  // Copy of Field to be gathered from

  Field_t field_m;
};

template <class FC, int Dim>
struct NGPScatterUseCache
{
  // Typedefs

  typedef NGPScatterUseCache<FC,Dim> This_t;
  typedef FC                         Field_t;
  typedef int                        PatchID_t;
  typedef int                        Size_t;

  // Constructor

  NGPScatterUseCache(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  NGPScatterUseCache(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~NGPScatterUseCache()
  {
  }

  // apply method implements NGP scatter using cached data

  template <class Patch1, class Patch2>
  void apply(const Patch1& attrib, const Patch2& cache, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = attrib.domain().size();
      if (n == 0)
        return;

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and scatter
      // particle attribute into the field elements.

      Size_t i;
      for (i=0; i<n; ++i)
        {
          // add the particle attribute value into the field using cached data
          
          fpatch(cache(i).ngp) += attrib(i);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class ValueT>
struct NGPScatterValueUseCache
{
  // Typedefs

  typedef NGPScatterValueUseCache<FC,Dim,ValueT> This_t;
  typedef FC                                     Field_t;
  typedef int                                    PatchID_t;
  typedef int                                    Size_t;
  typedef ValueT                                 Value_t;

  // Constructor

  NGPScatterValueUseCache(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  NGPScatterValueUseCache(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~NGPScatterValueUseCache()
  {
  }

  // apply method implements NGP scatter using cached data

  template <class Patch1>
  void apply(const Patch1& cache, PatchID_t pid) const
    {
      // Get the number of particles in this patch.
      // If there are no particles here, just return.

      Size_t n = cache.domain().size();
      if (n == 0)
        return;

      // Get this patch within the field we store.

      typedef typename Patch<Field_t>::Type_t Patch_t;
      Patch_t fpatch = field_m.patchLocal(pid);

      // Now loop through the particles on this patch and scatter
      // value into the field elements.

      Size_t i;
      for (i=0; i<n; ++i)
        {
          // add the value into the field using cached data
          
          fpatch(cache(i).ngp) += value_m;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTICLES_INTERPOLATORNGP_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: InterpolatorNGP.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
