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

#ifndef POOMA_PARTICLES_INTERPOLATORSUDS_H
#define POOMA_PARTICLES_INTERPOLATORSUDS_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * Specialization of Interpolator class template for subtracted dipole
 * scheme (SUDS) interpolation between Particle Attributes and Fields.
 * SUDS is a linear interpolation method that has an advantage over CIC
 * in three dimensions because it uses only 7 stencil points instead of 8.
 *
 * Specialization of Interpolator class template for 
 * cloud-in-cell interpolation.  Interpolation is performed
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
#include "Tiny/Vector.h"
#include "Utilities/PAssert.h"
#include "Utilities/ElementProperties.h"

#include <iostream>


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

// SUDS Functor classes for gather and scatter

template <class FC, int Dim, class T> struct SUDSGather;
template <class FC, int Dim, class T> struct SUDSScatter;
template <class FC, int Dim, class T, class ValueT>
struct SUDSScatterValue;

template <class FC, int Dim, class T> struct SUDSGatherFillCache;
template <class FC, int Dim, class T> struct SUDSScatterFillCache;
template <class FC, int Dim, class T, class ValueT>
struct SUDSScatterValueFillCache;

template <class FC, int Dim, class T> struct SUDSGatherUseCache;
template <class FC, int Dim, class T> struct SUDSScatterUseCache;
template <class FC, int Dim, class T, class ValueT>
struct SUDSScatterValueUseCache;


//-----------------------------------------------------------------------------
// Interpolator tag class for SUDS interpolation.
//-----------------------------------------------------------------------------

struct SUDS
{
};


//-----------------------------------------------------------------------------
// Struct to store cached interpolation data
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct SUDSCacheData
{
  // static variables and typedefs
  static const int dimensions = Dim;
  typedef T AxisType_t;

  Loc<Dim> ngp;
  Vector<Dim,T> dist;

  template <class Out>
  void print(Out& o) const
  {
    o << ngp;
    o << dist;
  }
};


//-----------------------------------------------------------------------------
// Global function for sending Interpolator cache data to an ostream
//-----------------------------------------------------------------------------

template <int Dim, class T>
std::ostream&
operator<<(std::ostream& o, const SUDSCacheData<Dim,T>& cache)
{
  cache.print(o);
  return o;
}

//-----------------------------------------------------------------------------
// Specialization of ElementProperties struct for SUDSCacheData.
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct ElementProperties< SUDSCacheData<Dim,T> > 
  : public TrivialElementProperties< SUDSCacheData<Dim,T> >
{ };

//-----------------------------------------------------------------------------
// Interpolator class template specialized for SUDS interpolation.
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct Interpolator<Dim,T,SUDS>
{
  // typedef for the CacheData struct

  typedef SUDSCacheData<Dim,T> Cache_t;

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

      // Check that the Field has adequate guard layers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }

      // Make sure the guard layers have been updated
      field.engine().fillGuards();

      // Create functor and give it the field to store
      SUDSGather<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSGather<FC,Dim,T>,
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

      // Check that the Field has adequate guard layers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }

      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field to store
      SUDSScatter<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSScatter<FC,Dim,T>,
                     PatchParticle2<false,false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute
      patchfun.block(attrib,pos);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field and value to store
      SUDSScatterValue<FC,Dim,T,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSScatterValue<FC,Dim,T,ValueT>,
                     PatchParticle1<false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute
      patchfun.block(pos);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Make sure the guard layers have been updated
      field.engine().fillGuards();

      // Create functor and give it the field to store
      SUDSGatherFillCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSGatherFillCache<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field to store
      SUDSScatterFillCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSScatterFillCache<FC,Dim,T>,
                     PatchParticle3<false,false,true> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute and caching the interpolation data
      patchfun.block(attrib,pos,cache);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field and value to store
      SUDSScatterValueFillCache<FC,Dim,T,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSScatterValueFillCache<FC,Dim,T,ValueT>,
                     PatchParticle2<false,true> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // particle position attribute and caching the interpolation data
      patchfun.block(pos,cache);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Make sure the guard layers have been updated
      field.engine().fillGuards();

      // Create functor and give it the field to store
      SUDSGatherUseCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSGatherUseCache<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field to store
      SUDSScatterUseCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSScatterUseCache<FC,Dim,T>,
                     PatchParticle2<false,false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // cached interpolation data
      patchfun.block(attrib,cache);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
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

      // Check that the Field has adequate GuardLayers for SUDS
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for SUDS!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field and value to store
      SUDSScatterValueUseCache<FC,Dim,T,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< SUDSScatterValueUseCache<FC,Dim,T,ValueT>,
                     PatchParticle1<false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // cached interpolation data
      patchfun.block(cache);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
    }
};


//-----------------------------------------------------------------------------
// Definitions for dimension-specific SUDSGatherFcn and
// SUDSScatterFcn functions
//-----------------------------------------------------------------------------

template <class T, class Patch, class AxisType>
inline
void SUDSGatherFcn(T& attrib, const Patch& field, const Loc<1>& index,
                   const Vector<1,AxisType>& delta)
{
  attrib = field(index) +
           0.5 * delta(0) *
           ( field(index + Loc<1>(1)) - field(index - Loc<1>(1)) );
}

template <class T, class Patch, class AxisType>
inline
void SUDSGatherFcn(T& attrib, const Patch& field, const Loc<2>& index,
                   const Vector<2,AxisType>& delta)
{
  attrib = field(index) +
           0.5 * delta(0) *
           ( field(index + Loc<2>(1,0)) - field(index - Loc<2>(1,0)) ) +
           0.5 * delta(1) *
           ( field(index + Loc<2>(0,1)) - field(index - Loc<2>(0,1)) );
}

template <class T, class Patch, class AxisType>
inline
void SUDSGatherFcn(T& attrib, const Patch& field, const Loc<3>& index,
                   const Vector<3,AxisType>& delta)
{
  attrib = field(index) +
           0.5 * delta(0) *
           ( field(index + Loc<3>(1,0,0)) - field(index - Loc<3>(1,0,0)) ) +
           0.5 * delta(1) *
           ( field(index + Loc<3>(0,1,0)) - field(index - Loc<3>(0,1,0)) ) +
           0.5 * delta(2) *
           ( field(index + Loc<3>(0,0,1)) - field(index - Loc<3>(0,0,1)) );
}

template <class T, class Patch, class AxisType>
inline
void SUDSScatterFcn(const T& value, const Patch& field,
                    const Loc<1>& index, const Vector<1,AxisType>& delta)
{
  field(index)             +=                  value;
  field(index + Loc<1>(1)) += 0.5 * delta(0) * value;
  field(index - Loc<1>(1)) -= 0.5 * delta(0) * value;
}

template <class T, class Patch, class AxisType>
inline
void SUDSScatterFcn(const T& value, const Patch& field,
                    const Loc<2>& index, const Vector<2,AxisType>& delta)
{
  field(index)               +=                  value;
  field(index + Loc<2>(1,0)) += 0.5 * delta(0) * value;
  field(index - Loc<2>(1,0)) -= 0.5 * delta(0) * value;
  field(index + Loc<2>(0,1)) += 0.5 * delta(1) * value;
  field(index - Loc<2>(0,1)) -= 0.5 * delta(1) * value;
}

template <class T, class Patch, class AxisType>
inline
void SUDSScatterFcn(const T& value, const Patch& field,
                    const Loc<3>& index, const Vector<3,AxisType>& delta)
{
  field(index)                 +=                  value;
  field(index + Loc<3>(1,0,0)) += 0.5 * delta(0) * value;
  field(index - Loc<3>(1,0,0)) -= 0.5 * delta(0) * value;
  field(index + Loc<3>(0,1,0)) += 0.5 * delta(1) * value;
  field(index - Loc<3>(0,1,0)) -= 0.5 * delta(1) * value;
  field(index + Loc<3>(0,0,1)) += 0.5 * delta(2) * value;
  field(index - Loc<3>(0,0,1)) -= 0.5 * delta(2) * value;
}


//-----------------------------------------------------------------------------
// Definitions for SUDSGather, SUDSScatter and SUDSScatterValue functors
//-----------------------------------------------------------------------------

template <class FC, int Dim, class T>
struct SUDSGather
{
  // Typedefs

  typedef SUDSGather<FC,Dim,T> This_t;
  typedef FC                   Field_t;
  typedef int                  PatchID_t;
  typedef int                  Size_t;

  // Constructor

  SUDSGather(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  SUDSGather(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~SUDSGather()
  {
  }

  // apply method implements SUDS gather

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
      typedef typename Mesh_t::PointType_t PointType_t;
      PointType_t gpos, delta, lpos, upos;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);

          // now compute position at the NGP and the normalized
          // distance between the particle and NGP positions

          gpos = mesh.vertexPosition(indx);
          lpos = mesh.vertexPosition(indx-1);
          upos = mesh.vertexPosition(indx+1);
          for (int idim=0; idim<Dim; ++idim)
            {
              if (pos(i)(idim) > gpos(idim))
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (upos(idim) - gpos(idim));
              else
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (gpos(idim) - lpos(idim));
            }

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the gather

          SUDSGatherFcn(attrib(i),fpatch,indx,delta);
        }
    }

  // Copy of Field to be gathered from.

  Field_t field_m;
};

template <class FC, int Dim, class T>
struct SUDSScatter
{
  // Typedefs

  typedef SUDSScatter<FC,Dim,T> This_t;
  typedef FC                    Field_t;
  typedef int                   PatchID_t;
  typedef int                   Size_t;

  // Constructor

  SUDSScatter(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  SUDSScatter(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~SUDSScatter()
  {
  }

  // apply method implements SUDS scatter

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
      typedef typename Mesh_t::PointType_t PointType_t;
      PointType_t gpos, delta, lpos, upos;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);

          // now compute position at the NGP and the normalized
          // distance between the particle and NGP positions

          gpos = mesh.vertexPosition(indx);
          lpos = mesh.vertexPosition(indx-1);
          upos = mesh.vertexPosition(indx+1);
          for (int idim=0; idim<Dim; ++idim)
            {
              if (pos(i)(idim) > gpos(idim))
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (upos(idim) - gpos(idim));
              else
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (gpos(idim) - lpos(idim));
            }

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          SUDSScatterFcn(attrib(i),fpatch,indx,delta);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class T, class ValueT>
struct SUDSScatterValue
{
  // Typedefs

  typedef SUDSScatterValue<FC,Dim,T,ValueT> This_t;
  typedef FC                                Field_t;
  typedef int                               PatchID_t;
  typedef int                               Size_t;
  typedef ValueT                            Value_t;

  // Constructor

  SUDSScatterValue(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  SUDSScatterValue(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~SUDSScatterValue()
  {
  }

  // apply method implements SUDS scatter

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
      typedef typename Mesh_t::PointType_t PointType_t;
      PointType_t gpos, delta, lpos, upos;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now compute position at the NGP and the normalized
          // distance between the particle and NGP positions

          gpos = mesh.vertexPosition(indx);
          lpos = mesh.vertexPosition(indx-1);
          upos = mesh.vertexPosition(indx+1);
          for (int idim=0; idim<Dim; ++idim)
            {
              if (pos(i)(idim) > gpos(idim))
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (upos(idim) - gpos(idim));
              else
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (gpos(idim) - lpos(idim));
            }

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          SUDSScatterFcn(value_m,fpatch,indx,delta);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};

//-----------------------------------------------------------------------------
// Definitions for SUDSGatherFillCache, SUDSScatterFillCache and
// SUDSScatterValueFillCache functors
//-----------------------------------------------------------------------------

template <class FC, int Dim, class T>
struct SUDSGatherFillCache
{
  // Typedefs

  typedef SUDSGatherFillCache<FC,Dim,T> This_t;
  typedef FC                            Field_t;
  typedef int                           PatchID_t;
  typedef int                           Size_t;

  // Constructor

  SUDSGatherFillCache(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  SUDSGatherFillCache(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~SUDSGatherFillCache()
  {
  }

  // apply method implements SUDS gather and fills cache

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
      typedef typename Mesh_t::PointType_t PointType_t;
      PointType_t gpos, delta, lpos, upos;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));

          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now compute position at the NGP and the normalized
          // distance between the particle and NGP positions

          gpos = mesh.vertexPosition(indx);
          lpos = mesh.vertexPosition(indx-1);
          upos = mesh.vertexPosition(indx+1);
          for (int idim=0; idim<Dim; ++idim)
            {
              if (pos(i)(idim) > gpos(idim))
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (upos(idim) - gpos(idim));
              else
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (gpos(idim) - lpos(idim));
            }

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the gather

          SUDSGatherFcn(attrib(i),fpatch,indx,delta);

          // Cache this interpolation data

          cache(i).ngp = indx;
          cache(i).dist = delta;
        }
    }

  // Copy of Field to be gathered from

  Field_t field_m;
};

template <class FC, int Dim, class T>
struct SUDSScatterFillCache
{
  // Typedefs

  typedef SUDSScatterFillCache<FC,Dim,T> This_t;
  typedef FC                             Field_t;
  typedef int                            PatchID_t;
  typedef int                            Size_t;

  // Constructor

  SUDSScatterFillCache(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  SUDSScatterFillCache(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~SUDSScatterFillCache()
  {
  }

  // apply method implements SUDS scatter and fills cache

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
      typedef typename Mesh_t::PointType_t PointType_t;
      PointType_t gpos, delta, lpos, upos;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now compute position at the NGP and the normalized
          // distance between the particle and NGP positions

          gpos = mesh.vertexPosition(indx);
          lpos = mesh.vertexPosition(indx-1);
          upos = mesh.vertexPosition(indx+1);
          for (int idim=0; idim<Dim; ++idim)
            {
              if (pos(i)(idim) > gpos(idim))
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (upos(idim) - gpos(idim));
              else
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (gpos(idim) - lpos(idim));
            }

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          SUDSScatterFcn(attrib(i),fpatch,indx,delta);

          // Cache this interpolation data

          cache(i).ngp = indx;
          cache(i).dist = delta;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class T, class ValueT>
struct SUDSScatterValueFillCache
{
  // Typedefs

  typedef SUDSScatterValueFillCache<FC,Dim,T,ValueT> This_t;
  typedef FC                                         Field_t;
  typedef int                                        PatchID_t;
  typedef int                                        Size_t;
  typedef ValueT                                     Value_t;

  // Constructor

  SUDSScatterValueFillCache(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  SUDSScatterValueFillCache(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~SUDSScatterValueFillCache()
  {
  }

  // apply method implements SUDS scatter and fills cache

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
      typedef typename Mesh_t::PointType_t PointType_t;
      PointType_t gpos, delta, lpos, upos;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // now compute position at the NGP and the normalized
          // distance between the particle and NGP positions

          gpos = mesh.vertexPosition(indx);
          lpos = mesh.vertexPosition(indx-1);
          upos = mesh.vertexPosition(indx+1);
          for (int idim=0; idim<Dim; ++idim)
            {
              if (pos(i)(idim) > gpos(idim))
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (upos(idim) - gpos(idim));
              else
                delta(idim) = (pos(i)(idim) - gpos(idim)) /
                              (gpos(idim) - lpos(idim));
            }

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          SUDSScatterFcn(value_m,fpatch,indx,delta);

          // Cache this interpolation data

          cache(i).ngp = indx;
          cache(i).dist = delta;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};


//-----------------------------------------------------------------------------
// Definitions for SUDSGatherUseCache, SUDSScatterUseCache and
// SUDSScatterValueUseCache functors
//-----------------------------------------------------------------------------

template <class FC, int Dim, class T>
struct SUDSGatherUseCache
{
  // Typedefs

  typedef SUDSGatherUseCache<FC,Dim,T> This_t;
  typedef FC                           Field_t;
  typedef int                          PatchID_t;
  typedef int                          Size_t;

  // Constructor

  SUDSGatherUseCache(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  SUDSGatherUseCache(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~SUDSGatherUseCache()
  {
  }

  // apply method implements SUDS gather using cached data

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
          // pass the necessary information to a dimension-specific 
          // function that actually does the gather

          SUDSGatherFcn(attrib(i),fpatch,cache(i).ngp,cache(i).dist);
        }
    }

  // Copy of Field to be gathered from

  Field_t field_m;
};

template <class FC, int Dim, class T>
struct SUDSScatterUseCache
{
  // Typedefs

  typedef SUDSScatterUseCache<FC,Dim,T> This_t;
  typedef FC                            Field_t;
  typedef int                           PatchID_t;
  typedef int                           Size_t;

  // Constructor

  SUDSScatterUseCache(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  SUDSScatterUseCache(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~SUDSScatterUseCache()
  {
  }

  // apply method implements SUDS scatter using cached data

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
          // pass the necessary information to a dimension-specific 
          // function that actually does the scatter

          SUDSScatterFcn(attrib(i),fpatch,cache(i).ngp,cache(i).dist);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class T, class ValueT>
struct SUDSScatterValueUseCache
{
  // Typedefs

  typedef SUDSScatterValueUseCache<FC,Dim,T,ValueT> This_t;
  typedef FC                                        Field_t;
  typedef int                                       PatchID_t;
  typedef int                                       Size_t;
  typedef ValueT                                    Value_t;

  // Constructor

  SUDSScatterValueUseCache(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  SUDSScatterValueUseCache(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~SUDSScatterValueUseCache()
  {
  }

  // apply method implements SUDS scatter using cached data

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
          // pass the necessary information to a dimension-specific 
          // function that actually does the scatter

          SUDSScatterFcn(value_m,fpatch,cache(i).ngp,cache(i).dist);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTICLES_INTERPOLATORSUDS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: InterpolatorSUDS.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
