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

#ifndef POOMA_PARTICLES_INTERPOLATORCIC_H
#define POOMA_PARTICLES_INTERPOLATORCIC_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * Specialization of Interpolator class template for cloud-in-cell
 * (CIC) interpolation between Particle Attributes and Fields.
 * CIC is also known as linear interpolation or volume weighting.
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

// CIC Functor classes for gather and scatter

template <class FC, int Dim, class T> struct CICGather;
template <class FC, int Dim, class T> struct CICScatter;
template <class FC, int Dim, class T, class ValueT>
struct CICScatterValue;

template <class FC, int Dim, class T> struct CICGatherFillCache;
template <class FC, int Dim, class T> struct CICScatterFillCache;
template <class FC, int Dim, class T, class ValueT>
struct CICScatterValueFillCache;

template <class FC, int Dim, class T> struct CICGatherUseCache;
template <class FC, int Dim, class T> struct CICScatterUseCache;
template <class FC, int Dim, class T, class ValueT>
struct CICScatterValueUseCache;


//-----------------------------------------------------------------------------
// Interpolator tag class for CIC interpolation.
//-----------------------------------------------------------------------------

struct CIC
{
};


//-----------------------------------------------------------------------------
// Struct to store cached interpolation data
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct CICCacheData
{
  // static variables and typedefs
  static const int dimensions = Dim;
  typedef T AxisType_t;

  Loc<Dim> lgp;
  Vector<Dim,T> dist;

  template <class Out>
  void print(Out& o) const
  {
    o << lgp;
    o << dist;
  }
};


//-----------------------------------------------------------------------------
// Global function for sending Interpolator cache data to an ostream
//-----------------------------------------------------------------------------

template <int Dim, class T>
std::ostream&
operator<<(std::ostream& o, const CICCacheData<Dim,T>& cache)
{
  cache.print(o);
  return o;
}

//-----------------------------------------------------------------------------
// Specialization of ElementProperties struct for CICCacheData.
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct ElementProperties< CICCacheData<Dim,T> > 
  : public TrivialElementProperties< CICCacheData<Dim,T> >
{ };

//-----------------------------------------------------------------------------
// Interpolator class template specialized for CIC interpolation.
//-----------------------------------------------------------------------------

template <int Dim, class T>
struct Interpolator<Dim,T,CIC>
{
  // typedef for the CacheData struct

  typedef CICCacheData<Dim,T> Cache_t;

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

      // Check that the Field has adequate guard layers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }

      // Make sure the guard layers have been updated
      field.engine().fillGuards();

      // Create functor and give it the field to store
      CICGather<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< CICGather<FC,Dim,T>,
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

      // Check that the Field has adequate guard layers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }

      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field to store
      CICScatter<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< CICScatter<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field and value to store
      CICScatterValue<FC,Dim,T,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< CICScatterValue<FC,Dim,T,ValueT>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Make sure the guard layers have been updated
      field.engine().fillGuards();

      // Create functor and give it the field to store
      CICGatherFillCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< CICGatherFillCache<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field to store
      CICScatterFillCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< CICScatterFillCache<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field and value to store
      CICScatterValueFillCache<FC,Dim,T,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< CICScatterValueFillCache<FC,Dim,T,ValueT>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Make sure the guard layers have been updated
      field.engine().fillGuards();

      // Create functor and give it the field to store
      CICGatherUseCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< CICGatherUseCache<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field to store
      CICScatterUseCache<FC,Dim,T> intfun(field);

      // Create a PatchFunction using this functor
      PatchFunction< CICScatterUseCache<FC,Dim,T>,
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

      // Check that the Field has adequate GuardLayers for CIC
      const GuardLayers<Dim>& gl = field.layout().internalGuards();
      for (int d=0; d<Dim; ++d)
        {
          PInsist(gl.lower(d)>=1 && gl.upper(d)>=1,
                  "Minimum GuardLayer width of 1 required for CIC!");
        }
              
      // Zero out the guard layers before scattering
      typename FC::Element_t zero(0);
      field.engine().setGuards(zero);
      setExternalGuards(field,zero);

      // Make sure setExternalGuards has completed before scattering
      Pooma::blockAndEvaluate();

      // Create functor and give it the field and value to store
      CICScatterValueUseCache<FC,Dim,T,ValueT> intfun(field,value);

      // Create a PatchFunction using this functor
      PatchFunction< CICScatterValueUseCache<FC,Dim,T,ValueT>,
                     PatchParticle1<false> > patchfun(intfun);
      
      // Apply the PatchFunction to the attribute using the
      // cached interpolation data
      patchfun.block(cache);

      // Accumulate values from guard layers into real elements
      field.engine().accumulateFromGuards();
    }
};


//-----------------------------------------------------------------------------
// Definitions for dimension-specific CICGatherFcn and CICScatterFcn functions
//-----------------------------------------------------------------------------

template <class T, class Patch, class AxisType>
inline void CICGatherFcn(T& attrib, const Patch& field, const Loc<1>& index,
                         const Vector<1,AxisType>& delta)
{
  attrib = (1.0 - delta(0)) * field(index) +
           delta(0)         * field(index + Loc<1>(1));
}

template <class T, class Patch, class AxisType>
inline void CICGatherFcn(T& attrib, const Patch& field, const Loc<2>& index,
                         const Vector<2,AxisType>& delta)
{
  attrib = (1.0 - delta(0)) * (1.0 - delta(1)) * field(index) +
           delta(0)         * (1.0 - delta(1)) * field(index + Loc<2>(1,0)) +
           (1.0 - delta(0)) * delta(1)         * field(index + Loc<2>(0,1)) +
           delta(0)         * delta(1)         * field(index + Loc<2>(1,1));
}

template <class T, class Patch, class AxisType>
inline void CICGatherFcn(T& attrib, const Patch& field, const Loc<3>& index,
                         const Vector<3,AxisType>& delta)
{
  attrib = (1.0 - delta(0)) * (1.0 - delta(1)) * (1.0 - delta(2)) *
             field(index) +
           delta(0)         * (1.0 - delta(1)) * (1.0 - delta(2)) *
             field(index + Loc<3>(1,0,0)) +
           (1.0 - delta(0)) * delta(1)         * (1.0 - delta(2)) *
             field(index + Loc<3>(0,1,0)) +
           delta(0)         * delta(1)         * (1.0 - delta(2)) *
             field(index + Loc<3>(1,1,0)) +
           (1.0 - delta(0)) * (1.0 - delta(1)) * delta(2) *
             field(index + Loc<3>(0,0,1)) +
           delta(0)         * (1.0 - delta(1)) * delta(2) *
             field(index + Loc<3>(1,0,1)) +
           (1.0 - delta(0)) * delta(1)         * delta(2) *
             field(index + Loc<3>(0,1,1)) +
           delta(0)         * delta(1)         * delta(2) *
             field(index + Loc<3>(1,1,1));
}

template <class T, class Patch, class AxisType>
inline void CICScatterFcn(const T& value, const Patch& field,
                          const Loc<1>& index, const Vector<1,AxisType>& delta)
{
  field(index)             += (1.0 - delta(0)) * value;
  field(index + Loc<1>(1)) += delta(0)         * value;
}

template <class T, class Patch, class AxisType>
inline void CICScatterFcn(const T& value, const Patch& field,
                          const Loc<2>& index, const Vector<2,AxisType>& delta)
{
  field(index)               += (1.0 - delta(0)) * (1.0 - delta(1)) * value;
  field(index + Loc<2>(1,0)) += delta(0)         * (1.0 - delta(1)) * value;
  field(index + Loc<2>(0,1)) += (1.0 - delta(0)) * delta(1)         * value;
  field(index + Loc<2>(1,1)) += delta(0)         * delta(1)         * value;
}

template <class T, class Patch, class AxisType>
inline void CICScatterFcn(const T& value, const Patch& field,
                          const Loc<3>& index, const Vector<3,AxisType>& delta)
{
  field(index)                 +=
    (1.0 - delta(0)) * (1.0 - delta(1)) * (1.0 - delta(2)) * value;
  field(index + Loc<3>(1,0,0)) +=
    delta(0)         * (1.0 - delta(1)) * (1.0 - delta(2)) * value;
  field(index + Loc<3>(0,1,0)) +=
    (1.0 - delta(0)) * delta(1)         * (1.0 - delta(2)) * value;
  field(index + Loc<3>(1,1,0)) +=
    delta(0)         * delta(1)         * (1.0 - delta(2)) * value;
  field(index + Loc<3>(0,0,1)) +=
    (1.0 - delta(0)) * (1.0 - delta(1)) * delta(2)         * value;
  field(index + Loc<3>(1,0,1)) +=
    delta(0)         * (1.0 - delta(1)) * delta(2)         * value;
  field(index + Loc<3>(0,1,1)) +=
    (1.0 - delta(0)) * delta(1)         * delta(2)         * value;
  field(index + Loc<3>(1,1,1)) +=
    delta(0)         * delta(1)         * delta(2)         * value;
}


//-----------------------------------------------------------------------------
// Definitions for CICGather, CICScatter and CICScatterValue functors
//-----------------------------------------------------------------------------

template <class FC, int Dim, class T>
struct CICGather
{
  // Typedefs

  typedef CICGather<FC,Dim,T> This_t;
  typedef FC                  Field_t;
  typedef int                 PatchID_t;
  typedef int                 Size_t;

  // Constructor

  CICGather(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  CICGather(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~CICGather()
  {
  }

  // apply method implements CIC gather

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
      PointType_t gpos, delta;
      int idim;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);

          // This is the nearest-grid-point, so convert this to the
          // lower-grid-point (LGP) by comparing grid-point position
          // with particle position and adjusting as needed.

          gpos = mesh.vertexPosition(indx);
          for (idim=0; idim<Dim; ++idim)
            if (gpos(idim)>pos(i)(idim))
              indx[idim] = indx[idim] - 1;

          // now compute position and spacings at the LGP

          gpos = mesh.vertexPosition(indx);
          delta = mesh.vertexPosition(indx+1) - gpos;

          // From this, we find the normalized distance between 
          // the particle and LGP positions.

          for (idim=0; idim<Dim; ++idim)
            delta(idim) = (pos(i)(idim) - gpos(idim)) / delta(idim);

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the gather

          CICGatherFcn(attrib(i),fpatch,indx,delta);
        }
    }

  // Copy of Field to be gathered from.

  Field_t field_m;
};

template <class FC, int Dim, class T>
struct CICScatter
{
  // Typedefs

  typedef CICScatter<FC,Dim,T> This_t;
  typedef FC                   Field_t;
  typedef int                  PatchID_t;
  typedef int                  Size_t;

  // Constructor

  CICScatter(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  CICScatter(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~CICScatter()
  {
  }

  // apply method implements CIC scatter

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
      PointType_t gpos, delta;
      int idim;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);

          // This is the nearest-grid-point, so convert this to the
          // lower-grid-point (LGP) by comparing grid-point position
          // with particle position and adjusting as needed.

          gpos = mesh.vertexPosition(indx);
          for (idim=0; idim<Dim; ++idim)
            if (gpos(idim)>pos(i)(idim))
              indx[idim] = indx[idim] - 1;

          // now compute position and spacings at the LGP

          gpos = mesh.vertexPosition(indx);
          delta = mesh.vertexPosition(indx+1) - gpos;

          // From this, we find the normalized distance between 
          // the particle and LGP positions.

          for (idim=0; idim<Dim; ++idim)
            delta(idim) = (pos(i)(idim) - gpos(idim)) / delta(idim);

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          CICScatterFcn(attrib(i),fpatch,indx,delta);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class T, class ValueT>
struct CICScatterValue
{
  // Typedefs

  typedef CICScatterValue<FC,Dim,T,ValueT> This_t;
  typedef FC                               Field_t;
  typedef int                              PatchID_t;
  typedef int                              Size_t;
  typedef ValueT                           Value_t;

  // Constructor

  CICScatterValue(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  CICScatterValue(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~CICScatterValue()
  {
  }

  // apply method implements CIC scatter

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
      PointType_t gpos, delta;
      int idim;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // This is the nearest-grid-point, so convert this to the
          // lower-grid-point (LGP) by comparing grid-point position
          // with particle position and adjusting as needed.

          gpos = mesh.vertexPosition(indx);
          for (idim=0; idim<Dim; ++idim)
            if (gpos(idim)>pos(i)(idim))
              indx[idim] = indx[idim] - 1;

          // now compute position and spacings at the LGP

          gpos = mesh.vertexPosition(indx);
          delta = mesh.vertexPosition(indx+1) - gpos;

          // From this, we find the normalized distance between 
          // the particle and LGP positions.

          for (idim=0; idim<Dim; ++idim)
            delta(idim) = (pos(i)(idim) - gpos(idim)) / delta(idim);

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          CICScatterFcn(value_m,fpatch,indx,delta);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};

//-----------------------------------------------------------------------------
// Definitions for CICGatherFillCache, CICScatterFillCache and
// CICScatterValueFillCache functors
//-----------------------------------------------------------------------------

template <class FC, int Dim, class T>
struct CICGatherFillCache
{
  // Typedefs

  typedef CICGatherFillCache<FC,Dim,T> This_t;
  typedef FC                           Field_t;
  typedef int                          PatchID_t;
  typedef int                          Size_t;

  // Constructor

  CICGatherFillCache(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  CICGatherFillCache(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~CICGatherFillCache()
  {
  }

  // apply method implements CIC gather and fills cache

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
      PointType_t gpos, delta;
      int idim;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));

          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // This is the nearest-grid-point, so convert this to the
          // lower-grid-point (LGP) by comparing grid-point position
          // with particle position and adjusting as needed.

          gpos = mesh.vertexPosition(indx);
          for (idim=0; idim<Dim; ++idim)
            if (gpos(idim)>pos(i)(idim))
              indx[idim] = indx[idim] - 1;

          // now compute position and spacings at the LGP

          gpos = mesh.vertexPosition(indx);
          delta = mesh.vertexPosition(indx+1) - gpos;

          // From this, we find the normalized distance between 
          // the particle and LGP positions.

          for (idim=0; idim<Dim; ++idim)
            delta(idim) = (pos(i)(idim) - gpos(idim)) / delta(idim);

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the gather

          CICGatherFcn(attrib(i),fpatch,indx,delta);

          // Cache this interpolation data

          cache(i).lgp = indx;
          cache(i).dist = delta;
        }
    }

  // Copy of Field to be gathered from

  Field_t field_m;
};

template <class FC, int Dim, class T>
struct CICScatterFillCache
{
  // Typedefs

  typedef CICScatterFillCache<FC,Dim,T> This_t;
  typedef FC                            Field_t;
  typedef int                           PatchID_t;
  typedef int                           Size_t;

  // Constructor

  CICScatterFillCache(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  CICScatterFillCache(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~CICScatterFillCache()
  {
  }

  // apply method implements CIC scatter and fills cache

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
      PointType_t gpos, delta;
      int idim;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // This is the nearest-grid-point, so convert this to the
          // lower-grid-point (LGP) by comparing grid-point position
          // with particle position and adjusting as needed.

          gpos = mesh.vertexPosition(indx);
          for (idim=0; idim<Dim; ++idim)
            if (gpos(idim)>pos(i)(idim))
              indx[idim] = indx[idim] - 1;

          // now compute position and spacings at the LGP

          gpos = mesh.vertexPosition(indx);
          delta = mesh.vertexPosition(indx+1) - gpos;

          // From this, we find the normalized distance between 
          // the particle and LGP positions.

          for (idim=0; idim<Dim; ++idim)
            delta(idim) = (pos(i)(idim) - gpos(idim)) / delta(idim);

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          CICScatterFcn(attrib(i),fpatch,indx,delta);

          // Cache this interpolation data

          cache(i).lgp = indx;
          cache(i).dist = delta;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class T, class ValueT>
struct CICScatterValueFillCache
{
  // Typedefs

  typedef CICScatterValueFillCache<FC,Dim,T,ValueT> This_t;
  typedef FC                                        Field_t;
  typedef int                                       PatchID_t;
  typedef int                                       Size_t;
  typedef ValueT                                    Value_t;

  // Constructor

  CICScatterValueFillCache(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  CICScatterValueFillCache(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~CICScatterValueFillCache()
  {
  }

  // apply method implements CIC scatter and fills cache

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
      PointType_t gpos, delta;
      int idim;
      for (i=0; i<n; ++i)
        {
          // Convert the particle position to an index into the Field's
          // domain using the Geometry.
          
          indx = mesh.cellContaining(pos(i));
          
          // check we are on the right patch
          
          PAssert(layout.globalID(indx) == gid);
          
          // This is the nearest-grid-point, so convert this to the
          // lower-grid-point (LGP) by comparing grid-point position
          // with particle position and adjusting as needed.

          gpos = mesh.vertexPosition(indx);
          for (idim=0; idim<Dim; ++idim)
            if (gpos(idim)>pos(i)(idim))
              indx[idim] = indx[idim] - 1;

          // now compute position and spacings at the LGP

          gpos = mesh.vertexPosition(indx);
          delta = mesh.vertexPosition(indx+1) - gpos;

          // From this, we find the normalized distance between 
          // the particle and LGP positions.

          for (idim=0; idim<Dim; ++idim)
            delta(idim) = (pos(i)(idim) - gpos(idim)) / delta(idim);

          // now pass the necessary information to a functor that is 
          // dimension-specific and will do the scatter

          CICScatterFcn(value_m,fpatch,indx,delta);

          // Cache this interpolation data

          cache(i).lgp = indx;
          cache(i).dist = delta;
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};


//-----------------------------------------------------------------------------
// Definitions for CICGatherUseCache, CICScatterUseCache and
// CICScatterValueUseCache functors
//-----------------------------------------------------------------------------

template <class FC, int Dim, class T>
struct CICGatherUseCache
{
  // Typedefs

  typedef CICGatherUseCache<FC,Dim,T> This_t;
  typedef FC                          Field_t;
  typedef int                         PatchID_t;
  typedef int                         Size_t;

  // Constructor

  CICGatherUseCache(const Field_t& field)
    : field_m(field)
  {
  }

  // Copy constructor

  CICGatherUseCache(const This_t& model)
    : field_m(model.field_m)
  {
  }

  // Destructor

  ~CICGatherUseCache()
  {
  }

  // apply method implements CIC gather using cached data

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

          CICGatherFcn(attrib(i),fpatch,cache(i).lgp,cache(i).dist);
        }
    }

  // Copy of Field to be gathered from

  Field_t field_m;
};

template <class FC, int Dim, class T>
struct CICScatterUseCache
{
  // Typedefs

  typedef CICScatterUseCache<FC,Dim,T> This_t;
  typedef FC                           Field_t;
  typedef int                          PatchID_t;
  typedef int                          Size_t;

  // Constructor

  CICScatterUseCache(const Field_t& field)
    : field_m(field)
    {
    }

  // Copy constructor

  CICScatterUseCache(const This_t& model)
    : field_m(model.field_m)
    {
    }

  // Destructor

  ~CICScatterUseCache()
  {
  }

  // apply method implements CIC scatter using cached data

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

          CICScatterFcn(attrib(i),fpatch,cache(i).lgp,cache(i).dist);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;
};

template <class FC, int Dim, class T, class ValueT>
struct CICScatterValueUseCache
{
  // Typedefs

  typedef CICScatterValueUseCache<FC,Dim,T,ValueT> This_t;
  typedef FC                                       Field_t;
  typedef int                                      PatchID_t;
  typedef int                                      Size_t;
  typedef ValueT                                   Value_t;

  // Constructor

  CICScatterValueUseCache(const Field_t& field, const Value_t& value)
    : field_m(field), value_m(value)
    {
    }

  // Copy constructor

  CICScatterValueUseCache(const This_t& model)
    : field_m(model.field_m), value_m(model.value_m)
    {
    }

  // Destructor

  ~CICScatterValueUseCache()
  {
  }

  // apply method implements CIC scatter using cached data

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

          CICScatterFcn(value_m,fpatch,cache(i).lgp,cache(i).dist);
        }
    }

  // Copy of Field to be scattered into

  Field_t field_m;

  // Value to be scattered

  Value_t value_m;
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTICLES_INTERPOLATORCIC_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: InterpolatorCIC.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
