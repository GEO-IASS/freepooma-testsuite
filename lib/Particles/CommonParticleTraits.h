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
// MPDynamicUniform, MPDynamicSpatial,
// MPRemoteDynamicUniform, MPRemoteDynamicSpatial
//-----------------------------------------------------------------------------

#ifndef POOMA_PARTICLES_COMMON_PARTICLE_TRAITS_H
#define POOMA_PARTICLES_COMMON_PARTICLE_TRAITS_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * These are some structs which can be used as a particle traits class
 * for defining the Particles' attribute engine tag type and layout
 * strategy type.  These are provided for the user's convenience, to
 * save their having to define it themselves.
 */

//-----------------------------------------------------------------------------
// Include files:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

struct DynamicTag;
struct Dynamic;
template <class Tag> struct Remote;
template <class LayoutTag, class PatchTag> struct MultiPatch;
class UniformLayout;
template <class Mesh, class FieldLayout> class SpatialLayout;


//-----------------------------------------------------------------------------
// MPDynamicUniform:  MP Dynamic engine for the particle attributes
// and a UniformLayout strategy.  Use MP Dynamic for single context only.
//-----------------------------------------------------------------------------

struct MPDynamicUniform
{
  typedef MultiPatch<DynamicTag,Dynamic> AttributeEngineTag_t;
  typedef UniformLayout                  ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// MPDynamicSpatial:  MP Dynamic engine for the particle attributes
// and a SpatialLayout strategy.  Use MP Dynamic for single context only.
// The template parameters are the mesh type and the field layout type.
//-----------------------------------------------------------------------------

template <class Mesh, class FL>
struct MPDynamicSpatial
{
  typedef MultiPatch<DynamicTag,Dynamic>       AttributeEngineTag_t;
  typedef Mesh                                 Mesh_t;
  typedef SpatialLayout<Mesh_t,FL>             ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// MPRemoteDynamicUniform:  MP Remote Dynamic engine for particle attributes
// and a UniformLayout strategy.  Use MP Remote Dynamic for multiple contexts.
//-----------------------------------------------------------------------------

struct MPRemoteDynamicUniform
{
  typedef MultiPatch< DynamicTag, Remote<Dynamic> > AttributeEngineTag_t;
  typedef UniformLayout                             ParticleLayout_t;
};


//-----------------------------------------------------------------------------
// MPRemoteDynamicSpatial:  MP Remote Dynamic engine for particle attributes
// and a SpatialLayout strategy.  Use MP Remote Dynamic for multiple contexts.
// The template parameters are the mesh type and the field layout type.
//-----------------------------------------------------------------------------

template <class Mesh, class FL>
struct MPRemoteDynamicSpatial
{
  typedef MultiPatch< DynamicTag, Remote<Dynamic> > AttributeEngineTag_t;
  typedef Mesh                                      Mesh_t;
  typedef SpatialLayout<Mesh_t,FL>                  ParticleLayout_t;
};



//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PARTICLES_COMMON_PARTICLE_TRAITS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CommonParticleTraits.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
