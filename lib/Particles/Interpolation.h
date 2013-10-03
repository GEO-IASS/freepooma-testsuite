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

#ifndef POOMA_PARTICLES_INTERPOLATION_H
#define POOMA_PARTICLES_INTERPOLATION_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Particles
 * @brief
 * General template for Interpolator class and
 * Gather/Scatter global function declarations for 
 * Particle-Field interpolation
 *
 * Global functions for gathering Field values
 * into a Particle Attribute and scattering Particle Attribute values
 * into a Field, using the particle positions and an interpolation 
 * stencil.  Gather and Scatter functions take as an argument an
 * Interpolation tag that indicates what type of stencil to use.
 * These functions will create the right type of Interpolator object
 * and ask it to do the gather or scatter operation.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Evaluator/PatchFunction.h"


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


///@name gather/scatter functions
/// using particle position attribute
//@{

template <class PA, class FC, class PPos, class InterpolatorTag>
void gather(const PA&, const FC&, const PPos&,
	    const InterpolatorTag&);

template <class PA, class FC, class PPos, class InterpolatorTag>
void scatter(const PA&, const FC&, const PPos&,
	     const InterpolatorTag&);

template <class T, class FC, class PPos, class InterpolatorTag>
void scatterValue(const T&, const FC&, const PPos&,
		  const InterpolatorTag&);
//@}

///@name gather/scatter functions
/// using particle position attribute and cached interpolation data
//@{

template <class PA, class FC, class PPos, class Cache, class InterpolatorTag>
void gatherCache(const PA&, const FC&, const PPos&, const Cache&,
		 const InterpolatorTag&);

template <class PA, class FC, class PPos, class Cache, class InterpolatorTag>
void scatterCache(const PA&, const FC&, const PPos&, const Cache&,
		  const InterpolatorTag&);

template <class T, class FC, class PPos, class Cache, class InterpolatorTag>
void scatterValueCache(const T&, const FC&, const PPos&, const Cache&,
		       const InterpolatorTag&);
//@}

///@name gather/scatter functions
/// using cached interpolation data
//@{

template <class PA, class FC, class Cache, class InterpolatorTag>
void gatherCache(const PA&, const FC&, const Cache&,
		 const InterpolatorTag&);

template <class PA, class FC, class Cache, class InterpolatorTag>
void scatterCache(const PA&, const FC&, const Cache&,
		  const InterpolatorTag&);

template <class T, class FC, class Cache, class InterpolatorTag>
void scatterValueCache(const T&, const FC&, const Cache&,
		       const InterpolatorTag&);
//@}

/////////////////////////////////////////////////////////////////////////

/** @class Interpolator
 * Interpolator is a general template for a class that performs
 * interpolation of data between arbitrary points in space
 * (e.g., particle positions) and field element positions.
 * Interpolator is templated on the dimensionality and axis type
 * of the spatial positions and a tag type that indicates what
 * sort of interpolation scheme to use.
 *
 * Each specialization of Interpolator should provide nine gather/scatter
 * methods corresponding to the nine global functions declared 
 * above.  Each method has the same interface as the corresponding
 * global function except that the InterpolatorTag argument is not
 * needed.  These methods actually implement the gather or scatter
 * using PatchParticle functors and dimension-specific specializations
 * of the actual interpolation computational kernel.
 *
 * In addition to these gather/scatter methods, the Interpolator
 * specialization should export a typedef called Cache_t that is
 * the type of an object suitable for caching all the necessary
 * interpolation data for this type of Interpolator.  The actual 
 * type that Cache_t refers to can be defined external to the 
 * Interpolator specialization, and an insertion operator for 
 * sending a Cache_t type to an ostream should be provided.
 *
 * Here is an example of what this might look like:
 * @code

// New interpolator tag
struct Foo { };

// declaration of interpolation data cache type
struct Bar;

// Definition of Interpolator class template specialization
template <int Dim, class T>
struct Interpolator<Dim,T,Foo>
{
  // typedef for the CacheData struct
  typedef Bar Cache_t;

  // gather/scatter using particle position attribute
  template <class PA, class FC, class PPos>
  static
  void gather(const PA&, const FC&, const PPos&) { ... }

  template <class PA, class FC, class PPos>
  static
  void scatter(const PA&, const FC&, const PPos&) { ... }

  template <class ValueT, class FC, class PPos>
  static
  void scatterValue(const ValueT&, const FC&, const PPos&) { ... }

  // gather/scatter using particle position attribute and cached
  // interpolation data
  template <class PA, class FC, class PPos, class EngineTag>
  static
  void gatherCache(const PA&, const FC&, const PPos&,
                   const DynamicArray<Cache_t,EngineTag>&) { ... }

  template <class PA, class FC, class PPos, class EngineTag>
  static
  void scatterCache(const PA&, const FC&, const PPos&,
                    const DynamicArray<Cache_t,EngineTag>&) { ... }

  template <class ValueT, class FC, class PPos, class EngineTag>
  static
  void scatterValueCache(const ValueT&, const FC&, const PPos&,
                         const DynamicArray<Cache_t,EngineTag>&) { ... }

  // gather/scatter using cached interpolation data
  template <class PA, class FC, class EngineTag>
  static
  void gatherCache(const PA&, const FC&,
                   const DynamicArray<Cache_t,EngineTag>&) { ... }

  template <class PA, class FC, class EngineTag>
  static
  void scatterCache(const PA&, const FC&,
                    const DynamicArray<Cache_t,EngineTag>&) { ... }

  template <class ValueT, class FC, class EngineTag>
  static
  void scatterValueCache(const ValueT&, const FC&,
                         const DynamicArray<Cache_t,EngineTag>&) { ... }
};

// definition of interpolation data cache type
struct Bar
{
  // Add data here for location of nearest grid point (a Loc<Dim>)
  // and distance to particle position (a Vector<Dim,AxisType_t>) as needed.
  // Also provide a templated print function that prints CacheData.
  // Publicly export Dim and AxisType_t.
};

// Send interpolation data cache type to an ostream
inline
std::ostream&
operator<<(std::ostream& o, const Bar&)
{
  // print out cache data
  return o;
}

 * @endcode
 *
 */

template <int Dim, class T, class InterpolatorTag>
struct Interpolator;


template <class Field>
void setExternalGuards(const Field&, typename Field::Element_t);


#include "Particles/Interpolation.cpp"


#endif     // POOMA_PARTICLES_INTERPOLATION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Interpolation.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
