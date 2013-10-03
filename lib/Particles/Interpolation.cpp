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
// Include files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <iostream>

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template <class T, class E> class DynamicArray;


//-----------------------------------------------------------------------------
// Gather/Scatter global function definitions for 
// Particle-Field interpolation
//-----------------------------------------------------------------------------

// gather/scatter using particle position attribute

template <class PA, class FC, class PPos, class InterpolatorTag>
void gather(const PA& attrib, const FC& field, const PPos& pos,
	    const InterpolatorTag&)
{
  // check that dimensions of field and position attribute match
  typedef typename PPos::Element_t PPosElement_t;
  CTAssert(int(FC::dimensions) == int(PPosElement_t::d1));

  // ask Interpolator class to perform the gather
  typedef typename PPosElement_t::Element_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::gather(attrib,field,pos);
}

template <class PA, class FC, class PPos, class InterpolatorTag>
void scatter(const PA& attrib, const FC& field, const PPos& pos,
	     const InterpolatorTag&)
{
  // check that dimensions of field and position attribute match
  typedef typename PPos::Element_t PPosElement_t;
  CTAssert(int(FC::dimensions) == int(PPosElement_t::d1));

  // ask Interpolator class to perform the scatter
  typedef typename PPosElement_t::Element_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::scatter(attrib,field,pos);
}

template <class T, class FC, class PPos, class InterpolatorTag>
void scatterValue(const T& value, const FC& field, const PPos& pos,
		  const InterpolatorTag&)
{
  // check that dimensions of field and position attribute match
  typedef typename PPos::Element_t PPosElement_t;
  CTAssert(int(FC::dimensions) == int(PPosElement_t::d1));

  // ask Interpolator class to perform the scatter
  typedef typename PPosElement_t::Element_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::scatterValue(value,field,pos);
}

// gather/scatter using particle position attribute and
// cache interpolation data

template <class PA, class FC, class PPos, class Cache, class InterpolatorTag>
void gatherCache(const PA& attrib, const FC& field, const PPos& pos,
                 const Cache& cache, const InterpolatorTag&)
{
  // check that dimensions of field and position attribute match
  typedef typename PPos::Element_t PPosElement_t;
  CTAssert((int)FC::dimensions == (int)PPosElement_t::d1);

  // ask Interpolator class to perform the gather and cache data
  typedef typename PPosElement_t::Element_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::gatherCache(attrib,field,pos,cache);
}

template <class PA, class FC, class PPos, class Cache, class InterpolatorTag>
void scatterCache(const PA& attrib, const FC& field, const PPos& pos,
                  const Cache& cache, const InterpolatorTag&)
{
  // check that dimensions of field and position attribute match
  typedef typename PPos::Element_t PPosElement_t;
  CTAssert((int)FC::dimensions == (int)PPosElement_t::d1);

  // ask Interpolator class to perform the scatter and cache data
  typedef typename PPosElement_t::Element_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::scatterCache(attrib,field,pos,cache);
}

template <class T, class FC, class PPos, class Cache, class InterpolatorTag>
void scatterValueCache(const T& value, const FC& field, const PPos& pos,
                       const Cache& cache, const InterpolatorTag&)
{
  // check that dimensions of field and position attribute match
  typedef typename PPos::Element_t PPosElement_t;
  CTAssert(int(FC::dimensions) == int(PPosElement_t::d1));

  // ask Interpolator class to perform the scatter and cache data
  typedef typename PPosElement_t::Element_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::scatterValueCache(value,field,pos,cache);
}

// gather/scatter using cached interpolation data

template <class PA, class FC, class Cache, class InterpolatorTag>
void gatherCache(const PA& attrib, const FC& field, const Cache& cache,
                 const InterpolatorTag&)
{
  // check that dimensions of field and cache data match
  typedef typename Cache::Element_t CacheData_t;
  CTAssert((int)FC::dimensions == (int)CacheData_t::dimensions);

  // ask Interpolator class to perform the gather using cached data
  typedef typename CacheData_t::AxisType_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::gatherCache(attrib,field,cache);
}

template <class PA, class FC, class Cache, class InterpolatorTag>
void scatterCache(const PA& attrib, const FC& field, const Cache& cache,
                  const InterpolatorTag&)
{
  // check that dimensions of field and cache data match
  typedef typename Cache::Element_t CacheData_t;
  CTAssert((int)FC::dimensions == (int)CacheData_t::dimensions);

  // ask Interpolator class to perform the scatter using cached data
  typedef typename CacheData_t::AxisType_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::scatterCache(attrib,field,cache);
}

template <class T, class FC, class Cache, class InterpolatorTag>
void scatterValueCache(const T& value, const FC& field, const Cache& cache,
                       const InterpolatorTag&)
{
  // check that dimensions of field and cache data match
  typedef typename Cache::Element_t CacheData_t;
  CTAssert((int)FC::dimensions == (int)CacheData_t::dimensions);

  // ask Interpolator class to perform the scatter using cached data
  typedef typename CacheData_t::AxisType_t AxisType_t;
  typedef Interpolator<FC::dimensions,AxisType_t,InterpolatorTag> Interp_t;
  Interp_t::scatterValueCache(value,field,cache);
}


template <class Field>
void setExternalGuards(const Field& f, typename Field::Element_t v)
{
  for (int i=0; i<Field::dimensions; ++i) {
    int d = f.layout().externalGuards().lower(i);
    if (d>0) {
      Interval<Field::dimensions> I(f.totalDomain());
      I[i] = Interval<1>(I[i].first(), I[i].first() + d-1);
      f(I) = v;
    }
    d = f.layout().externalGuards().upper(i);
    if (d>0) {
      Interval<Field::dimensions> I(f.totalDomain());
      I[i] = Interval<1>(I[i].last() - d+1, I[i].last());
      f(I) = v;
    }
  }
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Interpolation.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
