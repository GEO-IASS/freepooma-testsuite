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
// GuardLayers<Dim>
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_GUARDLAYERS_H
#define POOMA_LAYOUT_GUARDLAYERS_H

/** @file
 * @ingroup Layout
 * @brief
 * A simple container for a set of guard layer specifications.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Loc.h"

/**
 * This class is a simple container for two arrays of Dim integers, 
 * specifying the numbers of guard layers at the upper and lower extent
 * of each dimension. 
 */

template <int Dim>
class GuardLayers
{
public:

  //============================================================
  // Constructors
  //============================================================

  explicit GuardLayers(int gcs = 0)
  {
    PAssert(gcs >= 0);
    for (int i = 0; i < Dim; ++i)
      {
        lower_m[i] = gcs;
        upper_m[i] = gcs;
      }
  }
      
  GuardLayers(int lower[Dim], int upper[Dim])
  {
    for (int i = 0; i < Dim; ++i)
      {
        PAssert(lower[i] >= 0 && upper[i] >= 0);
        lower_m[i] = lower[i];
        upper_m[i] = upper[i];
      }
  }
      
  GuardLayers(const Loc<Dim> &lower, const Loc<Dim> &upper)
  {
    for (int i = 0; i < Dim; ++i)
      {
        PAssert(lower[i].first() >= 0 && upper[i].first() >= 0);
        lower_m[i] = lower[i].first();
        upper_m[i] = upper[i].first();
      }
  }
  
  //============================================================
  // Copy constructor and assignment
  //============================================================

  // Default is fine.
  
  //============================================================
  // Initialization functions.
  //============================================================
  
  void initialize(const Loc<Dim> &lower, const Loc<Dim> &upper)
  {
    for (int i = 0; i < Dim; ++i)
      {
        PAssert(lower[i].first() >= 0 && upper[i].first() >= 0);
        lower_m[i] = lower[i].first();
        upper_m[i] = upper[i].first();
      }
  }

  void initialize(const GuardLayers<Dim> &gl)
  {
    *this = gl;
  }
  
  //============================================================
  // Accessors
  //============================================================
  
  int lower(int i) const
  { 
#if POOMA_BOUNDS_CHECK
    PInsist(i<Dim&&i>=0," GuardLayers index out of range ");
#endif
    return lower_m[i]; 
  }
  int upper(int i) const 
  {   
#if POOMA_BOUNDS_CHECK
    PInsist(i<Dim&&i>=0," GuardLayers index out of range ");
#endif
    return upper_m[i]; 
  }
  
  //============================================================
  // Mutators
  //============================================================
  
  int &lower(int i) 
  {    
#if POOMA_BOUNDS_CHECK
    PInsist(i<Dim&&i>=0," GuardLayers index out of range ");
#endif
    return lower_m[i]; 
  }
  int &upper(int i) 
  {    
#if POOMA_BOUNDS_CHECK
    PInsist(i<Dim&&i>=0," GuardLayers index out of range ");
#endif
    return upper_m[i]; 
  }
  
  //============================================================
  // Operators
  //============================================================
  
  bool operator==(const GuardLayers<Dim> &gcs) const
  {
    bool result = true;
    for (int d = 0; d < Dim; ++d)
      {
        result = result && lower_m[d] == gcs.lower_m[d];
        result = result && upper_m[d] == gcs.upper_m[d];
      }
    return result;
  }
  
  bool operator==(int gcw) const
  {
    bool result = true;
    for (int d = 0; d < Dim; ++d)
      {
        result = result && lower_m[d] == gcw;
        result = result && upper_m[d] == gcw;
      }
    return result;
  }
  
  bool operator!=(const GuardLayers<Dim> &gcs) const
  {
    return !operator==(gcs);
  }
  
  bool operator!=(int gcw) const
  {
    return !operator==(gcw);
  }
  
  GuardLayers<Dim> operator-(const GuardLayers<Dim> &gcs)
  {
    GuardLayers<Dim> result;
    for (int d = 0; d < Dim; ++d)
      {
        result.lower(d) = lower_m[d] - gcs.lower_m[d];
        PAssert(result.lower(d) >= 0);
        result.upper(d) = upper_m[d] - gcs.upper_m[d];
        PAssert(result.upper(d) >= 0);
      }
    return result; 
  }
    
  GuardLayers<Dim> operator-(int dw)
  {
    GuardLayers<Dim> result;
    for (int d = 0; d < Dim; ++d)
      {
        result.lower(d) = lower_m[d] - dw;
        PAssert(result.lower(d) >= 0);
        result.upper(d) = upper_m[d] - dw;
        PAssert(result.upper(d) >= 0);
      }
    return result; 
  }
    
  //============================================================
  // Utility function
  //============================================================
  
  inline static void 
  addGuardLayers(Interval<Dim> &dom, const GuardLayers<Dim> &gcs)
  {
    for (int d = 0; d < Dim; ++d)
      {
        int a = dom[d].first() - gcs.lower(d);
        int b = dom[d].last()  + gcs.upper(d);
        dom[d] = Interval<1>(a,b);
      }
  }
  
  Interval<Dim> addGuardLayersToDomain(const Interval<Dim> &d) const
  {
    Interval<Dim> dom(d);
    addGuardLayers(dom, *this);
    
    return dom;
  }

  // --------------------------------------------------------------------------
  // Print a GuardLayers<Dim> on an output stream.
  
  template <class Ostream>
  void print(Ostream &ostr) const
  {
    ostr << "GuardLayers<" << Dim << "> [";
    for (int d = 0; d < Dim; ++d)
      {
        ostr << "l: " << lower_m[d] << ", u: " << upper_m[d];
        if (d != Dim - 1)
          ostr << "; ";
      }
    ostr << "]";
  }
  
private:

  int lower_m[Dim];
  int upper_m[Dim];
};

template<int Dim>
Interval<Dim> &growInPlace(Interval<Dim> &dom, const GuardLayers<Dim> &gcs)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() - gcs.lower(d);
      int b = dom[d].last()  + gcs.upper(d);
      dom[d] = Interval<1>(a,b);
    }
  return dom;
}

template<int Dim>
Interval<Dim> &shrinkInPlace(Interval<Dim> &dom, const GuardLayers<Dim> &gcs)
{
  for (int d = 0; d < Dim; ++d)
    {
      int a = dom[d].first() + gcs.lower(d);
      int b = dom[d].last()  - gcs.upper(d);
      dom[d] = Interval<1>(a,b);
    }
  return dom;
}

template<int Dim>
inline Interval<Dim> 
grow(const Interval<Dim> &dom, const GuardLayers<Dim> &gcs)
{
  Interval<Dim> ret(dom);
  return growInPlace(ret, gcs);
}

template<int Dim>
inline Interval<Dim> 
shrink(const Interval<Dim> &dom, const GuardLayers<Dim> &gcs)
{
  Interval<Dim> ret(dom);
  return shrinkInPlace(ret, gcs);
}

//-----------------------------------------------------------------------------
//
// ostream inserter for GuardLayers<Dim>.
//
//-----------------------------------------------------------------------------

template<int Dim>
std::ostream &operator<<(std::ostream &ostr, 
  const GuardLayers<Dim> &gl)
{
  gl.print(ostr);
  return ostr;
}

#endif     // POOMA_LAYOUT_GUARDLAYERS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: GuardLayers.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
