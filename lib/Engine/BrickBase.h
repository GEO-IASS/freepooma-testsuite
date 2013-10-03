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
//   BrickBase<Dim> 
//     - Base class for Brick-like engines
//   BrickViewBase<Dim> 
//     - Base class for BrickView-like engines
//   DoubleSliceHelper<Dim,Dim2,BaseDim>
//     - Helper functions for initializing slice-of-slice BrickViewBases
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_BRICKBASE_H
#define POOMA_ENGINE_BRICKBASE_H

/** @file
 * @ingroup Engine
 * @brief
 * Base classes for Brick- & BrickView-like engines. 
 * Encapsulates domain, stride, and subsetting operations.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Loc.h"
#include "Domain/SliceRange.h"
#include "Domain/SliceInterval.h"
#include "Layout/DomainLayout.h"
#include "Utilities/PAssert.h"

#if !defined(NOPassert)
#include <algorithm>
#endif

///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

namespace Pooma {

template <int Dim> class BrickViewBase;



/**
 * BrickBase<Dim> is the base-class for engines that have a Brick-like 
 * interpretation of a block of data; i.e. whose data is stored in a single
 * contiguous block of memory that is interpreted as a Dim-dimensional
 * "brick" using Fortran storage conventions.
 *
 * BrickBase caches a copy of the original strides and can zero and restore
 * the strides array at will. These operations are used by CompressibleBricks
 * and are only enabled if the BrickBase constructor specifies that the
 * compressible flag is true.
 *
 * The template parameters are:
 *   - Dim, the dimension of the Brick domain.
 *
 * BrickBase encapsulates the strides, domain, and the calculation
 * of indexing offsets brick-like engines.
 *
 * BrickViewBase<Dim>, which is defined below, is used to represent 
 * subsets of a data block described by a BrickBase. 
 * BrickViewBase serves as the base class for brick-view types of engines.
 */

template <int Dim>
class BrickBase 
{
public:

  //============================================================
  // Exported typedefs and constants
  //============================================================

  typedef Interval<Dim>  Domain_t;
  typedef DomainLayout<Dim> Layout_t;
  
  typedef BrickBase<Dim> This_t;
  
  enum { dimensions = Dim   };
  enum { brick      = true  };
  enum { zeroBased  = false };

  //============================================================
  // Constructors and Factory Methods
  //============================================================

  /// Default constructor. Creates an uninitialized BrickBase
  /// (i.e. empty domain, uninitialized strides, etc.)
  
  explicit BrickBase(bool compressible = false)
    : compressibleBase_m(compressible)
  { }

  /// Initialize the BrickBase with an Interval describing the Dim-
  /// dimensional domain to be indexed.
  
  explicit BrickBase(const Domain_t &domain, bool compressible = false);

  /// Initialize the BrickBase with a Layout describing the Dim-
  /// dimensional domain to be indexed.
  
  explicit BrickBase(const Node<Domain_t> &node, bool compressible = false);

  /// Initialize the BrickBase with a Layout describing the Dim-
  /// dimensional domain to be indexed.
  
  explicit BrickBase(const Layout_t &layout, bool compressible = false);

  // Copy constructor - Use compiler generated version.

  //============================================================
  // Destructor
  //============================================================

  ~BrickBase() {}

  //============================================================
  // Assignment operators - Use compiler generated version.
  //============================================================

  //============================================================
  // Accessor functions:
  //============================================================

  /// Return the domain.

  inline const Domain_t &domain() const { return layout_m.domain(); }

  /// Return the layout.

  inline const Layout_t &layout() const { return layout_m; }

  /// Return the strides array.
  
  inline const int *strides() const { return &strides_m[0]; }
  
  /// Return the true strides array.
  
  inline const int *originalStrides() const { return &ostrides_m[0]; }

  /// Used by constructors
  
  bool compressibleBase() const { return compressibleBase_m; }

// protected: // these should probably be protected, but we use them
              // publicly in testing.

  // Offset Calculations
  
  // These are provided for Domains and for sets of ints. 
  // Only BrickBase provides the offset0 & offsetC versions. 
  // offset0 treats the input domain as a zero-based domain giving 
  // offsets from the beginning in each dimension, rather than as 
  // points in the Brick's underlying domain. offsetC is similar, 
  // except that it multiplies the first offset by strides_m[0].
  // This is only for use in CompressibleBrick.
  
  template <class Domain>
  inline int offset(const Domain &dom) const
  { return off_m + offset0(dom); }
  
  template <class Domain>
  inline int offset0(const Domain &dom) const
  {
    CTAssert(Domain::dimensions == Dim);
    int offset = dom[0].first();
    for (int d = 1; d < Dim; ++d)
      offset += dom[d].first() * strides_m[d];
    return offset;
  }
  
  template <class Domain>
  inline int offsetC(const Domain &dom) const
  {
    CTAssert(Domain::dimensions == Dim);
    int offset = dom[0].first() * strides_m[0];
    for (int d = 1; d < Dim; ++d)
      offset += dom[d].first() * strides_m[d];
    return offset;
  }

  // The "int" versions do not assert that their dimensionality is
  // correct. The inheriting class should make these checks.
  // These are fairly short and there are a lot of them so I define
  // them in the class body for simplicity.

  inline int offset() const
  { return off_m; }
  
  inline int baseOffset() const
  { return off_m; }
  
  inline int offset(int i0) const
  { return off_m + i0; }
  
  inline int offset(int i0, int i1) const
  { return off_m + i0 + i1*strides_m[1]; }
  
  inline int offset(int i0, int i1, int i2) const
  { return off_m + i0 + i1*strides_m[1] + i2*strides_m[2]; }

  inline int offset(int i0, int i1, int i2, int i3) const
  { return off_m + i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]; }
                    
  inline int offset(int i0, int i1, int i2, int i3, int i4) const
  { return off_m + i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]
                      + i4*strides_m[4]; }
                    
  inline int offset(int i0, int i1, int i2, int i3, int i4, int i5) const
  { return off_m + i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]
                      + i4*strides_m[4] + i5*strides_m[5]; }
                    
  inline int offset(int i0, int i1, int i2, int i3, int i4, int i5, int i6) 
  const
  { return off_m + i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3] 
                      + i4*strides_m[4] + i5*strides_m[5] + i6*strides_m[6]; }
                    
  inline int offset0(int i0) const
  { return i0; }
  
  inline int offset0(int i0, int i1) const
  { return i0 + i1*strides_m[1]; }
  
  inline int offset0(int i0, int i1, int i2) const
  { return i0 + i1*strides_m[1] + i2*strides_m[2]; }
  
  inline int offset0(int i0, int i1, int i2, int i3) const
  { return i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]; }
  
  inline int offset0(int i0, int i1, int i2, int i3, int i4) const
  { return i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]
              + i4*strides_m[4]; }
  
  inline int offset0(int i0, int i1, int i2, int i3, int i4, int i5) const
  { return i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]
              + i4*strides_m[4] + i5*strides_m[5]; }
  
  inline int offset0(int i0, int i1, int i2, int i3, int i4, int i5, int i6) 
  const
  { return i0 + i1*strides_m[1] + i2*strides_m[2] + i3*strides_m[3]
              + i4*strides_m[4] + i5*strides_m[5] + i6*strides_m[6]; }
    
  inline int offsetC(int i0) const 
  { return i0*strides_m[0]; }
  
  inline int offsetC(int i0, int i1) const 
  { return i0*strides_m[0] + i1*strides_m[1]; }
  
  inline int offsetC(int i0, int i1, int i2) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2]; }

  inline int offsetC(int i0, int i1, int i2, int i3) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3]; }
                    
  inline int offsetC(int i0, int i1, int i2, int i3, int i4) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3] + i4*strides_m[4]; }
                    
  inline int offsetC(int i0, int i1, int i2, int i3, int i4, int i5) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3] + i4*strides_m[4] + i5*strides_m[5]; }
                    
  inline int offsetC(int i0, int i1, int i2, int i3, int i4, int i5, int i6) 
  const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3] + i4*strides_m[4] + i5*strides_m[5] 
         + i6*strides_m[6]; }
  
protected:

  //============================================================
  // Mutator functions:
  //============================================================

  // These are used to modify the strides array. They are only 
  // used by compressible engines and will assert if called
  // on a noncompressible engine.
  
  void zeroStrides() 
  { for (int d = 0; d < Dim; ++d) strides_m[d] = 0; }
  
  void restoreStrides() 
  { for (int d = 0; d < Dim; ++d) strides_m[d] = ostrides_m[d]; }
  
  //============================================================
  // Protected data
  //============================================================

  // Layout.

  Layout_t layout_m;

  // Strides through actual data block when stepping in different dimensions.
  // We keep two copies - strides_m is used by the offset calculations. 
  // If we are compressible, then when compressed, these will all be 
  // set to zero. 
  
  int strides_m[Dim];
  
  int ostrides_m[Dim];
    
  // Offset due to non-zero first elements.

  int off_m;
  
  // Flag indicating whether or not the stride compressions routines
  // are callable. This is only set by the constructor and cannot be
  // read or modified. It could be integrated into the type (as another
  // template parameter), but that would complicate certain operations,
  // and the benefit would probably not be worth the cost in complexity.
  
  bool compressibleBase_m;
};



/**
 * This class is used to implement the slice-of-slice constructors.
 * These have to be member functions, but we bounce them to functions
 * in a templated class as Metrowerks does not handle out-of-line
 * member templates very well. 
 * 
 * We implement the functionality as static members of a class to make 
 * pre-instantiation easier. Function preinstantiation requires listing 
 * the fully specialized prototype, while class preinstantiation simply 
 * requires a line like "template class DoubleSliceHelper<1,2,3>;". 
 *
 * DoubleSliceHelper should probably be in a BrickUtil namespace or something. 
 */ 

template <int Dim, int Dim2>
struct DoubleSliceHelper
{
  typedef Interval<Dim>  Domain_t;
  typedef DomainLayout<Dim> Layout_t;
  
  static void init(Domain_t &domain, 
                   int *strides,
                   int &baseOffset,
                   const BrickViewBase<Dim2> &bvbase,
                   const SliceInterval<Dim2,Dim> &dom);

  static void init(Domain_t &domain, 
                   int *strides,
                   int &baseOffset,
                   const BrickViewBase<Dim2> &bvbase,
                   const SliceRange<Dim2,Dim> &dom);
};


/**
 * BrickViewBase<Dim> is the base-class for engines that are "views"
 * into brick-like engines.
 *
 * The template parameters are:
 *   - Dim, the logical dimension of the view
 *   - BaseDim, the dimension of the underlying BrickBase being viewed
 *
 * If Dim < BaseDim, then the view is called "sliced".
 * 
 * BrickViewBase encapsulates the calculations of the strides, domains,
 * and indexing offsets for these complicated views. 
 *
 * The general template specifies the "sliced" view. We specialize below
 * for Dim == BaseDim, which avoids the complications of slicing.
 */

template <int Dim>
class BrickViewBase
{
public:
  
  //============================================================
  // Exported typedefs and constants
  //============================================================

  enum { dimensions = Dim  };
  enum { zeroBased  = true };
  
  typedef Interval<Dim>  Domain_t;
  typedef DomainLayout<Dim> Layout_t;
  
  //============================================================
  // Constructors
  //============================================================

  /// Default constructor - creates an uninitialized BrickViewBase 
  /// (which has an empty domain with uninitilazed
  /// strides, etc.).
  
  BrickViewBase() {}

  // Copy constructor - Use compiler generated one.

  typedef BrickViewBase<Dim> This_t;

  //@{

  /// This is a special copy constructor that can change the
  /// compressibility flag. Useful, for example, when constructing
  /// a brick-view of a compressible-brick-view.
  
  BrickViewBase(const This_t &, bool compressible);

  BrickViewBase(const BrickBase<Dim> &base, bool compressible)
  {
    *this = BrickViewBase<Dim>(base, base.domain());
    compressibleBase_m = compressible;
    if (!compressible) restoreStrides();
  }

  //@}

  //@{

  /// Subsetting Constructors.
  /// This class specialization is for strided, sliced, brick-shaped 
  /// views. We have to provide constructors that build such views from 
  /// Bricks and other views and from the appropriate types of brick-
  /// shaped subdomains. For slices, there are six possibilities:
  /// (NOTE: There is no use made of the difference between 
  /// SliceInterval and SliceRange, so why have SliceInterval??? 
  /// Getting rid of it would reduce the number of constructors to 4.)

  template<int BaseDim>
  BrickViewBase(const BrickBase<BaseDim> &bbase,
		const SliceRange<BaseDim,Dim> &dom)
    : domain_m(Pooma::NoInit()), 
      baseOffset_m(bbase.offset()),
      compressibleBase_m(bbase.compressibleBase())
  {    
    sliceInit(bbase.originalStrides(), dom);
  }


  template<int BaseDim>
  BrickViewBase(const BrickBase<BaseDim> &bbase,
		const SliceInterval<BaseDim,Dim> &dom)
    : domain_m(Pooma::NoInit()), 
      baseOffset_m(bbase.offset()),
      compressibleBase_m(bbase.compressibleBase())
  {
    sliceInit(bbase.originalStrides(), SliceRange<BaseDim,Dim>(dom));
  }

                
  BrickViewBase(const BrickBase<Dim> &, const Interval<Dim> &);
  BrickViewBase(const BrickBase<Dim> &, const Range<Dim> &);
  
  BrickViewBase(const This_t &, const Interval<Dim> &);
  BrickViewBase(const This_t &, const Range<Dim> &);
  
  template <int Dim2>
  BrickViewBase(const BrickViewBase<Dim2> &bvbase,
                const SliceRange<Dim2,Dim> &dom)
  : domain_m(Pooma::NoInit()),
    baseOffset_m(bvbase.baseOffset())
  {
    DoubleSliceHelper<Dim,Dim2>::
    init(domain_m, strides_m, baseOffset_m, 
         bvbase, dom);
    for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
  }
  
  template <int Dim2>
  BrickViewBase(const BrickViewBase<Dim2> &bvbase, 
                const SliceInterval<Dim2,Dim> &dom)
  : domain_m(Pooma::NoInit()),
    baseOffset_m(bvbase.baseOffset())
  {
    DoubleSliceHelper<Dim,Dim2>::
    init(domain_m, strides_m, baseOffset_m,
         bvbase, dom);
    for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
  }

  //@}

  //============================================================
  // Destructor
  //============================================================

  ~BrickViewBase() {}

  //============================================================
  // Assignment operators - Use compiler generated one.
  //============================================================

  //============================================================
  // Accessor functions:
  //============================================================

  /// Return our logical domain:

  inline const Domain_t &domain() const { return domain_m; }

  /// Return our layout:

  inline Layout_t layout() const { return Layout_t(domain_m); }

  /// Return the strides array.

  inline const int *strides() const { return &strides_m[0]; }
  
  /// Return the true strides array.
  
  inline const int *originalStrides() const { return &ostrides_m[0]; }

  /// Return the first index value for the specified dimension.
  /// (Always zero since views are zero-based).
  
  inline int first(int) const { return 0; }

  /// Used by constructors
  
  bool compressibleBase() const { return compressibleBase_m; }

  /// Base offset - offset from beginning of underlying Brick's beginning.
  
  int baseOffset() const { return baseOffset_m; }
  
  // Offset calculations
  
  template <class Domain>
  inline int offset(const Domain &dom) const
  {
    CTAssert(Domain::dimensions == Dim);
    int offset = dom[0].first() * strides_m[0];
    for (int d = 1; d < Dim; ++d)
      offset += dom[d].first() * strides_m[d];
    return offset;
  }
  
  inline int offset(int i0) const 
  { return i0*strides_m[0]; }
  
  inline int offset(int i0, int i1) const 
  { return i0*strides_m[0] + i1*strides_m[1]; }
  
  inline int offset(int i0, int i1, int i2) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2]; }

  inline int offset(int i0, int i1, int i2, int i3) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3]; }
                    
  inline int offset(int i0, int i1, int i2, int i3, int i4) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3] + i4*strides_m[4]; }
                    
  inline int offset(int i0, int i1, int i2, int i3, int i4, int i5) const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3] + i4*strides_m[4] + i5*strides_m[5]; }
                    
  inline int offset(int i0, int i1, int i2, int i3, int i4, int i5, int i6)
    const
  { return i0*strides_m[0] + i1*strides_m[1] + i2*strides_m[2] 
         + i3*strides_m[3] + i4*strides_m[4] + i5*strides_m[5] 
         + i6*strides_m[6]; }
  
protected:

  //============================================================
  // Mutator functions:
  //============================================================

  // These are used to modify the strides array. They are only 
  // used by compressible engines and will assert if called
  // on a noncompressible engine.
  
  void zeroStrides() 
  { for (int d = 0; d < Dim; ++d) strides_m[d] = 0; }
  
  void restoreStrides() 
  { for (int d = 0; d < Dim; ++d) strides_m[d] = ostrides_m[d]; }
  
  //============================================================
  // Utility functions
  //============================================================

  //---------------------------------------------------------------------------
  //
  // sliceInit(const int * &, 
  //           const SliceRange<BaseDim,Dim> &)
  //
  // Helper function used in taking a slice of a Brick or a Brick-like view.
  //---------------------------------------------------------------------------

  template<int BaseDim>
  void sliceInit(const int *baseStrides,
		 const SliceRange<BaseDim,Dim> &dom)
  {    
    typedef typename SliceRange<BaseDim,Dim>::TotalDomain_t TotalDomain_t;
    const TotalDomain_t &domain = dom.totalDomain();
    
    int dt, d;    
    for (dt = 0, d = 0; dt < BaseDim; ++dt)
    {
      if (!dom.ignorable(dt))
      {
	PAssert(d < Dim);
	domain_m[d]  = Interval<1>(domain[dt].length());
	strides_m[d] = baseStrides[dt] * domain[dt].stride();
	++d;
      }
        
      baseOffset_m += domain[dt].first() * baseStrides[dt];

    }

    PAssert(d == Dim);
  
    for (int d = 0; d < Dim; ++d) ostrides_m[d] = strides_m[d];
  }

  void sliceInit(const This_t &, const Range<Dim> &domain);

  void viewInit(const BrickBase<Dim> &bbase, const Range<Dim> &domain)
  {
    for (int d = 0; d < Dim; ++d)
    {
      domain_m[d]   = Interval<1>(domain[d].length());
      strides_m[d]  = bbase.originalStrides()[d] * domain[d].stride();
      ostrides_m[d] = strides_m[d];
      baseOffset_m  += domain[d].first() * bbase.originalStrides()[d];
    }
  }
  
  //============================================================
  // Data
  //============================================================

  // Domain for this engine:
  // (Somewhat wasteful since we are zero-based and only need
  // the length. However, we may want to return references to this
  // object for efficiency reasons, and thus we can't create it
  // on the fly when someone asks.)
  
  Domain_t domain_m;

  // Strides through actual data block when stepping in each dimension.
  // We keep two copies - strides_m is used by the offset calculations. 
  // If we are compressible, then when compressed, these will all be 
  // set to zero. 
  
  int strides_m[Dim];

  int ostrides_m[Dim];
    
  // Base offset - offset of beginning of view from underlying Brick's
  // beginning.
  
  int baseOffset_m;
  
  // Compressibility flag - see BrickBase comments.
  
  bool compressibleBase_m;
};



} // namespace Pooma

///////////////////////////////////////////////////////////////////////////////

// Include .cpp file to get out-of-line functions.

#include "Engine/BrickBase.cpp"

#endif // POOMA_ENGINE_BRICKBASE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BrickBase.h,v $   $Author: richi $
// $Revision: 1.19 $   $Date: 2004/11/29 16:21:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
