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
//   ViewIndexer           - A utility class for doing index calculations
//                           for view-type engines.
//-----------------------------------------------------------------------------

#ifndef POOMA_UTILITIES_VIEWINDEXER_H
#define POOMA_UTILITIES_VIEWINDEXER_H

/** @file
 * @ingroup Utilities
 * @brief
 * ViewIndexer translates a set of "local" indices for a view of some domain
 * into the "base" coordinates of the domain that ultimately spawned the view.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Domain/Range.h"
#include "Domain/SliceDomain.h"


/**
 * ViewIndexer translates indices in a Dim-dimensional domain back to the
 * corresponding indices in the original Dim2-dimensional domain, taking into
 * account things like slices.
 *
 * Exported typedefs:
 *  - BaseDomain_t - the type of base domain for this ViewIndexer.
 *  - Domain_t - the type of domain for this ViewIndexer.
 *  - Mask_t - the type of mask for this ViewIndexer.
 *  - This_t - this class.
 *
 * Constructors:
 *  - ViewIndexer() - default initialization creates empty domains and leaves
 *    the rest uninitialized. Initialization can be completed with operator=().
 *  - ViewIndexer(const SliceDomain<DT> &) - constructs a ViewIndexer for
 *     the first slice-type view of a domain.
 *  - ViewIndexer(const ViewIndexer<Dim, Dim2> &, const Domain<Dim, DT> &) -
 *     constructs a new ViewIndexer from another one by taking a non-slice
 *     view.
 *  - ViewIndexer(const ViewIndexer<OrigDim, Dim2> &, const SliceDomain<DT> &)
 *     constructs a new ViewIndexer from another one by taking a slice.
 *  - ViewIndexer(const This_t &) - copy constructor.
 *
 * Operators:
 *  - operator=(const This_t &) - copy assignment operator.
 *
 * Accessors:
 *  - baseDomain() - returns the base domain.
 *  - domain() - returns the domain.
 *  - indirection(int i ) - returns the ith indirection index.
 *  - mask() - returns the mask.
 *  - offset(int i ) - returns the ith offset.
 *  - stride(int i ) - returns the ith stride.
 *  - translate(int i0[, ... int i6], Loc<Dim2>) - translates 1-7 indices to
 *     a Loc in the base coordinate system.
 *
 * Utility functions:
 *   - localToBase - transforms a domain in the local coordinate system to
 *      the corresponding domain in base coordinates.
 *   - baseToLocal - transforms a domain in the base coordinate system to
 *      the corresponding domain in local coordinates.
 */

template<int Dim, int Dim2>
class ViewIndexer
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.

  // This class.
  
  typedef ViewIndexer<Dim, Dim2> This_t;
  
  // The domain type.
  
  typedef Interval<Dim> Domain_t;
  
  // The base domain type.
  
  typedef Range<Dim2> BaseDomain_t;
  
  // The mask type.
  
  typedef Loc<Dim2> Mask_t;

  //---------------------------------------------------------------------------
  // Constructors.

  // Sometimes it is convenient to use the default initializer and then
  // assign a new instance to an object, so we provide a default constructor
  // that initializes ViewIndexer with empty domains and leaves everything 
  // else uninitialized.
  
  ViewIndexer() { }
  
  template<class DT>
  ViewIndexer(const SliceDomain<DT> &dom)
  : domain_m(Pooma::NoInit()), baseDomain_m(dom.totalDomain())
  {
    // We are a slice and our dimensions must be consistent with us and the 
    // slice-domain we're being spawned by.

    CTAssert(Dim == DT::sliceDimensions);
    CTAssert(Dim2 == DT::dimensions);

    // Set up our offsets, domains, and strides. We are a slice so me must
    // do this using the non-ignorable domains. We set the mask according to
    // the ignorable domains

    int dt, d;
    const typename DT::TotalDomain_t &domain = dom.totalDomain();
    for (d = 0, dt = 0; dt < Dim2; ++dt)
      {
        if (!dom.ignorable(dt))
          {
            PAssert(d < Dim);
            offset_m[d] = domain[dt].first();
            stride_m[d] = domain[dt].stride();
            domain_m[d] = Interval<1>(domain[dt].length());
            ind_m[d] = dt;
            ++d;
          }
        else
          {
            PAssert(domain[dt].first() == domain[dt].last());
            mask_m[dt] = domain[dt].first();
          }
      }    
  }

  template<class DT>
  ViewIndexer(const ViewIndexer<Dim, Dim2> &orig, 
    const Domain<Dim, DT> &dom)
  : domain_m(Pooma::NoInit()), 
    baseDomain_m(Pooma::NoInit()), 
    mask_m(orig.mask())
  {
    // Fill in the base domain from the previous indexer. We'll overwrite
    // some dimensions below but we need to have the sliced out dimensions
    // in place.

    baseDomain_m = orig.baseDomain();

    // Set up our offsets, domains, and strides. This is easy since we're not
    // not being sliced by this domain meaning we don't need to fiddle with
    // the mask and indirection vectors.

    const typename DT::Domain_t &domain = dom.unwrap();
    for (int d = 0; d < Dim; ++d)
      {
        offset_m[d] = orig.offset(d) + orig.stride(d) * 
          domain[d].first();
        stride_m[d] = orig.stride(d) * domain[d].stride();
        domain_m[d] = Interval<1>(domain[d].length());
      
        ind_m[d] = orig.indirection(d);
      }

    localToBase(domain_m, baseDomain_m);
  }

  template<int OrigDim, class DT>
  ViewIndexer(const ViewIndexer<OrigDim, Dim2> &orig, 
    const SliceDomain<DT> &dom)
  : domain_m(Pooma::NoInit()), 
    baseDomain_m(orig.baseDomain()), mask_m(orig.mask())
  {  
    // Our dimensionality must be the same as the slice's reduced
    // dimensionality.
    
    CTAssert(DT::sliceDimensions == Dim);
    
    // The slice's dimensionality must match that of the previous view.
    
    CTAssert(DT::dimensions == OrigDim);
 
    // Set up our offsets, domains, and strides. We are a slice so me must
    // do this using the non-ignorable domains. We set the mask according to
    // the ignorable domains

    int dt, d;
    const typename DT::TotalDomain_t &domain = dom.totalDomain();
    for (d = 0, dt = 0; dt < OrigDim; ++dt)
      {
        if (!dom.ignorable(dt))
          {
            PAssert(d < Dim);
            offset_m[d] = orig.offset(dt) + orig.stride(dt) * 
              domain[dt].first();
            stride_m[d] = orig.stride(dt) * domain[dt].stride();
            domain_m[d] = Interval<1>(domain[dt].length());

            ind_m[d] = orig.indirection(dt);

            // Translate this part of the domain to base coordinates and
            // store.
            baseDomain_m[ind_m[d]] = Range<1>(
              offset_m[d],
              offset_m[d] + stride_m[d] * domain_m[d].last(), 
              stride_m[d]);

            ++d;
          }
        else
          {
            PAssert(domain[dt].first() == domain[dt].last());
            int m = orig.offset(dt) + orig.stride(dt) * 
              domain[dt].first();
            mask_m[orig.indirection(dt)] = m;
            baseDomain_m[orig.indirection(dt)] = 
              Range<1>(m, m, 1);
          }
      }
  }

  //---------------------------------------------------------------------------
  // Copy constructor.

  ViewIndexer(const This_t &model)
  : domain_m(model.domain()), baseDomain_m(model.baseDomain()),
    mask_m(model.mask())
  {
    for (int d = 0; d < Dim; d++)
      {
        ind_m[d]    = model.indirection(d);
        offset_m[d] = model.offset(d);
        stride_m[d] = model.stride(d);
      }
  }

  //---------------------------------------------------------------------------
  // Assignment operator.

  This_t &operator=(const This_t &rhs)
  {
    domain_m = rhs.domain();
    baseDomain_m = rhs.baseDomain();
    mask_m = rhs.mask();
    
    for (int d = 0; d < Dim; d++)
      {
        ind_m[d]    = rhs.indirection(d);
        offset_m[d] = rhs.offset(d);
        stride_m[d] = rhs.stride(d);
      }
    
    return *this;
  }
        
  //---------------------------------------------------------------------------
  // Accessors.

  const Domain_t &domain() const { return domain_m; }
  const BaseDomain_t &baseDomain() const { return baseDomain_m; }
  
  int indirection(int i) const { return ind_m[i]; }
  const Mask_t &mask() const { return mask_m; }
  int offset(int i) const { return offset_m[i]; }
  int stride(int i) const { return stride_m[i]; }

  //---------------------------------------------------------------------------
  // Translate the indices from the local coordinates to base coordinates.

  void translate(const Loc<Dim> &loc, Loc<Dim2> &oloc) const
  {
    oloc = mask_m;
  
    for (int d = 0; d < Dim; d++)
      oloc[ind_m[d]] = Loc<1>(offset_m[d] + stride_m[d] * loc[d].first());
  }

  void translate(int i0, Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
  }
  
  void translate(int i0, int i1, Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
    loc[ind_m[1]] = offset_m[1] + stride_m[1] * i1;
  }  
  
  void translate(int i0, int i1, int i2, 
    Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
    loc[ind_m[1]] = offset_m[1] + stride_m[1] * i1;
    loc[ind_m[2]] = offset_m[2] + stride_m[2] * i2;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
    loc[ind_m[1]] = offset_m[1] + stride_m[1] * i1;
    loc[ind_m[2]] = offset_m[2] + stride_m[2] * i2;
    loc[ind_m[3]] = offset_m[3] + stride_m[3] * i3;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    int i4, Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
    loc[ind_m[1]] = offset_m[1] + stride_m[1] * i1;
    loc[ind_m[2]] = offset_m[2] + stride_m[2] * i2;
    loc[ind_m[3]] = offset_m[3] + stride_m[3] * i3;
    loc[ind_m[4]] = offset_m[4] + stride_m[4] * i4;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    int i4, int i5, Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
    loc[ind_m[1]] = offset_m[1] + stride_m[1] * i1;
    loc[ind_m[2]] = offset_m[2] + stride_m[2] * i2;
    loc[ind_m[3]] = offset_m[3] + stride_m[3] * i3;
    loc[ind_m[4]] = offset_m[4] + stride_m[4] * i4;
    loc[ind_m[5]] = offset_m[5] + stride_m[5] * i5;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    int i4, int i5, int i6, Loc<Dim2> &loc) const
  {
    // Assign our masked loc, with the sliced out dimensions already filled in.
  
    loc = mask_m;
  
    // Compute the indices in the domain of the original view.
    // Use our indirection array to slot the indices into the current loc.
  
    loc[ind_m[0]] = offset_m[0] + stride_m[0] * i0;
    loc[ind_m[1]] = offset_m[1] + stride_m[1] * i1;
    loc[ind_m[2]] = offset_m[2] + stride_m[2] * i2;
    loc[ind_m[3]] = offset_m[3] + stride_m[3] * i3;
    loc[ind_m[4]] = offset_m[4] + stride_m[4] * i4;
    loc[ind_m[5]] = offset_m[5] + stride_m[5] * i5;
    loc[ind_m[6]] = offset_m[6] + stride_m[6] * i6;
  }  

  //---------------------------------------------------------------------------
  // Utility functions

  // Transforms a domain in the current coordinate system to one in base
  // coordinate system.

  template<class DT>
  BaseDomain_t &localToBase(const Domain<Dim, DT> &dlocal, 
    BaseDomain_t &base) const 
  {
    // The base domain contains the appropriate information for the
    //  sliced out dimensions.

    base = baseDomain_m;
  
    // We just need to fill in the non-sliced dimensions and transform
    // back to base coordinates.

    const typename DT::Domain_t &local = dlocal.unwrap();
    for (int d = 0; d < Dim; d++)
      {
        base[ind_m[d]] = Range<1>(
          offset_m[d] + stride_m[d] * local[d].first(),
          offset_m[d] + stride_m[d] * local[d].last(), 
          stride_m[d] * local[d].stride());
      }
    
    return base;
  }

  // Transforms a domain in the current coordinate system to one in base
  // coordinate system, returning a slice-range suitable for making a view.

  template<class DT>
  SliceRange<Dim2, Dim> &localToBase(const Domain<Dim, DT> &dlocal, 
    SliceRange<Dim2, Dim> &base) const 
  {
    // The base domain contains the appropriate information for the
    //  sliced out dimensions.

    base.totalDomain() = baseDomain_m;
  
    // We need to transform to base coordinates and fill in the appropriate
    // slots in the total domain and slice domain. We also need to label
    // the non-sliced dimensions as non-ignorable.

    const typename DT::Domain_t &local = dlocal.unwrap();
    for (int d = 0; d < Dim; d++)
      {
        Range<1> r(
          offset_m[d] + stride_m[d] * local[d].first(),
          offset_m[d] + stride_m[d] * local[d].last(), 
          stride_m[d] * local[d].stride());
        base.totalDomain()[ind_m[d]] = r;
        base.sliceDomain()[d] = r;
        base.cantIgnoreDomain(ind_m[d]);
      }
    
    return base;
  }

  // Transforms a domain in the base coordinate system to one in the
  // local coordinate system.

  Interval<Dim>  &baseToLocal(const BaseDomain_t &base, 
    Interval<Dim> &local) const
  {
    // We just need to strip out the sliced dimensions and transform back to
    // local coordinates. We check to make sure that stride really works out
    // to one.

    int j;
    for (int d = 0; d < Dim; d++)
      {
        j = ind_m[d];
        local[d] = Interval<1>(
          (base[j].first() - offset_m[d]) / stride_m[d],
          (base[j].last() - offset_m[d]) / stride_m[d]);
        PAssert(base[j].stride() / stride_m[d] == 1);
      }
    
    return local;
  }

  Range<Dim> &baseToLocal(const BaseDomain_t &base, 
    Range<Dim> &local) const
  {
    // We just need to strip out the sliced dimensions and transform back to
    // local coordinates.

    int j;
    for (int d = 0; d < Dim; d++)
      {
        j = ind_m[d];
        local[d] = Range<1>(
          (base[j].first() - offset_m[d]) / stride_m[d],
          (base[j].last() - offset_m[d]) / stride_m[d],
          base[j].stride() / stride_m[d]);
      }
    
    return local;
  }
  
  Interval<Dim>  &baseToLocalInterval(const Interval<Dim2> &base, 
				      Interval<Dim> &local) const
  {
    // We just need to strip out the sliced dimensions and transform back to
    // local coordinates.

    int j;
    for (int d = 0; d < Dim; d++)
    {
      j = ind_m[d];
      local[d] = Interval<1>((base[j].first() - offset_m[d]) / stride_m[d],
			     (base[j].last() - offset_m[d]) / stride_m[d]);

      PAssert(local[d].first() * stride_m[d] + offset_m[d] == base[j].first());
      PAssert(local[d].last() * stride_m[d] + offset_m[d] == base[j].last());
    }
    
    return local;
  }

private:

  // The current domain.

  Domain_t domain_m;
  
  // The base domain.

  BaseDomain_t baseDomain_m;
  
  // Strides and offsets.

  int stride_m[Dim], offset_m[Dim];
  
  // Mask loc and indirection vector.

  int ind_m[Dim];
  Mask_t mask_m;
};


/**
 * This is an extra-special version of View indexer that optimizes indexing for
 * the case where we have not taken a slice.
 *
 * Exported typedefs:
 *  - BaseDomain_t - the type of base domain for this ViewIndexer.
 *  - Domain_t - the type of domain for this ViewIndexer.
 *  - Mask_t - the type of mask for this ViewIndexer.
 *  - This_t - this class.
 *
 * Constructors:
 *  - ViewIndexer() - default initialization creates empty domains and leaves
 *     the rest uninitialized. Initialization can be completed with
 *     operator=().
 *  - ViewIndexer(const Domain<Dim, DT> &) - constructs a ViewIndexer for
 *     the first non-slice view of a domain.
 *  - ViewIndexer(const ViewIndexer<Dim, Dim> &, const Domain<Dim, DT> &) -
 *     constructs a new ViewIndexer from another one by taking a non-slice
 *     view.
 *  - ViewIndexer(const This_t &) - copy constructor.
 *
 * Operators:
 *  - operator=(const This_t &) - copy assignment operator.
 *
 * Accessors:
 *  - baseDomain() - returns the base domain.
 *  - domain() - returns the domain.
 *  - indirection(int i ) - returns the ith indirection index.
 *  - mask() - returns the mask.
 *  - offset(int i ) - returns the ith offset.
 *  - stride(int i ) - returns the ith stride.
 *  - translate(int i0[, ... int i6], Loc<Dim>) - translates 1-7 indices to
 *     a Loc in the base coordinate system.
 *
 * Utility functions:
 *   - localToBase - transforms a domain in the local coordinate system to
 *      the corresponding domain in base coordinates.
 *   - baseToLocal - transforms a domain in the base coordinate system to
 *      the corresponding domain in local coordinates.
 *   - baseToLocalInterval - transforms an interval in the base coordinate
 *      system to an interval in the local coordinates while ingnoring the
 *      fact that the stride doesn't work.  (The endpoints are transformed.)
 */

template<int Dim>
class ViewIndexer<Dim, Dim>
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.

  // This class.
  
  typedef ViewIndexer<Dim, Dim> This_t;
  
  // The domain type.
  
  typedef Interval<Dim> Domain_t;
  
  // The base domain type.
  
  typedef Range<Dim> BaseDomain_t;
  
  // The mask type.
  
  typedef Loc<Dim> Mask_t;

  //---------------------------------------------------------------------------
  // Constructors.
  
  ViewIndexer() { }
    
  template<class DT>
  ViewIndexer(const Domain<Dim, DT> &dom)
  : domain_m(Pooma::NoInit()), baseDomain_m(dom.unwrap())
  {
    // Set up our offsets, domains, and strides. This is easy since we're not
    // a slice. The mask is initialized to according to the default constructor
    // and we don't need to do anything. The indirection vector is a one-to-one
    // linear mapping.
    
    const typename DT::Domain_t &domain = dom.unwrap();
    for (int d = 0; d < Dim; ++d)
      {
        offset_m[d] = domain[d].first();
        stride_m[d] = domain[d].stride();
        domain_m[d] = Interval<1>(domain[d].length());
      }
  }

  // We can take a non-slice subset of an existing non-sliced view.
  
  template<class DT>
  ViewIndexer(const ViewIndexer<Dim, Dim> &orig, 
    const Domain<Dim, DT> &dom)
  : domain_m(Pooma::NoInit()), 
    baseDomain_m(dom.unwrap()), 
    mask_m(orig.mask())
  {
    // Set up our offsets, domains, and strides. This is easy since we're not
    // not being sliced by this domain meaning we don't need to fiddle with
    // the mask and indirection vectors.

    const typename DT::Domain_t &domain = dom.unwrap();
    for (int d = 0; d < Dim; ++d)
      {
        offset_m[d] = orig.offset(d) + orig.stride(d) * 
          domain[d].first();
        stride_m[d] = orig.stride(d) * domain[d].stride();
        domain_m[d] = Interval<1>(domain[d].length());
      }
      
    localToBase(domain_m, baseDomain_m);
}

  //---------------------------------------------------------------------------
  // Copy constructor.

  ViewIndexer(const This_t &model)
  : domain_m(model.domain()), baseDomain_m(model.baseDomain()),
    mask_m(model.mask())
  {
    for (int d = 0; d < Dim; d++)
      {
        offset_m[d] = model.offset(d);
        stride_m[d] = model.stride(d);
      }
  }

  //---------------------------------------------------------------------------
  // Assignment operator.

  This_t &operator=(const This_t &rhs)
  {
    domain_m = rhs.domain();
    baseDomain_m = rhs.baseDomain();
    mask_m = rhs.mask();
    
    for (int d = 0; d < Dim; d++)
      {
        offset_m[d] = rhs.offset(d);
        stride_m[d] = rhs.stride(d);
      }
    
    return *this;
  }
  
  //---------------------------------------------------------------------------
  // Accessors.

  const Domain_t &domain() const { return domain_m; }
  const BaseDomain_t &baseDomain() const { return baseDomain_m; }
  
  int indirection(int i) const { return i; }
  const Mask_t &mask() const { return mask_m; }
  int offset(int i) const { return offset_m[i]; }
  int stride(int i) const { return stride_m[i]; }

  //---------------------------------------------------------------------------
  // Translate the indices from the local coordinates to base coordinates.

  void translate(const Loc<Dim> &loc, Loc<Dim> &oloc) const
  {
    for (int d = 0; d < Dim; d++)
      oloc[d] = Loc<1>(offset_m[d] + stride_m[d] * loc[d].first());
  }

  void translate(int i0, Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
  }
  
  void translate(int i0, int i1, Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
    loc[1] = offset_m[1] + stride_m[1] * i1;
  }  
  
  void translate(int i0, int i1, int i2, 
    Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
    loc[1] = offset_m[1] + stride_m[1] * i1;
    loc[2] = offset_m[2] + stride_m[2] * i2;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
    loc[1] = offset_m[1] + stride_m[1] * i1;
    loc[2] = offset_m[2] + stride_m[2] * i2;
    loc[3] = offset_m[3] + stride_m[3] * i3;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    int i4, Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
    loc[1] = offset_m[1] + stride_m[1] * i1;
    loc[2] = offset_m[2] + stride_m[2] * i2;
    loc[3] = offset_m[3] + stride_m[3] * i3;
    loc[4] = offset_m[4] + stride_m[4] * i4;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    int i4, int i5, Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
    loc[1] = offset_m[1] + stride_m[1] * i1;
    loc[2] = offset_m[2] + stride_m[2] * i2;
    loc[3] = offset_m[3] + stride_m[3] * i3;
    loc[4] = offset_m[4] + stride_m[4] * i4;
    loc[5] = offset_m[5] + stride_m[5] * i5;
  }  

  void translate(int i0, int i1, int i2, int i3, 
    int i4, int i5, int i6, Loc<Dim> &loc) const
  {
    // Compute the indices in the domain of the original view.
  
    loc[0] = offset_m[0] + stride_m[0] * i0;
    loc[1] = offset_m[1] + stride_m[1] * i1;
    loc[2] = offset_m[2] + stride_m[2] * i2;
    loc[3] = offset_m[3] + stride_m[3] * i3;
    loc[4] = offset_m[4] + stride_m[4] * i4;
    loc[5] = offset_m[5] + stride_m[5] * i5;
    loc[6] = offset_m[6] + stride_m[6] * i6;
  }  

  //---------------------------------------------------------------------------
  // Utility functions

  // Transforms a domain in the current coordinate system to one in base
  // coordinate system.

  template<class DT>
  BaseDomain_t &localToBase(const Domain<Dim, DT> &dlocal, 
    BaseDomain_t &base) const 
  {  
    // We just transform back to base coordinates.

    const typename DT::Domain_t &local = dlocal.unwrap();
    for (int d = 0; d < Dim; d++)
      {
        base[d] = Range<1>(
          offset_m[d] + stride_m[d] * local[d].first(),
          offset_m[d] + stride_m[d] * local[d].last(), 
          stride_m[d] * local[d].stride());
      }
    
    return base;
  }

  // Transforms a domain in the current coordinate system to one in base
  // coordinate system, returning a slice-range suitable for making a view.

  // Do we even need this version? SliceRange<D,D> is never used, is it???
  
  template<class DT>
  SliceRange<Dim, Dim> &localToBase(const Domain<Dim, DT> &dlocal, 
    SliceRange<Dim, Dim> &base) const 
  {
    // The base domain contains the appropriate information for the
    //  sliced out dimensions.

    base.totalDomain() = baseDomain_m;
  
    // We need to transform to base coordinates and fill in the appropriate
    // slots in the total domain and slice domain. We also need to label
    // the non-sliced dimensions as non-ignorable.

    const typename DT::Domain_t &local = dlocal.unwrap();
    for (int d = 0; d < Dim; d++)
      {
        Range<1> r(
          offset_m[d] + stride_m[d] * local[d].first(),
          offset_m[d] + stride_m[d] * local[d].last(), 
          stride_m[d] * local[d].stride());
        base.totalDomain()[d] = r;
        base.sliceDomain()[d] = r;
        base.cantIgnoreDomain(d);
      }
    
    return base;
  }

  // Transforms a domain in the base coordinate system to one in the
  // local coordinate system.

  Interval<Dim> &baseToLocal(const BaseDomain_t &base, 
    Interval<Dim> &local) const
  {
    // We just transform back to local coordinates.
    // We check to make sure that stride really works out
    // to one.

    for (int d = 0; d < Dim; d++)
      {
        local[d] = Interval<1>(
          (base[d].first() - offset_m[d]) / stride_m[d],
          (base[d].last() - offset_m[d]) / stride_m[d]);
        PAssert(base[d].stride() / stride_m[d] == 1);
      }
    
    return local;
  }

  Range<Dim> &baseToLocal(const BaseDomain_t &base, 
    Range<Dim> &local) const
  {
    // We just transform back to local coordinates.

    for (int d = 0; d < Dim; d++)
      {
        local[d] = Range<1>(
          (base[d].first() - offset_m[d]) / stride_m[d],
          (base[d].last() - offset_m[d]) / stride_m[d],
          base[d].stride() / stride_m[d]);
      }
    
    return local;
  }
  
  // This version ignores strides.

  Interval<Dim> &baseToLocalInterval(const Interval<Dim> &base, 
				     Interval<Dim> &local) const
  {
    // We just transform back to local coordinates.

    for (int d = 0; d < Dim; d++)
    {
      local[d] = Interval<1>((base[d].first() - offset_m[d]) / stride_m[d],
			     (base[d].last() - offset_m[d]) / stride_m[d]);

      PAssert(local[d].first() * stride_m[d] + offset_m[d] == base[d].first());
      PAssert(local[d].last() * stride_m[d] + offset_m[d] == base[d].last());
    }
    
    return local;
  }

private:

  // The current domain.

  Domain_t domain_m;
  
  // The base domain.

  BaseDomain_t baseDomain_m;
  
  // Strides and offsets.

  int stride_m[Dim], offset_m[Dim];
  
  // Mask loc and indirection vector.

  Mask_t mask_m;
};

#endif // POOMA_UTILITIES_VIEWINDEXER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ViewIndexer.h,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo

