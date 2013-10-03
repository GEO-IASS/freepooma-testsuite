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
// Class:
// ScalarCodeInfo
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_SCALARCODEINFO_H
#define POOMA_EVALUATOR_SCALARCODEINFO_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * ScalarCodeInfo contains all the information necessary for evaluating a
 * piece of scalar code on several arguments.
 *
 *  You need to say which fields
 * are being written to, which ones use their guard layers, and the largest
 * extents of the computation (how far you go to look at neighbors).
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Layout/INode.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

class ScalarCodeInfo
{
public:
  typedef std::vector<int>  Extents_t;
  typedef std::vector<bool> BoolVector_t;

  ScalarCodeInfo()
  {
    arguments_m = 0;
    dimensions_m = 0;
  }

  /// The number of arguments of the ScalarCode functor. This needs to be
  /// called before any of write() and useGuards().
  void arguments(int n)
  {
    PAssert(n > 0);

    arguments_m = n;
    writers_m.resize(n);
    readers_m.resize(n);
    useGuards_m.resize(n);

    int i;
    for (i = 0; i < n; ++i)
    {
      writers_m[i]   = false;
      readers_m[i]   = true;
      useGuards_m[i] = true;
    }
    writers_m[0] = true;
    readers_m[0] = false;
  }

  /// The number of dimensions the arguments to the ScalarCode functor
  /// span. This needs to be specified before any of the lowerExtent()
  /// and upperExtent().
  void dimensions(int n)
  {
    PAssert(n > 0);

    dimensions_m = n;
    lower_m.resize(n);
    upper_m.resize(n);

    int i;
    for (i = 0; i < n; ++i)
    {
      lower_m[i] = 0;
      upper_m[i] = 0;
    }
  }

  /// Lower extent for dimension i. This specifies the (maximum) stencil size.
  /// Note that you need to make sure a view extending the physical domain by
  /// this amount can be taken of every argument to the functor.

  int &lowerExtent(int i)
  {
    return lower_m[i];
  }

  /// Upper extent for dimension i. This specifies the (maximum) stencil size.
  /// Note that you need to make sure a view extending the physical domain by
  /// this amount can be taken of every argument to the functor.

  int &upperExtent(int i)
  {
    return upper_m[i];
  }
 
  /// Specify whether argument i is written to, writing allows reading. This
  /// triggers notifying the engine after writing and possibly dirtying
  /// relations attached to a written field.
 
  void write(int i, bool f)
  {
    writers_m[i] = f;
    readers_m[i] = (!f);
  }

  BoolVector_t &writers()
  {
    return writers_m;
  }

  BoolVector_t &readers()
  {
    return readers_m;
  }

  void useGuards(int i, bool f)
  {
    useGuards_m[i] = f;
  }

  BoolVector_t &useGuards()
  {
    return useGuards_m;
  }

  /// The domain we take the view of before passing it to the functors
  /// operator() method.
  template<int D>
  inline Interval<D> extendDomain(const Interval<D> &domain)
  {
    Interval<D> ret;
    int d;
    for (d = 0; d < D; ++d)
    {
      ret[d] =
	        Interval<1>(
    		    domain[d].first() - lower_m[d],
    		    domain[d].last() + upper_m[d]
		    );
    }
    return ret;
  }

  /// The domain evaluation takes place on the viewed (zero-based!) engine.
  template<int D>
  inline Interval<D> evaluationDomain(const Interval<D> &domain)
  {
    Interval<D> ret;
    int d;
    for (d = 0; d < D; ++d)
    {
      ret[d] =
	Interval<1>(
		    lower_m[d],
		    domain[d].last() - domain[d].first()  + lower_m[d]
		    );
    }
    return ret;
  }

  template<int D>
  inline INode<D> extendDomain(const INode<D> &inode)
  {
    return INode<D>(inode, extendDomain(inode.domain()));
  }

private:

  int       arguments_m;
  int       dimensions_m;
  Extents_t upper_m;
  Extents_t lower_m;
  BoolVector_t useGuards_m;
  BoolVector_t writers_m;
  BoolVector_t readers_m;
};



//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_SCALARCODEINFO_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ScalarCodeInfo.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
