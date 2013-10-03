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
// Structs:
//   KernelTag1<Expr>
//   KernelTag<LHS,RHS>
//   CompressibleKernel<bool,bool>
// Tags:
//   InlineKernelTag
//   CompressibleViewKernelTag
//   CompressibleKernelTag
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_KERNELTAGS_H
#define POOMA_EVALUATOR_KERNELTAGS_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * Kernel Tags are used for picking the appropriate Kernel given the
 * engines in an expression.  Each Kernel tag represents a set of engines
 * that it is capable of dealing with.
 *
 * The external interface for KernelTags is the traits
 * - KernelTag1<Expr>::Tag_t
 * - KernelTag<LHSArray,RHSArray>::Kernel_t:
 *        this is the main function that produces an Kernel tag for the
 *        expression
 *
 * which yields an evaluator tag, given the Array type or the types for the
 * left hand side and right hand side.  To add new Kernels or new
 * engines, specialize the KernelForEngine struct to give the Kernel
 * tag for each engine, and specialize the NewKernel and KernelCombine
 * struct to determine how to chose a new Kernel given two Kernels.
 *
 * The implementation of this interface will probably change when other
 * kernels are added.  Currently we only have the basic inline kernel plus two
 * that deal with compression, so we pick the kernel tag by querying Compressible
 * about the left and right hand sides.
 *
 * Currently there are only three Kernels:
 * - InlineKernelTag: use the inline Kernel (simple loops, no patches)
 * - CompressibleViewKernelTag: for a compressible lhs, 
 *   takes a brickview of lhs then loops
 * - CompressibleKernelTag: checks if both sides are compressed to
 *   do compressed assign otherwise calls CVE
 *
 * The results for expressions with Bricks (B) and CompressibleBricks (C)
 * are:
 * - B = B+B;   InlineKernelTag
 * - B = C+B;   InlineKernelTag
 * - B = C+C;   InlineKernelTag
 * - C = B+B;   CompressibleViewKernelTag
 * - C = C+B;   CompressibleViewKernelTag
 * - C = C+C;   CompressibleKernelTag
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Evaluator/CompressibleEngines.h"
#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class Expr> struct CreateLeaf;


//-----------------------------------------------------------------------------
// Some Kernel tags (more tags can be defined elsewhere)
//-----------------------------------------------------------------------------

struct ErrorKernelTag 
{ 
  ErrorKernelTag(){}
  ~ErrorKernelTag(){}
};

struct InlineKernelTag 
{ 
  InlineKernelTag(){}
  ~InlineKernelTag(){}
};

struct CompressibleKernelTag 
{ 
  CompressibleKernelTag(){}
  ~CompressibleKernelTag(){}
};

struct CompressibleViewKernelTag 
{ 
  CompressibleViewKernelTag(){}
  ~CompressibleViewKernelTag(){}
};

//-----------------------------------------------------------------------------
// CompressibleKernel<bool,bool>
//
// Pick the appropriate kernel based on the compressibility of the left and
// right hand sides.
//-----------------------------------------------------------------------------

template<bool lhsComp,bool rhsComp>
struct CompressibleKernel
{
  CompressibleKernel(){}
  ~CompressibleKernel(){}
};

template<>
struct CompressibleKernel<false,false>
{
  CompressibleKernel(){}
  ~CompressibleKernel(){}
  typedef InlineKernelTag Kernel_t;
};

template<>
struct CompressibleKernel<false,true>
{
  CompressibleKernel(){}
  ~CompressibleKernel(){}
  typedef InlineKernelTag Kernel_t;
};

template<>
struct CompressibleKernel<true,false>
{
  CompressibleKernel(){}
  ~CompressibleKernel(){}
  typedef CompressibleViewKernelTag Kernel_t;
};

template<>
struct CompressibleKernel<true,true>
{
  CompressibleKernel(){}
  ~CompressibleKernel(){}
  typedef CompressibleKernelTag Kernel_t;
};


//-----------------------------------------------------------------------------
// KernelTag<LHS,RHS>
//
// Finally, this struct computes the Kernel tag for the whole expression
// given the types of the arrays on the left and right hand sides.
//-----------------------------------------------------------------------------

template<class Expr>
struct KernelTag1
{
  KernelTag1(){}
 ~KernelTag1(){}
  typedef typename Expr::Engine_t ExprEngine_t;
  typedef typename EngineFunctor<ExprEngine_t, Compressible>::Type_t Expr_t;
  enum { exprComp = Expr_t::val };
  typedef typename CompressibleKernel<exprComp,exprComp>::Kernel_t Kernel_t;
};


//-----------------------------------------------------------------------------
// KernelTag<LHS,RHS>
//
// Finally, this struct computes the Kernel tag for the whole expression
// given the types of the arrays on the left and right hand sides.
//-----------------------------------------------------------------------------

template<class LHS,class RHS>
struct KernelTag
{
  KernelTag(){}
 ~KernelTag(){}
  typedef typename LHS::Engine_t LHSEngine_t;
  typedef typename RHS::Engine_t RHSEngine_t;
  typedef typename EngineFunctor<LHSEngine_t,Compressible>::Type_t LHST_t;
  typedef typename EngineFunctor<RHSEngine_t,Compressible>::Type_t RHST_t;
  enum { lhsComp = LHST_t::val };
  enum { rhsComp = RHST_t::val };
  typedef typename CompressibleKernel<lhsComp,rhsComp>::Kernel_t Kernel_t;
};


#endif     // POOMA_EVALUATOR_KERNELTAGS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: KernelTags.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
