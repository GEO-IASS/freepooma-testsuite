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

/** @file
 * @ingroup Utilities
 * @brief
 * A tag for updating leafs in an expression.
 */

#ifndef POOMA_UTILITIES_PERFORMUPDATE_H
#define POOMA_UTILITIES_PERFORMUPDATE_H

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class T, class A> struct LeafFunctor;

/**
 * These leaf functor specializations are used to notify a field or expression
 * that it is going to be read and, therefore, needs to update itself. 
 *
 * The first LeafFunctor represents default behavior, which is to do nothing.
 *
 * Fields with engines that store internal fields AND don't copy those
 * fields' relations to its list must provide a specialization to get the 
 * internal fields to update.
 *
 * NOTE: we don't use the ExpressionApply machinery here because this really
 * operate on the engines.
 */

struct PerformUpdateTag {};

template<class Node>
struct LeafFunctor<Node, PerformUpdateTag>
{
  typedef int Type_t;

  inline static
  Type_t apply(const Node &, const PerformUpdateTag &)
    {
      return 0;
    }
};


#endif // POOMA_UTILITIES_PERFORMUPDATE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PerformUpdate.h,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
