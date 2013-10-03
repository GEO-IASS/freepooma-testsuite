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
 * A tag for telling terms they are about to be read.
 */

#ifndef POOMA_UTILITIES_NOTIFYPREREAD_H
#define POOMA_UTILITIES_NOTIFYPREREAD_H

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class T, class A> struct LeafFunctor;
template<class T> class Scalar;

//-----------------------------------------------------------------------------
// NotifyPreReadTag is a simple tag class used to indicate we're getting
// read to read a term.
//-----------------------------------------------------------------------------

struct NotifyPreReadTag { };

//----------------------------------------------------------------------
// Scalars don't need to do anything in order to be read.
//----------------------------------------------------------------------

template<class T>
struct LeafFunctor<Scalar<T>, NotifyPreReadTag>
{
  typedef bool Type_t;
  static Type_t apply(const Scalar<T> &, const NotifyPreReadTag &)
  {
    return true;
  }
};

#endif // POOMA_UTILITIES_NOTIFYPREREAD_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NotifyPreRead.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
