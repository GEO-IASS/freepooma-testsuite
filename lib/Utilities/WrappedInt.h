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
 * A tag class templated on an integer.  This class is intended to be
 * used to let you specialize a function on a compile time number.
 */

#ifndef POOMA_UTILITIES_WRAPPEDINT_H
#define POOMA_UTILITIES_WRAPPEDINT_H

//-----------------------------------------------------------------------------
// Class: WrappedInt
//-----------------------------------------------------------------------------

#include "Utilities/PurifyConstructors.h"


/** Helper class: WrappedInt<int> 
 *
 * A tag class templated on an integer.  This class is intended to be
 * used to let you specialize a function on a compile time number.
 *
 * For example, if you have an object of type T which you want to pass
 * to a subroutine foo, but you want to specialize that subroutine based on
 * the enum 'bar' defined in T, you could say:
 *
 * template<class T>  void foo(const T& t) 
 * { 
 *    foo(t,WrappedInt<T::bar>())
 * }
 *
 * You can then specialize foo(T,WrappedInt<x>) for different values
 * of x.
 *
 * With functor classes you can do this in a slightly slicker manner.
 * Define the general functor class like so:
 *
 * template< class T, class Switch = WrappedInt<T::bar> >
 * struct Foo
 * {
 *   static void foo()(const T& t);
 * }
 *
 * Then you can specialize foo for different values of Switch.  That
 * is you write things like Foo<T,WrappedInt<1>>::foo(const T&) to
 * specialize it for T::bar==1. 
 *
 * With this construction you don't have two function calls and the
 * WrappedInt is never constructed or passed.  Since it relies on Foo
 * being templated though, you can't use it nested in another class.
 */

template<int Integer> class WrappedInt
{
public:
  enum { val = Integer };
  POOMA_PURIFY_CONSTRUCTORS(WrappedInt)
};

#endif   // POOMA_UTILITIES_WRAPPEDINT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: WrappedInt.h,v $   $Author: richard $
// $Revision: 1.13 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
