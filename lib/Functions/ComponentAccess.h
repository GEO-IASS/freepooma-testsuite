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
// ComponentAccess
//-----------------------------------------------------------------------------

#ifndef POOMA_FUNCTIONS_COMPONENT_ACCESS_H
#define POOMA_FUNCTIONS_COMPONENT_ACCESS_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Utilities
 * @brief
 * ComponentAccess<Container, Comp> is a general functor class that users can
 * specialize to tell POOMA how to access components inside an object used as
 * an element in expressions.  (For example, Vector, Tensor etc.)
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Wrapper class to encapsulate a component specification.
 */

template<class Components>
class ComponentWrapper
{
public:
  explicit ComponentWrapper(const Components &c) : c_m(c) { }

  const Components &components() const { return c_m; }

private:
  const Components &c_m;
};


/**
 * ComponentAccess is an interface class that is used to provide an API for
 * accessing components of a composite type. The default version handles
 * scalars.
 */

template<class T, class Components>
struct ComponentAccess
{
  typedef T Element_t;
  typedef T &ElementRef_t;
  
  static inline ElementRef_t indexRef(T &v, const Components &)
  {
    return v;
  }
  
  static inline Element_t index(const T &v, const Components &)
  {
    return v;
  }
};

//-----------------------------------------------------------------------------
//
// POOMA_COMPONENT_ACCESS macro:
//
// This macro simplifies writing component accessors for user defined structs
// where you can safely just return references to the components.  For example:
//
// struct Bob
// {
//   double density;
//   Vector<2> velocity;
// }
//
// POOMA_COMPONENT_ACCESS(Bob,Density,double,density)
// POOMA_COMPONENT_ACCESS(Bob,Velocity,Vector<2>,velocity)
//
// Array<1,Bob> a;
//
// a.comp(Density()) = ...
// a.comp(Velocity()) = ...
//
// The four parameters to POOMA_COMPONENT_ACCESS are:
//   1) name of the user's struct or class
//   2) name of a tag class (which the macro defines for you)
//   3) type of the component
//   4) access method for the component.
//-----------------------------------------------------------------------------

#define POOMA_COMPONENT_ACCESS(IN, TAG, TYPE, MEMBER)                   \
                                                                        \
struct TAG                                                              \
{                                                                       \
  TAG() { }                                                             \
  TAG(const TAG &) { }                                                  \
};                                                                      \
                                                                        \
template<>                                                              \
struct ComponentAccess<IN, TAG>                                         \
{                                                                       \
  typedef TYPE Element_t;                                               \
  typedef TYPE &ElementRef_t;                                           \
                                                                        \
  static inline ElementRef_t indexRef(IN &in, const TAG &)              \
  {                                                                     \
    return in.MEMBER;                                                   \
  }                                                                     \
                                                                        \
  static inline Element_t index(const IN &in, const TAG &)              \
  {                                                                     \
    return in.MEMBER;                                                   \
  }                                                                     \
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_FUNCTIONS_COMPONENT_ACCESS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ComponentAccess.h,v $   $Author: richi $
// $Revision: 1.7 $   $Date: 2004/11/04 14:15:58 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
