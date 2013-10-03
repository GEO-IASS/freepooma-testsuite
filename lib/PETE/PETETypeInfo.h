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
// TypeInfo<> specializations for PETE things.
//-----------------------------------------------------------------------------

#ifndef POOMA_PETE_PETETYPEINFO_H
#define POOMA_PETE_PETETYPEINFO_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class T> struct TypeInfo;

struct OpAdd;
struct OpMultiply;
template<class Expr> struct Reference;
template<class Op, class Child> struct UnaryNode;
template<class Op, class Left, class Right> struct BinaryNode;
template<class T> struct Scalar;

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Compile-time TypeInfo for PETE things.
//-----------------------------------------------------------------------------

template<>
struct TypeInfo<OpAdd>
{
  static inline std::string name() { return "OpAdd"; }
};

template<>
struct TypeInfo<OpMultiply>
{
  static inline std::string name() { return "OpMultiply"; }
};

template<class T>
struct TypeInfo<Scalar<T> >
{
  static inline std::string name()
  {
    return "Scalar<" + TypeInfo<T>::name() + " >";
  }
};

template<class Expr>
struct TypeInfo<Reference<Expr> >
{
  static inline std::string name()
  {
    return "Reference<" + TypeInfo<Expr>::name() + " >";
  }
};

template<class Op, class Child>
struct TypeInfo<UnaryNode<Op, Child> >
{
  static inline std::string name()
  {
    return "UnaryNode<" + TypeInfo<Op>::name() + ","
      + TypeInfo<Child>::name()
      + " >";
  }
};

template<class Op, class Left, class Right>
struct TypeInfo<BinaryNode<Op, Left, Right> >
{
  static inline std::string name()
  {
    return "BinaryNode<" + TypeInfo<Op>::name() + ","
      + TypeInfo<Left>::name() + ","
      + TypeInfo<Right>::name() + " >";
  }
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_PETE_PETETYPEINFO_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PETETypeInfo.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:56 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
