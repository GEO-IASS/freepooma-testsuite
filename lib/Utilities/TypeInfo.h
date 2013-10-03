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
// TypeInfo
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * Undocumented.
 */

#ifndef POOMA_UTILITIES_TYPEINFO_H
#define POOMA_UTILITIES_TYPEINFO_H

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

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<class T>
struct TypeInfo
{
};

template<>
struct TypeInfo<double>
{
  static inline std::string name() { return "double"; }
};

template<>
struct TypeInfo<int>
{
  static inline std::string name() { return "int"; }
};

template<int D>
struct TypeInfoInt
{
};

template<>
struct TypeInfoInt<1>
{
  static inline std::string name() { return "1"; }
};

template<>
struct TypeInfoInt<2>
{
  static inline std::string name() { return "2"; }
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_UTILITIES_TYPEINFO_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TypeInfo.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
