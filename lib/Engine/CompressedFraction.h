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
// functions:
//   compressed()
//   compressedFraction()
//   elementsCompressed()
//   compress()
//   uncompress()
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_COMPRESSEDFRACTION_H
#define POOMA_ENGINE_COMPRESSEDFRACTION_H

/** @file
 * @ingroup Engine
 * @brief
 * External functions that can be applied to any engine or container to answer
 * questions about compression and to force compression or uncompression.
 */

//-----------------------------------------------------------------------------
//
// bool compressed()
//
/// Return whether or not something is currently compressed.
//
//-----------------------------------------------------------------------------

template <class Expr>
inline bool compressed(const Expr &expr)
{
  return false;
}

//-----------------------------------------------------------------------------
//
// double compressedFraction()
//
/// Compute the fraction of the total domain that is currently compressed.
//
//-----------------------------------------------------------------------------

template <class Expr>
double compressedFraction(const Expr &expr)
{
  return static_cast<double>(elementsCompressed(expr)) / expr.domain().size();
}


//-----------------------------------------------------------------------------
//
// long elementsCompressed()
//
/// Compute the number of the elements that are currently compressed.
//
//-----------------------------------------------------------------------------

template <class Expr>
inline long elementsCompressed(const Expr &expr)
{
  return 0L;
}

//-----------------------------------------------------------------------------
//
// void compress()
//
/// (Try to) compress an engine or container.
//
//-----------------------------------------------------------------------------

template <class Expr>
inline void compress(Expr &)
{ }


//-----------------------------------------------------------------------------
//
// void uncompress()
//
/// Manually uncompress an engine or container.
//
//-----------------------------------------------------------------------------

template <class Expr>
inline void uncompress(Expr &)
{ }

#endif     // POOMA_ENGINE_COMPRESSEDFRACTION_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CompressedFraction.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
