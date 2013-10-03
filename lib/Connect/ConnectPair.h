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

#ifndef POOMA_CONNECT_CONNECT_PAIR_H
#define POOMA_CONNECT_CONNECT_PAIR_H

//-----------------------------------------------------------------------------
// Classes:
// ConnectPair<T1,T2>
//-----------------------------------------------------------------------------

/** @file
 * @group Connect
 * @brief
 * ConnectPair<T1,T2> stores two items of types T1 and T2, basically the
 * same as std::pair from the STL.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * ConnectPair<T1,T2> stores two items of types T1 and T2, basically the
 * same as std::pair from the STL.
 *
 * The items are stored as copies.  The
 * public interface is either methods first() and second(), or the two
 * public data objects first_m and second_m.  This class is used for
 * storing two items in a single object so that Connector can be specialized
 * on the pair.
 */

template<class T1, class T2>
class ConnectPair
{
public:
  //============================================================
  // ConnectPair Constructor
  //============================================================

  /// Initialize from the two objects, and store copies

  ConnectPair(const T1 &a, const T2 &b)
    : first_m(a), second_m(b)
  {
  }

  /// Destructor does nothing special

  ~ConnectPair()
  {
  }

  //============================================================
  // ConnectPair accessors
  //============================================================

  ///@name Return references to the two elements, first and second
  //@{

  T1 &first()              { return first_m; }
  const T1 &first() const  { return first_m; }

  T2 &second()             { return second_m; }
  const T2 &second() const { return second_m; }

  //@}

  //============================================================
  // ConnectPair public data
  //============================================================

  /// The two elements in the pair - we make them public to match
  /// the philosophy (but not the exact interface) of std::pair

  T1 first_m;
  T2 second_m;
};


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_CONNECT_CONNECT_PAIR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ConnectPair.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:16 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
