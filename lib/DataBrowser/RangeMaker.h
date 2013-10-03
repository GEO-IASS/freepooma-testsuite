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
//   RangeMaker
//-----------------------------------------------------------------------------

#ifndef POOMA_DATABROWSER_RANGEMAKER_H
#define POOMA_DATABROWSER_RANGEMAKER_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
// 
// Classes:
//
// RangeMaker : ???
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------
#include "Domain/Range.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//
// Full Description:
//
// Classes:
//
// RangeMaker:
// 
// RangeMaker is a functor class ....
// 
//-----------------------------------------------------------------------------

// See DataBrowser.h, in further partial specializations of dbprint<>(), for
// why this RangeMaker class exists:

template<int Dim, int NumIntArgs>
class RangeMaker;

// Run through all the meaningful combinations of Dim and NumIntArgs and make
// partial specializations of RangeMaker:

// For Dim=NumIntArgs, interpret as requesting a single element:


// --------------------------------------------------------
// Dim = 1

template<>
class RangeMaker<1,1>
{
public:
  RangeMaker() { }
  Range<1> operator()(const int& i0)
  {
    return Range<1>(i0,i0,1);
  }
};

template<>
class RangeMaker<1,2>
{
public:
  RangeMaker() { }
  Range<1> operator()(const int& i0, const int& i1)
  {
    return Range<1>(i0, i1, 1);
  }
};

template<>
class RangeMaker<1,3>
{
public:
  RangeMaker() { }
  Range<1> operator()(const int& i0, const int& i1, const int& i2)
  {
    return Range<1>(Range<1>(i0, i1, i2));
  }
};


// --------------------------------------------------------
// Dim = 2

template<>
class RangeMaker<2,1>
{
public:
  RangeMaker() { }
  Range<2> operator()(const int& i0)
  {
    return Range<2>(Range<1>(i0), Range<1>(i0));
  }
};

template<>
class RangeMaker<2,2>
{
public:
  RangeMaker() { }
  Range<2> operator()(const int& i0, const int& i1)
  {
    return Range<2>(Range<1>(i0, i0, 1), Range<1>(i1, i1, 1));
  }
};

template<>
class RangeMaker<2,4>
{
public:
  RangeMaker() { }
  Range<2> operator()(const int& i0, const int& i1, const int& i2, const int& i3)
  {
    return Range<2>(Range<1>(i0, i1, 1), Range<1>(i2, i3, 1));
  }
};

template<>
class RangeMaker<2,6>
{
public:
  RangeMaker() { }
  Range<2> operator()(const int& i0, const int& i1, const int& i2, const int& i3, 
                      const int& i4, const int& i5)
  {
    return Range<2>(Range<1>(i0, i1, i2), Range<1>(i3, i4, i5));
  }
};


// --------------------------------------------------------
// Dim = 3

template<>
class RangeMaker<3,1>
{
public:
  RangeMaker() { }
  Range<3> operator()(const int& i0)
  {
    return Range<3>(Range<1>(i0), Range<1>(i0), Range<1>(i0));
  }
};

template<>
class RangeMaker<3,3>
{
public:
  RangeMaker() { }
  Range<3> operator()(const int& i0, const int& i1, const int& i2)
  {
    return Range<3>(Range<1>(i0, i0, 1), Range<1>(i1, i1, 1), 
                    Range<1>(i2, i2, 1));
  }
};

template<>
class RangeMaker<3,6>
{
public:
  RangeMaker() { }
  Range<3> operator()(const int& i0, const int& i1, const int& i2, 
                      const int& i3, const int& i4, const int& i5)
  {
    return Range<3>(Range<1>(i0, i1, 1), Range<1>(i2, i3, 1), 
                    Range<1>(i4, i5, 1));
  }
};

template<>
class RangeMaker<3,9>
{
public:
  RangeMaker() { }
  Range<3> operator()(const int& i0, const int& i1, const int& i2, 
                      const int& i3, const int& i4, const int& i5, 
                      const int& i6, const int& i7, const int& i8)
  {
    return Range<3>(Range<1>(i0, i1, i2), Range<1>(i3, i4, i5), 
                    Range<1>(i6, i7, i8));
  }
};


// --------------------------------------------------------
// Dim = 4

template<>
class RangeMaker<4,1>
{
public:
  RangeMaker() { }
  Range<4> operator()(const int& i0)
  {
    return Range<4>(Range<1>(i0), Range<1>(i0), Range<1>(i0), Range<1>(i0));
  }
};

template<>
class RangeMaker<4,4>
{
public:
  RangeMaker() { }
  Range<4> operator()(const int& i0, const int& i1, const int& i2, 
                      const int& i3)
  {
    return Range<4>(Range<1>(i0, i0, 1), Range<1>(i1, i1, 1), 
                    Range<1>(i2, i2, 1), Range<1>(i3, i3, 1));
  }
};

template<>
class RangeMaker<4,6>
{
public:
  RangeMaker() { }
  Range<4> operator()(const int& i0, const int& i1, const int& i2, 
                      const int& i3, const int& i4, const int& i5)
  {
    return Range<4>(Range<1>(i0, i1, 1), Range<1>(i2, i3, 1), 
                    Range<1>(i4, i5, 1));
  }
};

template<>
class RangeMaker<4,12>
{
public:
  RangeMaker() { }
  Range<4> operator()(const int& i0, const int& i1, const int& i2, 
                      const int& i3, const int& i4, const int& i5, 
                      const int& i6, const int& i7, const int& i8, 
                      const int& i9, const int& i10, const int& i11)
  {
    return Range<4>(Range<1>(i0, i1, i2), Range<1>(i3, i4, i5), 
                    Range<1>(i6, i7, i8), Range<1>(i9, i10, i11));
  }
};

#endif     // POOMA_DATABROWSER_RANGEMAKER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RangeMaker.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:27 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
