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

#ifndef POOMA_DOMAIN_ARITH_OPS_TRAITS_H
#define POOMA_DOMAIN_ARITH_OPS_TRAITS_H


//-----------------------------------------------------------------------------
// Class:
// DomainArithOpsTraits
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * DomainArithOpsTraits is intended to be used to select the return type of 
 * arithmetic operations between domains and pseudo-domains.
 *
 * DomainArithOpsTraits comprises the following list:
 *
 * Loc<1> Loc<N> Interval<N> Range<N> IndirectionList<int> Grid<N>
 *
 * Valid combinations are
 *  - Loc<1> + - * / Loc<1>               => Loc<1>
 *  - Loc<1> + - * / Loc<N>               => Loc<N>
 *  - Loc<1> + - * / Interval<N>          => Interval<N>
 *  - Loc<1> + - * / Range<N>             => Range<N>
 *  - Loc<1> + - * / IndirectionList<int> => IndirectionList<int>
 *  - Loc<1> + - * / Grid<N>              => Grid<N>
 *  - Loc<N> + - * / Loc<N>               => Loc<N>
 *  - Loc<N> + - * / Interval<N>          => Interval<N>
 *  - Loc<N> + - * / Range<N>             => Range<N>
 *  - Loc<N> + - * / Grid<N>              => Grid<N>
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------


template<int Dim>
class Loc;

template<int Dim>
class Interval;

template<int Dim>
class Range;

template<class T> 
class IndirectionList;

template<int Dim>
class Grid;


template<class T1,class T2 >
struct DomainArithOpsTraits
{
};

template <>
struct DomainArithOpsTraits< Loc<1> , Loc<1> >
{
  typedef Loc<1> AddResult_t;
  typedef Loc<1> SubResult_t;
  typedef Loc<1> MultResult_t;
};

template<int Dim>
struct DomainArithOpsTraits<Loc<1> , Loc<Dim> >
{
  typedef Loc<Dim> AddResult_t;
  typedef Loc<Dim> SubResult_t;
  typedef Loc<Dim> MultResult_t;
};

template<int Dim>
struct DomainArithOpsTraits<Loc<Dim> , Loc<1> >
{
  typedef  Loc<Dim> AddResult_t;
  typedef  Loc<Dim> SubResult_t;
  typedef  Loc<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Loc<Dim> , Loc<Dim2> >
{
  typedef  Loc<Dim> AddResult_t;
  typedef  Loc<Dim> SubResult_t;
  typedef  Loc<Dim> MultResult_t;
};


//interval

template<int Dim>
struct DomainArithOpsTraits<Loc<1> , Interval<Dim> >
{
  typedef  Interval<Dim> AddResult_t;
  typedef  Range<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

template<int Dim>
struct DomainArithOpsTraits<Interval<Dim>, Loc<1>  >
{
  typedef  Interval<Dim> AddResult_t;
  typedef  Interval<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Loc<Dim> , Interval<Dim2> >
{
  typedef  Interval<Dim> AddResult_t;
  typedef  Range<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Interval<Dim>, Loc<Dim2>  >
{
  typedef  Interval<Dim> AddResult_t;
  typedef  Interval<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

//Range

template<int Dim>
struct DomainArithOpsTraits<Loc<1> , Range<Dim> >
{
  typedef  Range<Dim> AddResult_t;
  typedef  Range<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

template<int Dim>
struct DomainArithOpsTraits<Range<Dim>, Loc<1>  >
{
  typedef  Range<Dim> AddResult_t;
  typedef  Range<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Loc<Dim> , Range<Dim2> >
{
  typedef  Range<Dim> AddResult_t;
  typedef  Range<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Range<Dim>, Loc<Dim2>  >
{
  typedef  Range<Dim> AddResult_t;
  typedef  Range<Dim> SubResult_t;
  typedef  Range<Dim> MultResult_t;
};


//IndirectionList<int>


template< >
struct DomainArithOpsTraits<Loc<1> , IndirectionList<int> >
{
  typedef  IndirectionList<int> AddResult_t;
  typedef  IndirectionList<int> SubResult_t;
  typedef  IndirectionList<int> MultResult_t;
};

template< >
struct DomainArithOpsTraits<IndirectionList<int>, Loc<1>  >
{
  typedef  IndirectionList<int> AddResult_t;
  typedef  IndirectionList<int> SubResult_t;
  typedef  IndirectionList<int> MultResult_t;
};

//Grid

 
template<int Dim>
struct DomainArithOpsTraits<Loc<1> , Grid<Dim> >
{
  typedef  Grid<Dim> AddResult_t;
  typedef  Grid<Dim> SubResult_t;
  typedef  Grid<Dim> MultResult_t;
};

template<int Dim>
struct DomainArithOpsTraits<Grid<Dim>, Loc<1>  >
{
  typedef  Grid<Dim> AddResult_t;
  typedef  Grid<Dim> SubResult_t;
  typedef  Grid<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Loc<Dim> , Grid<Dim2> >
{
  typedef  Grid<Dim> AddResult_t;
  typedef  Grid<Dim> SubResult_t;
  typedef  Grid<Dim> MultResult_t;
};

template<int Dim,int Dim2>
struct DomainArithOpsTraits<Grid<Dim>, Loc<Dim2>  >
{
  typedef  Grid<Dim> AddResult_t;
  typedef  Grid<Dim> SubResult_t;
  typedef  Grid<Dim> MultResult_t;
};


#endif
