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
// SparseTileLayout test: Create and use SparseTileLayout objects
//-----------------------------------------------------------------------------
#include "Pooma/Pooma.h"
#include "Pooma/Domains.h"
#include "Utilities/Tester.h"
#include <iostream>
#include <iterator>

/** @file
 * @ingroup Domain
 * @brief
 * TBD.
 */

template <int Dim> 
std::vector<Interval<Dim> >
DomainRemoveOverlap(const Interval<Dim> & s,const Interval<Dim> &r)
{
  typedef Interval<Dim>         Domain_t;
  typedef std::vector<Domain_t> DomainList_t;

  DomainList_t result,temp;

  result.push_back(s);

  for (int i=0;i<Dim;++i)
    {
      typename DomainList_t::iterator start = result.begin();
      typename DomainList_t::iterator end = result.end();
      for ( ; start!=end; ++start)
	{ 
	  if (touches( (*start)[i], Loc<1>(r[i].min())))
	    {
	      Domain_t lower=*start,upper=*start;
	      if(r[i].min()-1>=lower[i].min())
		{
		  lower[i] = Interval<1>(lower[i].min(),r[i].min()-1);
		  temp.push_back(lower);
		}	      upper[i] = Interval<1>(r[i].min(),upper[i].max());
	
	      temp.push_back(upper);
	    }
	  else 
	    temp.push_back(*start);
	}
      result = temp;
      temp.clear();
      
      start = result.begin();
      end = result.end();
      for ( ; start!=end; ++start)
	{ 
	  if (touches( (*start)[i], Loc<1>(r[i].max())))
	    {
	      Domain_t lower=*start,upper=*start;
	      lower[i] = Interval<1>(lower[i].min(),r[i].max());
	      temp.push_back(lower);
	      if( r[i].max()+1 <= upper[i].max())
		{
		  upper[i] = Interval<1>(r[i].max()+1,upper[i].max());
		  temp.push_back(upper);
		}
	    }
	  else
	    temp.push_back(*start);
	}
      result=temp;
      temp.clear();
    }
  
  typename DomainList_t::iterator start = result.begin();
  typename DomainList_t::iterator end = result.end();
  for ( ; start!=end ; ++start)
    {
      if (!touches(*start,r ) ) 
	temp.push_back(*start);
    }
  return temp;
}
