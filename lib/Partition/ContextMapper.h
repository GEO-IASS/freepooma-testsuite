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

#ifndef POOMA_CONTEXTMAPPER_H
#define POOMA_CONTEXTMAPPER_H

//-----------------------------------------------------------------------------
// Classes:
// ContextMapper
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Partition
 * @brief
 * ContextMapper is the base class for classes used to map a node list
 * to contexts.
 *
 * Available mappers are LocalMapper for non-distributed use and
 * DistributedMapper for distributed use.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Domain/Interval.h"
#include "Utilities/PAssert.h"
#include "Layout/Node.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<int Dim>
class ContextMapper
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef Interval<Dim>                       Domain_t;
  typedef Node<Domain_t>                      Value_t;
  typedef std::vector<Value_t *>              List_t;

  //============================================================
  // Constructors
  //============================================================

  ContextMapper(){};

  virtual ~ContextMapper(){};

  virtual void map(const List_t & templist) const = 0;
  
  void setAffinity(const List_t & templist) const;

};

template<int Dim>
void ContextMapper<Dim>::setAffinity(const List_t & templist) const
{
  int affinityMax = Smarts::concurrency();
  int idMax = 0;

  typename List_t::const_iterator start = templist.begin();
  typename List_t::const_iterator end = templist.end();

  for ( ; start != end ; ++start)
    if((*start)->context()==Pooma::context())
      {
	(*start)->localID()=idMax;
	++idMax;
      }

  start = templist.begin();
  for ( ; start != end ; ++start)
    { 
      if((*start)->context()==Pooma::context())
	(*start)->affinity() = static_cast<int>
	  ( affinityMax * ( (*start)->localID()
			    / static_cast<double>(idMax) ) );
    }

  return;
} 


template<int Dim>
class LocalMapper
  : public ContextMapper<Dim>
{ 
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef Interval<Dim>                       Domain_t;
  typedef Node<Domain_t>                      Value_t;
  typedef std::vector<Value_t *>              List_t;

  template<class Partitioner>
  LocalMapper(const Partitioner &)
  {}
  
  LocalMapper()
  {}
  
  void map(const List_t & templist) const;

};

template<int Dim>
void LocalMapper<Dim>::map(const List_t & templist) const
{
  int idMax = templist.size();
  int naff = Smarts::concurrency();
  for (int i = 0; i< templist.size(); ++i)
    {
      templist[i]->context() = -1;
      templist[i]->localID() = i;
      templist[i]->affinity() = static_cast<int>( ( naff * ( i / 
				static_cast<double>(idMax) ) ) );
    }
}


#endif     // POOMA_CONTEXTMAPPER_H
