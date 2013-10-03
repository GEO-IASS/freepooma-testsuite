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
// DistributedMapper<Dim>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Partition
 * @brief
 * DistributedMapper is a ContextMapper implementation.
 *
 * It dispatches to one of ContiguousMapper, UniformMapper
 * and BisectionMapper based on the size and structure of contexts
 * and nodes (patches).
 */

#ifndef POOMA_DISTRIBUTEDMAPPER_H
#define POOMA_DISTRIBUTEDMAPPER_H

#include "Partition/ContextMapper.h"
#include "Partition/ContiguousMapper.h"
#include "Partition/BisectionMapper.h"
#include "Partition/UniformMapper.h"
#include "Utilities/WrappedInt.h"



/**
 * DistributedMapper is a ContextMapper implementation.
 *
 * It dispatches to one of ContiguousMapper, UniformMapper
 * and BisectionMapper based on the size and structure of contexts
 * and nodes (patches).
 */

template<int Dim>
class DistributedMapper
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
  DistributedMapper(const Partitioner & gp) 
    : blocks_m(gp.blocks())
  {
  }


  void map(const List_t & templist) const;

  void uniformMap(const Loc<1> &blocks,
		  const List_t &templist,
		  const WrappedInt<1>&) const
  {
    UniformMapper(blocks).map(templist);
  }

  template <int D>
  void uniformMap(const Loc<D> &,
		  const List_t &,
		  const WrappedInt<D>&) const
  {
    PAssert(false); // uniformMapper is 1 D only 
  }

  // member data
private:

  Loc<Dim> blocks_m;

};

template<int Dim>
void DistributedMapper<Dim>::map(const List_t & templist) const
{
  int ncontexts = Pooma::contexts();
  int npc = templist.size()/ncontexts;
  // If there are more contexts than patches, assign one
  // patch per context for as many patches as there are. 
  if(ncontexts> templist.size())
    {
      // we should probably alert the user here!!
      npc = 1;
      ncontexts = templist.size();
    }

  if (Dim == 1)
    {
      // work around, since UniformMapper is 1-dim
      uniformMap(blocks_m,templist,WrappedInt<Dim>());
    }
  else if(npc<3)
    {
      ContiguousMapper<Dim>(blocks_m).map(templist);
    }
  else
    {
      BisectionMapper<Dim>(blocks_m).map(templist);
    }
  return;      
}


#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DistributedMapper.h,v $   $Author: richi $
// $Revision: 1.12 $   $Date: 2004/11/10 22:17:08 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
