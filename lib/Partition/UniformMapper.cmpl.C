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

// ----------------------------------------------------------------------
// Member function map() of class UniformMapper
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Include files
// ----------------------------------------------------------------------

#include "Partition/UniformMapper.h"
#include "Pooma/Pooma.h"

// namespace Pooma {

void UniformMapper::map(const List_t& templist) const
{
  // assign the Nodes in the list to contexts in a uniform manner

  int contexts = Pooma::contexts();
  int npc = blocks_m.first() / contexts;
  int remainder = blocks_m.first() % contexts;

  int index = 0;
  for (int c=0; c<contexts; ++c) 
  {
    // put an equal number of Nodes on each context
    for (int i=0; i<npc; ++i)
      templist[index++]->context() = c;
    // handle any remainder
    if (c < remainder)
      templist[index++]->context() = c;
  }

  // assign local ID's and affinity values to local Nodes
  // defer this to ContextMapper base class

  setAffinity(templist);
}

// } // namespace Pooma

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformMapper.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
