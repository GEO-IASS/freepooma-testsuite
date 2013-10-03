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

#ifndef POOMA_PARTITION_UNIFORM_MAPPER_H
#define POOMA_PARTITION_UNIFORM_MAPPER_H

//-----------------------------------------------------------------------------
// Classes:
//   UniformMapper
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Partition/ContextMapper.h"
#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Layout/Node.h"

#include <vector>

/** @file
 * @ingroup Partition
 * @brief
 * UniformMapper is a ContextMapper implementation
 * specifically for 1D patches.
 */

// namespace Pooma {

/**
 * UniformMapper is a ContextMapper specifically for 1D patches.
 * All it does is put a roughly equal number of patches on each context.
 */

class UniformMapper
  : public ContextMapper<1>
{ 
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef Interval<1>                         Domain_t;
  typedef Node<Domain_t>                      Value_t;
  typedef std::vector<Value_t *>              List_t;

  template <class Partitioner>
  inline
  UniformMapper(const Partitioner& gp) 
    : blocks_m(gp.blocks())
  {
  }

  inline
  UniformMapper(const Loc<1>& blocks)
    : blocks_m(blocks)
  {
  }

  inline
  UniformMapper(int blocks = 1)
    : blocks_m(blocks)
  {
  }

  virtual ~UniformMapper(){}

  void map(const List_t&) const;

  // member data
private:

  Loc<1> blocks_m;
};


// } // namespace Pooma

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTITION_UNIFORM_MAPPER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: UniformMapper.h,v $   $Author: richi $
// $Revision: 1.7 $   $Date: 2004/11/10 22:17:08 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
