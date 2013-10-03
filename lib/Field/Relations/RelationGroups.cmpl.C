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
// Functions supporting relation groups.
//-----------------------------------------------------------------------------

#include "Field/Relations/RelationGroups.h"
#include "Pooma/Pooma.h"

namespace Pooma {

  namespace {
    
    unsigned int activeGroups_g = 1;
    unsigned int numGroups_g = 1;
      
  } // anonymous namespace
      
  unsigned int activeRelationGroups()
    {
      return activeGroups_g;
    }
      
  bool isRelationGroupActive(unsigned int groups)
    {
      return (groups & activeGroups_g) != 0;
    }
    
  void activateRelationGroup(unsigned int group)
    {
      blockAndEvaluate();
      activeGroups_g |= group;
    }
      
  void deactivateRelationGroup(unsigned int group)
    {
      blockAndEvaluate();
      activeGroups_g &= ~group;
    }
      
  unsigned int newRelationGroup()
    {
      unsigned int n = (1 << numGroups_g++);
      activateRelationGroup(n);
      return n;
    }

}  // namespace Pooma

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RelationGroups.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
