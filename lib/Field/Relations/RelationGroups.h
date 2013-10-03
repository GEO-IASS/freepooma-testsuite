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

/** @file
 * @ingroup Relations
 * @brief
 * Functions supporting relation groups (undocumented).
 */

#ifndef POOMA_FIELD_RELATIONS_RELATIONGROUPS_H
#define POOMA_FIELD_RELATIONS_RELATIONGROUPS_H

namespace Pooma {

  unsigned int activeRelationGroups();
  bool isRelationGroupActive(unsigned int groups);
    
  void activateRelationGroup(unsigned int group);
  void deactivateRelationGroup(unsigned int group);
  unsigned int newRelationGroup();

}  // namespace Pooma

#endif // POOMA_FIELD_RELATIONS_RELATIONGROUPS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RelationGroups.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
