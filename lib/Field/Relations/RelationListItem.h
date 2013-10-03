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
// RelationListItem: the ultimate base class for all relation objects.
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Relations
 * @brief
 * the ultimate base class for all relation objects.
 */

#ifndef POOMA_FIELD_RELATIONS_RELATIONLISTITEM_H
#define POOMA_FIELD_RELATIONS_RELATIONLISTITEM_H

#if defined(POOMA_DEBUG_RELATIONS)
#include <iostream>
#endif

#include "Field/Relations/RelationGroups.h"

/**
 * RelationListItem is an abstract base class for all relations:
 *
 * <PRE>
 *   RelationListItem
 *          |
 *   RelationRetargetBase<Target>
 *          |
 *   RelationBase<Target, Functor>
 *          |
 *   Relation<Target, Functor>
 * </PRE>
 *
 * The hierarchy exists to incrementally supply services and expose the Target
 * and Functor template parameters.
 *
 * RelationListItem provides an interface for the most basic services: responding
 * to pre-read and post-write events, setting/clearing the dirty flag, and
 * managing the priority/group membership.
 *
 * This class is not templated because we want to store these puppies in
 * RelationList containers, which are also not templated. 
 */

class RelationListItem {
public:

  //---------------------------------------------------------------------------
  // Default and copy constructors.

  RelationListItem() 
  : priority_m(0), 
    groups_m(Pooma::activeRelationGroups()),
    dirty_m(true) 
    { }
  
  RelationListItem(const RelationListItem &model)
  : priority_m(model.priority_m),
    groups_m(model.groups_m),
    dirty_m(model.dirty_m) 
    { }
  
  //---------------------------------------------------------------------------
  // Trivial destructor, but virtual since we will be subclassing.

  virtual ~RelationListItem() { }
  
  //---------------------------------------------------------------------------
  // Relation apply function. Subclasses must override to make this relation
  // actually do something.
  
  virtual void apply() = 0; 

  //---------------------------------------------------------------------------
  // Notification functions. 

  virtual void notifyPreRead()
    {
      // This function is called if somebody is getting ready to read somewhere 
      // and we need to update. By default, we check the dirty flag and update
      // if it is set. In either case, the dirty flag should be clear on exit.
    
      if (Pooma::isRelationGroupActive(groups_m) && dirty_m)
        {
          apply();
          clearDirty();
        }
    }
  
  virtual void notifyPostWrite()
    {
      // This function is called if somebody has already written somewhere.
      // By default, we simply set the dirty flag.
    
      setDirty();
    }

  //---------------------------------------------------------------------------
  // Accessors for the dirty flag and priority respectively.

  inline bool dirty() const { return dirty_m; }
  inline unsigned int priority() const { return priority_m; }
    
  //---------------------------------------------------------------------------
  // Modifiers used to set/clear the dirty flag.
  
  virtual void setDirty()
    {
#if defined(POOMA_DEBUG_RELATIONS)
      std::cout << "Setting dirty flag for " << (void *) this 
                << std::endl;
#endif
      dirty_m = true;
    }
     
  virtual void clearDirty()
    {
#if defined(POOMA_DEBUG_RELATIONS)
      std::cout << "Clearing dirty flag for " << (void *) this 
                << std::endl;
#endif
      dirty_m = false;
    }  
 
  //---------------------------------------------------------------------------
  // Modifier for the priority.

  inline void setPriority(unsigned int p)
    {
      priority_m = p;
    }
        
private:

  unsigned int priority_m;
  unsigned int groups_m;
  bool dirty_m;
  
};

#endif // POOMA_FIELD_RELATIONS_RELATIONLISTITEM_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RelationListItem.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
