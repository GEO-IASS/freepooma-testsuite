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
// RelationRetargetBase & RelationBase: the base classes for all relations.
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Relations
 * @brief
 * RelationRetargetBase & RelationBase: the base classes for all relations.
 */

#ifndef POOMA_FIELD_UPDATER_RELATIONBASES_H
#define POOMA_FIELD_UPDATER_RELATIONBASES_H

#include "Field/Relations/RelationListItem.h"


namespace Pooma {

  /**
   * This tag is used to tell constructors not to copy relations.
   */

  struct DontCopyRelations { };

}


/**
 * RelationRetargetBase is an abstract base class for all relations:
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
 * It is sometimes necessary to be able to replace the target of an
 * relation with another object. RelationListItem is not templated, so the
 * idea is to do a dynamic_cast to an RelationRetargetBase and then call the
 * virtual retarget function, which needs to be implemented in an RelationBase
 * subclass.
 *
 * RelationRetargetBase provides an interface for storing and accessing the 
 * target; that is, the object that is going to be updated by this relation.
 *
 * Subclasses must define the retarget function, which is used to create a
 * clone of this Relation, but with a new target.
 */

template<class Target>
class RelationRetargetBase : public RelationListItem {
public:

  //---------------------------------------------------------------------------
  // Constructors. 

  // Initialize the target.
  
  RelationRetargetBase(const Target &target)
  : target_m(target, Pooma::DontCopyRelations()) 
    {  }

  // Copy constructor. 

  RelationRetargetBase(const RelationRetargetBase<Target> &model)
  : RelationListItem(model),
    target_m(model.target_m) 
    { }
      
  //---------------------------------------------------------------------------
  // Trivial destructor, but virtual since we will be subclassing.

  virtual ~RelationRetargetBase() { }

  //---------------------------------------------------------------------------
  // Target accessor.
  
  Target &target() { return target_m; }
  const Target &target() const { return target_m; }

  //---------------------------------------------------------------------------
  // Retarget function. Clones this relation for a new target.
 
  virtual RelationListItem *retarget(const Target &target) const = 0;

protected:

  Target target_m;
    
};


/** RelationBase is an abstract base class for all automatic relations.
 * It is sometimes necessary to be able to replace the target of an
 * relation with another object. RelationListItem is not templated so the
 * idea is to do a dynamic_cast to an RelationRetargetBase and then call the
 * virtual retarget function, which needs to be implemented in RelationBase
 * subclasses.
 */

template<class Target, class Functor>
class RelationBase : public RelationRetargetBase<Target> {
public:

  //---------------------------------------------------------------------------
  // Constructors. 

  // Initialize the target and functor.
  
  RelationBase(const Target &t, const Functor &f)
  : RelationRetargetBase<Target>(t),
    functor_m(f, t)
    { }

  // Copy constructor. 

  RelationBase(const RelationBase<Target, Functor> &model)
  : RelationRetargetBase<Target>(model),
    functor_m(model.functor_m) 
    { }
     
  //---------------------------------------------------------------------------
  // Trivial destructor, but virtual since we will be subclassing.

  virtual ~RelationBase() { }
  
  //---------------------------------------------------------------------------
  // Functor accessor.
  
  Functor &functor() { return functor_m; }
  const Functor &functor() const { return functor_m; }
    
protected:

  Functor functor_m;
    
};

#endif // POOMA_FIELD_UPDATER_RELATIONBASES_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RelationBases.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
