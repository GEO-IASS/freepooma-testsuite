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
// RelationList: manages a list of relations.
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Relations
 * @brief
 * manages a list of relations.
 */

#ifndef POOMA_FIELD_RELATIONS_RELATIONLIST_H
#define POOMA_FIELD_RELATIONS_RELATIONLIST_H

#include "Field/Relations/RelationListItem.h"
#include "Field/Relations/RelationBases.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/RefCounted.h"

#include <vector>

#if defined(POOMA_DEBUG_RELATIONS)
#include <iostream>
#endif

/**
 * RelationListData is a ref-counted object that holds the data for an
 * relation list.
 */

class RelationListData : public RefCounted
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  // The default constructor implicitly calls the default constructor for the
  // contained vector.
  
  RelationListData()
    { }
    
  // Copy constructor.

  RelationListData(const RelationListData &model)
  : data_m(model.data_m)
    { }
  
  
  //---------------------------------------------------------------------------
  // Copy assignment operator.
  
  RelationListData &operator=(const RelationListData &rhs)
    {
      data_m = rhs.data_m;
    
      return *this;
    }
  

  //---------------------------------------------------------------------------
  // We need to delete the actual RelationListItems.

  ~RelationListData() 
    {
      typedef List_t::size_type size_type;
      for (size_type i = 0; i < data_m.size(); i++)
        delete data_m[i]; 
    }


  //---------------------------------------------------------------------------
  // Return the number of RelationListItems.
  
  inline int size() const
    {
      return data_m.size();
    }

  
  //---------------------------------------------------------------------------
  // Returns the ith element.
  
  inline RelationListItem *elem(int i) const
    {
      return data_m[i];
    }
  
  inline RelationListItem* &elem(int i)
    {
      return data_m[i];
    }


  //---------------------------------------------------------------------------
  // Pushes an item on the list.

  void add(RelationListItem *item)
    {
      // Add the item.
    
      data_m.push_back(item);
    
      // Move it to its position. Make use of the fact that item is at 
      // data_m[i] each iteration through the loop.
    
      int i = size() - 1;
      while (i > 0)
        {
          if (data_m[i]->priority() <= data_m[i]->priority())
            break;
          data_m[i] = data_m[i - 1];
          data_m[i - 1] = item;
          --i;
        }
    }
   
private:

  typedef std::vector<RelationListItem *> List_t;

  List_t data_m;
};


/**
 * RelationList is a container that dispatches events to the list of boundary 
 * conditions it contains. 
 */

class RelationList {
public:

  //---------------------------------------------------------------------------
  // The default constructor simply makes an empty list. 
  
  RelationList()
  : list_m(new RelationListData)
    { }


  //---------------------------------------------------------------------------
  // The copy constructor makes a shallow copy of the data.
  
  RelationList(const RelationList &model)
  : list_m(model.list_m)
    { }


  //---------------------------------------------------------------------------
  // Destructor: trivial.
  
  ~RelationList() { }


  //---------------------------------------------------------------------------
  // Replaces the current list with a private copy of itself.
  
  template<class Target>
  void makeOwnCopy(const Target &t)
    {
      // Make a copy of the list. 

      list_m = new RelationListData(*list_m);
       
      // Now, we need to replace the individual RelationListItems
      // with cloned versions with this Target as the subject. 
      
      for (int i = 0; i < list_m->size(); i++)
        {
          RelationRetargetBase<Target> *u = 
            dynamic_cast<RelationRetargetBase<Target> *>(list_m->elem(i));
          if (u != NULL)
            list_m->elem(i) = u->retarget(t);
        }
    }


  //---------------------------------------------------------------------------
  // Replaces the current list with an empty list.
  
  inline void erase()
    {
     list_m = new RelationListData;
    }


  //---------------------------------------------------------------------------
  // Adds an relation to the end of our list.

  void addRelation(RelationListItem *item)
    {
      list_m->add(item);

#if defined(POOMA_DEBUG_RELATIONS)
      std::cout << "Adding Relation " << (void *) item << std::endl;
#endif
    }

  
  //---------------------------------------------------------------------------
  // Notify the relations about pre-read/post-write events.

  void notifyPreRead() const
    {
      for (int i = 0; i < list_m->size(); ++i)
        list_m->elem(i)->notifyPreRead();
    }

  void notifyPostWrite() const
    {
      for (int i = 0; i < list_m->size(); ++i)
        list_m->elem(i)->notifyPostWrite();
    }

  
  //---------------------------------------------------------------------------
  //@name dirty flag handling, deprecated
  //@{

  /// Set the dirty flags for all relations.

  void setDirty() const
    {
      for (int i = 0; i < list_m->size(); ++i)
        list_m->elem(i)->setDirty();
    }

  /// Set the dirty flags for all relations.

  void clearDirty() const
    {
      for (int i = 0; i < list_m->size(); ++i)
        list_m->elem(i)->clearDirty();
    }
  
  /// Query if any of the relations is dirty.

  bool dirty() const
    {
      for (int i = 0; i < list_m->size(); ++i)
        if (list_m->elem(i)->dirty())
	  return true;
      return false;
    }

  //@}
  
  //---------------------------------------------------------------------------
  // Give access to a specific relation.

  inline RelationListItem *operator()(int i) const
    {
#if POOMA_BOUNDS_CHECK
      PInsist2(i >= 0 && i < list_m->size(),
        "RelationList bounds error: index = %d, size = %d.",
        i, list_m->size());
#endif

      return list_m->elem(i);
    }

  inline RelationListItem *operator()(int i)
    {
#if POOMA_BOUNDS_CHECK
      PInsist2(i >= 0 && i < list_m->size(),
        "RelationList bounds error: index = %d, size = %d.",
        i, list_m->size());
#endif

      return list_m->elem(i);
    }


  //---------------------------------------------------------------------------
  // Return the number of relations.
  
  inline int size() const { return list_m->size(); }
  
  
private:

  RefCountedPtr<RelationListData> list_m;
};


#endif // POOMA_FIELD_RELATIONS_RELATIONLIST_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RelationList.h,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/10 22:06:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
