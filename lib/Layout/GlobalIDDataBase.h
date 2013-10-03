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
// GlobalIDDataBase
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_GLOBALIDDATABASE_H
#define POOMA_LAYOUT_GLOBALIDDATABASE_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Layout
 * @brief
 * GlobalIDDataBase stores global patch ID's for INodes.
 *
 * Since an INode can
 * come from intersecting several different layouts, the global ID's are
 * stored with a unique ID that comes from the layout.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <vector>
#include <map>
#include <utility>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * GlobalIDDataBase stores global patch ID's for INodes.  Since an INode can
 * come from intersecting several different layouts, the global ID's are
 * stored with a layout ID that comes from the layout.
 *
 * Basically this database stores a map: (layoutID,nodeKey)->globalID
 *
 * The nodeKey is a unique key generated for each INode that is pushed
 * into this database.  The database also tracks where a node came from
 * in the intersection process.  When you perform touches on a given INode
 * to produce new INodes, the original INode's key is used as the parent
 * key.  That way we can trace back through the parents to find the global
 * ID for layouts that were previously intersected.
 *
 * The interface GlobalIDDataBase provides is:
 *
 * GlobalIDDataBase gidStore;
 *
 * nodeKey = gidStore.push(layoutID,globalID);
 *  - adds a new node to the database for the layout with ID layoutID.
 *
 * nodeKey = gidStore.push(layoutID,globalID,parentKey);
 *  - adds a new node with parent node with key parentKey.
 *
 * globalID = gidStore.globalID(layoutID,myKey);
 *  - get the global ID given a layout ID and a node key.
 *
 * gidStore.shared(layoutID1,layoutID2);
 *  - tells the database that we didn't intersect layout layoutID1 because it
 *    is the same as layoutID2 (they're identical views).
 */

class GlobalIDDataBase
{
public:

  typedef int LayoutID_t;
  typedef int GlobalID_t;
  typedef int NodeKey_t;  // don't change this to unsigned!

  typedef std::map<LayoutID_t, GlobalID_t> Shared_t;

  //----------------------------------------------------------------------
  // Simple constructor.  Since this data base is an internal pooma object
  // that gets used during expression evaluation, it should never be
  // copied or assigned.

  GlobalIDDataBase() { }

  //----------------------------------------------------------------------
  // nullNodeKey() is a nodeKey that acts like end()

  inline static
  NodeKey_t nullNodeKey()
  {
    return -1;
  }

  //----------------------------------------------------------------------
  // push() puts a new node in the database.  You provide a layoutID
  // and the globalID and it returns a nodeKey to access the database
  // later.

  NodeKey_t push(LayoutID_t layoutID, int context, GlobalID_t globalID);

  //----------------------------------------------------------------------
  // This version of push() takes a nodeKey that represents the parent
  // node.  (Parent nodes are the result of previous intersections that
  // contain a given INode.)

  NodeKey_t push(LayoutID_t layoutID,
		 int context,
		 GlobalID_t globalID,
		 NodeKey_t parentNode);

  //----------------------------------------------------------------------
  // Inform the database that a given layout was not intersected because
  // it was identical to another layout.

  void shared(LayoutID_t idNew, LayoutID_t idOld);

  //----------------------------------------------------------------------
  // Access the globalID for a given layoutID and nodeKey.  We search
  // through a node and all its parents for the right layout id.

  GlobalID_t globalID(LayoutID_t layoutID, NodeKey_t key) const;

  //----------------------------------------------------------------------
  // Access the context for a given layoutID and nodeKey.  We search
  // through a node and all its parents for the right layout id.

  int context(LayoutID_t layoutID, NodeKey_t key) const;

  //----------------------------------------------------------------------
  // Access the most common context for a given nodeKey.
  // The current version just picks the first context which is only a bad
  // choice when we're performing reductions on expressions with multiple
  // unaligned arrays.

  int context(NodeKey_t key) const;

  //----------------------------------------------------------------------
  // Asks the question, does a given context participate in the
  // intersections used to create a given inode?
  
  bool contextParticipates(int context, NodeKey_t key) const;

  //----------------------------------------------------------------------
  // print method for debugging purposes

  template<class OSTR>
  inline void print(OSTR &ostr)
  {
    typedef std::vector<Pack> Store_t;
    typedef typename Store_t::const_iterator Iterator_t;
    Iterator_t p = data_m.begin();
    for (; p != data_m.end(); ++p)
    {
      ostr << "(" << (*p).layoutID() << ","
	   << (*p).globalID() << ","
	   << (*p).context() << ","
	   << (*p).parent() << ")";
    }
  }

private:

  // Pack is a utility structure containing the database records
  // (layout ID, global ID, parent node)

  struct Pack
  {
    Pack()
      : layoutID_m(0), context_m(0), globalID_m(0), parent_m(0)
    { }

    inline
    Pack(LayoutID_t layoutID, int context, GlobalID_t globalID,
	 NodeKey_t parent)
      : layoutID_m(layoutID),
	context_m(context),
	globalID_m(globalID),
	parent_m(parent)
    { }

    inline LayoutID_t layoutID() const { return layoutID_m; }
    inline int        context() const { return context_m; }
    inline GlobalID_t globalID() const { return globalID_m; }
    inline NodeKey_t  parent() const { return parent_m; }

    LayoutID_t layoutID_m;
    int        context_m;
    GlobalID_t globalID_m;
    NodeKey_t  parent_m;
  };

  // The actual database is stored in this vector.

  std::vector<Pack> data_m;

  // Information on layouts which share entries:

  Shared_t shared_m;
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_LAYOUT_GLOBALIDDATABASE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: GlobalIDDataBase.h,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
