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

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//
// GlobalIDDataBase stores global patch ID's for INodes.  Since an INode can
// come from intersecting several different layouts, the global ID's are
// stored with a unique ID that comes from the layout.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Layout/GlobalIDDataBase.h"

//----------------------------------------------------------------------
// push() puts a new node in the database.  You provide a layoutID
// and the globalID and it returns a nodeKey to access the database
// later.

GlobalIDDataBase::NodeKey_t
GlobalIDDataBase::push(LayoutID_t layoutID, int context, GlobalID_t globalID)
{
  int ret = data_m.size();
  data_m.push_back(Pack(layoutID, context, globalID, nullNodeKey()));
  return ret;
}

//----------------------------------------------------------------------
// This version of push() takes a nodeKey that represents the parent
// node.  (Parent nodes are the result of previous intersections that
// contain a given INode.)

GlobalIDDataBase::NodeKey_t
GlobalIDDataBase::push(LayoutID_t layoutID,
		       int context,
		       GlobalID_t globalID,
		       NodeKey_t parentNode)
{
  int ret = data_m.size();
  data_m.push_back(Pack(layoutID, context, globalID, parentNode));
  return ret;
}

//----------------------------------------------------------------------
// Inform the database that a given layout was not intersected because
// it was identical to another layout.

void
GlobalIDDataBase::shared(LayoutID_t idNew, LayoutID_t idOld)
{
  // if the old layout was also not stored because it's has the
  // same intersections as a previous one, then we want to point
  // to the previous layout

  Shared_t::const_iterator p = shared_m.find(idOld);
  if (p != shared_m.end())
  {
    idOld = p->second;
  }

  Shared_t::value_type i(idNew, idOld);
  shared_m.insert(i);
}

//----------------------------------------------------------------------
// Access the globalID for a given layoutID and nodeKey.  We search
// through a node and all its parents for the right layout id.

GlobalIDDataBase::GlobalID_t
GlobalIDDataBase::globalID(LayoutID_t layoutID, NodeKey_t key) const
{
  // First check if the layout is in the list of layouts that were
  // bypassed:

  Shared_t::const_iterator p = shared_m.find(layoutID);
  if (p != shared_m.end())
  {
    layoutID = p->second;
  }

  while( key != nullNodeKey() )
  {
    if (data_m[key].layoutID() == layoutID)
    {
      return data_m[key].globalID();
    }
    else
    {
      key = data_m[key].parent();
    }
  }

  // If we reach this point, then the database has been corrupted
  // or we're asking for the id from a layout that wasn't intersected.

  PAssert(false);
  return -1;
}

//----------------------------------------------------------------------
// Access the context for a given layoutID and nodeKey.  We search
// through a node and all its parents for the right layout id.

int
GlobalIDDataBase::context(LayoutID_t layoutID, NodeKey_t key) const
{
  // First check if the layout is in the list of layouts that were
  // bypassed:

  Shared_t::const_iterator p = shared_m.find(layoutID);
  if (p != shared_m.end())
  {
    layoutID = p->second;
  }

  while( key != nullNodeKey() )
  {
    if (data_m[key].layoutID() == layoutID)
    {
      return data_m[key].context();
    }
    else
    {
      key = data_m[key].parent();
    }
  }

  // If we reach this point, then the database has been corrupted
  // or we're asking for the id from a layout that wasn't intersected.

  PAssert(false);
  return -1;
}

//----------------------------------------------------------------------
// Access the most common context for a given nodeKey.
// The current version just picks the first context which is only a bad
// choice when we're performing reductions on expressions with multiple
// unaligned arrays.

int
GlobalIDDataBase::context(NodeKey_t key) const
{
  PAssert(key != nullNodeKey());
  
  return data_m[key].context();
}

//----------------------------------------------------------------------
// Asks the question, does a given context participate in the
// intersections used to create a given inode?
  
bool
GlobalIDDataBase::contextParticipates(int context, NodeKey_t key) const
{
  while( key != nullNodeKey() )
  {
    if (data_m[key].context() == context)
    {
      return true;
    }
    else
    {
      key = data_m[key].parent();
    }
  }

  return false;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: GlobalIDDataBase.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
