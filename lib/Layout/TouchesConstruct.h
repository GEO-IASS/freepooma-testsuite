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
// Classes:
//   TouchesConstructNodePtr
//   TouchesConstructNodeObj
//   TouchesConstructINode
//
// Functions:
//   touchesConstruct
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_TOUCHESCONSTRUCT_H
#define POOMA_LAYOUT_TOUCHESCONSTRUCT_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Layout
 * @brief
 * touchesConstruct is a factory method that is used to build Nodes and INodes
 * by various layout touches methods.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/Unique.h"
#include "Layout/Node.h"
#include "Layout/GlobalIDDataBase.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
// TouchesConstructNodePtr, TouchesConstructNodeObj, and
// TouchesConstructINode
//   - tags for selecting specializations of...
//
// touchesConstruct
//   - overloaded factory method. Uses above tags to select version.
//     The appropriate version returns the object named in the tag.
//
//-----------------------------------------------------------------------------

struct TouchesConstructNodePtr {
  TouchesConstructNodePtr(){};
  ~TouchesConstructNodePtr(){};
};
struct TouchesConstructNodeObj {
  TouchesConstructNodeObj(){};
  ~TouchesConstructNodeObj(){};
};

// Build Node and return pointer; caller is responsible for deleting

template<class Domain>
inline Node<Domain> *
touchesConstruct(const Domain &owned, 
		 int affinity, int c, int gid, int lid,
		 const TouchesConstructNodePtr &)
{
  return new Node<Domain>(affinity, owned, c, gid, lid);
}

template<class Domain, class AllocatedDomain>
inline Node<Domain,AllocatedDomain> *
touchesConstruct(const Domain &owned, const AllocatedDomain &allocated,
		 int affinity, int c, int gid, int lid,
		 const TouchesConstructNodePtr &)
{
  return new Node<Domain,AllocatedDomain>
             (affinity, owned, allocated, c, gid, lid);
}

// Build Node object

template<class Domain>
inline Node<Domain>
touchesConstruct(const Domain &owned, 
		 int affinity, int c, int gid, int lid, 
		  const TouchesConstructNodeObj &)
{
  return Node<Domain>(affinity, owned, c, gid, lid);
}

template<class Domain, class AllocatedDomain>
inline Node<Domain,AllocatedDomain>
touchesConstruct(const Domain &owned, const AllocatedDomain &allocated,
		 int affinity, int c, int gid, int lid, 
		 const TouchesConstructNodeObj &)
{
  return Node<Domain,AllocatedDomain>(affinity, owned, allocated, c, gid, lid);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_LAYOUT_TOUCHESCONSTRUCT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TouchesConstruct.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
