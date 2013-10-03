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
//   INode           
//     - a node-like class for communicating information regarding
//       intersections of layouts. Currently very simple.
//
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_INODE_H
#define POOMA_LAYOUT_INODE_H

/** @file
 * @ingroup Layout
 * @brief
 *   a node-like class for communicating information regarding
 *   intersections of layouts. Currently very simple.
 */

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "Domain/Contains.h"
#include "Domain/Interval.h"
#include "Domain/Loc.h"
#include "Layout/Node.h"
#include "Layout/GlobalIDDataBase.h"
#include "Layout/TouchesConstruct.h"
#include <iosfwd>


template<int Dim> class INode;

/**
 * TouchesConstructINode is used to construct an INode during the touches()
 * operation in layouts.
 */

template<int Dim>
struct TouchesConstructINode
{
  typedef GlobalIDDataBase::NodeKey_t NodeKey_t;
  typedef Unique::Value_t             LayoutID_t;

  TouchesConstructINode(LayoutID_t layoutID,
			NodeKey_t parent,
			GlobalIDDataBase *globalIDDataBase)
    : layoutID_m(layoutID), parent_m(parent),
      globalIDDataBase_m(globalIDDataBase)
  { }

  inline LayoutID_t layoutID() const { return layoutID_m; }
  inline NodeKey_t  parent() const { return parent_m; }
  inline GlobalIDDataBase *globalIDDataBase() const
  {
    return globalIDDataBase_m;
  }

  LayoutID_t        layoutID_m;
  NodeKey_t         parent_m;
  GlobalIDDataBase *globalIDDataBase_m;
};

/**
 * INode is a class for communicating information regarding
 * intersections of layouts. It is used by engine classes to generate
 * efficient views. Currently, this class is very simple.
 */

template <int Dim>
class INode
{
public:
  
  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef INode<Dim>                                    This_t;
  typedef Interval<Dim>                                 Domain_t;
  typedef GlobalIDDataBase::LayoutID_t                  LayoutID_t;
  typedef GlobalIDDataBase::GlobalID_t                  GlobalID_t;
  typedef GlobalIDDataBase::NodeKey_t                   NodeKey_t;
  
  enum { dimensions = Dim };
  
  
  //===========================================================================
  // Constructors
  //===========================================================================

  // Default constructor leaves our domain unitialized.
  
  inline INode() : domain_m(Pooma::NoInit()) { }
  
  // Copy constructor does a memberwise copy.
  
  inline INode(const INode<Dim> &model)
    : domain_m(model.domain_m),
      globalIDDataBase_m(model.globalIDDataBase_m),
      key_m(model.key_m)
  { }
  
  // This version lets you modify the domain of a given INode while
  // maintaining the global ID information.  (used by stencils, where
  // we're guaranteed to have enough guard cells to expand the domain.
  // In general this is a dangerous operation.)
  
  template<int D2, class Dom>  
  inline INode(const INode<D2> &model, const Dom &dom)
    : domain_m(dom),
      globalIDDataBase_m(model.globalIDDataBase()),
      key_m(model.key())
  { }
  
  // Constructors that take a domain, layout id, global id,
  // global id database, and optionally a parent key.
  // These versions store the global id information in the
  // database for fast lookup later.
  
  inline
  INode(const Interval<Dim> &dom, LayoutID_t layoutID, int context,
	GlobalID_t globalID,
	GlobalIDDataBase *globalIDDataBase, NodeKey_t parent = -1)
    : domain_m(dom),
      globalIDDataBase_m(globalIDDataBase)
  {
    key_m = globalIDDataBase_m->push(layoutID, context, globalID, parent);
  }

  template<class Alloc>
  inline
  INode(const Node<Interval<Dim>, Alloc> &node, LayoutID_t layoutID,
	GlobalIDDataBase *globalIDDataBase)
    : domain_m(node.domain()),
      globalIDDataBase_m(globalIDDataBase)
  {
    key_m = globalIDDataBase_m->push(layoutID, node.context(),
				     node.globalID());
  }
  
  // TouchesConstruct versions all take (domain, context, globalID, tcin)

  inline
  INode(const Interval<Dim> &dom, int context, GlobalID_t globalID,
	const TouchesConstructINode<Dim> &tcin)
    : domain_m(dom),
      globalIDDataBase_m(tcin.globalIDDataBase())
  {
    key_m = globalIDDataBase_m->push(tcin.layoutID(), context, globalID,
				     tcin.parent());
  }
  
  inline
  INode(const INode<Dim> &inode, int context, GlobalID_t globalID,
	const TouchesConstructINode<Dim> &tcin)
    : domain_m(inode.domain()),
      globalIDDataBase_m(tcin.globalIDDataBase())
  {
    key_m = globalIDDataBase_m->push(tcin.layoutID(), context, globalID,
				     tcin.parent());
  }

  template<class Alloc>
  inline
  INode(const Node<Interval<Dim>, Alloc> &node, int context,
	GlobalID_t globalID,
	const TouchesConstructINode<Dim> &tcin)
    : domain_m(node.domain()),
      globalIDDataBase_m(tcin.globalIDDataBase())
  {
    key_m = globalIDDataBase_m->push(tcin.layoutID(), context, globalID,
				     tcin.parent());
  }
  
  // Range versions store an interval that spans the range.

  inline INode(const Range<Dim> &range, int context, int globalID,
	       const TouchesConstructINode<Dim> &tcin)
    : globalIDDataBase_m(tcin.globalIDDataBase())
  {
    key_m = globalIDDataBase_m->push(tcin.layoutID(), context, globalID,
				     tcin.parent());
    int i;
    for (i = 0; i < Dim; ++i)
    {
      domain_m[i] = Interval<1>(range[i].first(), range[i].last());
    }
  }

  //===========================================================================
  // Assignment operator
  //===========================================================================
  
  // Assignment operator does a memberwise assign.
  
  inline INode<Dim> &operator=(const INode<Dim> &rhs)
  {
    // Check for self-assignment.
    
    if (&rhs != this)
    {
      domain_m   = rhs.domain();
      globalIDDataBase_m = rhs.globalIDDataBase();
      key_m      = rhs.key();
    }
      
    return *this;
  }
  
  //===========================================================================
  // Destructor
  //===========================================================================

  // Trival since members manage their own data.
  
  inline ~INode() { }
  

  //===========================================================================
  // Accessors
  //===========================================================================

  // Return the domain.

  inline const Domain_t &domain() const { return domain_m; }

  // Find the global ID for this node and a given layout ID in the
  // database.

  inline
  GlobalID_t globalID(LayoutID_t id) const
  {
    PAssert(globalIDDataBase_m != NULL);
    return globalIDDataBase_m->globalID(id, key_m);
  }

  inline int context() const
  {
    PAssert(globalIDDataBase_m != NULL);
    return globalIDDataBase_m->context(key_m);
  }

  inline int context(LayoutID_t id) const
  {
    PAssert(globalIDDataBase_m != NULL);
    return globalIDDataBase_m->context(id, key_m);
  }

  inline bool contextParticipates(int context) const
  {
    PAssert(globalIDDataBase_m != NULL);
    return globalIDDataBase_m->contextParticipates(context, key_m);
  }

  // Return the global ID database pointer.

  inline
  GlobalIDDataBase *globalIDDataBase() const { return globalIDDataBase_m; }

  // Return the node key for the database.

  inline NodeKey_t key() const { return key_m; }

  //===========================================================================
  // Factories
  //===========================================================================

  // Make a TouchesConstructINode<Dim> object for use in intersections.
  // (This function helps to hide details about the GlobalIDDataBase from
  // higher level functions.)

  inline
  TouchesConstructINode<Dim> touchesConstructINode(LayoutID_t layoutID)
  {
    return TouchesConstructINode<Dim>(layoutID, key_m, globalIDDataBase_m);
  }

  template<int Dim2>
  inline static
  TouchesConstructINode<Dim> touchesConstructINode(LayoutID_t layoutID,
						   const INode<Dim2> &inode)
  {
    return TouchesConstructINode<Dim>(layoutID, inode.key(),
				      inode.globalIDDataBase());
  }

  //============================================================
  // I/O
  //============================================================

  // output an INode to an output stream, in the format
  //   {domain: key=value}

  template<class Out>
  void print(Out &o) const 
  {
    o << "{" << domain();
    o << ": key=" << key();
    o << "}";
  }

private:

  //===========================================================================
  // Data members
  //===========================================================================
  
  Domain_t          domain_m;

  GlobalIDDataBase *globalIDDataBase_m;
  NodeKey_t         key_m;
};


//-----------------------------------------------------------------------------
//
// Add a Loc to an INode, returning a new one with the adjusted domain.
//
//-----------------------------------------------------------------------------

template<int Dim>
inline INode<Dim> operator+(const INode<Dim> &inode, const Loc<Dim> &loc)
{
  return INode<Dim>(inode, inode.domain() + loc);
}


//-----------------------------------------------------------------------------
//
// A specialization of the Inform traits used to say that node has
// a print method.
//
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &o, const INode<Dim> &inode)
{
  inode.print(o);
  return o;
}


//-----------------------------------------------------------------------------
//
// We want to use an INode in some places as a domain, even though it isn't.
// Define a version of contains and a poor man's DomainTraits to support this.
// Same with TemporaryNewDomain1.
//
//-----------------------------------------------------------------------------

template<int Dim>
inline bool contains(const Interval<Dim> &i, const INode<Dim> &n)
{
  return contains(i, n.domain());
}

template<int Dim>
struct DomainTraits<INode<Dim> >
{
  enum { singleValued = 0 };
};

template<class Domain, class Sub>
struct TemporaryNewDomain1;

template<class Domain, int N>
struct TemporaryNewDomain1<Domain, INode<N> > 
{
  typedef INode<N> SliceType_t;
  static inline
  const SliceType_t &combineSlice(const Domain &, const INode<N> &i)
  {
    return i;
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

// Build INode object

template<class Domain, int Dim>
inline INode<Dim>
touchesConstruct(const Domain &d,
		 int, int context, int gid, int,
		 const TouchesConstructINode<Dim> &tcin)
{
  return INode<Dim>(d, context, gid, tcin);
}

template<class Domain, class AllocatedDomain, int Dim>
inline INode<Dim>
touchesConstruct(const Domain &d, const AllocatedDomain &,
		 int, int context, int gid, int,
		 const TouchesConstructINode<Dim> & tcin)
{
  return INode<Dim>(d, context, gid, tcin);
}



#endif // POOMA_LAYOUT_INODE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: INode.h,v $   $Author: richard $
// $Revision: 1.25 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
