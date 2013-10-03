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
//   DomainLayout<int>
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_DOMAIN_LAYOUT_H
#define POOMA_LAYOUT_DOMAIN_LAYOUT_H

/** @file
 * @ingroup Layout
 * @brief
 * DomainLayout<int>:
 * Layout class that just wraps around an Interval.  It provides the std
 * layout interface, but is not shareable or have a dynamic interface.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Intersect.h"
#include "Layout/DynamicEvents.h"
#include "Layout/GuardLayers.h"
#include "Layout/Node.h"
#include "Layout/TouchesConstruct.h"
#include "Pooma/Pooma.h"
#include "Utilities/DerefIterator.h"
#include "Utilities/PAssert.h"

#include <iosfwd>

// namespace Pooma {


/**
 * DomainTag class.
 */

struct DomainTag { };


/**
 * DomainLayout is used to act as a layout object for all engines that
 * just need a simple layout for a single domain of dimension Dim.  It is
 * not shareable, and it does not provide a dynamic create/destroy interface.
 * It is possible to change the domain of this layout, by calling
 * initialize.
 *
 * DomainLayout provides the same interface for begin/end iterators over
 * Nodes as all other layouts.  It will appear to have just one local
 * patch, and no remote patches.
 */

template<int Dim>
class DomainLayout
{
public:

  //============================================================
  // Typedefs and enumerations
  //============================================================

  // General public typedefs.

  typedef DomainLayout<Dim>                  This_t;
  typedef Interval<Dim>                      Domain_t;
  typedef Node<Domain_t>                     Value_t;
  typedef DynamicEvents::PatchID_t           PatchID_t;
  typedef DynamicEvents::CreateSize_t        CreateSize_t;
  typedef GuardLayers<Dim>                   GuardLayers_t;

  // Iterator through our one node. Just a pointer.  
  
  typedef Value_t                            *iterator;
  typedef const Value_t                      *const_iterator;

  // General public enums.

  enum { dimensions = Dim };
  enum { dynamic = false };


  //============================================================
  // Constructor
  //============================================================

  // Default constructor: stores an empty domain.

  inline DomainLayout()
    {
    }

  // Construct a layout with a global domain.  This will make a new
  // data storage object, and make this object the first user.

  explicit DomainLayout(const Domain_t &dom)
    : node_m(0, dom, Pooma::context(), 0, 0)
    {
    }

  // Construct a layout with a global domain and guard layers.

  DomainLayout(const Domain_t &dom, const GuardLayers_t &g)
    : node_m(0, dom, grow(dom, g), Pooma::context(), 0, 0)
    {
    }

  // Construct a layout with a global domain and guard layers.

  explicit DomainLayout(const Value_t &node)
    : node_m(node)
    {
    }

  // Copy: Construct a layout object from another layout object.  This layout
  // object will share usage of the other layout's data object.

  DomainLayout(const This_t &layout)
    : node_m(layout.node_m)
    {
    }

  // Initialize this object with a new domain.

  void initialize(const Domain_t &dom)
    {
      node_m = Value_t(0, dom, Pooma::context(), 0, 0);
    }

  // Initialize this object with a new domain and guard layers.

  void initialize(const Domain_t &dom, const GuardLayers_t &g)
    {
      node_m = Value_t(0, dom,  grow(dom, g), Pooma::context(), 0, 0);
    }

  // Initialize this object with the settings from another layout.

  void initialize(const This_t &layout)
  {
    node_m = layout.node_m;
  }


  //============================================================
  // Destructor
  //============================================================

  // Nothing to see here, move along

  inline ~DomainLayout()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  // Return whether or not this layout has been initialized.

  inline bool initialized() const
    {
      return domain().initialized();
    }

  // d'th component of the lower left of the inner domain.

  inline int first(int d) const { return innerDomain()[d].first(); }

  // A reference to our node object

  inline Value_t &node()
    {
      return node_m;
    }

  inline const Value_t &node() const
    {
      return node_m;
    }

  // Number of blocks in each dimension.

  inline Loc<Dim> blocks() const { return Loc<Dim>(1); }

  // Return the global domain.

  inline const Domain_t &domain() const
    {
      return node_m.allocated();
    }

  // Return the global domain less the external guard layers.
  
  inline const Domain_t &innerDomain() const
  {
    return node_m.domain();
  }

  // Return the global allocated domain.

  inline const Domain_t &allocated() const
    {
      return node_m.allocated();
    }

  inline GuardLayers_t internalGuards() const
  {
    return GuardLayers_t(0);
  }

  inline GuardLayers_t externalGuards() const
  {
    GuardLayers_t gl;
    for (int i = 0; i < Dim; i++)
      {
        gl.lower(i) = node_m.domain()[i].min() - node_m.allocated()[i].min();
        gl.upper(i) = node_m.allocated()[i].max() - node_m.domain()[i].max();
      }
        
    return gl;
  }

  inline const Domain_t &domain(int i) const
  {
    PAssert(i==0);
    return node_m.allocated();
  }
  
 inline const Domain_t &ownedDomain(int i) const
  {
    PAssert(i==0);
    return node_m.domain();
  }

  // Return the allocated domain

  inline const Domain_t &allocatedDomain(int i) const
    {
      PAssert(i==0);
      return node_m.allocated();
    }

  // Compare to another Layout.  The layouts are the same if:
  //   1. They have the same base domain.

  template<class L>
  inline bool operator==(const L &layout) const
    {
      return (domain() == layout.domain());
    }

  // Compare for inequality.

  template<class L>
  inline bool operator!=(const L &layout) const
    {
      return !(*this == layout);
    }

  //============================================================
  // Iterators
  //============================================================

  // Return begin and end iterators for the list of all subdomains

  inline iterator begin()
    { 
      return &node_m;
    }
  inline iterator end()
    { 
      return &node_m + 1;
    }
  inline const_iterator begin() const
    { 
      return &node_m;
    }
  inline const_iterator end() const
    { 
      return &node_m + 1;
    }

  // Return the total number of nodes (patches)

  inline long size() const
    { 
      return 1;
    }

  // Return begin and end iterators for the list of all local subdomains

  inline iterator beginLocal()
    { 
      return begin();
    }
  inline iterator endLocal()
    { 
      return end();
    }
  inline const_iterator beginLocal() const 
    { 
      return begin();
    }
  inline const_iterator endLocal() const
    { 
      return end();
    }
  inline long sizeLocal() const
    {
      return size();
    }

  // Return begin and end iterators for the list of all global subdomains

  inline iterator beginGlobal()
    { 
      return begin();
    }
  inline iterator endGlobal()
    { 
      return end();
    }
  inline const_iterator beginGlobal() const 
    { 
      return begin();
    }
  inline const_iterator endGlobal() const
    { 
      return end();
    }
  inline long sizeGlobal() const
    {
      return size();
    }

  // Return begin and end iterators for the list of all remote subdomains

  inline iterator beginRemote()
    {
      return iterator(0);
    }
  inline iterator endRemote()
    {
      return iterator(0);
    }
  inline const_iterator beginRemote() const
    {
      return const_iterator(0);
    }
  inline const_iterator endRemote() const
    {
      return const_iterator(0);
    }
  inline long sizeRemote() const
    {
      return 0;
    }

  // Accessors for getting the global ID of the patch containing a
  // particular element.  Here, all points should be in patch 0, if
  // they are not its an error.

  inline int globalID(const Loc<Dim> &loc) const
    {
      PAssert(contains(domain(), loc));
      return 0;
    }

  inline int globalID(int i1) const
    {
      return globalID(Loc<1>(i1));
    }

  inline int globalID(int i1, int i2) const
    {
      return globalID(Loc<2>(i1, i2));
    }

  inline int globalID(int i1, int i2, int i3) const
    {
      return globalID(Loc<3>(i1, i2, i3));
    }

  inline int globalID(int i1, int i2, int i3, int i4) const
    {
      return globalID(Loc<4>(i1, i2, i3, i4));
    }

  inline int globalID(int i1, int i2, int i3, int i4, int i5) const
    {
      return globalID(Loc<5>(i1, i2, i3, i4, i5));
    }

  inline int globalID(int i1, int i2, int i3, int i4, int i5,
		      int i6) const
    {
      return globalID(Loc<6>(i1, i2, i3, i4, i5, i6));
    }

  inline int globalID(int i1, int i2, int i3, int i4, int i5,
		      int i6, int i7) const
    {
      return globalID(Loc<7>(i1, i2, i3, i4, i5, i6, i7));
    }

  //============================================================
  // Touch methods
  //============================================================

  // Find all subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.  Return
  // the number of touching elements. This version of touches can build
  // either pointers or objects.

  template<class OtherDomain, class OutIter, class ConstructTag>
  int touches(const OtherDomain &d, OutIter o, ConstructTag ctag) const;

  // Find local subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.  Return
  // the number of touching elements. This version of touches can build
  // either pointers or objects.

  template<class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesLocal(const OtherDomain &d, OutIter o, 
                          const ConstructTag &ctag) const
    {
      return touches(d, o, ctag);
    }

  // Find remote subdomains that touch on a given domain, and insert the
  // intersection of these subdomains into the given output iterator.  Return
  // the number of touching elements. This version of touches can build
  // either pointers or objects.

  template<class OtherDomain, class OutIter, class ConstructTag>
  inline int touchesRemote(const OtherDomain &, OutIter,
			   const ConstructTag &) const
    {
      return 0;
    }
    
  // Find all/local/remote subdomains that touch on a given domain, 
  // and insert the intersection of these subdomains into the 
  // given output iterator.  Return the number of touching elements.
  // These versions of touches can build only objects (not pointers).

  template<class OtherDomain, class OutIter>
  inline int touches(const OtherDomain &d, OutIter o) const
    {
      return touches(d, o, TouchesConstructNodeObj());
    }

  template<class OtherDomain, class OutIter>
  inline int touchesLocal(const OtherDomain &d, OutIter o) const
    {
      return touchesLocal(d, o, TouchesConstructNodeObj());
    }

  template<class OtherDomain, class OutIter>
  inline int touchesRemote(const OtherDomain &d, OutIter o) const
    {
      return touchesRemote(d, o, TouchesConstructNodeObj());
    }

  //============================================================
  // I/O
  //============================================================

  template<class Out>
  void print(Out &o) const 
    {
      o << "DomainLayout: Node = " << node_m;
    }

private:

  //============================================================
  // Data
  //============================================================

  // The domain that we store

  Value_t node_m;
};

template<int Dim>
template<class OtherDomain, class OutIter, class ConstructTag>
int DomainLayout<Dim>::touches(const OtherDomain &d, OutIter o,
			       ConstructTag ctag) const
{
  int i, count = 0;

  // type of output elements

  typedef typename IntersectReturnType<Domain_t,OtherDomain>::Type_t 
    OutDomain_t;
  typedef Node<OutDomain_t> OutNode_t;

  // find the intersection of our domain and the given one

  OutDomain_t outDomain = intersect(d, domain());

  // add in touching domain if there is anything that intersects

  if (!outDomain.empty())
    {
      ++count;
      *o = touchesConstruct(outDomain,
			    node().affinity(),
			    node().context(),
			    node().globalID(),
			    node().localID(),
			    ctag);
    }

  // return the number of non-empty domains we found

  return count;
}


template <int Dim>
std::ostream &operator<<(std::ostream &o, const DomainLayout<Dim> &layout)
{
  layout.print(o);
  return o;
}


/** 
 * NewDomain1 traits classes for DomainLayout
 *
 * This is so an array can be initialized with a DomainLayout.
 */

template<int Dim>
struct NewDomain1< DomainLayout<Dim> >
{
  typedef DomainLayout<Dim> &Type_t;

  inline static Type_t combine(const DomainLayout<Dim> &a)
    {
      return const_cast<Type_t>(a);
    }
};

// } // namespace POOMA

#endif // POOMA_LAYOUT_DOMAIN_LAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainLayout.h,v $   $Author: richard $
// $Revision: 1.31 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
