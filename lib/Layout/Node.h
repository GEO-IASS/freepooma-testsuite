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
// Node<Dom>
// ReturnNodePtr and ReturnNodeObj tags
//-----------------------------------------------------------------------------

#ifndef POOMA_LAYOUT_NODE_H
#define POOMA_LAYOUT_NODE_H

/** @file
 * @ingroup Layout
 * @brief
 * A simple class which stores the following information:
 *   -# Two domains (the class is templated on the domain type)
 *      specifying the owned and allocated domains for the patch.
 *   -# A context to which the domain has been assigned.
 *   -# A global ID value for the node.
 *   -# A local ID value for the node.
 *   -# A memory affinity value.
 * Layout objects store lists of nodes, and they are used elsewhere to
 * refer to the information about a subdomain block of a larger domain.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Domain/Contains.h"
#include "Domain/DomainTraits.h"
#include "Domain/Interval.h"
#include "Utilities/Pooled.h"
#include "Utilities/PAssert.h"
#include <iosfwd>

/**
 * Node is a quite simple object which stores information about where a
 * subdomain should be located.  It is templated on the type of domains it
 * is locating.  It stores domains describing the owned and allocated
 * domains of the patch, the latter being in reference to the actual 
 * underlying layout domain, the context where the patch should be stored,
 * and two ID values:
 *   -# a global ID, which, for a set of N subdomains partitioned across ALL
 *      the contexts, should be a unique value from 0 ... N-1
 *   -# a local ID, which, for a set of M subdomains all located on a single
 *      context, should be a unique value from 0 ... M-1.
 *
 * The entities that create Node objects must assign the ID values to the
 * Nodes.  This object is pooled for faster creation/destruction time, since
 * many Node objects are created and destroyed.
 */

template<class Dom, class OrigDom = Dom>
class Node : public Pooled<Node<Dom, OrigDom> >
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // The types for the domains, the context value, and ID values.
  
  typedef Dom      Domain_t;
  typedef OrigDom  AllocatedDomain_t;
  typedef int      Context_t;
  typedef int      ID_t;
  
  // For convenience:
  
  typedef Node<Dom,OrigDom> This_t;
  
  //============================================================
  // Constructors
  //============================================================

  // The default constructor initializing to empty domains.
  
  Node()
    : local_m(-1), global_m(0), context_m(0), affinity_m(-1)
  {
  }

  // The constructor for Node expects a domain, context, global ID, and
  // optionally the local ID. Versions are available with and without guard
  // cells.  If no local ID is given, -1 will be used, which means "this is
  //  not local".

  Node(const Domain_t &owned, const AllocatedDomain_t &allocated, 
       Context_t c, ID_t gid, ID_t lid = (-1))
    : domain_m(owned), allocated_m(allocated), 
      local_m(lid), global_m(gid),
      context_m(c), affinity_m(-1)
  {
    PAssert(owned.size() >= 0);
    PAssert(allocated.size() >= 0);
    PAssert(gid >= 0);
  }

  Node(const Domain_t &d, Context_t c, ID_t gid, ID_t lid = (-1))
    : domain_m(d), allocated_m(d),
      local_m(lid), global_m(gid),
      context_m(c), affinity_m(-1)
  {
    PAssert(d.size() >= 0);
    PAssert(gid >= 0);
  }

  // construct with affinity
  
  Node(int affinity, const Domain_t &owned, const AllocatedDomain_t &allocated, 
       Context_t c, ID_t gid, ID_t lid = (-1))
    : domain_m(owned), allocated_m(allocated), 
      local_m(lid), global_m(gid), context_m(c),
      affinity_m(affinity)
  {
    PAssert(owned.size() >= 0);
    PAssert(allocated.size() >= 0);
    PAssert(gid >= 0);
  }

  Node(int affinity, const Domain_t &d,
       Context_t c, ID_t gid, ID_t lid = (-1))
    : domain_m(d), allocated_m(d),
      local_m(lid), global_m(gid),
      context_m(c), affinity_m(affinity)
  {
    PAssert(d.size() >= 0);
    PAssert(gid >= 0);
  }

  // The copy constructor
  
  Node(const This_t &n)
    : domain_m(n.domain_m), allocated_m(n.allocated_m),
      local_m(n.local_m), global_m(n.global_m),
      context_m(n.context_m), affinity_m(n.affinity_m)
  {
  }

  // The copy constructor (from a different sort of Node)
  
  template<class ODom, class OAlloc>
  Node(const Node<ODom,OAlloc> &n)
    : domain_m(n.domain()), allocated_m(n.allocated()),
      local_m(n.localID()), global_m(n.globalID()),
      context_m(n.context()), affinity_m(n.affinity())
  {
  }

  // Initialize the node to the given set of values
  
  void initialize(const Domain_t &owned, const AllocatedDomain_t &allocated, 
                  Context_t c, ID_t gid, ID_t lid = (-1))
  {
    PAssert(owned.size() >= 0);
    PAssert(gid >= 0);
    domain_m = owned;
    allocated_m = allocated;
    context_m = c;
    local_m = lid;
    global_m = gid;
  }

  void initialize(const Domain_t &d, Context_t c, ID_t gid, ID_t lid = (-1))
  {
    PAssert(d.size() >= 0);
    PAssert(gid >= 0);
    domain_m = d;
    allocated_m = d;
    context_m = c;
    local_m = lid;
    global_m = gid;
  }

  //============================================================
  // Destructor
  //============================================================

  // The destructor here does nothing
  
  ~Node()
  {
  }


  //============================================================
  // Accessors
  //============================================================

  // Return the owned domain of this node
  
  const Domain_t &domain() const { return domain_m; }
  
  // Return the allocated domain for this node
  
  const AllocatedDomain_t &allocated() const { return allocated_m; }

  // Return the context value
  
  Context_t context() const { return context_m; }

  // Return the local and global ID values
  
  ID_t localID() const { return local_m; }
  ID_t globalID() const { return global_m; }

  // Return true if this node is local, false otherwise.  If it is local,
  // then localID() returns a value >= 0.
  
  bool isLocal() const { return (local_m >= 0); }

  // Return the memory affinity
  
  int affinity() const { return affinity_m; }

  //============================================================
  // Mutators
  //============================================================

  // Modify the affinity setting.
  
  int& affinity() { return affinity_m; }

  // Modify the context setting.
  
  int& context() { return context_m; }

  // Modify the local ID .

  int& localID() { return local_m; }

  // Modify the owned domain.

  void setDomain(const Domain_t &dom) { domain_m = dom; }
  Domain_t &domain() { return domain_m; }

  // Modify the allocated domain.

  void setAllocated(const AllocatedDomain_t &dom) { allocated_m = dom; }
  AllocatedDomain_t &allocated() { return allocated_m; }

  //============================================================
  // Operators
  //============================================================

  // assignment operator
  
  This_t &operator=(const This_t &n) 
  {
    domain_m = n.domain();
    allocated_m = n.allocated();
    context_m = n.context();
    local_m = n.localID();
    global_m = n.globalID();
    affinity_m = n.affinity();
    return *this;
  }

  // assignment operator (from another sort of node)
  
  template<class ODom, class OAlloc>
  This_t &operator=(const Node<ODom,OAlloc> &n) 
  {
    domain_m = n.domain();
    allocated_m = n.allocated();
    context_m = n.context();
    local_m = n.localID();
    global_m = n.globalID();
    affinity_m = n.affinity();
    return *this;
  }

  // assignment operator (from just a domain).  Does not affect
  // anything else. (Should this also set allocated_m??? Where 
  // is this used???)
  
  This_t &operator=(const Domain_t &d) 
  {
    domain_m = d;
    return *this;
  }

  //============================================================
  // I/O
  //============================================================

  // output a Node to an output stream, in the format
  //   {domain: allocated=dom, con=#, aff=#, gid=#, lid=#}

  template<class Out>
  void print(Out &o) const 
  {
    o << "{" << domain();
    o << ": allocated=" << allocated();
    o << ", con=" << context();
    o << ", aff=" << affinity();
    o << ", gid=" << globalID();
    o << ", lid=" << localID();
    o << "}";
  }

private:
  
  enum { dim = DomainTraits<Dom>::dimensions };
  enum { origDim = DomainTraits<OrigDom>::dimensions };
  
  // The owned domain for this node
  
  Domain_t domain_m;

  // The allocated domain for this node.
  
  AllocatedDomain_t allocated_m;
  
  // The local and global ID values
  
  ID_t local_m;
  ID_t global_m;

  // The context for this node
  
  Context_t context_m;

  // The memory affinity value.
  
  int affinity_m;
};


//-----------------------------------------------------------------------------
//
// A specialization of the Inform traits used to say that node has
// a print method.
//
//-----------------------------------------------------------------------------

template <class D, class A>
std::ostream &operator<<(std::ostream &o, const Node<D,A> &node)
{
  node.print(o);
  return o;
}


//-----------------------------------------------------------------------------
//
// We want to use an INode in some places as a domain, even though it isn't.
// Define a version of contains along with a poor man's DomainTraits to 
// support this. Same with TemporaryNewDomain1.
//
//-----------------------------------------------------------------------------

template<int Dim, class Dom, class OrigDom>
inline bool contains(const Interval<Dim> &i, const Node<Dom, OrigDom> &n)
{
  return contains(i, n.domain());
}

template<class Dom, class OrigDom>
struct DomainTraits<Node<Dom, OrigDom> >
{
  enum { singleValued = 0 };
};

template<class Domain, class Sub>
struct TemporaryNewDomain1;

template<class Domain, class OwnedDomain, class AllocatedDomain>
struct TemporaryNewDomain1<Domain, Node<OwnedDomain, AllocatedDomain> > 
{
  typedef Node<OwnedDomain,AllocatedDomain> SliceType_t;
  static inline
  const SliceType_t &combineSlice(const Domain &, 
    const Node<OwnedDomain,AllocatedDomain> &n)
  {
    return n;
  }
};


#endif     // POOMA_LAYOUT_NODE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Node.h,v $   $Author: richard $
// $Revision: 1.39 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
