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

#ifndef POOMA_DOMAIN_DOMAIN_MAP_H
#define POOMA_DOMAIN_DOMAIN_MAP_H

//-----------------------------------------------------------------------------
// Classes:
// DomainMap<Domain,Data>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Domain
 * @brief
 * DomainMap<Domain,Data> stores a list of N domains (of type Domain), each
 * with an associated piece of data (of type Data). 
 *
 * The domains are stored
 * in a tree-like structure which makes it efficient to find all subdomains
 * which touch a given domain.  Using a tree, the touch operation can be done
 * in O(log(N)) time instead of O(N), since the domains are sorted.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Domain/Split.h"
#include "Domain/Contains.h"
#include "Domain/Touches.h"
#include "Utilities/Pooled.h"
#include "Utilities/PAssert.h"
#include <utility>
#include <list>
#include <iosfwd>


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {


/**
 * DomainMapNode is a node in a tree, where each Node has a list of domains
 * and a left and right branch.  This class is pooled since it will be
 * created and deleted often.  It contains methods for finding the left
 * and right nodes in its leaves.
 */

template<class Dom, class T>
class DomainMapNode : public Pooled<DomainMapNode<Dom,T> >
{
public:
  // Typedefs
  typedef Dom                        Domain_t;
  typedef T                          Data_t;
  typedef std::pair<Domain_t,Data_t> Value_t;
  typedef std::list<Value_t>         List_t;
  typedef DomainMapNode<Dom,T>       Node_t;
  typedef typename List_t::iterator  iterator;

  // Constructor: tell this node its domain and parent
  DomainMapNode(const Domain_t &d, Node_t *p = 0)
    : domain_m(d), left_m(0), right_m(0), parent_m(p) {
  }

  // Destructor: delete the left and right leaves as well
  ~DomainMapNode() {
    if (left_m != 0)
      delete left_m;
    if (right_m != 0)
      delete right_m;
  }

  // Return the domain of this node
  const Domain_t &domain() const { return domain_m; }

  // Return begin/end iterators for our list of values
  iterator begin() { return list_m.begin(); }
  iterator end() { return list_m.end(); }

  // Insert a Value_t object into the node's list, or into its leaves' list.
  void insert(const Value_t &v) {
    
    // We must have it that the inserted domain is contained in our
    // current domain.

    PAssert(contains(domain_m, v.first));

    // Make sure we have left and right branches, even though they'll
    // initially be empty.  Do this by splitting our current domain.

    if (left_m == 0) 
      {
	Domain_t leftdom, rightdom;
	split(domain_m, leftdom, rightdom);
	left_m  = new Node_t(leftdom, this);
	right_m = new Node_t(rightdom, this);
      }

    // Now figure out which one this goes in, otherwise keep it here if it
    // straddles both left and right or if our size has gone down to 1.

    if (contains(v.first, domain_m))
      {
	list_m.push_back(v);
      }
    else if (contains(left_m->domain_m, v.first))
      {
	left_m->insert(v);
      }
    else if (contains(right_m->domain_m, v.first))
      {
	right_m->insert(v);
      }
    else
      {
	list_m.push_back(v);
      }
  }

  // Get the next node which is to the right of this one

  Node_t *nextRightNode() {
    Node_t *y, *p = this;
    if (p->right_m != 0) {
      // a right node is available ... go there, and then all the way left
      p = p->right_m;
      while (p->left_m != 0)
        p = p->left_m;
    } else {
      // there is no right, so go up until we can go right
      for (y=p->parent_m; y != 0 && p == y->right_m; y=y->parent_m)
        p = y;
      p = y;
    }

    // return the node we've found, which may be a null pointer
    return p;
  }

  // Get the next node which is to the right of this one, touching the domain

  Node_t *nextRightTouchNode(const Domain_t &d) {
    Node_t *p = this;
    Node_t *y = right_m;
    if (y != 0 && touches(d, y->domain_m)) {
      // The right side exists and touches this domain, so try it
      p = y;
      for (y=y->left_m; y != 0 && touches(d, y->domain_m); y=y->left_m)
        p = y;
    } else {
      // There is no right, so go up until we can go right
      // No need to test for touching on the way because we wouldn't be here
      // if the parent didn't touch
      for (y=p->parent_m; y != 0 && p == y->right_m; y=y->parent_m)
        p = y;
      p = y;
    }

    // return the node we've found, which may be a null pointer
    return p;
  }

  // Get the leftmost non-empty node

  Node_t *findLeftNode() {
    Node_t *p = this;

    // first, go as far left as we can
    while (p->left_m != 0)
      p = p->left_m;

    // Then, check if it is empty.  If it is, move to the right until
    // we find a non-empty node.
    while (p != 0 && p->begin() == p->end())
      p = p->nextRightNode();

    // return the node we've found, which may be a null pointer
    return p;
  }

  // Get the leftmost non-empty node which touches the given domain

  Node_t *findLeftTouchNode(const Domain_t &d) {
    Node_t *y, *p = this;
    for (y=p->left_m; y != 0 && touches(d, y->domain_m); y=y->left_m)
      p = y;
    return p;
  }

private:
  // This node's domain
  Domain_t domain_m;

  // left, right, and parent nodes
  Node_t *left_m, *right_m, *parent_m;

  // The list of values in this node.
  List_t list_m;
};


/**
 * An iterator for a DomainMap.  This has forward-iterator semantics.  It
 * is initially given a starting node and location in that node's list of
 * elements; it will iterate through the elements in the node, and then move
 * on to the next node until there are no nodes left.  When it reaches the
 * end, it sets node pointer to 0.
 */

template<class Dom, class T>
class DomainMapIterator
{
public:
  // Typedefs
  typedef Dom                       Domain_t;
  typedef DomainMapNode<Dom,T>      Node_t;
  typedef typename Node_t::Data_t   Value_t;
  typedef typename Node_t::iterator NodeIter_t;

  // Default constructor: initialize pointers to zero
  DomainMapIterator() : node_m(0) { }

  // Constructor: initialize with the node to start iterating from,
  // and where in the node's list to start
  DomainMapIterator(Node_t *n, NodeIter_t i) : node_m(n), iter_m(i) { }

  // Destructor: does nothing
  ~DomainMapIterator() { }

  // Return the node and iterator values
  Node_t *getNode() const { return node_m; }
  NodeIter_t getIter() const { return iter_m; }

  // Comparison operators
  bool operator==(const DomainMapIterator<Dom,T> &dmi) {
    return (node_m == dmi.node_m && (node_m == 0 || iter_m == dmi.iter_m));
  }
  bool operator!=(const DomainMapIterator<Dom,T> &dmi) {
    return !(*this == dmi);
  }

  // Dereference the iterator.  Return by reference since this is a
  // non-const iterator.
  Value_t &operator*() {
    PAssert(node_m != 0);
    return (*iter_m).second;
  }

  // Return the domain of the current iterator position.  Return by ref
  // since this is a non-const iterator.
  Domain_t &domain() {
    PAssert(node_m != 0);
    return (*iter_m).first;
  }

  // Increment the iterator.
  DomainMapIterator<Dom,T> &operator++() {
    PAssert(node_m != 0);
 
    // try to increment while in the current node; if this hits the
    // end, move to next available node
    if ((++iter_m) == node_m->end()) {
      do {
        node_m = node_m->nextRightNode();
      } while (node_m != 0 && (iter_m=node_m->begin()) == node_m->end());
    }

    // op++ must return a ref to the iterator
    return *this;
  }

private:
  // The current node.
  Node_t *node_m;

  // Where in the current node we're pointing.
  NodeIter_t iter_m;
};


/**
 * An iterator for a DomainMap.  This has forward-iterator semantics.  It
 * is initially given a starting node and location in that node's list of
 * elements; it will iterate through the elements in the node, and then move
 * on to the next node until there are no nodes left.  When it reaches the
 * end, it sets node pointer to 0.
 *
 * This is the const version of the iterator, so that the deref operator
 * returns a copy of instead of a reference to the data.
 */

template<class Dom, class T>
class DomainMapConstIterator
{
public:
  // Typedefs
  typedef Dom                       Domain_t;
  typedef DomainMapNode<Dom,T>      Node_t;
  typedef typename Node_t::Data_t   Value_t;
  typedef typename Node_t::iterator NodeIter_t;

  // Default constructor: initialize pointers to zero
  DomainMapConstIterator() : node_m(0) { }

  // Constructor: initialize with the node to start iterating from,
  // and where in the node's list to start
  DomainMapConstIterator(Node_t *n, NodeIter_t i)
    : node_m(n), iter_m(i) { }

  // Constructor: initialize from a non-const iterator
  DomainMapConstIterator(const DomainMapIterator<Dom,T> &dmi)
    : node_m(dmi.getNode()), iter_m(dmi.getIter()) { }

  // Destructor: does nothing
  ~DomainMapConstIterator() { }

  // Comparison operators
  bool operator==(const DomainMapConstIterator<Dom,T> &dmi) {
    return (node_m == dmi.node_m && (node_m == 0 || iter_m == dmi.iter_m));
  }
  bool operator!=(const DomainMapConstIterator<Dom,T> &dmi) {
    return !(*this == dmi);
  }

  // Dereference the iterator.  Return by value since this is a const iter.
  Value_t operator*() {
    PAssert(node_m != 0);
    return (*iter_m).second;
  }

  // Return the domain of the current iterator position.  Return by value
  // since this is a const iterator.
  Domain_t domain() {
    PAssert(node_m != 0);
    return (*iter_m).first;
  }

  // Increment the iterator.
  DomainMapConstIterator<Dom,T> &operator++() {
    PAssert(node_m != 0);
 
    // try to increment while in the current node; if this hits the
    // end, move to next available node
    if ((++iter_m) == node_m->end()) {
      do {
        node_m = node_m->nextRightNode();
      } while (node_m != 0 && (iter_m=node_m->begin()) == node_m->end());
    }

    // op++ must return a ref to the iterator
    return *this;
  }

private:
  // The current node.
  Node_t *node_m;

  // Where in the current node we're pointing.
  NodeIter_t iter_m;
};


/**
 * The touch iterator for a DomainMap.  This has forward-iterator semantics.
 * This is similar to the regular DomainMapIterator, except that it only
 * returns domains which touch a given domain.  There is no const version
 * of this class.
 */

template<class Dom, class T>
class DomainMapTouchIterator
{
public:
  // Typedefs
  typedef Dom                       Domain_t;
  typedef DomainMapNode<Dom,T>      Node_t;
  typedef typename Node_t::Data_t   Value_t;
  typedef typename Node_t::iterator NodeIter_t;

  // Default constructor: initialize pointers to zero
  DomainMapTouchIterator() : node_m(0) { }

  // Constructor: initialize with the node to start iterating from,
  // and where in the node's list to start; also need the touch domain
  DomainMapTouchIterator(Node_t *n, NodeIter_t i, const Domain_t &d)
    : node_m(n), iter_m(i), domain_m(d) { }

  // Destructor: does nothing
  ~DomainMapTouchIterator() { }

  // Comparison operators
  bool operator==(const DomainMapTouchIterator<Dom,T> &dmi) {
    return (node_m == dmi.node_m && (node_m == 0 || iter_m == dmi.iter_m));
  }
  bool operator!=(const DomainMapTouchIterator<Dom,T> &dmi) {
    return !(*this == dmi);
  }

  // Dereference the iterator.  Return by reference since this is a
  // non-const iterator.
  Value_t &operator*() {
    PAssert(node_m != 0);
    return (*iter_m).second;
  }

  // Return the domain of the current iterator position.  Return by ref
  // since this is a non-const iterator.
  Domain_t &domain() {
    PAssert(node_m != 0);
    return (*iter_m).first;
  }

  // Increment the iterator.
  DomainMapTouchIterator<Dom,T> &operator++() {
    PAssert(node_m != 0);

    // try to increment while in the current node; if this hits the
    // end, move to next available node
    while ((++iter_m) != node_m->end()) {
      if (touches(domain_m, (*iter_m).first))
        return *this;
    }

    // we reached the end of the current node ... try to find the next one
    do {
      if ((node_m = node_m->nextRightTouchNode(domain_m)) != 0)
        for (iter_m = node_m->begin(); iter_m != node_m->end(); ++iter_m)
          if (touches(domain_m, (*iter_m).first))
            return *this;
    } while (node_m != 0);
 
    return *this;
  }

private:
  // The current node.
  Node_t *node_m;

  // Where in the current node we're pointing.
  NodeIter_t iter_m;

  // The touch domain we are checking
  Domain_t domain_m;
};


/**
 * DomainMap<Domain,Data> is templated on the type of domains it is storing,
 * and the Data type it stores for each Domain.  The purpose of DomainMap
 * is to store a set of N domains in a way that is very fast for 'touches'
 * operations.  This operation is done when you want to find out which
 * subdomains in a list happen to touch a given domain.  This domain is
 * typically the extent for an expression, but it could be anything.
 *
 * DomainMap maintains a binary tree of domains, where each node in the tree
 * is of type 'DomainMapNode' and stores the following information:
 *   -# The domain for that node.  This is a section of the total domain,
 *      which is obtained by splitting the domain of the parent node.  The
 *      root node has a domain equal to the total domain.  Under this are
 *      two nodes with the parent domain split in two, and so on.
 *   -# A list of domains which are part of that node.  When a subdomain
 *      is inserted into the DomainMap, it is inserted into the root node,
 *      which checks to see if the subdomain is contained by the left or
 *      right split domains.  If it is, the subdomain is inserted in the
 *      left or right.  But if it spans both left and right, it is inserted
 *      in the current node's list.
 *
 * A DomainMap is constructed either with a default constructor, or with
 * a global domain which should represent the "bounding box" of the DomainMap.
 * Subsequent insertions of subdomains into the DomainMap should be for
 * subdomains contained within the bounding box domain.  DomainMap keeps
 * the root node for its tree of DomainMapNode's, and a count of how many
 * domains have been inserted.  If the default constructor is used, the
 * 'initialize(domain)' method must be called before the DomainMap can be
 * used in any other way.
 *
 * For each subdomain inserted into the DomainMap, there is a data element
 * of type 'Data', where Data is the second template parameter for DomainMap.
 * For example, for Layout objects, Data is an int that stores the context
 * number.  The elements stored in the lists in each DomainMapNode are actually
 * of type pair<Domain,Data>; this typedef'd to DomainMap<Domain,Data>::Value_t
 * and it is elements of type Value_t which are inserted:
 *   DomainMap<...> dmap;
 *   dmap.insert(DomainMap<...>::Value_t(domain, context));
 *
 * After a number of elements have been inserted, the user should call
 * DomainMap.update(), which resets an internal pointer in DomainMap to point
 * to the leftmost node.  If update() is not called after an insertion, then
 * the 'touch' method will not function properly.  However, you can peform
 * multiple insert() operations between calls to update without a problem.
 * The typical method is to create a DomainMap, insert all the subdomains in
 * some kind of loop, and then call update() after everything has been
 * inserted.  update() can be called any number of times.  DomainMap could
 * have done an implicit call to update() at the end of each insert(), but
 * this is inefficient and in most cases unnecessary.
 *
 * The elements of DomainMap can be iterated over, using begin() and end()
 * methods.  DomainMap has a size() method as well.  The iterators are
 * of type DomainMap::iterator and DomainMap::const_iterator.  These
 * iterators have forward-iterator semantics only; dereferencing an iterator
 * returns an item of type DomainMap::Value_t, which is a pair<Domain,Data>.
 * Elements will be iterated over from "left node" to "right node".  The
 * DomainMapIterator and DomainMapConstIterator classes are used to implement
 * these iterators.
 *
 * Finally, the key use of DomainMap is to perform a 'touch' operation.  The
 * touch(domain) method returns a pair of iterators, as
 *   pair<DomainMap::touch_iterator,DomainMap::touch_iterator>
 * which is typedef'd as DomainMap::Touch_t.  This pair of iterators is
 * a begin/end pair which can be used to iterate through all subdomains which
 * touch the domain given to the 'touch' method.  touch_iterator has
 * forward-iterator semantics, and dereferencing returns a Value_t pair.
 */

template<class Dom, class T>
class DomainMap
{
public:
  //============================================================
  // Typedefs
  //============================================================

  typedef Dom                                      Domain_t;
  typedef Dom                                      key_type;
  typedef T                                        Data_t;
  typedef T                                        mapped_type;
  typedef std::pair<Domain_t,Data_t>               Value_t;
  typedef std::pair<Domain_t,Data_t>               value_type;
  typedef DomainMapIterator<Domain_t,Data_t>       iterator;
  typedef DomainMapConstIterator<Domain_t,Data_t>  const_iterator;
  typedef DomainMapTouchIterator<Domain_t,Data_t>  touch_iterator;
  typedef std::pair<touch_iterator,touch_iterator> Touch_t;
  typedef std::pair<touch_iterator,touch_iterator> touch_type;
  typedef long                                     Size_t;
  typedef long                                     size_type;

  //============================================================
  // Constructors
  //============================================================

  // Default constructor for DomainMap: if this is used, the initialize()
  // method should be used before using this DomainMap in any way.
  DomainMap() : size_m(0), root_m(0) { }

  // Create a DomainMap with a bounding-box domain.  All domains inserted
  // into this DomainMap should be contained within this bounding box.
  DomainMap(const Domain_t &d) : size_m(0), root_m(0) {
    initialize(d);
  }

  // Perform initialization, which creates a root node with the given
  // bounding box.  The user still needs to call update() after inserting
  // values.
  void initialize(const Domain_t &d) {
    PAssert(root_m == 0 && size_m == 0);
    root_m = new Node_t(d);
  }

  //============================================================
  // Destructor
  //============================================================

  // Delete the root node, which takes out all the others.
  ~DomainMap() {
    if (root_m != 0)
      delete root_m;
  }

  //============================================================
  // Accessors
  //============================================================

  // return iterators to the domains
  iterator begin() { return left_m; }
  iterator end() { return iterator(); }
  const_iterator begin() const { return const_iterator(left_m); }
  const_iterator end() const { return const_iterator(); }

  // return our size
  Size_t size() const { return size_m; }

  // return a pair of iterators which describe those domains which touch
  // the given domain
  Touch_t touch(const Domain_t &d) const {
    Node_t *p = root_m;

    if (p != 0)
    {
      // First dive left, checking touches
      p = p->findLeftTouchNode(d);

      // Now look for a node which has an element
      do
      {
	// check current node for a touching element
	for (NodeIter_t a = p->begin(); a != p->end(); ++a)
	  if (touches(d, (*a).first))
	    return Touch_t(touch_iterator(p, a, d), touch_iterator());

	// if none found, move on to next node
	p = p->nextRightTouchNode(d);
      } while (p != 0);
    }

    // if we're here, we didn't find anything, so return empty iterators
    return Touch_t(touch_iterator(), touch_iterator());
  }

  //============================================================
  // Modifiers
  //============================================================

  // insert a new element into the DomainMap.  This does NOT update the
  // pointer to the leftmost element; the user is responsible for doing that
  // by calling update() when the insertions are complete.
  void insert(const Value_t &v) {
    PAssert(root_m != 0);
    root_m->insert(v);
    size_m++;
  }

  // update this DomainMap's leftmost-element pointer.  If this is not
  // done between when a domain is inserted and when a touch() operation
  // is performed, the results can be inaccurate.  This could be done
  // automatically after each insert() operation, but that would be a little
  // inefficient.
  void update() {
    Node_t *leftnode;
    if (root_m != 0 && size_m > 0 && (leftnode=root_m->findLeftNode()) != 0)
      left_m = iterator(leftnode, leftnode->begin());
    else
      left_m = iterator();
  }

  // clear out our current domain list; just leave the empty root
  // node intact
  void clear() {
    if (root_m != 0) {
      Node_t *newroot = new Node_t(root_m->domain());
      delete root_m;
      root_m = newroot;
      size_m = 0;
      update();
    }
  }

  // when you want to start over ... resets everything to the initial
  // state when constructed.

  void zap() {
    if (root_m != 0)
      delete root_m;
    size_m = 0;
    root_m = 0;
  }

  //============================================================
  // I/O
  //============================================================

  // output a DomainMap to an output stream, by doing a touches
  // operation on the entire global domain.

  template<class Out>
  void print(Out &o) const {
    if (size() < 1 || root_m == 0)
      {
	o << "DomainMap: empty.";
      }
    else
      {
	o << "DomainMap: Total domain = " << root_m->domain();
	o << ", touching domains:\n";
	Touch_t touchiters = touch(root_m->domain());
	while (touchiters.first != touchiters.second)
	  {
	    o << "  " << touchiters.first.domain() << " ==> ";
	    o << *(touchiters.first) << "\n";
	    ++(touchiters.first);
	  }
      }
  }

  void print() const {
    if (size() < 1 || root_m == 0)
      {
	std::cout << "DomainMap: empty.";
      }
    else
      {
	std::cout << "DomainMap: Total domain = " << root_m->domain();
	std::cout << ", touching domains:\n";
	Touch_t touchiters = touch(root_m->domain());
	while (touchiters.first != touchiters.second)
	  {
	    std::cout << "  " << touchiters.first.domain() << " ==> ";
	    std::cout << *(touchiters.first) << "\n";
	    ++(touchiters.first);
	  }
      }
  }



private:
  // A typedef for a node in our tree, and an iterator in that node
  typedef DomainMapNode<Domain_t,Data_t> Node_t;
  typedef typename Node_t::iterator      NodeIter_t;

  // The number of elements in our list
  Size_t size_m;

  // The root node for our tree
  Node_t *root_m;

  // An iterator pointing to the leftmost node in the tree
  iterator left_m;

  // Disable the copy constructor and operator= for now
  DomainMap(const DomainMap<Domain_t,Data_t> &);
  DomainMap<Domain_t,Data_t> &operator=(const DomainMap<Domain_t,Data_t> &);
};


/// A specialization of the std::ostream output operator to say that
/// DomainMap has a print method.

template <class Dom, class T>
std::ostream &operator<<(std::ostream &o, const DomainMap<Dom, T> &dmap)
{
  dmap.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_MAP_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainMap.h,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
