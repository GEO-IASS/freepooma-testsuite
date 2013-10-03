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

#ifndef POOMA_PARTICLES_ATTRIBUTE_LIST_H
#define POOMA_PARTICLES_ATTRIBUTE_LIST_H

//-----------------------------------------------------------------------------
// Classes:
//   AttributeList
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Particles
 * @brief
 * AttributeList stores a list of Attribute pointers, that are used to
 * store a heterogenous collection of Attributes.
 *
 * When an AttributeList
 * is destroyed, it will delete all the Attributes it contains.  It
 * provides the same interface as Attribute, but it will loop over all
 * the Attributes it stores and perform these operations on each one in
 * turn.  You can add new Attributes to an AttributeList or delete
 * Attributes from a list by referring to the Attributes index.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/AttributeWrapper.h"
#include "Utilities/PAssert.h"
#include <vector>
#include <iosfwd>


//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//
// Full Description of AttributeList:
//
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

class AttributeList
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  typedef AttributeList              This_t;
  typedef std::vector<Attribute *>   List_t;
  typedef List_t::size_type Size_t;


  //============================================================
  // Constructors
  //============================================================

  // AttributeList just has a default constructor; it initially has
  // no attributes at all.

  AttributeList();


  //============================================================
  // Destructor
  //============================================================

  // AttributeList will delete all attributes that it owns when
  // it is deleted.

  ~AttributeList();


  //============================================================
  // Accessors
  //============================================================

  // Return the number of attributes

  inline Size_t size() const
    {
      return list_m.size();
    }

  // Return the Nth attribute

  inline Attribute *attribute(Size_t n)
    {
      PAssert(n < size());
      return list_m[n];
    }

  inline const Attribute *attribute(Size_t n) const
    {
      PAssert(n < size());
      return list_m[n];
    }


  //============================================================
  // Modifiers
  //============================================================

  // Add a new attribute to our list.  This will make an AttributeWrapper
  // instance that will wrap the item, and add it to our list.
  // Return the index of the newly added attribute.

  template<class T>
  inline Size_t add(T &item)
    {
      list_m.push_back(new AttributeWrapper<T>(item));
      return (size() - 1);
    }

  // Remove (and delete) the Nth attrib from our list.  This does not
  // delete the underlying wrapped object, just the Attribute wrapper/
  // container.
  // Return success.

  bool remove(Size_t);


  //============================================================
  // Methods to operate on all elements of the list.
  //============================================================

  // Print the contents of each Attribute to the given ostream.

  void print(std::ostream&) const;

private:

  // The list of Attribute pointers

  List_t list_m;

  // Make copy constructor and operator= private and undefined since
  // they should not be used

  AttributeList(const AttributeList &)
    {
      PInsist(false, "Called AttributeList copy constructor.");
    }

  const This_t &operator=(const This_t &)
    {
      PInsist(false, "Called AttributeList operator=.");
      return *this;
    }
};


// Declaration of operator<< for ostream and AttributeList

std::ostream& operator<<(std::ostream&, const AttributeList&);


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_ATTRIBUTE_LIST_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: AttributeList.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
