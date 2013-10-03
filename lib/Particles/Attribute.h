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

#ifndef POOMA_PARTICLES_ATTRIBUTE_H
#define POOMA_PARTICLES_ATTRIBUTE_H

//-----------------------------------------------------------------------------
// Classes:
//   Attribute
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Particles
 * @brief
 * Attribute is a non-templated base class used to provide an interface
 * to DynamicArray objects used as attributes in Particle classes.
 *
 * The dynamic operations such as create, destroy, copy, etc., are actually
 * performed via requests to a layout object that each DynamicArray will
 * use, but this class defines a small set of virtual/nonvirtual functions to:
 * print the DynamicArray contents to a stream (for debugging)
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include <iosfwd>

//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------

template <class T> class IndirectionList;
template <int Dim> class Interval;

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * Attribute is an abstract base class used to implement an external
 * polymorphism interface to DynamicArray objects.  The AttributeWrapper
 * subclass is templated on the type of Array (really, DynamicArray) that
 * the user wants to provide an abstract interface to.  This is used to
 * let users create heterogenous collections of DynamicArray's, and to
 * perform common tasks on all of them such as print.
 */

class Attribute
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  //============================================================
  // Constructors
  //============================================================

  /// Default constructor.

  Attribute()
    {
    }

  /// Copy constructor.

  Attribute(const Attribute &)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  /// Attribute needs a virtual destructor, since we will be deleting
  /// Attribute's from a base class pointer.

  virtual ~Attribute()
    {
    }


  //============================================================
  // Public virtual interface
  //============================================================

  /// Print the contents of the Array to the given stream.

  virtual void print(std::ostream &) const = 0;

  /* Omit this until we have serialize/deserialize functions for DynamicArray.

  /// serialize/deserialize the Array using the given stream.

  virtual int serialize(std::ostream &) const = 0;
  virtual int serialize(std::fstream &) const = 0;
  virtual int deserialize(std::ostream &) = 0;
  virtual int deserialize(std::fstream &) = 0;

  */

#if POOMA_MESSAGING

  /// packSize, pack and unpack function interface for particle swapping

  virtual int packSize(int) const = 0;
  virtual int pack(int, const IndirectionList<int> &, char *) const = 0;
  virtual int unpack(int, const Interval<1> &, char *) = 0;

#endif // POOMA_MESSAGING

};


//-----------------------------------------------------------------------------
//
/// A specialization of the Inform traits used to say that Attribute has
/// a print method.
//
//-----------------------------------------------------------------------------

inline
std::ostream&
operator<<(std::ostream& o, const Attribute& attrib)
{
  attrib.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_ATTRIBUTE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Attribute.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
