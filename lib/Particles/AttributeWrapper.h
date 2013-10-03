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

#ifndef POOMA_PARTICLES_ATTRIBUTE_WRAPPER_H
#define POOMA_PARTICLES_ATTRIBUTE_WRAPPER_H

//-----------------------------------------------------------------------------
// Classes:
//   AttributeWrapper<T>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Particles
 * @brief
 * AttributeWrapper<T> is a subclass of Attribute that implements
 * the basic Attribute interface by passing on the operations in the
 * interface to an object of type T that AttributeWrapper wraps.
 *
 * This is basically a standard external polymorphism mechanism for objects of
 * various types, for example for wrapping DynamicArray objects.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/Attribute.h"
#include "Utilities/Inform.h"
#include "Utilities/PAssert.h"

#if POOMA_MESSAGING
#include "Tulip/Messaging.h"
#endif

#include <iostream>
#include <fstream>

//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * AttributeWrapper<T> is a subclass of the abstract base class Attribute.
 * It is templated on a class T, where T should be some form of Array or
 * DynamicArray that supports a dynamic data structure interface.
 *
 * AttributeWrapper is meant to be used as an external polymorphism derived
 * class.  You create an AttributeWrapper and give it an object to wrap around;
 * the abstract base class is used to provide an abstract interface to a
 * heterogenous collection of AttributeWrappers from some other user, for
 * example a Particles class.  Particles actually uses an AttributeList object
 * to hold a collection of Attributes.
 */

template<class T>
class AttributeWrapper : public Attribute
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  typedef Attribute            Base_t;
  typedef AttributeWrapper<T>  This_t;
  typedef T                    Wrapped_t;


  //============================================================
  // Constructors
  //============================================================

  // AttributeWrapper should be initialized with an object of type T
  // that it will wrap.  It will store a reference to this object.

  AttributeWrapper(Wrapped_t &array)
    : wrapped_m(array)
    {
    }


  //============================================================
  // Destructor
  //============================================================

  // AttributeWrapper does not need to do anything in the destructor,
  // since it just wraps a reference to an object.
  virtual ~AttributeWrapper()
    {
    }


  //============================================================
  // Accessors
  //============================================================

  // Return a reference to our wrapped object.

  Wrapped_t &array()
    {
      return wrapped_m;
    }

  const Wrapped_t &array() const
    {
      return wrapped_m;
    }


  //============================================================
  // Public virtual interface
  //============================================================

  // Print the contents to the given Inform stream.

  virtual void print(std::ostream &o) const
    {
      o << array() << std::endl;
    }

  /* Omit this until we have serialize/deserialize functions for DynamicArray.

  // serialize/deserialize the Array using the given stream.

  virtual int serialize(std::ostream &o) const
  {
    return serialize(o,array());
  }

  virtual int serialize(std::fstream &f) const
  {
    return serialize(f,array());
  }

  virtual int deserialize(std::ostream &o)
  {
    return deserialize(array(),o);
  }

  virtual int deserialize(std::fstream &f)
  {
    return deserialize(array(),f);
  }

  */

#if POOMA_MESSAGING

  // packSize, pack and unpack functions for particle swapping

  virtual int packSize(int elems) const
  {
    typedef typename Wrapped_t::Element_t Element_t;
    return ( elems * Cheetah::Serialize<Cheetah::CHEETAH, Element_t>::
	     size(Element_t()) );
  }

  virtual int pack(int pid, const IndirectionList<int> &list,
                   char *buffer) const
  {
    return array().engine().localPatch(pid).pack(list,buffer);
  }

  virtual int unpack(int pid, const Interval<1>& dom, char *buffer)
  {
    return array().engine().localPatch(pid).unpack(dom,buffer);
  }

#endif // POOMA_MESSAGING

private:
  // The object that we're wrapping

  Wrapped_t &wrapped_m;

  // Make copy consructor, default constructor, and operator= private
  // and undefined since they should not be used

  AttributeWrapper();

  AttributeWrapper(const This_t &);

  const This_t &operator=(const This_t &);
};


//-----------------------------------------------------------------------------
//
// A specialization of the Inform traits used to say that AttributeWrapper has
// a print method.
//
//-----------------------------------------------------------------------------

template<class T>
std::ostream &operator<<(std::ostream &o, const AttributeWrapper<T> &attrib)
{
  attrib.print(o);
  return o;
}


// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_PARTICLES_ATTRIBUTE_WRAPPER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: AttributeWrapper.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
