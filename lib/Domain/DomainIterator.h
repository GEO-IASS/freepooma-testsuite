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

#ifndef POOMA_DOMAIN_DOMAIN_ITERATOR_H
#define POOMA_DOMAIN_DOMAIN_ITERATOR_H

//-----------------------------------------------------------------------------
// Classes: 
//   DomainIterator<Dom>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Domain
 * @brief
 *   DomainIterator<Dom> - Iterates through domain data (of type Dom)
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Domain/DomainTraits.h"

#include <stddef.h> // ptrdiff_t
#include <iterator>

//-----------------------------------------------------------------------------
// Open POOMA namespace:
//-----------------------------------------------------------------------------

// namespace POOMA {

/**
 * A simple iterator class to iterate through all of the points
 * in a given domain of type Dom. The individual points are returned
 * as Loc<Dim>s or Region<Dim,T>'s when the iterator is dereferenced.
 *
 * This is an input-iterator, in the STL sense.  It only defines deref,
 * ->, and ++ operators.
 */

template <class Dom>
class DomainIterator
{
public:

  //============================================================
  // Typedefs and static data
  //============================================================

  typedef DomainIterator<Dom>            This_t;
  typedef Dom                            Domain_t;
  typedef typename Domain_t::AskDomain_t Value_t;

  typedef std::input_iterator_tag        iterator_category;
  typedef Value_t                        value_type;
  typedef ptrdiff_t                      difference_type;
  typedef const Value_t*                 pointer;
  typedef const Value_t&                 reference;
  
  enum { dimensions = DomainTraits<Dom>::dimensions };

  //============================================================
  // Constructors
  //============================================================

  // The main DomainIterator stores the given domain, and sets all its 1D
  // iterators to the start. This constructor sets up a "begin" iterator.

  DomainIterator(const Dom &d, int size = 0)
    : domain_m(d), loc_m(d.firsts()), index_m(size)
    {
      PAssert(index_m >= 0 && index_m <= domain_m.size());
      for (int i=0; i < dimensions; ++i)
	current_m[i] = 0;
    }

  // Copy constructor.

  DomainIterator(const This_t &model)
    : domain_m(model.domain_m), loc_m(model.loc_m), index_m(model.index_m)
    {
      PAssert(index_m >= 0 && index_m <= domain_m.size());
      for (int i=0; i < dimensions; ++i)
	current_m[i] = model.current_m[i];
    }

  // The default constructor constructs an end iterator for an empty
  // domain.

  DomainIterator()
    : index_m(0)
    {
      for (int i=0; i < dimensions; ++i)
	current_m[i] = 0;
    }

  //============================================================
  // Accessors
  //============================================================

  // Dereference operator. Returns const ref to internal Loc.

  inline const Value_t &operator*() const
    {
      PAssert(!done());
      return loc_m;
    }

  // Member selection operator. Allows calling const Loc
  // member functions. Not too useful, but it is part of
  // the required input iterator interface. 

  inline const Value_t *operator->() const
    {
      PAssert(!done());
      return &loc_m;
    }

  // Equality tests.
  // Note that any two iterators that are both marked
  // as being at the end of iteration will compare equal.

  inline bool operator==(const This_t &rhs) const
    {
      return (index_m == rhs.index_m);
    }

  inline bool operator!=(const This_t &rhs) const
    {
      return (index_m != rhs.index_m);
    }

  // At-end (false) test.
  // Returns true if this iterator is at-end.

  bool done() const
    {
      return (index_m >= domain_m.size());
    }

  //============================================================
  // Mutators
  //============================================================

  // Assignment operator.

  This_t &operator=(const This_t &model)
    {
      index_m   = model.index_m;
      domain_m  = model.domain_m;
      loc_m     = model.loc_m;
      for (int i=0; i < dimensions; ++i)
	current_m[i] = model.current_m[i];
      return *this;
    }

  // Pre-increment operator. 
  // Takes us to the next point in the Interval<Dim> space of
  // points. This is done in Fortran (column-major) order.

  This_t &operator++()
    { 
      increment();
      return *this;
    }

  // Post-increment operator.
  // This has to make a copy, so prefer the above if possible.

  This_t operator++(int)
    {
      DomainIterator<Dom> save(*this);
      increment();
      return save;
    }

private:

  //============================================================
  // Data members.
  //============================================================

  // The domain we're iterating over.

  Domain_t domain_m;

  // Our current value, stored as a point domain.

  Value_t loc_m;

  // Our current position in each dimension.

  int current_m[dimensions];

  // Our current total index.

  int index_m;

  //============================================================
  // Implementation functions
  //============================================================

  // Increment iterator.

  void increment()
    {
      PAssert(!done());

      for (int i = 0; i < dimensions; ++i)
	{
	  if (++(current_m[i]) >= domain_m[i].length())
	    {
	      if (i < dimensions-1)
		{
		  current_m[i] = 0;
		  loc_m[i] = domain_m[i].first();
		}
	    }
	  else
	    {
	      loc_m[i] = domain_m[i](current_m[i]);
	      break;
	    }
	}

      // increase our total index

      ++index_m;
    }
};


// } // namespace POOMA

#endif // POOMA_DOMAIN_DOMAIN_ITERATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainIterator.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
