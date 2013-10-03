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

#ifndef POOMA_DOMAIN_DOMAIN_BLOCK_ITERATOR_H
#define POOMA_DOMAIN_DOMAIN_BLOCK_ITERATOR_H

//-----------------------------------------------------------------------------
// Classes: 
//   DomainBlockIterator<Dom>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Domain
 * @brief
 *   DomainBlockIterator<Dom> - Iterates through domain data (of type Dom),
 *   and returns block domains (Interval or Region)
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Domain/DomainTraits.h"

//-----------------------------------------------------------------------------
// Open POOMA namespace:
//-----------------------------------------------------------------------------

// namespace POOMA {

/**
 * A simple iterator class to iterate through all of the points
 * in a given domain of type Dom.  This iterator returns Interval or
 * Region objects that define the blocks formed as the "cells" between
 * the "vertices" that are the domain points.
 *
 * Dereferencing a DomainBlockIterator returns an Interval or Region with
 * the current cell.  You can also call the following methods on a
 * DomainBlockIterator:
 *
 *  Loc<Dim> point() - returns the block-index of the current block.
 *                     This is a Loc<Dim> from 0 ... # blocks in each
 *                     dimension.
 *
 *  int index() - returns an index for the current cell.  The index values
 *                start at zero and increment by one each time you move
 *                to a new cell.
 *
 *  bool done() - returns true if the iterator is done, that is, if it
 *                at the end and would compare equal to an "end" iterator.
 *
 * This is an input-iterator, in the STL sense.  It only defines deref,
 * ->, and ++ operators.
 */

template <class Dom>
class DomainBlockIterator
{
public:

  //============================================================
  // Typedefs and static data
  //============================================================

  typedef DomainBlockIterator<Dom>          This_t;
  typedef Dom                               Domain_t;
  typedef typename Domain_t::OneDomain_t    OneDomain_t;
  typedef typename Domain_t::AskDomain_t    Value_t;
  typedef typename Domain_t::BlockDomain_t  Block_t;
  typedef typename Block_t::OneDomain_t     OneBlock_t;
  typedef typename OneDomain_t::iterator    iterator;
  typedef typename Domain_t::Element_t      Element_t;

  enum { dimensions = DomainTraits<Dom>::dimensions };

  //============================================================
  // Constructors
  //============================================================

  // Default constructor constructs an "end" iterator.

  DomainBlockIterator()
    : index_m(-1)
    {
    }

  // The main DomainBlockIterator stores the given domain, and sets all its 1D
  // iterators to the start. This constructor sets up a "begin" iterator.

  DomainBlockIterator(const Dom &d);

  // Copy constructor.

  DomainBlockIterator(const This_t &model);

  //============================================================
  // Destructor
  //============================================================

  ~DomainBlockIterator()
    {
    }

  //============================================================
  // Accessors
  //============================================================

  // Dereference operator. Returns const ref to internal block.

  inline const Block_t &operator*() const
    {
      PAssert(!done());
      return block_m;
    }

  // Member selection operator. Allows calling const block
  // member functions. Not too useful, but it is part of
  // the required input iterator interface. 

  inline const Block_t *operator->() const
    {
      PAssert(!done());
      return block_m;
    }

  // Return the upper-left corner of the current block; this is just a
  // single point, not a whole block.

  inline const Value_t &point() const
    {
      PAssert(!done());
      return loc_m;
    }

  // Return the current block index
  
  inline int index() const
    {
      PAssert(!done());
      return index_m;
    }

  // Equality tests.
  // Note that any two iterators that are both marked
  // as being at the end of iteration will compare equal.

  inline bool operator==(const This_t &rhs) const
    {
      if (done() || rhs.done())
	return (done() && rhs.done());
      else
	return (block_m == *rhs);
    }

  inline bool operator!=(const This_t &rhs) const
    {
      return !(operator==(rhs));
    }

  // At-end (false) test.
  // Returns true if this iterator is at-end.

  bool done() const
    {
      return (index_m < 0);
    }

  //============================================================
  // Mutators
  //============================================================

  // Assignment operator.

  This_t &operator=(const This_t &model)
    {
      domain_m  = model.domain_m;
      loc_m     = model.loc_m;
      block_m   = model.block_m;
      index_m   = model.index_m;
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
      This_t save(*this);
      increment();
      return save;
    }

private:

  //============================================================
  // Data members.
  //============================================================

  // The domain we're iterating over

  Domain_t domain_m;

  // Our current left and right 

  iterator current_m[dimensions];

  // Our current corner point, stored as a point domain

  Value_t loc_m;

  // Our current block

  Block_t block_m;

  // The current block index.  If this is < 0, we're at the end.

  int index_m;


  //============================================================
  // Implementation functions
  //============================================================

  // Set our done flag to true

  void setDone()
    {
      index_m = (-1);
    }

  // Increment iterator.

  void increment();
};


//-----------------------------------------------------------------------------
//
// DomainBlockIterator out-of-line method definitions
//
//-----------------------------------------------------------------------------


/// The main DomainBlockIterator stores the given domain, and sets all its 1D
/// iterators to the start. This constructor sets up a "begin" iterator.

template<class Dom>
DomainBlockIterator<Dom>::DomainBlockIterator(const Dom &d)
  : domain_m(d), loc_m(0), index_m(0)
{
  for (int i=0; i < dimensions; ++i)
    {
      if (d[i].begin() == d[i].end())
	{
	  setDone();
	  break;
	}
      else
	{
	  current_m[i] = d[i].begin();
	  iterator next = current_m[i];
	  Element_t a, b;
	  if (++next == d[i].end())
	    {
	      a = b = (*(current_m[i])).first();
	    }
	  else
	    {
	      a = (*(current_m[i])).first();
	      b = (*next).first();
	    }
	  if (b < a)
	    block_m[i] = OneBlock_t(b + 1, a);
	  else if (b == a)
	    block_m[i] = OneBlock_t(a, a);
	  else
	    block_m[i] = OneBlock_t(a, b - 1);
	}
    }
}


//-----------------------------------------------------------------------------
// Copy constructor.
//-----------------------------------------------------------------------------

template<class Dom>
DomainBlockIterator<Dom>::DomainBlockIterator(const This_t &model)
  : domain_m(model.domain_m), loc_m(model.loc_m),
    block_m(model.block_m), index_m(model.index_m)
{
  for (int i=0; i < dimensions; ++i)
    current_m[i] = model.current_m[i];
}


//-----------------------------------------------------------------------------
// Increment iterator.
//-----------------------------------------------------------------------------

template<class Dom>
void DomainBlockIterator<Dom>::increment()
{
  PAssert(!done());

  for (int i = 0; i < dimensions; ++i)
    {
      if (++(current_m[i]) == domain_m[i].end())
	{
	  // if we go past the very end, this dimension must
	  // have only one point, so go on to the next one.  We don't
	  // need to change block_m or loc_m, since they were set in the
	  // constructor

	  if (i >= dimensions-1)
	    {
	      setDone();
	      break;
	    }
	}
      else
	{
	  // we haven't gone past the end, so get the next point

	  iterator next = current_m[i];
	  if (++next == domain_m[i].end())
	    {
	      // We've reached the end of this dimension, so move to
	      // the next.  If there isn't a next, we're done.
	      
	      if (i < dimensions-1)
		{
		  current_m[i] = domain_m[i].begin();
		  next = current_m[i];
		  Element_t a, b;
		  if (++next == domain_m[i].end())
		    {
		      a = b = (*(current_m[i])).first();
		    }
		  else
		    {
		      a = (*(current_m[i])).first();
		      b = (*next).first();
		    }
		  if (b < a)
		    block_m[i] = OneBlock_t(b + 1, a);
		  else if (b == a)
		    block_m[i] = OneBlock_t(a, a);
		  else
		    block_m[i] = OneBlock_t(a, b - 1);
		  loc_m[i] = 0;
		}
	      else
		{
		  setDone();
		  break;
		}
	    }
	  else
	    {
	      // We're still going strong in this dimension, go get this
	      // value and quit.

	      Element_t a = (*(current_m[i])).first();
	      Element_t b = (*next).first();
	      if (b < a)
		block_m[i] = OneBlock_t(b + 1, a);
	      else
		block_m[i] = OneBlock_t(a, b - 1);
	      loc_m[i] = loc_m[i].first() + 1;
	      index_m++;
	      break;
	    }
	}
    }
}

// } // namespace POOMA

#endif // POOMA_DOMAIN_DOMAIN_BLOCK_ITERATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainBlockIterator.h,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
