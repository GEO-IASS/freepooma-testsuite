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

#ifndef POOMA_DOMAIN_INDIRECTION_LIST_H
#define POOMA_DOMAIN_INDIRECTION_LIST_H

//-----------------------------------------------------------------------------
// Class:
// IndirectionList
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * IndirectionList<T> is a class that takes any object with an
 * Array-like interface and stores the information as a list of elements
 * of type T.
 *
 * It is meant to store lists of indices for indirection-list
 * operations.  T can be an int, or multidimensional items like Loc<Dim>.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/DataBlockPtr.h"
#include "Utilities/PAssert.h"
#include <iosfwd>


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * IndirectionList is a domain representing a N dimension set of integer
 * sequences.  It is always a 1D list, but can store data of multiple
 * dimensions (if, say, the type T is Loc<Dim>)
 *
 * IndirectionList can be constructed using any of the Array class objects 
 * (Array, DynamicArray) . The constructor accepts only one such
 * object.  It is a templated constructor, however, and can work with any
 * object that has an Array-like interface.
 *
 * The default IndirectionList constructor initializes the IndirectionList
 * to be empty, that is, to have a length() of zero.  In that case, the
 * endpoints are undefined, as is any operation involving the IndirectionList.
 *
 * In addition to the constructors, IndirectionList has the following public
 * interface, similar to all 1D domain objects.
 *
 * IndirectionList<T> interface:
 *  - long size() - return the 'volume' of the domain, which is the product
 *      of the lengths of the N 1D IndirectionLists
 *  - bool empty() - return if any of the IndirectionList objects have
 *      length == 0
 *  - IndirectionList operator[](int N) - return the Nth IndirectionList
 *      in a multidimensional IndirectionList<M>. For IndirectionList
 *      objects, this just returns the object back.
 *  - comparison operators: <, >, !=, ==, <=, >= : compare a IndirectionList
 *      to another domain object.  The compared domains must have the same
 *      number of dimensions.
 *  - arithmetic accumulation operators +=, -=, *=, /= : add or subtract in a
 *      given domain.  The added domain must have the same number of
 *      dimensions, or a dimension of 1 (in which case, the same value
 *      is used for all dimensions), and be known to be single-valued (which
 *      is true for Loc and int's).  Note that for IndirectionList, *= and /=
 *      ARE allowed *= and /= result in scaling of each element in the
 *      IndirectionList which leaves
 *      the length (and size) the same.  += and -= shift the beginning
 *      endpoints by the given values, also leaving the length and size the
 *      same.  Negation of a IndirectionList negates the elements.
 *      Be aware that negation
 *      of an IndirecitonList will result in an object that is malformed for
 *      the DynamicArray.destroy() operators. 
 *  - binary arithmethic operators +, -, *, / : for + and -, adding a
 *      IndirectionList to another Loc or int returns a new IndirectionList.
 *      For * and /, scaling
 *      by a Loc or int also returns a IndirectionList object
 *  - increment/decrement operator ++, -- : only prefix versions of ++ and --
 *      are provided; they act just like += 1 and -= 1 operations.
 *  - int length() - number of elements (including endpoints) of the domain. 
 *  - long size() - same as above. 
 *  - int first() - the beginning endpoint. 
 *  - int last() - the ending endpoint.  
 *  - int min(), int max() - min or max of the endpoints.
 *  - IndirectionList::iterator begin() and end() - return iterators for
 *    the 1D domain.  These act like (at least) forward iterators.
 */

template<class T>
class IndirectionList
{
public:
  //
  // Typedefs and static data.
  //

  typedef IndirectionList<T>  This_t;
  typedef T                   Element_t;
  typedef IndirectionList<T>  Domain_t;
  typedef DataBlockPtr<T>     Storage_t;

  typedef long                Size_t;

  enum { dimensions      = 1 };
  enum { loopAware       = false };
  enum { singleValued    = false };
  enum { unitStride      = false };

  //
  // Constructors.
  //

  // Default constructor : initialize to an empty domain

  IndirectionList()
    : size_m(0)
    {
    }

  // Copy constructor

  IndirectionList(const This_t &a)
    : iList_m(a.iList_m), size_m(a.size_m)
    {
    }

  // Construct with just a number, indicating the number of elements
  // in the list, but do not give the list any values yet.

  IndirectionList(int num)    : iList_m(num), size_m(num) { PAssert(num >= 0); }
  IndirectionList(long num)   : iList_m(num), size_m(num) { PAssert(num >= 0); }
  IndirectionList(unsigned int num)  : iList_m(num), size_m(num) { }
  IndirectionList(unsigned long num) : iList_m(num), size_m(num) { }

  // Construct with a first, stride, num triple.  A list is constructed
  // that stores the list of points.  It must be possible to convert
  // from type T1 to type T, and type T1 must support += operations.

  template<class T1>
  IndirectionList(const T1 &first, const T1 &stride, Size_t num)
    : iList_m(num), size_m(num)
    {
      PAssert(num >= 0);

      T1 val = first;
      for (Size_t i=0; i < num; ++i)
	{
	  iList_m[i] = val;
	  val += stride;
	}
    }

  // General templated constructor.  This will work if the object contains an
  // engine that has a dataBlock() method (and internal object) and a
  // domain() method.

  template<class T1>
  explicit IndirectionList(const T1 &a)
    : iList_m(a.engine().dataBlock()), size_m(a.domain().size())
    {
    }

  //
  // Destructor.  DataBlockPtr destructor will take care of clean-up.
  //

  ~IndirectionList()
    {
    }

  // Copy assignment.

  This_t &operator=(const This_t &newilist) 
    {
      iList_m = newilist.iList_m;
      size_m = newilist.size_m;
      return *this;
    }

  // Templated assignment.

  template<class T1>
  This_t &operator=(const T1 &a)
    {
      iList_m = a.engine().dataBlock();
      size_m = a.domain().size();
      return *this;
    }

  // Basic domain operations.

  This_t &operator[](int i)
    {
      PAssert(i == 0);
      return *this;
    }

  const This_t &operator[](int i) const
    {
      PAssert(i == 0);
      return *this;
    }

  const T &operator()(Size_t i1) const
    {
      return iList_m[i1];
    }

  T &operator()(Size_t i1)
    {
      return iList_m[i1];
    }

  Size_t length() const 
    {
      return size_m;
    }

  T first() const 
    {
      return iList_m[0];
    }

  T last() const
    {
      return iList_m[size_m-1];
    }

  T stride() const
    {
      return T(0);
    }

  T min() const
    {
      // let's find the smallest item, assuming T has an operator<
      T result = iList_m[0];
      for (Size_t i=1; i<size_m; ++i)
	result = (iList_m[i] < result) ? iList_m[i] : result;
      return result;
    }

  T max() const
    {
      // let's find the largest item, assuming T has an operator<
      T result = iList_m[0];
      for (Size_t i=1; i<size_m; ++i)
	result = (result < iList_m[i]) ? iList_m[i] : result;
      return result;
    }

  Size_t size() const
    {
      return size_m;
    }

  bool empty() const
    {
      return (size() == 0);
    }

  bool initialized() const
    {
      return (!empty());
    }

  // arithmetic operators.  Just perform operation on each element.
  // These will make sure that we have our own copy of the data first,
  // so that we do not modify other user's lists.

  template<class T1>
  This_t &operator+=(const T1 &val)
    {
      Size_t n = size();
      if (n > 0)
	{
	  iList_m.makeOwnCopy();
	  for (Size_t i=0; i < n; ++i)
	    iList_m[i] += val;
	}
      return *this;
    }

  template<class T1>
  This_t &operator-=(const T1 &val)
    {
      Size_t n = size();
      if (n > 0)
	{
	  iList_m.makeOwnCopy();
	  for (Size_t i=0; i < n; ++i)
	    iList_m[i] -= val;
	}
      return *this;
    }

  template<class T1>
  This_t &operator*=(const T1 &val)
    {
      Size_t n = size();
      if (n > 0)
	{
	  iList_m.makeOwnCopy();
	  for (Size_t i=0; i < n; ++i)
	    iList_m[i] *= val;
	}
      return *this;
    }

  template<class T1>
  This_t &operator/=(const T1 &val)
    {
      Size_t n = size();
      if (n > 0)
	{
	  iList_m.makeOwnCopy();
	  for (Size_t i=0; i < n; ++i)
	    iList_m[i] /= val;
	}
      return *this;
    }

  //
  // I/O
  //

  // print an IndirectionList to a stream, in the format
  //   "[" val1, val2, ... , valN "]"

  template<class Out>
  void print(Out &o) const
    {
      o << "[";
      for (Size_t i=0; i < size(); ++i)
	{
	  o << (*this)(i);
	  if (i < (size() - 1))
	    o << ",";
	}
      o << "]";
    }

private:

  // Data storage, as a reference-counted data block

  Storage_t iList_m;

  // Number of elements

  Size_t size_m;
};


/// print an IndirectionList to a stream, in the format
///   "[" val1, val2, ... , valN "]"

template<class T>
std::ostream& operator<<(std::ostream &o, const IndirectionList<T> &list)
{
  list.print(o);
  return o;
}


// } // namespace Pooma
//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_INDIRECTION_LIST_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IndirectionList.h,v $   $Author: richard $
// $Revision: 1.20 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
