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
//   FieldOffset
//   FieldOffsetList
// Functions:
//   accumulate
//   sum
//   av
//   min
//   max
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_OFFSET_H
#define POOMA_FIELD_OFFSET_H

/** @file
 * @ingroup Field
 * @brief
 * specifies a relative cell offset and subfield number
 *
 * FieldOffset
 *   - specifies a relative cell offset and subfield number
 * FieldOffsetList
 *   - is a sequence of FieldOffset's
 * FieldOffsetList reductions
 *   - computations using the entries in a FieldOffsetList
 */

// FIXME: Perhaps we wish to store pointers to input and output
// FIXME: centerings in a FieldOffset.  Storing the input centering
// FIXME: will permit assertion checking when a FieldOffset is used to
// FIXME: reference a field.  The stored output centering is used in
// FIXME: data-parallel uses.

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include <iostream>
#include <vector>
#include <functional>
#include "Field/Field.h"

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

template <int Dim>
class FieldOffset;

template <int Dim>
class FieldOffsetList;


//-----------------------------------------------------------------------------
// FieldOffset.
//-----------------------------------------------------------------------------

/**
 * Given a field f, a Loc loc, and a field offset (offset,num), a
 * field value can be obtained.  Since each value specified by the
 * field's centering is stored in a separate subfield, the offset is
 * used to specify the appropriate subfield.
 */

template <int Dim>
class FieldOffset {
public:

  //---------------------------------------------------------------------------
  /// User-callable constructors. These ctors are meant to be called by users.

  FieldOffset(const Loc<Dim> &loc, const int subFieldNumber = 0)
    : cell_offset_m(loc), subfield_number_m(subFieldNumber)
  {
#if POOMA_BOUNDS_CHECK
    PInsist(subfield_number_m >= 0, "Erroneous FieldOffset subfield number.");
#endif
    return;
  }

  //---------------------------------------------------------------------------
  /// Internal POOMA constructor.  These operations are used internally
  /// by POOMA.  They are not really meant to be called by users.

  FieldOffset()
    : cell_offset_m(Loc<Dim>()), subfield_number_m(0)
  {}

  inline void setSubFieldNumber(int subFieldNumber)
    {
      subfield_number_m = subFieldNumber;
      return;
    }

  inline Loc<Dim> &modifyCellOffset()
    {
      return cell_offset_m;
    }

  //---------------------------------------------------------------------------
  // Accessors.

  /// Get the cell offset.

  inline const Loc<Dim> &cellOffset() const
    {
      return cell_offset_m;
    }

  /// Get the subfield number.

  inline int subFieldNumber() const
    {
      return subfield_number_m;
    }

private:

  /// The cell offset.
  Loc<Dim> cell_offset_m;

  /// The subfield number, if appropriate.
  int subfield_number_m;
};


//-----------------------------------------------------------------------------
/// Overload the << operator to print a FieldOffset to a stream.
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &o,
			 const FieldOffset<Dim> &fieldOffset)
{
  return o << "FieldOffset: (" << fieldOffset.cellOffset()
	   << ", " << fieldOffset.subFieldNumber() << ")";
}

//-----------------------------------------------------------------------------
/// Define equality and inequality operators for FieldOffsets/
//-----------------------------------------------------------------------------

template <int Dim>
inline bool
operator==(const FieldOffset<Dim> &fieldOffset1,
	   const FieldOffset<Dim> &fieldOffset2)
{
  return
    fieldOffset1.cellOffset() == fieldOffset2.cellOffset() &&
    fieldOffset1.subFieldNumber() == fieldOffset2.subFieldNumber();
}

template <int Dim>
inline bool
operator!=(const FieldOffset<Dim> &fieldOffset1,
	   const FieldOffset<Dim> &fieldOffset2)
{
  return !(fieldOffset1 == fieldOffset2);
}


//-----------------------------------------------------------------------------
// FieldOffsetList.
//-----------------------------------------------------------------------------

/**
 * A set of FieldOffset's can be stored in a FieldOffsetList.  The
 * list has a fixed size.  The following member operations are
 * supported:
 * - size_type size() const - return the number of FieldOffset's.
 * - const_reference operator[](size_type n) - return the nth
 *   FieldOffset.
 */

template <int Dim>
class FieldOffsetList
{
public:

  // Exported typedefs.
  
  typedef size_t size_type;
  typedef FieldOffset<Dim>& reference;
  typedef const FieldOffset<Dim>& const_reference;

  /// Return the number of FieldOffset's.

  size_type size() const
  {
    return v_m.size();
  }

  /// Return a FieldOffset.

  const_reference operator[](const size_type n) const
  {
#if POOMA_BOUNDS_CHECK
    PInsist(n < size(), "Erroneous FieldOffsetList index.");
#endif
    return v_m[n];
  }

  //---------------------------------------------------------------------------
  // Internal POOMA operators.  These operations are used internally
  // by POOMA.  They are not really meant to be called by users.

  /// Create an empty list.  This is used for arrays or std::vectors.

  FieldOffsetList() {}

  /// Create a list that can hold the specified number of entries.

  FieldOffsetList(const size_type sz)
  {
#if POOMA_BOUNDS_CHECK
    PInsist(sz > 0, "Erroneous FieldOffsetList size.");
#endif
    v_m.reserve(sz);
  }

  /// Construct from a vector.

  FieldOffsetList(const std::vector<FieldOffset<Dim> > &v) {
    v_m.resize(v.size());
    std::copy(v.begin(), v.end(), v_m.begin());
  }

  /// Copy a vector's entries to this FieldOffsetList.

  FieldOffsetList &operator=(const std::vector<FieldOffset<Dim> > &v)
  {
    v_m.resize(v.size());
    std::copy(v.begin(), v.end(), v_m.begin());
    return *this;
  }

  /// Permit adding the specified entry.

  reference operator[](const size_type n)
  {
#if POOMA_BOUNDS_CHECK
    PInsist(n < size(), "Erroneous FieldOffsetList index.");
#endif
    return v_m[n];
  }

private:

  std::vector<FieldOffset<Dim> > v_m;
};

//-----------------------------------------------------------------------------
/// Overload the << operator to print a FieldOffsetList to a stream.
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &o,
			 const FieldOffsetList<Dim> &fieldOffsetList)
{
  o << "FieldOffsetList:\n";
  for (int index = 0; index < fieldOffsetList.size(); ++index)
    o << fieldOffsetList[index] << std::endl;
  return o;
}


//-----------------------------------------------------------------------------
// FieldOffsetList Reductions.
//-----------------------------------------------------------------------------

///@name FieldOffsetList<Dim> Reductions
///
/// FIXME: What do we do for data-parallel statements?
///
/// Many computations using `FieldOffsetList's perform similar
/// computations.  The provided reduction functions support these
/// computations.  All take a field, a FieldOffsetList, and a Loc.  The
/// Loc specifies a cell within the given field.  Together with a
/// FieldOffsetList and a field, a field value is specified.  The
/// reduction function combines all the specified field values
/// corresponding to `FieldOffset's in the list.  The list must not be
/// empty.
///
/// accumulate:	Sequentially applies the given binary function to all
///		field offsets in the list.
/// sum:		Adds the values indicated by the field offset list.
/// av:		Averages the values indicated by the field offset
///		list.  Note the division is computed using the element
///		type, e.g., integral or floating point division.
/// min:		Returns the smallest of the indicated values.
/// max:		Returns the largest of the indicated values.

//@{

/// Accumulate all the specified field locations using the supplied STL
/// binary function.  Then, for each FieldOffset fo in the list, result
/// = binary_op(result, fv), where fv is the corresponding field value.
/// FIXME: Add data-parallel code.

template<class GeometryTag, class T, class Expr, int Dim,
         class BinaryFunction>
inline
T
accumulate(BinaryFunction binary_op,
	   const Field<GeometryTag, T, Expr>& field, 
	   const FieldOffsetList<Dim> &lst, 
	   const Loc<Dim> &loc)
{
  typedef typename Field<GeometryTag, T, Expr>::T_t T_t;
  typedef typename FieldOffsetList<Dim>::size_type size_type;
  CTAssert((Field<GeometryTag, T, Expr>::dimensions == Dim));

  const size_type lstLength = lst.size();
  PInsist(lstLength > 0, "accumulate must be given a nonempty list.");
  T_t init = field(lst[0], loc);
  for (size_type i = 1; i < lstLength ; ++i)
    init = binary_op(init, field(lst[i], loc));
  return init;
}


/// Sum all the values at the field locations.

template<class GeometryTag, class T, class Expr, int Dim>
inline
T
sum(const Field<GeometryTag, T, Expr>& field, 
    const FieldOffsetList<Dim> &lst, 
    const Loc<Dim> &loc)
{
  typedef typename Field<GeometryTag, T, Expr>::T_t T_t;
  CTAssert((Field<GeometryTag, T, Expr>::dimensions == Dim));
  return accumulate(std::plus<T_t>(), field, lst, loc);
}


/// Average all the values at the field locations.  Note the return
/// value has the same type as the field types so integer division may
/// be used.

template<class GeometryTag, class T, class Expr, int Dim>
inline
T
av(const Field<GeometryTag, T, Expr>& field, 
   const FieldOffsetList<Dim> &lst, 
   const Loc<Dim> &loc)
{
  typedef typename Field<GeometryTag, T, Expr>::T_t T_t;
  CTAssert((Field<GeometryTag, T, Expr>::dimensions == Dim));
  return sum(field, lst, loc) / lst.size();
}


/// Return the minimum value of the field locations.

template <class T>
struct fomin : public std::binary_function<T, T, T>
{
  T operator()(const T &op1, const T &op2) const {
    return std::min(op1, op2);
  }
};

template<class GeometryTag, class T, class Expr, int Dim>
inline
T
min(const Field<GeometryTag, T, Expr>& field, 
    const FieldOffsetList<Dim> &lst, 
    const Loc<Dim> &loc)
{
  typedef typename Field<GeometryTag, T, Expr>::T_t T_t;
  CTAssert((Field<GeometryTag, T, Expr>::dimensions == Dim));
  return accumulate(fomin<T_t>(), field, lst, loc);
}


/// Return the maximum value of the field locations.

template <class T>
struct fomax : public std::binary_function<T, T, T>
{
  T operator()(const T &op1, const T &op2) const {
    return std::max(op1, op2);
  }
};

template<class GeometryTag, class T, class Expr, int Dim>
inline
T
max(const Field<GeometryTag, T, Expr>& field, 
    const FieldOffsetList<Dim> &lst, 
    const Loc<Dim> &loc)
{
  typedef typename Field<GeometryTag, T, Expr>::T_t T_t;
  CTAssert((Field<GeometryTag, T, Expr>::dimensions == Dim));
  return accumulate(fomax<T_t>(), field, lst, loc);
}

//@}

//-----------------------------------------------------------------------------
// replicate.
//-----------------------------------------------------------------------------

/// Copy field values to the specified locations.  The first field
/// parameter specifies the field supplying the values to replicate.
/// The second std::vector<FieldOffsetList> parameter specifies, for
/// each value in the returned field, which input field value to use.
/// The vector's length must match the number of values in each output
/// field's cell.  For example, the output field's first value is
/// copied from the location specified by the vector's first list.  The
/// third parameter indicates the returned field's centering.

template <class GeometryTag, class T, class Expr, int Dim>
inline
typename
View2<Field<GeometryTag, T, Expr>, std::vector<FieldOffset<Dim> >,
      Centering<Dim> >::Type_t
replicate(const Field<GeometryTag, T, Expr>& field,
          const std::vector<FieldOffsetList<Dim> > &vec,
          const Centering<Dim> &centering)
{
  CTAssert((Field<GeometryTag, T, Expr>::dimensions == Dim));
  typedef typename std::vector<FieldOffsetList<Dim> >::size_type vsize_type;
  PInsist(vec.size() > 0, "Cannot replicate no values.");
  PInsist(vec.size() == centering.size(),
	  "Vector and output centering sizes must match.");

  std::vector<FieldOffset<Dim> > vecFO(vec.size());
  for (vsize_type i = 0; i < vec.size(); ++i) {
    PInsist(vec[i].size() == 1, "Can replicate only one value.");
    vecFO[i] = vec[i][0];
  }

  return field(vecFO, centering);
}


#endif // POOMA_FIELD_OFFSET_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldOffset.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:43 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
