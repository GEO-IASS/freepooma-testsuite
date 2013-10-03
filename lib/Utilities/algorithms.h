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
// Function templates: 
//   Algorithms::delete_backfill
//   Algorithms::delete_shiftup
//   Algorithms::copy
//-----------------------------------------------------------------------------

#ifndef POOMA_UTILITIES_ALGORITHMS_H
#define POOMA_UTILITIES_ALGORITHMS_H

/** @file
 * @ingroup Utilities
 * @brief
 * copy, delete_backfill and delete_shiftup algorithms.
 *
 *   copy(start, end, dest)
 *    - An STL-like copy algorithm. This one first checks, using 
 *      ElementProperties, for whether or not the element is "concrete".
 *      If so, it performs the copy with std::copy or with std::memmove. 
 *      If not, it does the copy explicitly using ElementProperties::construct.
 *
 *   delete_backfill(data_begin, data_end, kill_begin, kill_end)
 *    - Deletes the elements described by the range [kill_begin,kill_end)
 *      from the data described by the range [data_begin,data_end), filling
 *      the deleted slots from the end of the list. All iterators are assumed
 *      to be random access iterators. Returns the number of deleted elements.
 *
 *   delete_shiftup(data_begin, data_end, kill_begin, kill_end)
 *    - Deletes the elements described by the range [kill_begin,kill_end)
 *      from the data described by the range [data_begin,data_end), filling
 *      the deleted slots from the end of the list. All iterators are assumed
 *      to be random access iterators. Returns the number of deleted elements.
*/

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/ElementProperties.h"
#include "Utilities/PAssert.h"

#include <iterator> // iterator_traits, reverse_iterator
#include <string.h> // memmove
#include <algorithm>// copy

// Open namespaces:

namespace Pooma {
namespace Algorithms {


/**
 * Tag class used to choose the proper specialization of copy.
 */

#define POOMA_TAG_DUMMY_MEMBERS(tagclassname)                \
  inline tagclassname() { };                                 \
  inline tagclassname(const tagclassname &) { };             \
  inline tagclassname &operator=(const tagclassname &) { return *this; }
  
template <bool type>
struct IsConcrete 
{ 
POOMA_TAG_DUMMY_MEMBERS(IsConcrete);
};


/// void copy(Iterator begin, Iterator end, Iterator2 dest)
///
/// General version of Pooma's copy routine. Uses ElementProperties to
/// determine whether elements are concrete or not and calls the appropriate
/// specialization.
///
/// This will work with overlapping regions if the destination is below
/// the source (std::copy has the same limitation).

template <class It, class It2>
inline It2 copy(It begin, It end, It2 dest)
{
  PAssert(begin <= end);
  PAssert(dest <= begin);
  typedef typename std::iterator_traits<It>::value_type Value_t;
  typedef ElementProperties<Value_t> EP_t;
  typedef IsConcrete<EP_t::concrete> ConcreteType_t;
  
  return copy_special(begin, end, dest, ConcreteType_t());
}


/// void copy_special(Iterator begin, Iterator end, 
///                   Iterator2 dest, IsConcrete<true>)
///
/// Version of Pooma's copy routine for "concrete" types. Uses either the stl
/// copy algorithm or memmove to copy the data, depending on the data size.
/// (Some stl copy algorithms may already make this optimization in which case
/// this could be simplified - probably should have a cpp macro that comments
/// out the memmove version for platform where copy is just as fast.)

template <class It, class It2>
inline It2 copy_special(It begin, It end, It2 dest, IsConcrete<true>)
{
  typedef std::iterator_traits<It>                DataTraits_t;
#if POOMA_NONSTANDARD_ITERATOR
  typedef typename DataTraits_t::distance_type    Diff_t;
#else
  typedef typename DataTraits_t::difference_type  Diff_t;
#endif
  typedef typename DataTraits_t::value_type       Value_t;

  // Comment out this test for long data block and just use std::copy
  // because std::advance does not work correctly with MS header files
  // when tmp is not a pointer or a std::iterator type.
#if 0
  const Diff_t len = end - begin;

  if (len < 100)
    return std::copy(begin, end, dest);
  else
    {
      memmove(&(*dest), &(*begin), len*sizeof(Value_t));
      It2 tmp(dest);
      std::advance(tmp,len);
      return tmp;  
    }
#else
  return std::copy(begin, end, dest);
#endif
}


/// void copy_special(Iterator begin, Iterator end, 
///                   Iterator2 dest, IsConcrete<false>)
///
/// Version of Pooma's copy routine for non-"concrete" types. This uses
/// ElementProperties::construct to do the copy.

template <class It, class It2>
inline It2 copy_special(It begin, It end, It2 dest, IsConcrete<false>)
{
  typedef typename std::iterator_traits<It>::value_type Value_t;
  while (begin < end)
    {
      ElementProperties<Value_t>::construct(dest++, *begin++);
    }
  return dest;
}


/// difference_type delete_backfill(DataIterator, DataIterator, 
///                                 KillIterator, KillIterator, difference_type)
///
/// Loop through the data and destroy the desired elements, replacing
/// them with elements from the end of the list.  The KillIterators
/// refer to a list of indices of elements to be killed.  The final
/// optional argument provides an offset for the kill list index values,
/// for cases in which they are not zero-based.
/// Return the number of elements killed.
///
/// NOTE: This algorithm assumes that the killList is sorted
/// and that all iterators are random-access iterators. 

template <class DataIterator, class KillIterator>
inline
#if POOMA_NONSTANDARD_ITERATOR
typename std::iterator_traits<DataIterator>::distance_type
#else
typename std::iterator_traits<DataIterator>::difference_type
#endif
delete_backfill(DataIterator data_begin, DataIterator data_end, 
  const KillIterator kill_begin, const KillIterator kill_end,
#if POOMA_NONSTANDARD_ITERATOR
  typename std::iterator_traits<DataIterator>::distance_type offset = 0)
#else
  typename std::iterator_traits<DataIterator>::difference_type offset = 0)
#endif
{
  PAssert(data_end >= data_begin);
  PAssert(kill_end >= kill_begin);
  
  // No data has to be moved if we're destroying values that are at
  // the end of the sequence. Thus we look for these first.

#if POOMA_NONSTANDARD_ITERATOR
  std::reverse_iterator<KillIterator, const int> rk_pos(kill_end);
  std::reverse_iterator<KillIterator, const int> rk_end(kill_begin);
#else
  std::reverse_iterator<KillIterator>            rk_pos(kill_end);
  std::reverse_iterator<KillIterator>            rk_end(kill_begin);
#endif
  
  typedef std::iterator_traits<DataIterator>     DataTraits_t;
#if POOMA_NONSTANDARD_ITERATOR
  typedef typename DataTraits_t::distance_type   Diff_t;
#else
  typedef typename DataTraits_t::difference_type Diff_t;
#endif

  Diff_t last = data_end - data_begin - 1;
    
  while (!(rk_pos == rk_end))
    {
    // If we're no longer deleting the last element, we break and
    // finish up below. (Can't do post inc/dec here as it would put
    // the count off by one when we break.)
        
      if ((*rk_pos - offset) != last) break;

    // We're still on the last element. Deleting it is simply a matter
    // of moving our last pointer and continuing.
              
      --last;
      ++rk_pos;
    }

  // Now we're deleting non-last elements. This only executes if
  // the above loop was terminated early by the break.
  // You have to go through the data backward to ensure that you 
  // don't fill with an element you're going to later want to delete.
  
  DataIterator last_pos = data_begin + last;
  
  while (!(rk_pos == rk_end))
    {
    // For BackFill, we overwrite the deleted element with
    // the current last element and then decrement the "last"
    // index since that element has effectively been moved.
          
      *(data_begin + (*rk_pos++ - offset)) = *last_pos--;
    }

  return kill_end - kill_begin;
}


/// difference_type delete_shiftup(DataIterator, DataIterator, 
///                                KillIterator, KillIterator, difference_type)
///
/// Loop through the data and destroy the desired elements, shifting 
/// remaining elements forward in the list so as to maintain their order.
/// The KillIterators refer to a list of indices of elements to be killed. 
/// The final optional argument provides an offset for the kill list index
/// values, for cases in which they are not zero-based.
/// Return the number of elements killed.
///
/// NOTE: This algorithm assumes that the killList is sorted
/// and that all iterators are random-access iterators. 

template <class DataIterator, class KillIterator>
inline
#if POOMA_NONSTANDARD_ITERATOR
typename std::iterator_traits<DataIterator>::distance_type
#else
typename std::iterator_traits<DataIterator>::difference_type
#endif
delete_shiftup(DataIterator data_begin, DataIterator data_end,
  KillIterator kill_begin, KillIterator kill_end,
#if POOMA_NONSTANDARD_ITERATOR
  typename std::iterator_traits<DataIterator>::distance_type offset = 0)
#else
  typename std::iterator_traits<DataIterator>::difference_type offset = 0)
#endif
{
  DataIterator insert_pos = data_begin + (*kill_begin - offset);
  KillIterator kill_pos   = kill_begin;
  
  typedef std::iterator_traits<DataIterator>          DataTraits_t;
  typedef typename DataTraits_t::value_type           Value_t;
#if POOMA_NONSTANDARD_ITERATOR
  typedef typename DataTraits_t::distance_type        Diff_t;
#else
  typedef typename DataTraits_t::difference_type      Diff_t;
#endif

  while (kill_pos < kill_end)
  {
    Diff_t copy_index = *kill_pos + 1;
    while ( (kill_pos + 1) < kill_end && copy_index == *(kill_pos + 1) )
    {
      ++copy_index;
      ++kill_pos;
    }

    DataIterator copy_begin = data_begin + (copy_index - offset);
    DataIterator copy_end;
	  
    if (copy_begin < data_end)
    {
      if (kill_pos + 1 < kill_end)
        copy_end = data_begin + (*(kill_pos + 1) - offset);
      else
        copy_end = data_end;
      
      const Diff_t length = copy_end - copy_begin;

      Pooma::Algorithms::copy(copy_begin, copy_end, insert_pos);
      
      insert_pos += length;
    }
    ++kill_pos;
  }

  return kill_end - kill_begin;
}


/// DataIterator find_most_common(DataIterator, DataIterator)
///
/// Loop through the data, counting each distinct value and 
/// return an iterator pointing to the most common one (or the end of
/// the sequence if not found).
///
/// NOTE: This algorithm assumes that the data is sorted.

template<class DataIterator>
inline DataIterator
find_most_common(DataIterator dataBegin, DataIterator dataEnd)
{
  DataIterator checkValue, mostCommonValue = dataEnd;
  int checkCount, count = 0;

  while (dataBegin != dataEnd)
    {
      checkValue = dataBegin++;
      checkCount = 1;
      while (dataBegin != dataEnd && *dataBegin == *checkValue)
	{
	  ++checkCount;
	  ++dataBegin;
	}

      if (checkCount > count)
        {
	  mostCommonValue = checkValue;
	  count = checkCount;
	}
    }

  return mostCommonValue;
}

} // namespace Algorithms
} // namespace Pooma

#endif // POOMA_UTILITIES_ALGORITHMS
