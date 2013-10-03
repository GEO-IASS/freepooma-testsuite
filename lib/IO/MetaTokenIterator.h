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
// Class:
//   Pooma::MetaTokenIterator
//-----------------------------------------------------------------------------

/** @file
 * @ingroup IO
 * @brief
 * MetaTokenIterator is a helper class designed to efficiently parse
 * a line from the "DiscField" .meta file.
 */

#ifndef POOMA_IO_METATOKENITERATOR_H
#define POOMA_IO_METATOKENITERATOR_H

#include <string>

namespace Pooma {

  /**
   * MetaTokenIterator is a helper class designed to efficiently parse
   * a line from the "DiscField" .meta file. MetaTokenIterator views
   * each line as having the form:
   *
   *     word0 [=] word1 word2 word3 [#comment]
   * 
   * where the words are separated by whitespace. The iterator returns
   * the sequence of words only, ignoring the '=' and any comments. 
   */

  // I don't inherit from std::iterator here because VC++'s iterator
  // class is badly screwed up. Instead I just typedef things
  // directly. 

  class MetaTokenIterator
  {
  public:

    //-------------------------------------------------------------------------
    // Typedefs
    //-------------------------------------------------------------------------

    // Iterator requirements

    typedef std::input_iterator_tag iterator_category;
    typedef std::string             value_type;
    typedef long                    difference_type;
    typedef std::string*            pointer;
    typedef std::string&            reference;

    // Convenience

    typedef std::string::size_type Size_t;

    //-------------------------------------------------------------------------
    // Constructors 
    //-------------------------------------------------------------------------

    MetaTokenIterator()
      : line_m(std::string("")), begIdx_m(std::string::npos)
    { }
      
    MetaTokenIterator(const std::string &line)
      : line_m(line), begIdx_m(0), endIdx_m(0), firstWord_m(true)
    {
      token_m.reserve(32);

      // Advance to the first word.

      next();
    }

    MetaTokenIterator(const MetaTokenIterator &model)
      : line_m(model.line_m), 
        begIdx_m(model.begIdx_m), 
        endIdx_m(model.endIdx_m),
        firstWord_m(model.firstWord_m)
    { 
      token_m.reserve(32);
    }
      
    //-------------------------------------------------------------------------
    // Iterator interface
    //-------------------------------------------------------------------------

    inline const std::string &operator*() const
    {
      token_m = line_m.substr(begIdx_m, endIdx_m - begIdx_m);
      return token_m;
    }
    
    inline const std::string *operator->() const
    {
      token_m = line_m.substr(begIdx_m, endIdx_m - begIdx_m);
      return &token_m;
    }

    inline MetaTokenIterator &operator++()
    {
      if (firstWord_m) skipEquals();
      firstWord_m = false;
      next();
      return *this;
    }

    inline MetaTokenIterator operator++(int)
    {
      MetaTokenIterator tmp(*this);
      if (firstWord_m) skipEquals();
      firstWord_m = false;
      next();
      return tmp;
    }

    inline bool operator==(const MetaTokenIterator &iter) const
    {
      // We don't actually check that we're in the same string as this
      // would complicate comparison to the end iterator. 

      if (iter.begIdx_m == std::string::npos) 
        return begIdx_m == std::string::npos;
      else 
        return (begIdx_m == iter.begIdx_m && endIdx_m == iter.endIdx_m);
    }

    inline bool operator!=(const MetaTokenIterator &iter) const
    { 
      return !operator==(iter);
    }

  private:

    //-------------------------------------------------------------------------
    // Utility functions
    //-------------------------------------------------------------------------

    void skipEquals();
    void next();

    //-------------------------------------------------------------------------
    // Data
    //-------------------------------------------------------------------------

    const std::string &line_m;

    Size_t begIdx_m;
    Size_t endIdx_m;

    // String used to store token - needed to implement operator->()
    // and operator*().  Mutable since these are const member
    // functions.

    mutable std::string token_m;

    // Are we in the first word? Used to handle the optional '=' after
    // the "keyword".

    bool firstWord_m;
  };

} // namespace Pooma

#endif // POOMA_IO_METATOKENITERATOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MetaTokenIterator.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
