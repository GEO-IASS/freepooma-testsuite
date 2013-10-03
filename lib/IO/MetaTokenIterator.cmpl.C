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
// Include files
//-----------------------------------------------------------------------------

#include "IO/MetaTokenIterator.h"

namespace Pooma {

  //---------------------------------------------------------------------------
  // Handy constants
  //---------------------------------------------------------------------------

  namespace {
    const std::string::size_type npos = std::string::npos;
    const char *const delim = " \t\n";
  }

  //---------------------------------------------------------------------------
  // MetaTokenIterator member function implementations
  //---------------------------------------------------------------------------

  void MetaTokenIterator::skipEquals()
  {
    // If the next character is an '=', advance the end iterator
    // past it.

    if (endIdx_m < line_m.length() && line_m[endIdx_m] == '=')
      {
        ++endIdx_m;
      }
    else
      {
        // If the next "word" starts with an '=', advance the
        // endIdx past it so that the "next word" search will skip
        // over the '='.

        Size_t idx = line_m.find_first_not_of(delim, endIdx_m);

        if (idx != npos && line_m[idx] == '=') endIdx_m = idx + 1;
      }
  }

  void MetaTokenIterator::next()
  {
    // Find the first character that is not white-space

    // If we're in the first word and the next character is '=',
    // move the endIdx past the = sign. 

    begIdx_m = line_m.find_first_not_of(delim, endIdx_m);

    // If it is not found, then there are no more words - return.

    if (begIdx_m == npos) return;

    // If the next word starts with '#', we can ignore the rest of the
    // line.

    if (line_m[begIdx_m] == '#') 
      {
        begIdx_m = npos;
        return;
      }

    // Find the next character that is white-space

    endIdx_m = line_m.find_first_of(delim, begIdx_m);

    // If it is npos, then the word ends the line - sent endIdx_m to be
    // the length of the line (one past the last character in the
    // word).

    if (endIdx_m == npos) endIdx_m = line_m.length();

    // If the first word ends with '=', back off a character.

    if (firstWord_m && line_m[endIdx_m-1] == '=') --endIdx_m;
  }

} // namespace Pooma

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MetaTokenIterator.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
