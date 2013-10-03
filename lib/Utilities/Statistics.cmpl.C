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

// Include files

#include "Utilities/Inform.h"
#include "Utilities/Statistics.h"

#include <iostream>
#include <iomanip>
#include <string.h>

namespace Pooma {
  

///////////////////////////////////////////////////////////////////////////////
//
// Constructor. Simple, but doesn't need to be inline.
//
///////////////////////////////////////////////////////////////////////////////

Statistics::Statistics()
{
}


///////////////////////////////////////////////////////////////////////////////
//
// Destructor: delete all of the StatData_t's.
//
///////////////////////////////////////////////////////////////////////////////

Statistics::~Statistics() 
{
  for (int i = 0; i < statList_m.size(); ++i) 
    {
      delete statList_m[i];
    }
}


///////////////////////////////////////////////////////////////////////////////
//
// Default Filtering function (does nothing).
//
///////////////////////////////////////////////////////////////////////////////

long Statistics::defaultFilter(long val)
{
  return val;
}


///////////////////////////////////////////////////////////////////////////////
//
// Print out the statistics to the given Inform object. The filter function
// exists in case we want to process the value before printing it out. One
// interesting applications is if we are in a multi-context program. In that
// case we want to reduce over contexts and print out on one of them.
//
///////////////////////////////////////////////////////////////////////////////

void Statistics::print(Inform &o, long (*filter)(long)) 
{
  int i, j;

  // If we have no stats, just return.
  
  if (statList_m.size() == 0)
    return;

  // For each statistic, print out the description, a set of ...'s, and
  // the stat, right-justified to 10 places.
  
  o << "Runtime statistics summary:" << std::endl;
  for (i = 0; i < statList_m.size(); ++i) 
    {
      o << statList_m[i]->description() << " ";

      int numperiods = 53 - strlen(statList_m[i]->description().c_str());
      if (numperiods < 2)
        numperiods = 2;
    
      for (j = 0; j < numperiods; ++j)
        o << ".";
      
      o << " " << std::setw(12) << filter(statList_m[i]->value()) << std::endl;
    }
}

} // namespace Pooma

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Statistics.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
