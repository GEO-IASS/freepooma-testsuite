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

//-----------------------------------------------------------------------------
// AttributeList method implementations.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/AttributeList.h"

//-----------------------------------------------------------------------------
// Constructor.
// AttributeList just has a default constructor; it initially has
// no attributes at all.
//-----------------------------------------------------------------------------

AttributeList::AttributeList()
{
}


//-----------------------------------------------------------------------------
// Destructor.
// AttributeList will delete all attributes that it owns when
// it is deleted.
//-----------------------------------------------------------------------------

AttributeList::~AttributeList()
{
  while (size() > 0)
    remove(size() - 1);
}


//-----------------------------------------------------------------------------
// Remove the Nth attrib from our list, deleting it in the process.
// Return success.
//-----------------------------------------------------------------------------

bool AttributeList::remove(Size_t n)
{
  if (n < size())
    {
      delete attribute(n);
      list_m.erase(list_m.begin() + n);
      return true;
    }

  return false;
}


//-----------------------------------------------------------------------------
// Print out AttributeList to an ostream.
//-----------------------------------------------------------------------------

void AttributeList::print(std::ostream& o) const
{
  Size_t n = size();
  for (Size_t i = 0; i < n; ++i)
		{
			attribute(i)->print(o);
		}
	return;
}

//-----------------------------------------------------------------------------
//
// When AttributeList is passed to an ostream, call the print method.
//
//-----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& o, const AttributeList& alist)
{
  alist.print(o);
  return o;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: AttributeList.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
