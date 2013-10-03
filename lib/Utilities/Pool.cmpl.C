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

// include files
#include "Utilities/Pool.h"
#include "Utilities/PAssert.h"


//----------------------------------------------------------------------
//
// Make a new pool with a given size block.
//
//----------------------------------------------------------------------

Pool::Pool(size_t sz)
  :
  // The first one. Start out with nothing there.
  head_m(0),
  // The number of outstanding allocs.
  outstandingAllocs_m(0),
  // The size of each block
  bsize_m(roundToAlign(sz)),     
  // Number of blocks
  nblock_m(blocksInPage(bsize_m))
{
}

//----------------------------------------------------------------------
//
// Make in invalid pool using the null ctor.
//
//----------------------------------------------------------------------

Pool::Pool()
  :
  // The first one. Start out with nothing there.
  head_m(0),
  // The number of outstanding allocs.
  outstandingAllocs_m(0),
  // The size of each block
  bsize_m(0),
  // Number of blocks
  nblock_m(0)
{
}

//----------------------------------------------------------------------
//
// Delete a pool.
//
//----------------------------------------------------------------------

Pool::~Pool()
{
  PInsist(outstandingAllocs_m==0,"Not all of the pooled memory was freed!");

  // Loop over the allocated chunks.
  for (std::vector<char*>::iterator p=chunks_m.begin(); p!=chunks_m.end(); ++p)
    // Delete it.
    delete [] *p;
}

//----------------------------------------------------------------------
//
// Grow a Pool.
//
//----------------------------------------------------------------------

void Pool::grow()
{
  // Calculate the size to allocate.
  size_t alloc_this;

  // If we need more than a page, allocate that
  if ( bsize_m>page )
    alloc_this = bsize_m;
  else
    // Otherwise allocate a page.
    alloc_this = page;

  // Allocate a chunk.
  char *start = new char[alloc_this]; 

  // Put it in the list of things to delete.
  chunks_m.push_back(start);

  // Get a pointer to the last one in the chunk.
  char *last  = start + (nblock_m-1)*bsize_m; 

  // For all but the last one
  for (char *p=start; p!=last; p+=bsize_m)  
    // point to the next.
    ((Link*)p)->next_m = (Link*)(p+bsize_m);  

  // The last points to the current head of the list.
  ((Link*)last)->next_m = head_m;		  

  // Reset the head to the first in this chunk.
  head_m = (Link*)start;			  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Pool.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
