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
// Pool
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * A class for maintaining a large chunks of memory and handing
 * out small blocks very quickly.
 *
 * Intended to be used in new and delete operators of small classes.
 */

#ifndef POOMA_UTILITIES_POOL_H
#define POOMA_UTILITIES_POOL_H

///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include <stddef.h>
#include <string.h>
#include <vector>


/**
 * A Pool maintains a set of page-sized chunks of memory, and hands out 
 * small blocks very quickly.  It does this by considering the large 
 * chunks to be a set of small blocks and connecting the blocks in a
 * singly linked list.  When asked to hand over a block, it returns
 * first one in the list.  When a block is handed back, it goes at
 * the front of the list.
 * 
 * The intent is that a user class will have a static Pool and implement 
 * new and delete operators that use the pool to get memory.  For example,
 * if the user is building a lightweight Node class, the relevant parts
 * would look like:
 *
 * class Node
 * {
 * public:
 *
 *   // Get memory from the pool.
 *   // Don't need to pass the size_t to the pool because
 *   // that is in the ctor for the pool.
 *   void *operator new(size_t) { return pool_s.alloc() }
 *
 *   // Return memory to the pool.
 *   // Once again, the size_t is already in the pool.
 *   void *operator delete(void *p, size_t) { pool_s.free(p); }
 *
 * private:
 *
 *   // The static memory pool.
 *   static Pool pool_s;
 * };
 * 
 * Then in the Node.cpp file you would have:
 *
 * // Initialize the memory pool with the size of the blocks it will manage.
 * Pool Node::pool_s(sizeof(Node));
 */ 

class Pool 
{
public: 

  // Make a new pool with a given block size.
  Pool(size_t sz);

  // Make an invalid Pool.  Don't try to use it, but you
  // can construct a new one on top of it.
  Pool();

  // Delete the pool and the chunks it has allocated.
  ~Pool();

  // Allocate a block from the pool.
  inline void* alloc()
    {
      // Record an allocation.
      outstandingAllocs_m += 1;

      // If the free list is empty, get more memory.
      if ( head_m==0 )
	grow();

      // Get the first block.  We'll return this.
      Link *p = head_m;

      // Make the next one the new head of the list.
      // We can't do head_m = p->next_m since p will soon be treated
      // as something other than a Link.  By doing this assignment
      // with memcpy, we ensure that p->next_m will be read before
      // it is clobbered. 
      memcpy(&head_m, &p->next_m, sizeof(head_m));
      
      // Return the requested block.
      return p;
    }

  // Release a block to the pool.
  inline void free(void *b)
    {
      // Record a free.
      outstandingAllocs_m -= 1;

      // Cast the pointer to the right type.
      Link *p = (Link*)b;

      // Make it point to the current head of the free list.
      p->next_m = head_m;

      // Make the next one the new head of the list.
      // We can't do head_m = p->next_m since p will soon be treated
      // as something other than a Link.  By doing this assignment
      // with memcpy, we ensure that p->next_m will be read before
      // it is clobbered. 
      memcpy(&p->next_m, &head_m, sizeof(head_m));

      // Make it the head of the free list.
      head_m = p;
    }

private:

  // The Pool builds a linked list through each allocated block.
  // The links are of type Link.
  struct Link { Link *next_m; };

  // Some enums for calculating sizes of things.

  // page: The size of the large chunks to allocate.
  // This number is chosen to let the malloc fit in a single
  // page on most machines.
  enum { page=4096-8 };

  // align: The number of bytes to align the block on.
  // 8 means align on double words.
  // This must be a power of two.
  enum { align = 8 };

  // alignMask: masks of the bits that aren't aligned.
  enum { alignMask = align-1 };

  // Calculate the number of blocks in a page.
  static inline int blocksInPage(size_t sz)
    {
      return (page>sz)?(page/sz):1;
    }

  // Given a size, round up to an aligned size.
  static inline size_t roundToAlign(size_t s)
    {
      if (s)
	s = (s & ~alignMask) + ((s&alignMask)?align:0);
      else
	s = align;
      return s;
    }

  // Allocate another chunk and put its blocks in the free list.
  void grow();

  // The first one.
  Link *head_m;

  // The number of blocks in the user's hands.
  int outstandingAllocs_m;
  
  // How big is each block.
  size_t bsize_m;			

  // How many to allocate at once.
  size_t nblock_m;		

  // The currently allocated chunks.
  std::vector<char*> chunks_m;
};

//////////////////////////////////////////////////////////////////////

#endif // POOL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Pool.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
