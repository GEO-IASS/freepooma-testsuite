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
// Classes
//  PatchSizeSyncer
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Tulip
 * @brief
 * PatchSizeSyncer is used to synchronize a set of Grid objects that are
 * used to represent a set of contiguous patches (e.g. in DynamicLayout).
 */

#ifndef POOMA_CHEETAH_PATCHSIZESYNCER_H
#define POOMA_CHEETAH_PATCHSIZESYNCER_H

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tulip/Messaging.h"
#include "Domain/Grid.h"

#include <vector>
#include <utility>

namespace Pooma {


/**
 * This class encapsulates the communication that must occur when
 * synchronizing the patch domains for a DynamicLayout.
 * PatchSizeSyncer is used by instantiating a version with a Grid<1>
 * object representing one's local patches and then calling
 * calcGlobalGrid with a grid object that will be filled with the
 * redistributed patches for the global Grid. The DynamicLayout can
 * then re-label the domains in its Nodes to match the new patch
 * distribution.
 */

class PatchSizeSyncer
{
public:
  //============================================================
  // Typedefs
  //============================================================

  typedef Grid<1> Grid_t;

  //============================================================
  // Constructor & Destructor
  //============================================================

  // This sets up the local data for the calculation.

  PatchSizeSyncer(int contextKey, Grid_t &localGrid);

  // Cleans up the grid list on context 0.

  ~PatchSizeSyncer();

  //============================================================
  // Mutators
  //============================================================

  // This gathers the local grids on context 0, renormalizes
  // the local domains, constructs a new global grid, broadcasts
  // it to all contexts, and returns it to the caller.

 void calcGlobalGrid(Grid_t &globalGrid);

private:

  //============================================================
  // Data
  //============================================================

  int myContext_m;

  int numContexts_m;

  int localKey_m;
  
  Grid_t localGrid_m;
  
  typedef std::pair<int,Grid_t *> Elem_t;
  
  std::vector<Elem_t> gridList_m;

  // PatchSizeSyncers should not be copied or assigned.

  PatchSizeSyncer(const PatchSizeSyncer &);
  PatchSizeSyncer &operator=(const PatchSizeSyncer &);

  // This is the tag used for message volleys. It is shared
  // by all PatchSizeSyncer objects.

  static int tag_s;

};

} // namespace Pooma


#if POOMA_MESSAGING

namespace Cheetah {

/**
 * Serialize specialization for std::pair<int,Grid<1> >
 */

template<>
class Serialize<CHEETAH, std::pair<int,Grid<1> > >
{
public:

  typedef Grid<1>                    Grid_t;
  typedef std::pair<int,Grid_t>      Pair_t;
  typedef Serialize<CHEETAH,Grid_t>  GridPacker_t;
  
  static inline int
  size(const Pair_t &a)
  {
    return sizeof(int) + GridPacker_t::size(a.second);
  }

  static inline int
  pack(const Pair_t &a, char *buffer)
  {
    *reinterpret_cast<int *>(buffer) = a.first;
    
    int nBytes = sizeof(int);
    nBytes += GridPacker_t::pack(a.second,buffer+nBytes);
        
    return nBytes;
  }

  static inline int
  unpack(Pair_t* &a, char *buffer)
  {
    int key = *reinterpret_cast<int *>(buffer);

    int nBytes = sizeof(int);
    
    Grid_t *pg;
    nBytes += GridPacker_t::unpack(pg, buffer+nBytes);
    
    a = new Pair_t(key,*pg);
    
    return nBytes;
  }

  static inline void
  cleanup(Pair_t* a)
  {
    delete a;
  }
};

} // namespace Cheetah

#endif // POOMA_MESSAGING

#endif // POOMA_CHEETAH_PATCHSIZESYNCER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PatchSizeSyncer.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
