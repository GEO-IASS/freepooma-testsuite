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
// non-inline definitions for PatchSizeSyncer
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include "Tulip/Messaging.h"
#include "Pooma/Pooma.h"
#include "Tulip/PatchSizeSyncer.h"
#include "Tulip/RemoteProxy.h"
#include "Tulip/CollectFromContexts.h"

#include <utility>
#include <algorithm>

namespace Pooma {

//-----------------------------------------------------------------------------
// PatchSize constructor & destructor...
//-----------------------------------------------------------------------------

PatchSizeSyncer::PatchSizeSyncer(int contextKey, Grid_t &localGrid)
  : myContext_m(Pooma::context()),
    numContexts_m(Pooma::contexts()),
    localKey_m(contextKey),
    localGrid_m(localGrid)
{ 
  if (myContext_m == 0) gridList_m.reserve(numContexts_m);
}

PatchSizeSyncer::~PatchSizeSyncer()
{
  PAssert(myContext_m == 0 || gridList_m.size() == 0);
  for (int i = 0; i < gridList_m.size(); ++i)
    delete gridList_m[i].second;
}

//-----------------------------------------------------------------------------
// PatchSizeSyncer::calcGlobalGrid
//
// Does a reduction of the grids, sending all the local grids to
// context 0, re-normalizing them to form a global grid, and broadcasting
// this global grid back to all contexts.
//-----------------------------------------------------------------------------

namespace {

// Functor used to sort the grid list after context 0 receives
// the remote grids.

struct ElemCompare
{
  typedef std::pair<int,Grid<1>*> Elem_t;
  inline bool operator()(const Elem_t &l, const Elem_t &r)
  {
    return l.first < r.first;
  }
};

}

void PatchSizeSyncer::calcGlobalGrid(Grid_t &globalGrid)
{
#if POOMA_MESSAGING

  Grid<1> result;

  CollectFromContexts<std::pair<int, Grid_t> > collection
	(std::make_pair(localKey_m,localGrid_m));
  if (myContext_m == 0)
    {
      // The grid list is full. We sort it and then renormalize the
      // domains to make them globally consistent.  The
      // renormalization is done by looking through the list and
      // correcting each grid by adding to it the number of elements
      // have been added on the previous grids. We simultaneously
      // calculate the total number of points, needed to form the
      // global result.

      for (int j = 0; j < numContexts_m; ++j)
        gridList_m.push_back(Elem_t(collection[j].first,
				    new Grid_t(collection[j].second)));

      std::sort(gridList_m.begin(),gridList_m.end(),ElemCompare());

      int total_points = gridList_m[0].second->size() - 1;

      for (int j = 1; j < numContexts_m; ++j)
	{
	  Grid<1> &l = *gridList_m[j-1].second;
	  Grid<1> &r = *gridList_m[j].second;
	  r += l.last() - r.first();
	  total_points += r.size() - 1;
	}

      ++total_points; // for the final endpoint.

      // Finally, construct a composite Grid representing the global layout.

      IndirectionList<int> glist(total_points);

      int k = 0;
      for (int j = 0; j < numContexts_m; ++j)
	for (int i = 0; i < gridList_m[j].second->size() - 1; ++i)
	  glist(k++) = (*gridList_m[j].second)(i);

      glist(k) = gridList_m[numContexts_m-1].second->last();
      result = Grid<1>(glist);
    }

  // Broadcast the result.

  RemoteProxy<Grid<1> > broadcast(result,0);
  globalGrid = Grid<1>(broadcast.value());

#else  // !POOMA_MESSAGING

  globalGrid = localGrid_m;

#endif // POOMA_MESSAGING
}


} // namespace Pooma


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PatchSizeSyncer.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:15 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
