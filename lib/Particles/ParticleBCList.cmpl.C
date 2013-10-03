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
// ParticleBCList method implementations.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Particles/ParticleBCList.h"

//-----------------------------------------------------------------------------
// Constructor.
// ParticleBCList just has a default constructor; it initially has
// no ParticleBCs at all.
//-----------------------------------------------------------------------------

ParticleBCList::ParticleBCList()
{
}


//-----------------------------------------------------------------------------
// Copy Constructor.
// Copy contents of an existing ParticleBCList.
//-----------------------------------------------------------------------------

ParticleBCList::ParticleBCList(const ParticleBCList& model)
{
  // copy the BC list.
  Size_t i, n = model.bcList_m.size();
  bcList_m.reserve(n);
  for (i=0; i<n; ++i)
    bcList_m.push_back(model.bcList_m[i]);
}

//-----------------------------------------------------------------------------
// Destructor.
// ParticleBCList will delete all ParticleBCs that it owns when
// it is deleted.
//-----------------------------------------------------------------------------

ParticleBCList::~ParticleBCList()
{
  while (size() > 0)
    removeBC(size() - 1);
}


//-----------------------------------------------------------------------------
// Remove the ith ParticleBC from our list, deleting it in the process.
//-----------------------------------------------------------------------------

void ParticleBCList::removeBC(Size_t i)
{
  PAssert(i < bcList_m.size());
  delete bcList_m[i];
  bcList_m.erase(bcList_m.begin() + i);
  return;
}


//-----------------------------------------------------------------------------
// Print the ParticleBCList to the given ostream.
//-----------------------------------------------------------------------------

void ParticleBCList::print(std::ostream& o) const
{
  Size_t n = size();

  for (Size_t i = 0; i < n; ++i)
    bcList_m[i]->print(o);
	return;
}


// operator<< for ostream and ParticleBCList

std::ostream&
operator<<(std::ostream& o, const ParticleBCList& bcList)
{
  bcList.print(o);
  return o;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ParticleBCList.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
