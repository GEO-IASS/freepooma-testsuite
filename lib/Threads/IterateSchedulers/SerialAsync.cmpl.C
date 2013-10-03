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
// This program was prepared by the Regents of the University of California
// at Los Alamos National Laboratory (the University) under Contract No. 
// W-7405-ENG-36 with the U.S. Department of Energy (DOE). The University has 
// certain rights in the program pursuant to the contract and the program 
// should not be copied or distributed outside your organization. All rights
// in the program are reserved by the DOE and the University. Neither the U.S.
// Government nor the University makes any warranty, express or implied, or
// assumes any liability or responsibility for the use of this software
//-----------------------------------------------------------------------------
// Class:
// IterateScheduler<SerialAsync>
// Iterate<SerialAsync>
// DataObject<SerialAsync>
//-----------------------------------------------------------------------------

/*
LIBRARY:
        SerialAsync

CLASSES: IterateScheduler

CLASSES: DataObject

CLASSES: Iterate

OVERVIEW 
        SerialAsync IterateScheduler is a policy template to create a
        dependence graphs and executes the graph respecting the
        dependencies without using threads. There is no parallelism,
        but Iterates may be executed out-of-order with respect to the
        program text.

-----------------------------------------------------------------------------*/

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
// Smarts classes for times when you want no threads but you do want
// dataflow evaluation.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/IterateSchedulers/SerialAsync.h"

namespace Smarts {

std::list<RunnablePtr_t> SystemContext::workQueueMessages_m;
std::list<RunnablePtr_t> SystemContext::workQueue_m;
#if POOMA_MPI
  const int SystemContext::max_requests;
  MPI_Request SystemContext::requests_m[SystemContext::max_requests];
  std::map<int, SystemContext::IteratePtr_t> SystemContext::allocated_requests_m;
  std::set<int> SystemContext::free_requests_m;
#endif
std::stack<int> IterateScheduler<SerialAsync>::generationStack_m;

}

//////////////////////////////////////////////////////////////////////

/***************************************************************************
 * $RCSfile: SerialAsync.cmpl.cpp,v $   $Author: richard $
 * $Revision: 1.5 $   $Date: 2004/11/01 18:17:08 $
 ***************************************************************************/
