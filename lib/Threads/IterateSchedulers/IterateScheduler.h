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
// SMARTS Runtime System: No Thread Iterate Scheduler
//
// This program was prepared by the Regents of the University of California
// at Los Alamos National Laboratory (the University) under Contract No. 
// W-7405-ENG-36 with the U.S. Department of Energy (DOE). The University has 
// certain rights in the program pursuant to the contract and the program 
// should not be copied or distributed outside your organization. All rights
// in the program are reserved by the DOE and the University. Neither the U.S.
// Government nor the University makes any warranty, express or implied, or
// assumes any liability or responsibility for the use of this software
//
#ifndef ITERATE_SCHEDULER_H
#define ITERATE_SCHEDULER_H

/** @file
 * @ingroup IterateSchedulers
 * @brief
 * Templates for the scheduler classes.
 *
 * This is sort of like abstract base classes, since it doesn't 
 * implement anything and you can't build one of these directly.
 * They are implemented by specializing for a tag class like
 * Stub or SerialAsync.
 */ 

namespace Smarts {

template<class T>
class IterateScheduler;

template<class T>
class Iterate;

template<class T>
class DataObject;

} // close namespace Smarts

#endif // ITERATE_SCHEDULER_H
