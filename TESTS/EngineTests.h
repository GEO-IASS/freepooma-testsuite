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
// Functions:
// typesExist(tag)
// typesExist(engine)
// checkStorage(engine)
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_TESTS_ENGINETESTS_H
#define POOMA_ENGINE_TESTS_ENGINETESTS_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Engine/Engine.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//
// typesExist(tag)
// typesExist(engine)
//
// This function just performs a compile time check that the required typedefs
// and enums exist for a given engine.
//
//-----------------------------------------------------------------------------

template<class Tag>
inline void
typesExist(const Tag &)
{
  typedef Engine<1, double, Tag> E1_t;
  typedef Engine<2, int, Tag>    E2_t;
  typedef Engine<5, bool, Tag>   E5_t;

  // First check that Tag_t is Tag:

  typedef typename E1_t::Tag_t E1Tag_t;
  const bool tagPass = SameType<Tag, E1Tag_t>::same;
  CTAssert(tagPass);

  // Is there a domain?

  typedef typename E2_t::Domain_t E2Domain_t;
  CTAssert(E2Domain_t::dimensions == 2);

  //
  // It's not clear that Layout_t is or should be a requirement,
  // so this section is commented out:
  //
  // Is there a layout?

  //  typedef E5_t::Layout_t E5Layout_t;
  //  CTAssert(E5Layout_t::dimensions == 5);

  // Is Element_t defined and is it the same as the template parameter
  // T?  (We decided that it must be.  ElementRef_t is allowed to be
  // different.)

  typedef typename E5_t::Element_t E5Element_t;
  const bool elemPass = SameType<bool, E5Element_t>::same;
  CTAssert(elemPass);

  // enums:

  const int e5dim = E5_t::dimensions;
  CTAssert(e5dim == 5);
  const bool e5dataObject = E5_t::hasDataObject;
  const bool e5multiPatch = E5_t::multiPatch;
  const bool e5zeroBased  = E5_t::zeroBased;
  const bool e5dynamic    = E5_t::dynamic;
}

template<int Dim, class T, class Tag>
inline int
typesExist(const Engine<Dim, T, Tag> &)
{
  typedef Engine<Dim, T, Tag> Engine_t;

  // First check that Tag_t is Tag:

  typedef typename Engine_t::Tag_t EngineTag_t;
  const bool tagPass = SameType<Tag, EngineTag_t>::same;
  CTAssert(tagPass);

  // Is there a domain?

  typedef typename Engine_t::Domain_t EngineDomain_t;
  CTAssert(EngineDomain_t::dimensions == Dim);

  // Is Element_t defined and is it the same as the template parameter
  // T?  (We decided that it must be.  ElementRef_t is allowed to be
  // different.)

  typedef typename Engine_t::Element_t EngineElement_t;
  const bool elemPass = SameType<T, EngineElement_t>::same;
  CTAssert(elemPass);

  // enums:

  const int enginedim = Engine_t::dimensions;
  CTAssert(enginedim == Dim);
  const bool enginedataObject = Engine_t::hasDataObject;
  const bool enginemultiPatch = Engine_t::multiPatch;
  const bool enginezeroBased  = Engine_t::zeroBased;
  const bool enginedynamic    = Engine_t::dynamic;

  // use the above values to avoid warnings.

  return (enginedataObject || enginemultiPatch || enginezeroBased ||
	  enginedynamic) ? enginedim : 1;
}

//-----------------------------------------------------------------------------
//
// checkStorage(engine)
//
// This function performs some runtime checks on a read-write engine to make
// sure that it can store things.  It takes an engine that stores ints, or
// things that are equivalent to ints.  The engine is modified by this
// function.
//
//-----------------------------------------------------------------------------

template<class T, class Tester>
void
checkIntStore(T &check, Tester & tester)
{
  int checkI;

  check = -17;
  checkI = check;

  bool storeInt = (checkI == -17);
  tester.check("can engine store ints", storeInt);
}

template<class T, class Tag, class Tester>
inline void
checkStorage(Engine<1, T, Tag> &engine, Tester &tester)
{
  T store;
  checkIntStore(store, tester);

  typedef typename Engine<1, T, Tag>::Domain_t Domain_t;
  CTAssert(Domain_t::dimensions == 1);

  Domain_t domain = engine.domain();
  int i0;
  for (i0 = domain[0].first(); i0 < domain[0].last(); ++i0)
  {
    engine(i0) = 5 * i0;
  }

  bool passed = true;

  for (i0 = domain[0].first(); i0 < domain[0].last(); ++i0)
  {
    if (engine.read(i0) != 5 * i0)
    {
      passed = false;
      tester.out() << "storage failure at (" << i0 << ")"
		   << std::endl;
    }
  }

  tester.check("engine<1> storage test", passed);
}

template<class T, class Tag, class Tester>
inline void
checkStorage(Engine<2, T, Tag> &engine, Tester &tester)
{
  T store;
  checkIntStore(store, tester);

  typedef typename Engine<2, T, Tag>::Domain_t Domain_t;
  CTAssert(Domain_t::dimensions == 2);

  Domain_t domain = engine.domain();
  int i0, i1;
  for (i0 = domain[0].first(); i0 < domain[0].last(); ++i0)
  {
    for (i1 = domain[1].first(); i1 < domain[1].last(); ++i1)
    {
      engine(i0, i1) = 3 * i0 + 7 * i1;
    }
  }

  bool passed = true;

  for (i0 = domain[0].first(); i0 < domain[0].last(); ++i0)
  {
    for (i1 = domain[1].first(); i1 < domain[1].last(); ++i1)
    {
      if (engine.read(i0, i1) != 3 * i0 + 7 * i1)
      {
	passed = false;
	tester.out() << "storage failure at (" << i0 << "," << i1 << ")"
		     << std::endl;
      }
    }
  }

  tester.check("engine<2> storage test", passed);
}

template<class T, class Tag, class Tester>
inline void
checkViews(Engine<2, T, Tag> &engine, Tester &tester)
{
  typedef Engine<2, T, Tag> Engine_t;

  typedef typename Engine<2, T, Tag>::Domain_t Domain_t;
  CTAssert(Domain_t::dimensions == 2);

  Domain_t domain = engine.domain();

  int i0f = domain[0].first();
  int i0l = domain[0].last();
  int i1f = domain[1].first();
  int i1l = domain[1].last();

  Interval<1> i1 = Interval<1>(i0f + 1, i0l - 1);
  Interval<1> i2 = Interval<1>(i1f + 1, i1l - 1);

  Interval<2> sub1(i1, i2);

  typedef typename NewEngine<Engine_t, Interval<2> >::Type_t EngineSub1_t;

  EngineSub1_t engineSub1(NewEngineEngine<Engine_t,
			  Interval<2> >::apply(engine, sub1),
			  NewEngineDomain<Engine_t,
			  Interval<2> >::apply(engine, sub1));

  typesExist(engineSub1);
  checkStorage(engineSub1, tester);

  engineSub1(0, 0) = 42;
  tester.check("interval view aligned", engine(i0f + 1, i1f + 1) == 42);
}

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_TESTS_ENGINETESTS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EngineTests.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
