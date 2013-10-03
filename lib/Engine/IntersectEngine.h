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
// Class:
// IntersectEngine
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_INTERSECTENGINE_H
#define POOMA_ENGINE_INTERSECTENGINE_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Engine
 * @brief
 * IntersectEngine provides a common interface for applying the intersector
 * object to various engines.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Engine/EngineFunctor.h"
#include "PETE/PETE.h"


/**
 * ExpressionApply<IntersectorTag<Intersector> >
 *
 * IntersectEngine is simple wrapper that is used to apply intersector
 * objects to engines.  It contains a reference to the intersector object
 * and for engines with multiple patches it should hand the engine back
 * to the intersector.  Typical use would look something like:
 *
 * IntersectEngine<Intersector> ie(intersector);
 * engineFunctor(eng, ie);
 *
 * This level of indirection allows us to short-circut intersection for
 * trivial engines and scalars, and the use of engineFunctor automatically
 * deals with expression engines.
 *
 * The return value for intersection is a boolean that is currently unused.
 * (The result of the intersection is stored in the intersector object.)
 */

template<class Inter>
struct IntersectorTag
{
  inline IntersectorTag(Inter &i) : intersector_m(i) { }

  Inter &intersector_m;
};


/**
 * The default behaviour for IntersectEngine is to simply return true.
 * We assert that the engine is not multi-patch.
 */

template<class Eng, class Intersect>
struct DefaultExpressionApply<Eng, IntersectorTag<Intersect> >
{
  typedef int Type_t;

  inline static
  Type_t apply(const Eng &,
	       const ExpressionApply<IntersectorTag<Intersect> > &)
  {
    // Engines that are multipatch must specialize this functor
    // to perform the correct intersection.
    CTAssert(!(Eng::multiPatch));
    return true;
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_ENGINE_INTERSECTENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IntersectEngine.h,v $   $Author: richard $
// $Revision: 1.15 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
