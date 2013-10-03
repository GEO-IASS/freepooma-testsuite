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

#ifndef POOMA_POOMA_PARTICLES_H
#define POOMA_POOMA_PARTICLES_H

/** @file
 * @ingroup Pooma
 * @brief
 * A one-stop-shopping header file that sets up everything one needs to use
 * all Pooma II particle objects and layouts.
 */

// Include files for the basic particle base class

#include "Particles/Particles.h"

// Include file for DynamicArray, used for particle attributes:

#include "Pooma/DynamicArrays.h" 

// Include file for canned particles traits structs

#include "Particles/CommonParticleTraits.h"

// Include files for the particle layout classes

#include "Particles/SpatialLayout.h"
#include "Particles/UniformLayout.h"

// Include files for particle boundary conditions

#include "Particles/AbsorbBC.h"
#include "Particles/KillBC.h"
#include "Particles/PeriodicBC.h"
#include "Particles/ReflectBC.h"
#include "Particles/ReverseBC.h"

// Interpolators

#include "Particles/InterpolatorNGP.h"
#include "Particles/InterpolatorCIC.h"
#include "Particles/InterpolatorSUDS.h"

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Particles.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
