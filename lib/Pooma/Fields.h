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

#ifndef POOMA_POOMA_FIELDS_H
#define POOMA_POOMA_FIELDS_H

/** @file
 * @ingroup Pooma
 * @brief
 * A one-stop-shopping header file that sets up everything one needs to use
 * all POOMA Fields (and, by inclusion of Arrays.h, Arrays).
 */

// Include files

// Arrays:

#include "Pooma/Arrays.h"

// Field class:

#include "Field/Field.h"

// FieldEngines:

#include "Field/FieldEngine/FieldEngine.h"
#include "Field/FieldEngine/FieldEngine.ExprEngine.h"

// Meshes:

#include "Field/Mesh/NoMesh.h"
#include "Field/Mesh/UniformRectilinearMesh.h"
#include "Field/Mesh/RectilinearMesh.h"
#include "Field/Mesh/MeshFunctions.h"
#include "Field/Mesh/PositionFunctions.h"

// Relations:

#include "Field/Relations/ConstantFaceBC.h"
#include "Field/Relations/PosReflectFaceBC.h"
#include "Field/Relations/PeriodicFaceBC.h"

// Differential Operators:

#include "Field/DiffOps/Div.UR.h"
#include "Field/DiffOps/Div.h"

// Other stuff:

#include "Functions/Reductions.h"
#include "Field/PrintField.h"
#include "Field/FieldOperatorSpecializations.h"

// Field offsets:

#include "Field/FieldOffset.h"
#include "Field/DiffOps/FieldShiftEngine.h"
#include "Field/DiffOps/FieldOffsetReduction.h"

#endif // POOMA_POOMA_FIELDS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Fields.h,v $   $Author: richard $
// $Revision: 1.18 $   $Date: 2004/11/01 18:17:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
