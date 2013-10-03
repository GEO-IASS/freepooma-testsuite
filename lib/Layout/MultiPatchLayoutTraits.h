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

#ifndef POOMA_LAYOUT_MULTIPATCHLAYOUTTRAITS_H
#define POOMA_LAYOUT_MULTIPATCHLAYOUTTRAITS_H

//-----------------------------------------------------------------------------
// Classes: 
//   MultiPatchLayoutTraits<LayoutTag,Dim>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Layout
 * @brief
 *     a traits class specifying the layout and layout view types for
 *     a particular layout tag.
 */

// namespace Pooma {

/**
 * General version of the MultiPatchLayoutTraits class:
 * All layouts that are to be used by the MultiPatch-engine need to define 
 * a layout tag and specialize this class to define their type and their view
 * type.
 */

template <class LayoutTag, int Dim>
struct MultiPatchLayoutTraits
{ };

// }; // namespace Pooma

#endif // POOMA_LAYOUT_MULTIPATCHLAYOUTTRAITS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MultiPatchLayoutTraits.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:54 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
