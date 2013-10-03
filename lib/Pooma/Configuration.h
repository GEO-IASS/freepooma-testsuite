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

#ifndef POOMA_POOMA_CONFIGURATION_H
#define POOMA_POOMA_CONFIGURATION_H

/*-----------------------------------------------------------------------------
// Classes:
// (none; this file is used to properly include configuration-option files)
//---------------------------------------------------------------------------*/

/** @file
 * @ingroup Pooma
 * @brief
 * Main configuration file.
 *
 * Config settings are found in a file 'lib/$(SUITE)/PoomaConfiguration.h';
 * but this file does not exist for all compiler environments.  This file
 * is used to determine if that other configuration file should be included.
 * If that file should be used, POOMA should also be compiled with a -I
 * statement to look for include files in 'lib/$(SUITE)'.
 *
 * On systems where there is no configure-generated include file, define
 * the symbol POOMA_NO_POOMA_CONFIGURATION_FILE, and put the actual settings
 * for various symbols in another location.
 */

#ifndef POOMA_NO_POOMA_CONFIGURATION_FILE
# include "PoomaConfiguration.h"
#endif


#endif     /* POOMA_POOMA_CONFIGURATION_H */

/* ACL:rcsinfo */
/* ----------------------------------------------------------------------
 * $RCSfile: Configuration.h,v $   $Author: richard $
 * $Revision: 1.8 $   $Date: 2004/11/01 18:17:03 $
 * ----------------------------------------------------------------------
 */
/* ACL:rcsinfo */
