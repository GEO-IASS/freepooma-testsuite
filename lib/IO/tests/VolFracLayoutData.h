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
// VolFracLayoutData.h - array definition containing a character dump
// of the VolFrac.layout file from Jean Marshall. 
//-----------------------------------------------------------------------------

#ifndef POOMA_IO_TEST_VOLFRACLAYOUTDATA_H
#define POOMA_IO_TEST_VOLFRACLAYOUTDATA_H

//-----------------------------------------------------------------------------
// Test data
//-----------------------------------------------------------------------------

// This array was taken from a dump of VolFrac.layout, supplied by
// Jean Marshall, and written on the SGI, which is BigEndian.

// Character dump of file VolFrac.layout

char VolFrac_layout_dump[] = {
  0x00,0x00,0x00,0x04, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x06, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x03, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x06, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x06, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x02, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x03, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x01, 
  0x00,0x00,0x00,0x06, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00
};

// This is the same data, but byte-reversed by hand.

char VolFrac_layout_dump_reversed[] = {
  0x04,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x06,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x03,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x06,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x06,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x02,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x03,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00, 0x01,0x00,0x00,0x00, 
  0x06,0x00,0x00,0x00, 0x00,0x00,0x00,0x00, 
  0x00,0x00,0x00,0x00
};

#endif // POOMA_IO_TEST_VOLFRACLAYOUTDATA_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: VolFracLayoutData.h,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
