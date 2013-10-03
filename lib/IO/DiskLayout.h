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
//   DiskLayout<Dim>
//-----------------------------------------------------------------------------

/** @file
 * @ingroup IO
 * @brief
 * DiskLayout<Dim> manages the reading of the .layout file in a
 * DiskField fileset and the communication between contexts
 * necessary for each reading process to have full "layout"
 * information.
 *
 * The information in the .layout file is actually
 * redundant and thus writing is handled by the same code that
 * writes the .offset file.
 *
 * NOTE: The current version assumes a single file-set being read
 * from context 0. It is mostly coded for the general case, but we
 * do not do the communication of the localNodes in order to
 * calculate the allNodes set and we also need to broadcast some
 * bools when things fail. 
 *
 * The file format is the same as that used by the POOMA r1
 * DiscField .layout file, designed by Bill Humphrey.
 */

#ifndef POOMA_IO_DISKLAYOUT_H
#define POOMA_IO_DISKLAYOUT_H

#include "Domain/Interval.h"

#include <fstream>  // file I/O
#include <vector>   // node lists
#include <string>

/**
 * Simple struct containing a context + domain, and a constructor 
 * that can construct these from an array of 6*Dim ints.
 */

template <int Dim>
struct DiskNode
{
  typedef Interval<Dim> Domain_t;
  typedef int DomainData_t[6*Dim];

  DiskNode();
  DiskNode(int context, const DomainData_t &dd);

  int context_m;
  Domain_t domain_m;
};

/**
 * Class encapsulating the reading of a "DiscField" .layout file. 
 * Handles byte-ordering correction, unlike the r1 DiscField
 */

template <int Dim>
class DiskLayout
{
public:

  //---------------------------------------------------------------------------
  // Typedefs and private types
  //---------------------------------------------------------------------------

  typedef Interval<Dim>       Domain_t;
  typedef DiskNode<Dim>       Node_t;
  typedef std::vector<Node_t> NodeList_t;

  //---------------------------------------------------------------------------
  // Constructors and destructor
  //---------------------------------------------------------------------------

  DiskLayout(const char *filesetname);

  ~DiskLayout();

  //---------------------------------------------------------------------------
  // DiskLayout interface
  //---------------------------------------------------------------------------

  // Open the layout file.

  bool open();

  // Read the next layout in the file.

  bool read();

  //---------------------------------------------------------------------------
  // Accessors
  //---------------------------------------------------------------------------

  // allNodes will be broadcast to all contexts.

  const NodeList_t &allNodes() const { return allNodes_m; }

  // localNodes will be empty on all but IO contexts.

  const NodeList_t &localNodes() const { return localNodes_m; }

  // The domain will be broadcast to all contexts.

  const Domain_t &domain() const { return domain_m; }

  // bytesReversed is only valid on the IO contexts.

  bool bytesReversed() const { return bytesReversed_m; }

private:
  
  //---------------------------------------------------------------------------
  // Utility functions
  //---------------------------------------------------------------------------

  // This is the function that actually reads data from the file. This
  // is only called on IO contexts, internally.

  bool readLocal();

  // Check that the set of domains given by allNodes_m completely
  // covers the total domain with no overlaps; i.e. that it is a valid
  // partitioning of the global domain.

  bool checkLayout();

  //---------------------------------------------------------------------------
  // Data
  //---------------------------------------------------------------------------

  // Name of the layout file

  std::string filename_m; 

  // File stream for reading the .layout file

  std::ifstream fin_m;

  // Layout information - local nodes, global nodes, and the total domain

  NodeList_t localNodes_m;
  NodeList_t allNodes_m;

  Domain_t domain_m;

  // The average number of blocks in each direction - used by
  // checkLayout, which builds a local multi-patch Array in order to
  // check layout consistency.

  int avgblocks_m[Dim];

  // Local context

  int mycontext_m;

  // IO context

  int iocontext_m;

  // Do we need to reverse the bytes? 

  bool bytesReversed_m;

  //---------------------------------------------------------------------------
  // Disabled functions
  //---------------------------------------------------------------------------

  DiskLayout(const DiskLayout &);
  DiskLayout &operator=(const DiskLayout &);
};

#endif // POOMA_IO_DISKLAYOUT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DiskLayout.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
