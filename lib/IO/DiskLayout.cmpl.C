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
// Implementation of DiskLayout member functions.
//
// The DiskLayout template is preinstantiated for Dim = 1, 2, 3. If
// higher dimension are needed, modify the instantiation statements
// near the end of the file.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

#include "IO/DiskLayout.h"      // Class definition
#include "Pooma/Pooma.h"        // Pooma::context()
#include "Utilities/PAssert.h"  // PInsist
#include "IO/ByteOrder.h"       // reverseBytes()

// For checkLayout:

#include "Array/Array.h"
#include "Engine/CompressibleBrick.h"
#include "Engine/MultiPatchEngine.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/UniformGridLayout.h"
#include "Tulip/RemoteProxy.h"

#if POOMA_NO_IOS_HEADER
#include <iostream>
#else
#include <ios>
#endif

//-----------------------------------------------------------------------------
// BARRIER macro (for testing)
//-----------------------------------------------------------------------------

#define BARRIER

#ifndef BARRIER
# if POOMA_CHEETAH
#  define BARRIER Pooma::controller()->barrier()
# else
#  define BARRIER
# endif
#endif

//-----------------------------------------------------------------------------
// DiskNode implementation
//-----------------------------------------------------------------------------

// Default constructor

template <int Dim>
inline DiskNode<Dim>::DiskNode()
{ }

// This constructor is only used in this file, so inline it.

template <int Dim>
inline DiskNode<Dim>::DiskNode(int context, const DomainData_t &dd)
  : context_m(context)
{ 
  for (int d = 0; d < Dim; ++d)
    {
      int first  = dd[6*d+1];
      int stride = dd[6*d+2];
      int length = dd[6*d+3];
      PInsist(stride == 1, "Attempt to read non-unit-stride file");
      domain_m[d] = Interval<1>(first,first+(length-1));
    }
}

//-----------------------------------------------------------------------------
// DiskLayout implementation
//-----------------------------------------------------------------------------

template <int Dim>
DiskLayout<Dim>::DiskLayout(const char *fileset)
  : mycontext_m(Pooma::context()), iocontext_m(0), // for now!!!
    bytesReversed_m(false)
{
  if (mycontext_m == iocontext_m)
    {
      filename_m = std::string(fileset) + ".layout";
    }
}

template <int Dim>
bool DiskLayout<Dim>::open()
{
  bool success;

  if (mycontext_m == iocontext_m)
    {
      fin_m.open(filename_m.c_str(), std::ios::binary);
      if (!fin_m) 
        {
          success = false;
        }
      else
        {
          // Read the first domain from the file, check the first stride,
          // and figure out whether we need to reverse bytes or not.

          struct TestData
          {
            int n;
            int domainData[6];
          };

          TestData td;
          fin_m.read((char*)&td, sizeof(td));

          fin_m.seekg(0); // Reset the file to the beginning

          // Is the stride byte reversed?

          int teststride = td.domainData[2];
          char *pc = reinterpret_cast<char*>(&teststride);

          int lastbyte = sizeof(int)-1;

          PInsist(pc[0] == 1 || pc[lastbyte] == 1, 
                  "Domain has non-unit stride");

          if (teststride != 1) bytesReversed_m = true;

          success = true;
        }
    }
  return RemoteProxy<bool>(success);
}

template <int Dim>
DiskLayout<Dim>::~DiskLayout()
{
  if (mycontext_m == iocontext_m) fin_m.close();
}

template <int Dim>
bool DiskLayout<Dim>::readLocal()
{
  PAssert(mycontext_m == iocontext_m);

  localNodes_m.clear();
  allNodes_m.clear();
  domain_m = Domain_t();

  int numNodes;
  fin_m.read((char*)&numNodes, sizeof(numNodes));
  if (fin_m.gcount() != sizeof(int) || numNodes <= 0) return false;
  if (bytesReversed_m) reverseBytes(numNodes);

  localNodes_m.reserve(numNodes);
  allNodes_m.reserve(numNodes);

  int domainData[6*Dim];

  for (int i = 0; i < numNodes; ++i)
    {
      fin_m.read((char*)&domainData[0], sizeof(domainData));
      if (fin_m.gcount() != sizeof(domainData)) return false;

      if (bytesReversed_m)
        {
          for (int i = 0; i < 6*Dim; ++i)
            reverseBytes(domainData[i]);
        }

      localNodes_m.push_back(DiskNode<Dim>(mycontext_m, domainData));
    }

  // For now there is only one fileset, so allNodes == localNodes
  // THIS NEEDS FIXED FOR MULTIPLE FILESETS!!!

  allNodes_m = localNodes_m;

  // Calculate the total domain and check that the patches span the
  // domain. Also, calculate the average number of patches in each
  // direction.

  int imin[Dim];
  int imax[Dim];
  int iavg[Dim];

  int d;
  for (d = 0; d < Dim; ++d)
    {
      imax[d] = allNodes_m[0].domain_m[d].last();
      imin[d] = allNodes_m[0].domain_m[d].first();
      iavg[d] = imax[d] - imin[d];
    }

  for (d = 0; d < Dim; ++d)
    {
      for (unsigned int i = 1; i < allNodes_m.size(); ++i)
        {
          int dl = allNodes_m[i].domain_m[d].last();
          int df = allNodes_m[i].domain_m[d].first();
          if (dl > imax[d]) imax[d] = dl;
          if (df < imin[d]) imin[d] = df;
          iavg[d] += (dl - df);
        }
      iavg[d] /= allNodes_m.size();
    }

  for (d = 0; d < Dim; ++d)
    {
      domain_m[d]    = Interval<1>(imin[d], imax[d]);
      avgblocks_m[d] = domain_m[d].size() / iavg[d];
    }

  // Check that the nodes cover the domain with no overlap.

  return checkLayout();
}


template <int Dim>
bool DiskLayout<Dim>::read()
{
  typedef unsigned int uint;

  bool success = true;

  // If we're on the IO context, read the data.

  if (mycontext_m == iocontext_m) success = readLocal();

  // Broadcast the value of success, and if the read failed, everyone
  // returns.

  bool iosuccess = RemoteProxy<bool>(success);
  if (!iosuccess) return false;

  // Read was successfull on the IO context.

  // We need to broadcast allNodes and the domain to the other contexts. 
  // (Once we allow multiple filesets, we'll first have to construct
  // allNodes from the assorted localNodes calculations!!!)

  // Since iocontext_m == 0, for now we just broadcast the domain
  // list. We're probably making more copies than we need here, but
  // this is simple. If profiling shows this is significant, we'll
  // optimize it. Anyway, first we must construct the domain list:

  typedef std::vector<Domain_t> DomainList_t;
  DomainList_t domainList;

  if (mycontext_m == iocontext_m)
    {
      const int sz = allNodes_m.size();
      domainList.resize(sz);
      for (int ip = 0; ip < sz; ++ip) 
        {
          PAssert(allNodes_m[ip].context_m == 0); // for now
          domainList[ip] = allNodes_m[ip].domain_m;
        }
    }
    
  // Broadcast the data

  RemoteProxy<DomainList_t> domainListProxy(domainList);
  RemoteProxy<Domain_t> domainProxy(domain_m);

  // Unpack the broadcast data and initialize domain_m and allNodes_m.

  if (mycontext_m != iocontext_m)
    {
      // The domain is easy

      domain_m = domainProxy.value();

      // Get the domainList

      DomainList_t domainList = domainListProxy.value();
      const int sz = domainList.size();

      // Reset the allNodes array and initialize from domainList

      allNodes_m.clear();
      allNodes_m.resize(sz);
      for (int ip = 0; ip < sz; ++ip)
        {
          allNodes_m[ip].context_m = 0; // for now!!!
          allNodes_m[ip].domain_m = domainList[ip];
        }
    }

  return true;
}

// Check that the NodeList completely covers the global domain with no
// overlaps. This can be done cleverly using an Array<int> - we just
// initialize to 1 and then subtract 1 for each node-view. If any
// patches are missing, we'll get 1's in the final arrays, and if
// there are any overlaps, we'll get negative values in the final
// array. This does take some memory, though. In order to minimize
// memory I create a local multipatch array roughly the same number of
// domains in each direction as the input layout, and I use
// compressible patches. At the beginning and end the entire array
// should be compressed. If the domains don't overlap exactly, then
// some patches will uncompress during the calculation, but should
// recompress. There are obviously cases where most of the array could
// be uncompressed - for example, if the layout nodes form a
// checkerboard in the first half. Hopefully that won't happen too
// often. 

template <int Dim>
bool DiskLayout<Dim>::checkLayout()
{
  typedef CompressibleBrick                 PatchTag_t;
  typedef MultiPatch<UniformTag,PatchTag_t> ETag_t;
  typedef Array<Dim, signed char, ETag_t>   Array_t;

  bool success;

  if (mycontext_m == iocontext_m)
    {
      // Create a uniform grid layout with roughly the same
      // partitioning in each direction.

      Loc<Dim> blocks;
      for (int d = 0; d < Dim; ++d)
        {
          blocks[d] = avgblocks_m[d];
        }

      UniformGridPartition<Dim> partition(blocks);   
      UniformGridLayout<Dim> layout(domain_m, partition, ReplicatedTag());
      
      // Create our local array.

      Array_t domaincheck(layout);

      // Set it to 1 and then subtract 1 from each domain view.

      domaincheck = 1;

      typedef unsigned int uint;
      for (uint i = 0; i < allNodes_m.size(); ++i)
        {
          DiskNode<Dim> &dnode = allNodes_m[i];
          domaincheck(dnode.domain_m) -= 1;
        }

      // Check that all elements of domaincheck are zero

      success = all(domaincheck == 0);
    }

  // We should probably broadcast this result.

  return success;
}

//-----------------------------------------------------------------------------
// Template preinstantiations
//-----------------------------------------------------------------------------

template class DiskLayout<1>;
template class DiskLayout<2>;
template class DiskLayout<3>;

// Uncomment this if you need to read DiskFields with Dim > 3.

#if 0
template class DiskLayout<4>;
template class DiskLayout<5>;
template class DiskLayout<6>;
template class DiskLayout<7>;
#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DiskLayout.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.11 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
