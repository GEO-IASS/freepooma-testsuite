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
// Includes
//-----------------------------------------------------------------------------

// Pooma

#include "IO/FileSetReader.h"    // Class definition
#include "IO/ByteOrder.h"        // reverseBytes functions
#include "Utilities/PAssert.h"   // PInsist

#include "Array/Array.h"
#include "Engine/RemoteEngine.h"
#include "Engine/CompressibleBrick.h"

// System

#include <iostream>

#if !POOMA_NO_IOS_HEADER
# include <ios>
#endif

//-----------------------------------------------------------------------------
// OffsetData member function
//-----------------------------------------------------------------------------

template <int Dim>
template <class T>
inline void FileSetReader<Dim>::OffsetData<T>::reverseBytes()
{
  for (int i = 0; i < 6*Dim; ++i)
    ::reverseBytes(nodedata[i]);
  ::reverseBytes(offset);
  ::reverseBytes(compressedValue);
}

//-----------------------------------------------------------------------------
// FileSetReader template member functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// readOffsetData
// 
// Encapsulate the reading of the OffsetData struct from the .offset
// file. This allows us to reuse the same code amongst all
// instantiations of read having the same Dim, T, simplifying the code
// and avoiding code bloat.
//-----------------------------------------------------------------------------

template <int Dim>
template <class T>
bool FileSetReader<Dim>::
readOffsetData(OffsetData<T> &odata, Interval<Dim> &dom)
{
  PInsist(mycontext_m == iocontext_m, 
          "This function should only be called on an IO context");

  // Read the entry from the offset file for this field and
  // calculate the domain of the patch. Return the (byte-corrected)
  // OffsetData struct and the domain.

  const int sz = sizeof(OffsetData<T>);
  foffset_m.read((char*)&odata, sz);

  if (foffset_m.gcount() != sz) 
    {
      errorMsg_m = "Incomplete .offset record";
      return false;
    }

  // Fix the byte ordering, if necessary

  if (bytesReversed_m) odata.reverseBytes();

  // Calculate the domain of the patch

  int first, stride, length;
  for (int d = 0; d < Dim; ++d)
    {
      first  = odata.nodedata[6*d+1];
      stride = odata.nodedata[6*d+2];
      length = odata.nodedata[6*d+3];
      dom[d] = Interval<1>(first,first+(length-1));
    }

  return true;
}

//-----------------------------------------------------------------------------
// readFieldData
// 
// template function to encapsulate the actual reading of data into
// the local engine. This code is only run on the IO context (and only
// for nodes that belong to the IO context), and it reads directly
// into the compressible brick engine.
//
// By encapsulating into a separate template function, we avoid code
// bloat by having a different version for all the various engine tags
// for which read is instantiated.
//-----------------------------------------------------------------------------

template <int Dim>
template <class T>
bool FileSetReader<Dim>::
readFieldData(const DiskNode<Dim> &dnode,
              const Engine<Dim, T, CompressibleBrick> &e, 
              OffsetData<T> &odata)
{
  PInsist(mycontext_m == iocontext_m && mycontext_m == dnode.context_m,
          "This should only be called for nodes on the IO context.");

  if (odata.u.isCompressed)
    {
      // If compressed, simply set the compressed value. 

      e.compressedReadWrite() = odata.compressedValue;
    }
  else
    {
      // The data needs to be read from the ".data" file.

      // Get a pointer to the beginning of the uncompressed buffer
      // (dataBlock() causes an uncompress). (It is important to keep
      // a copy of the datablock on the stack, rather than using
      // dataBlock().beginPointer() directly, as the latter produces a
      // temporary that recompresses as soon as the statement is
      // finished, invalidating the buffer pointer.)

      DataBlockPtr<T> dbptr = e.dataBlock();
      T *buffptr = dbptr.beginPointer();

      // Find the specified position in the file and read the data.
      // NOTE: the offset in the OffsetData struct is a number of
      // elements (i.e. a pointer offset), not the number of bytes.

      fdata_m.seekg(odata.offset*sizeof(T));
      const int sz = dnode.domain_m.size();
      fdata_m.read((char*)buffptr, sz*sizeof(T));

      if (fdata_m.gcount() != sz*sizeof(T)) 
        {
          errorMsg_m = "Incomplete .data record";
          return false;
        }

      // Fix the byte order if necessary. 

      if (bytesReversed_m)
        {
          for (T *pt = buffptr; pt < buffptr+sz; ++pt)
            {
              reverseBytes(*pt);
            }
        }

    } // !compressed
  
  return true;
}

//-----------------------------------------------------------------------------
// read(Array_t &)
//
// The general read interface for arrays
//-----------------------------------------------------------------------------

template <int Dim>
template <class T, class ETag>
bool FileSetReader<Dim>::read(Array<Dim,T,ETag> &a)
{
  // If we've read everything, return false

  if (currentRecord_m == numRecords_m) return false;

  // If this is the first "field" of a record, read the layout.

  bool success;

  if (currentField_m == 0)
    {
      success = diskLayout_m.read(); 
      if (!success) return false; // diskLayout broadcasts its return
    }

  // Read the field ID from the file and check that it agrees with our
  // field position. Only check the value on the IO context. 

  int fieldID;

  if (mycontext_m == iocontext_m)
    {
      fieldID = readFieldID();
    }

  // Broadcast the value and check for errors

  PAssert(RemoteProxy<int>(fieldID).value() == currentField_m);
  
  // Get the number of patches from the layout and check that it
  // agrees with the .meta file.

  int numpatches = diskLayout_m.allNodes().size();

  if (mycontext_m == iocontext_m)
    {
      success = 
        (numpatches == pmetafile_m->patchesPerRecord()[currentRecord_m]);
    }

  PAssert(RemoteProxy<bool>(success).value());

  // Next, loop through all of the patches, reading the offset data,
  // performing some error checking, reading the data patch, and
  // assigning the patch to the input array.

  // NOTE: We only do patch-by-patch error checking on the various
  // reads with PAssert, so bad data files will probably just crash a
  // production code. Unfortunately, well behaved error-checking
  // requires either making copies of data or sending several message
  // for each patch.

  for (int i = 0; i < numpatches; ++i)
    {
      // Get the next "node"

      const DiskNode<Dim> &dnode = diskLayout_m.allNodes()[i];

      // Get the OffsetData and the domain of the patch.  For arrays,
      // this should match the domain in the current dnode. This is
      // only done on the reading contexts.

      OffsetData<T> odata;
      Interval<Dim> dom;

      if (mycontext_m == iocontext_m && dnode.context_m == iocontext_m)      
        {
          success = readOffsetData(odata, dom);
        }

      PAssert(RemoteProxy<bool>(success).value());

      if (mycontext_m == iocontext_m && dnode.context_m == iocontext_m)      
        {
          // Check the domain consistency...

          // For fields, the domain in the .layout file might be
          // slightly bigger than the domain of the patch. This
          // shouldn't be true for arrays, so for now we assert that
          // they're equal. 

          // QUESTION: Do we want to allow users to read Field data
          // into an Array. If so, we'd need to broadcast the actual
          // domain and use its size for the remote buffer.

          success = (dom == dnode.domain_m);
        }

      PAssert(RemoteProxy<bool>(success).value());

      // Now that we've checked on the domain, construct a Remote
      // engine. The RemoteEngine requires special construction since
      // the owning context must be set, so we do this in two steps. 
      // (Since we're not using the domain, we could put this before
      // reading the offset file, but I wanted the structure to mirror
      // that of the Field version, which has to communicate the
      // actual .offset domain in some instances.)

      typedef Remote<CompressibleBrick> PatchTag_t;
      typedef Array<Dim, T, PatchTag_t> RemotePatchBuff_t;
      typedef typename RemotePatchBuff_t::Engine_t RemoteEngine_t;

      RemotePatchBuff_t rbuff;
      rbuff.engine() = RemoteEngine_t(iocontext_m, dnode.domain_m);

      // Next, initialize the rbuff data from either the compressed
      // value or the data from the .data file.

      if (mycontext_m == iocontext_m && dnode.context_m == iocontext_m)
        {
          // Read the data into the local engine. 

          success = readFieldData(dnode, rbuff.engine().localEngine(), odata);
        }

      PAssert(RemoteProxy<bool>(success).value());

      // Finally, we execute an array assignment into a view of the
      // input array. This statement must be executed on all
      // contexts. The evaluator will handle all communication.

      a(dnode.domain_m) = rbuff;

    } // for each patch...

  // Increment the field count. 

  ++currentField_m;

  // If that was the last field of the record, reset field count
  // and increment the record count.

  if (currentField_m == fieldsPerRecord_m) 
    {
      currentField_m = 0;
      ++currentRecord_m;
    }
  
  return true;
}

//-----------------------------------------------------------------------------
// readSubField(Field_t &)
//
// Separate the logic for reading a subfield into a separate function
// to make logic a little cleaner. This should only be called by
// read(Field_t &) as it requires some setup (like reading the next
// layout, etc.).
//
// The logic is pretty much identical to that for Array above, except
// for enforcing domain consistency between the patch and the layout.
//-----------------------------------------------------------------------------

template <int Dim>
template <class Mesh, class T, class ETag>
bool FileSetReader<Dim>::readSubField(Field<Mesh,T,ETag> &sf)
{
  // Read the field ID from the file and check that it agrees with our
  // field position. Only check the value on the IO context. 

  int fieldID;
  if (mycontext_m == iocontext_m)
    {
      fieldID = readFieldID();
    }

  // Broadcast the value and check for errors

  PAssert(RemoteProxy<int>(fieldID).value() == currentField_m);
  
  // Get the number of patches from the layout and check that it
  // agrees with the .meta file.

  int numpatches = diskLayout_m.allNodes().size();

  bool success;
  if (mycontext_m == iocontext_m)
    {
      success = 
        (numpatches == pmetafile_m->patchesPerRecord()[currentRecord_m]);
    }

  PAssert(RemoteProxy<bool>(success).value());

  // Loop over all the local patches, getting offset information and
  // checking that it agrees with the layout.

  for (int i = 0; i < numpatches; ++i)
    {
      // Get the next "node"

      const DiskNode<Dim> &dnode = diskLayout_m.allNodes()[i];

      // SubField nodes may have a different domain from the layout,
      // which is vert-centered. We'll use the following node to store
      // the correct domain information.

      DiskNode<Dim> fnode = dnode; 

      // Get the OffsetData and the domain of the patch.  This is only
      // done on the reading contexts.

      OffsetData<T> odata;
      Interval<Dim> dom;

      if (mycontext_m == iocontext_m && dnode.context_m == iocontext_m)      
        {
          readOffsetData(odata, dom);
        }

      PAssert(RemoteProxy<bool>(success).value());

      if (mycontext_m == iocontext_m && dnode.context_m == iocontext_m)      
        {
          // For fields, the domain in the .layout file might be
          // slightly bigger than the domain of the patch since the
          // layout stores vert domains, while the patches may have
          // data on smaller domains.

          if (dom != dnode.domain_m)
            {
              success = contains(dnode.domain_m,dom);
              fnode.domain_m = dom;
            }
        }

      PAssert(RemoteProxy<bool>(success).value());

      // Now that we've got the right domain, construct a Remote
      // engine. The RemoteEngine requires special construction since
      // the owning context must be set, so we do this in two steps.

      typedef Remote<CompressibleBrick> PatchTag_t;
      typedef Array<Dim, T, PatchTag_t> RemotePatchBuff_t;
      typedef typename RemotePatchBuff_t::Engine_t RemoteEngine_t;

      RemotePatchBuff_t rbuff;
      rbuff.engine() = RemoteEngine_t(iocontext_m, fnode.domain_m);

      // Next, initialize the rbuff data from either the compressed
      // value or the data from the .data file.

      if (mycontext_m == iocontext_m && dnode.context_m == iocontext_m)
        {
          // Read the data into the local engine. 

          success = readFieldData(fnode, rbuff.engine().localEngine(), odata);
        }

      PAssert(RemoteProxy<bool>(success).value());

      // Finally, we execute field assignment into a view of the
      // subfield.  This statement must be executed on all contexts.
      // The evaluator will handle all communication.

      sf(fnode.domain_m) = rbuff;

    } // for each patch ...

  return true;
}

template <int Dim>
template <class Mesh, class T, class ETag>
bool FileSetReader<Dim>::read(Field<Mesh,T,ETag> &f)
{
  // If we've read everything, return false

  if (currentRecord_m == numRecords_m) return false;

  // Find out the number of materials and centering points.

  int nMaterials = numMaterials(f);
  int nCentering = centeringSize(f);

  // Fields with subfields take one record each, so if we're not at
  // the beginning of a record, it is an error.

  int numSubFields = nMaterials * nCentering;
  if (currentField_m + numSubFields > fieldsPerRecord_m) return false;

  // Read the layout for the next record

  bool success;
  if (currentField_m == 0)
    {
      success = diskLayout_m.read();
      if (!success) return false; // diskLayout broadcasts its return
    }

  // Loop through each material and centering

  for (int m = 0; m < nMaterials; ++m)
    {
      for (int c = 0; c < nCentering; ++c)
        {
          // Create a subfield view of the field we're reading. 

          Field<Mesh,T,ETag> sf = subField(f,m,c);

          // Read the data for the subfield.

          readSubField(sf);

          // Advance to next "field" of the record.

          ++currentField_m;
        }
    }

  // Update the record position

  if (currentField_m == fieldsPerRecord_m)
    {
      currentField_m = 0;
      ++currentRecord_m;
    }

  return true;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetReader.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
