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
// functions: 
//   pack
//   unpack
//-----------------------------------------------------------------------------

#ifndef POOMA_FUNCTIONS_PACK_UNPACK_H
#define POOMA_FUNCTIONS_PACK_UNPACK_H

/** @file
 * @ingroup Utilites
 * @brief
 * pack(field) and unpack(field, rcbp) are used to provide the user with a
 * local 1D view of all the data in a field not including guard layers that
 * belongs to the local processor.
 *
 * The local data is returned from pack() in a RefCountedBlockPtr which
 * can provide raw pointers to the data via the beginPointer member function.
 * Currently pack and unpack copy the data to and from a separate block of
 * memory.  If we need to later, we can perform various optimizations under
 * the hood.  For example, if there is only one local patch, and the domain
 * is correct (no guards and it doesn't have unused points because of the
 * centering), we could return the RefCountedBlockPtr to the underlying brick
 * data.  Given that this function will typically be applied to compressed
 * data, which we would need to uncompress anyway, the copy is probably not
 * going to be expensive.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/RefCountedBlockPtr.h"
#include "Engine/RemoteEngine.h"
#include "Pooma/Pooma.h"

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Full Description of pack/unpack.
//
// RefCountedBlockPtr<T> pack(Field<Mesh, T, Eng> f)
//   - copy all the data from f (excluding guard layers) that belongs to this
//     processor into a block and return a ref-counted pointer to the data.
// 
// unpack(Field<Mesh, T, Eng> f, RefCountedBlockPtr<T>)
//   - copy the packed data from the ref-counted block back into the field.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// pack
//
// PackLocalPatches
//   - functor that gets handed to LoopApplyEvaluator to copy the data
//
// pack()
//   - loops through local patches, computes the total size, allocates a
//     ref-counted block and loops again to copy the local patches into the
//     block.
//-----------------------------------------------------------------------------

template<class InputField>
struct PackLocalPatches
{
  typedef typename InputField::Element_t Element_t;

  PackLocalPatches(RefCountedBlockPtr<Element_t> block)
    : block_m(block)
  {
  }

  inline void operator()(const Element_t &t)
  {
    *block_m = t;
    ++block_m;
  }

  RefCountedBlockPtr<Element_t> block_m;
  int total_m;
};

template<class InputField>
RefCountedBlockPtr<typename InputField::Element_t>
pack(const InputField &field)
{
  Pooma::blockAndEvaluate();

  typedef typename InputField::Element_t Element_t;

  int size, i;
  size = 0;
  for (i = 0; i < field.numPatchesLocal(); ++i)
  {
    size += field.patchLocal(i).domain().size();
  }

  RefCountedBlockPtr<Element_t> ret(size);
  RefCountedBlockPtr<Element_t> current = ret;

  for (i = 0; i < field.numPatchesLocal(); ++i)
  {
    typedef typename Patch<InputField>::Type_t PatchField_t;
    PatchField_t patch = field.patchLocal(i);
    PackLocalPatches<PatchField_t> packFunctor(current);
    EngineBlockSerialize::apply(packFunctor, patch, patch.domain());
    current += patch.domain().size();
  }

  return ret;
}

//-----------------------------------------------------------------------------
// unpack
//
// opposite of pack.
//-----------------------------------------------------------------------------

template<class InputField>
struct UnPackLocalPatches
{
  typedef typename InputField::Element_t Element_t;

  UnPackLocalPatches(RefCountedBlockPtr<Element_t> block)
    : block_m(block)
  {
  }

  inline void operator()(Element_t &t)
  {
    t = *block_m;
    ++block_m;
  }

  RefCountedBlockPtr<Element_t> block_m;
  int total_m;
};

template<class InputField, class T>
void
unpack(const InputField &field, RefCountedBlockPtr<T> block)
{
  Pooma::blockAndEvaluate();

  int i;

  RefCountedBlockPtr<T> current = block;

  for (i = 0; i < field.numPatchesLocal(); ++i)
  {
    typedef typename Patch<InputField>::Type_t PatchField_t;
    PatchField_t patch = field.patchLocal(i);
    UnPackLocalPatches<PatchField_t> unpackFunctor(current);
    EngineBlockSerialize::apply(unpackFunctor, patch, patch.physicalDomain());
    current += patch.physicalDomain().size();
  }
}


#endif // POOMA_FUNCTIONS_PACK_UNPACK_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PackUnpack.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:49 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
