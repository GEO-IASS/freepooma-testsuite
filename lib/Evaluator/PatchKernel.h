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
// PatchKernel
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_PATCHKERNEL_H
#define POOMA_EVALUATOR_PATCHKERNEL_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Evaluator
 * @brief
 * A PatchKernel encapsulates performing operations on a patch of an
 * expression.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Threads/PoomaSmarts.h"
#include "Evaluator/InlineEvaluator.h"
#include "Evaluator/EvaluatorTags.h"
#include "Evaluator/RequestLocks.h"
#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Pooma/Configuration.h"
#include "Threads/PoomaCSem.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<class A1,class Function>
class PatchKernel
  : public Pooma::Iterate_t
{
public:
  PatchKernel(const A1& a1, const Function& function, bool write = true)
    : Pooma::Iterate_t(Pooma::scheduler()),
      write_m(write), a1_m(a1), function_m(function)
  {
    DataObjectRequest<BlockAffinity> getAffinity;
    hintAffinity(engineFunctor(a1_m.engine(),getAffinity));

    // Request locks
    // currently we ignore write_m, because I'm not sure if iterates
    // will run if they don't have a write lock requested.

    DataObjectRequest<WriteRequest> writeReq(*this);
    engineFunctor(a1_m.engine(),writeReq);

  }	      

  virtual ~PatchKernel()
  {

    DataObjectRequest<WriteRelease> writeReq;
    engineFunctor(a1_m.engine(),writeReq);
  }

  virtual void run()
  {
    function_m.apply(a1_m);
  }

private:

  bool write_m;
  A1 a1_m;
  Function function_m;
};

template<class A1, class A2, class Function>
class PatchKernel2
  : public Pooma::Iterate_t
{
public:
  PatchKernel2(const A1 &a1, const A2 &a2,
	      const Function &function)
    : Pooma::Iterate_t(Pooma::scheduler()),
      a1_m(a1), a2_m(a2), function_m(function)
  {
    DataObjectRequest<BlockAffinity> getAffinity;
    hintAffinity(engineFunctor(a1_m.engine(),getAffinity));

    // Request locks

    // First make the write request.
    // The write request tag remembers the data block
    // for the left hand side so we can check if there is a stencil
    // on the right.

    DataObjectRequest<WriteRequest> writeReq(*this);
    engineFunctor(a1_m.engine(),writeReq);

    // Now make the read request.
    // Use the remembered write request block to check and see
    // if that block is used on the right.  If so, do a notify to the
    // iterated instead of a request of the data object. 

    DataObjectRequest<ReadRequest> readReq(writeReq);
    engineFunctor(a2_m.engine(),readReq);
  }	      

  virtual ~PatchKernel2()
  {

    // The write request remembers the data block it sees on the left
    // so that it can check and see if it finds it on the right.

    DataObjectRequest<WriteRelease> writeReq;
    engineFunctor(a1_m.engine(),writeReq);

    // The read request checks to see if the data object for the left
    // appears on the right.  If it does, it doesn't do a release for it
    // since a request wasn't generated above.

    DataObjectRequest<ReadRelease> readReq(writeReq);
    engineFunctor(a2_m.engine(),readReq);
  }

  virtual void run()
  {
    function_m.apply(a1_m, a2_m);
  }

private:

  A1 a1_m;
  A2 a2_m;
  Function function_m;
};

template<class A1, class A2, class A3, class Function>
class PatchKernel3 : public Pooma::Iterate_t
{
public:
  PatchKernel3(const A1 &a1, const A2 &a2, const A3 &a3,
	       const Function &function)
    : Pooma::Iterate_t(Pooma::scheduler()),
      a1_m(a1), a2_m(a2), a3_m(a3), function_m(function)
  {
    DataObjectRequest<BlockAffinity> getAffinity;
    hintAffinity(engineFunctor(a1_m.engine(),getAffinity));

    // Request locks

    // First make the write request.
    // The write request tag remembers the data block
    // for the left hand side so we can check if there is a stencil
    // on the right.

    DataObjectRequest<WriteRequest> writeReq(*this);
    engineFunctor(a1_m.engine(),writeReq);

    // Now make the read request.
    // Use the remembered write request block to check and see
    // if that block is used on the right.  If so, do a notify to the
    // iterated instead of a request of the data object. 

    DataObjectRequest<ReadRequest> readReq(writeReq);
    engineFunctor(a2_m.engine(),readReq);
    engineFunctor(a3_m.engine(),readReq);
  }	      

  virtual ~PatchKernel3()
  {
    // The write request remembers the data block it sees on the left
    // so that it can check and see if it finds it on the right.

    DataObjectRequest<WriteRelease> writeReq;
    engineFunctor(a1_m.engine(),writeReq);

    // The read request checks to see if the data object for the left
    // appears on the right.  If it does, it doesn't do a release for it
    // since a request wasn't generated above.

    DataObjectRequest<ReadRelease> readReq(writeReq);
    engineFunctor(a2_m.engine(),readReq);
    engineFunctor(a3_m.engine(),readReq);
  }

  virtual void run()
  {
    function_m.apply(a1_m, a2_m, a3_m);
  }

private:
  A1 a1_m;
  A2 a2_m;
  A3 a3_m;
  Function function_m;
};

//-----------------------------------------------------------------------------
// ParticleKernels
//-----------------------------------------------------------------------------

template<class Array, class Function>
class ParticleKernel
  : public Pooma::Iterate_t
{
public:
  ParticleKernel(const Array& array, const Function& function, int patchID,
		  bool write1)
    : Pooma::Iterate_t(Pooma::scheduler()),
      write1_m(write1), array_m(array), function_m(function),
      patchID_m(patchID)
  {
    hintAffinity(engineFunctor(array_m.engine(),
			       DataObjectRequest<BlockAffinity>()));

    // Request locks

    DataObjectRequest<WriteRequest> writeReq(*this);
    DataObjectRequest<ReadRequest>  readReq(writeReq);

    if (write1_m)
    {
      engineFunctor(array_m.engine(),writeReq);
    }
    else
    {
      engineFunctor(array_m.engine(),readReq);
    }
  }	      

  virtual ~ParticleKernel()
  {
    DataObjectRequest<WriteRelease> writeReq;
    DataObjectRequest<ReadRelease>  readReq(writeReq);

    if (write1_m)
    {
      engineFunctor(array_m.engine(),writeReq);
    }
    else
    {
      engineFunctor(array_m.engine(),readReq);
    }
  }

  virtual void run()
  {
    function_m.apply(array_m,patchID_m);
  }

private:
  bool write1_m;
  Array array_m;
  Function function_m;
  int patchID_m;
};

template<class Array, class Function>
class ParticleKernelBlock
  : public Pooma::Iterate_t
{
public:

  ParticleKernelBlock(const Array& array, const Function& function,
		      int patchID, bool write1, Pooma::CountingSemaphore *csem)
    : Pooma::Iterate_t(Pooma::scheduler()),
      array_m(array), function_m(function), patchID_m(patchID),
      write1_m(write1), csem_m(csem)
  {
    hintAffinity(engineFunctor(array_m.engine(),
			       DataObjectRequest<BlockAffinity>()));

    // Request locks

    DataObjectRequest<WriteRequest> writeReq(*this);
    DataObjectRequest<ReadRequest>  readReq(writeReq);

    if (write1_m)
    {
      engineFunctor(array_m.engine(),writeReq);
    }
    else
    {
      engineFunctor(array_m.engine(),readReq);
    }
  }	      

  virtual ~ParticleKernelBlock()
  {
    DataObjectRequest<WriteRelease> writeReq;
    DataObjectRequest<ReadRelease>  readReq(writeReq);

    if (write1_m)
    {
      engineFunctor(array_m.engine(),writeReq);
    }
    else
    {
      engineFunctor(array_m.engine(),readReq);
    }
  }

  virtual void run()
  {
    function_m.apply(array_m,patchID_m);
    csem_m->incr();
  }

private:
  bool write1_m;
  Array array_m;
  Function function_m;
  int patchID_m;
  Pooma::CountingSemaphore *csem_m;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<class Array1, class Array2, class Function>
class ParticleKernel2
  : public Pooma::Iterate_t
{
public:
  ParticleKernel2(const Array1 &array1,
		   const Array2 &array2,
		   const Function &function,
		   int id, bool write1, bool write2)
    : Pooma::Iterate_t(Pooma::scheduler()),
      array1_m(array1), array2_m(array2), function_m(function),
      id_m(id), write1_m(write1), write2_m(write2)
  {
    hintAffinity(engineFunctor(array1_m.engine(),
			       DataObjectRequest<BlockAffinity>()));

    // Request locks

    DataObjectRequest<WriteRequest> writeReq(*this);

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    DataObjectRequest<ReadRequest> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);

  }	      

  virtual ~ParticleKernel2()
  {
    DataObjectRequest<WriteRelease> writeReq;

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    DataObjectRequest<ReadRelease> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);
  }

  virtual void run()
  {
    function_m.apply(array1_m, array2_m, id_m);
  }

private:
  Array1 array1_m;
  Array2 array2_m;
  Function function_m;
  int id_m;
  bool write1_m;
  bool write2_m;
};


template<class Array1, class Array2, class Function>
class ParticleKernel2Block
  : public Pooma::Iterate_t
{
public:
  ParticleKernel2Block(const Array1 &array1,
		       const Array2 &array2,
		       const Function &function,
		       int id, bool write1, bool write2,
		       Pooma::CountingSemaphore *csem)
    : Pooma::Iterate_t(Pooma::scheduler()),
      array1_m(array1), array2_m(array2), function_m(function),
      id_m(id), write1_m(write1), write2_m(write2), csem_m(csem)
  {
    hintAffinity(engineFunctor(array1_m.engine(),DataObjectRequest<BlockAffinity>()));

    // Request locks

    DataObjectRequest<WriteRequest> writeReq(*this);

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    DataObjectRequest<ReadRequest> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);
  }	      

  virtual ~ParticleKernel2Block()
  {
    DataObjectRequest<WriteRelease> writeReq;

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    DataObjectRequest<ReadRelease> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);
  }

  virtual void run()
  {
    function_m.apply(array1_m, array2_m, id_m);
    csem_m->incr();
  }

private:
  Array1 array1_m;
  Array2 array2_m;
  Function function_m;
  int id_m;
  bool write1_m;
  bool write2_m;
  Pooma::CountingSemaphore *csem_m;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<class Array1, class Array2, class Array3, class Function>
class ParticleKernel3
  : public Pooma::Iterate_t
{
public:
  ParticleKernel3(const Array1 &array1,
		   const Array2 &array2,
		   const Array3 &array3,
		   const Function &function, int id,
		   bool write1, bool write2, bool write3)
    : Pooma::Iterate_t(Pooma::scheduler()),
      array1_m(array1), array2_m(array2), array3_m(array3),
      function_m(function), id_m(id),
      write1_m(write1), write2_m(write2), write3_m(write3)
  {
    hintAffinity(engineFunctor(array1_m.engine(),
			       DataObjectRequest<BlockAffinity>()));

    // Request locks

    DataObjectRequest<WriteRequest> writeReq(*this);

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    if (write3_m)
      engineFunctor(array3_m.engine(),writeReq);

    DataObjectRequest<ReadRequest> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);

    if (!write3_m)
      engineFunctor(array3_m.engine(),readReq);
  }	      

  virtual ~ParticleKernel3()
  {
    // Request locks

    DataObjectRequest<WriteRelease> writeReq;

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    if (write3_m)
      engineFunctor(array3_m.engine(),writeReq);

    DataObjectRequest<ReadRelease> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);

    if (!write3_m)
      engineFunctor(array3_m.engine(),readReq);
  }

  virtual void run()
  {
    function_m.apply(array1_m, array2_m, array3_m, id_m);
  }

private:
  Array1 array1_m;
  Array2 array2_m;
  Array3 array3_m;
  Function function_m;
  int id_m;
  bool write1_m;
  bool write2_m;
  bool write3_m;
};


template<class Array1, class Array2, class Array3, class Function>
class ParticleKernel3Block
  : public Pooma::Iterate_t
{
public:
  ParticleKernel3Block(const Array1 &array1,
		       const Array2 &array2,
		       const Array3 &array3,
		       const Function &function, int id,
		       bool write1, bool write2, bool write3,
		       Pooma::CountingSemaphore *csem)
    : Pooma::Iterate_t(Pooma::scheduler()),
      array1_m(array1), array2_m(array2), array3_m(array3),
      function_m(function), id_m(id),
      write1_m(write1), write2_m(write2), write3_m(write3), csem_m(csem)
  {
    hintAffinity(engineFunctor(array1_m.engine(),
			       DataObjectRequest<BlockAffinity>()));

    // Request locks

    DataObjectRequest<WriteRequest> writeReq(*this);

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    if (write3_m)
      engineFunctor(array3_m.engine(),writeReq);

    DataObjectRequest<ReadRequest> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);

    if (!write3_m)
      engineFunctor(array3_m.engine(),readReq);
  }	      

  virtual ~ParticleKernel3Block()
  {
    // Request locks

    DataObjectRequest<WriteRelease> writeReq;

    if (write1_m)
      engineFunctor(array1_m.engine(),writeReq);

    if (write2_m)
      engineFunctor(array2_m.engine(),writeReq);

    if (write3_m)
      engineFunctor(array3_m.engine(),writeReq);

    DataObjectRequest<ReadRelease> readReq(writeReq);

    if (!write1_m)
      engineFunctor(array1_m.engine(),readReq);

    if (!write2_m)
      engineFunctor(array2_m.engine(),readReq);

    if (!write3_m)
      engineFunctor(array3_m.engine(),readReq);
  }

  virtual void run()
  {
    function_m.apply(array1_m, array2_m, array3_m, id_m);
    csem_m->incr();
  }

private:
  Array1 array1_m;
  Array2 array2_m;
  Array3 array3_m;
  Function function_m;
  int id_m;
  bool write1_m;
  bool write2_m;
  bool write3_m;
  Pooma::CountingSemaphore *csem_m;
};

//////////////////////////////////////////////////////////////////////

#endif     // POOMA_EVALUATOR_PATCHKERNEL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PatchKernel.h,v $   $Author: richard $
// $Revision: 1.22 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
