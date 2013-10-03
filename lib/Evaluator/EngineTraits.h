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
// Tags:
//   EvaluatorTypeTag
//   EvaluatorCombineTag
//   MainEvaluatorTag
//   MultiPatchEvaluatorTag
//   SinglePatchEvaluatorTag
//   RemoteMultiPatchEvaluatorTag
//   RemoteSinglePatchEvaluatorTag
//   ScalarEngineTag
//   DistributionTraits
// Structs:
//   EvaluatorEngineTraits
//-----------------------------------------------------------------------------

#ifndef POOMA_EVALUATOR_ENGINETRAITS_H
#define POOMA_EVALUATOR_ENGINETRAITS_H

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

struct Brick;
struct BrickView;

struct CompressibleBrick;
struct CompressibleBrickView;

struct Dynamic;
struct DynamicView;

template<class LayoutTag, class PatchTag>
struct MultiPatch;
template<class LayoutTag, class PatchTag, int Dim2>
struct MultiPatchView;

template<class Functor>
struct IndexFunction;
template<int Dim2, class Functor>
struct IndexFunctionView;

template<class Eng, class Components>
struct CompFwd;

template<class A1,class A2>
struct IndirectionTag;

struct ConstantFunction;

template<class Expr> 
struct ExpressionTag;

template<class Tag>
struct Remote;

struct DistributedTag;
struct ReplicatedTag;

/** @file
 * @ingroup Evaluator
 * @brief
 * EvaluatorEngineTraits<EngineTag> associates evaluator types with engine 
 * tags.
 *
 * This struct must be specialized whenever a new engine-type is added.
 *
 * The type of evaluator is given by the required Evaluator_t typedef and,
 * currently, must be one of the following choices:
 *   - SinglePatchEvaluatorTag
 *   - MultiPatchEvaluatorTag
 *   - RemoteSinglePatchEvaluatorTag
 *   - RemoteMultiPatchEvaluatorTag
 *
 * depending on whether or not the engine consists of single or multiple
 * patches or whether or not it involves remote objects.
 */

//-----------------------------------------------------------------------------
// Special tag to represent a scalar.
//-----------------------------------------------------------------------------

struct ScalarEngineTag 
{
  ScalarEngineTag() {}
  ~ScalarEngineTag() {}

};


//-----------------------------------------------------------------------------
// The evaluator tags. 
//-----------------------------------------------------------------------------

// The most general evaluator.

struct MainEvaluatorTag 
{
  MainEvaluatorTag() {}
  ~MainEvaluatorTag() {}
};

// The evaluator for single-patch expressions involving no remote objects.

struct SinglePatchEvaluatorTag 
{ 
  SinglePatchEvaluatorTag() {}
  ~SinglePatchEvaluatorTag() {}
};

// The evaluator for multi-patch expressions involving no remote objects.

struct MultiPatchEvaluatorTag 
{ 
  MultiPatchEvaluatorTag() {}
  ~MultiPatchEvaluatorTag() {}
};

// The evaluator for single-patch expressions involving remote objects.

struct RemoteSinglePatchEvaluatorTag 
{ 
  RemoteSinglePatchEvaluatorTag() {}
  ~RemoteSinglePatchEvaluatorTag() {}
};

// The evaluator for multi-patch expressions involving remote objects.

struct RemoteMultiPatchEvaluatorTag 
{ 
  RemoteMultiPatchEvaluatorTag() {}
  ~RemoteMultiPatchEvaluatorTag() {}
};


//-----------------------------------------------------------------------------
// Functor tags to interface with PETE. EvaluatorTypeTag is used to discover
// the type of evaluator at the leaf of a parse tree. EvaluatorCombineTag is
// used to combine evaluators from the left and right sides of the expression
// to produce a single evaluator for the expression. 
//-----------------------------------------------------------------------------

struct EvaluatorTypeTag 
{
  EvaluatorTypeTag() {}
  ~EvaluatorTypeTag() {}
};

struct EvaluatorCombineTag 
{
  EvaluatorCombineTag() {}
  ~EvaluatorCombineTag() {}
};


//-----------------------------------------------------------------------------
// Base template and specializations of EvaluatorEngineTraits<EngineTag>. These
// associate evaluator types with engine types.
//-----------------------------------------------------------------------------

template<class EngineTag>
struct EvaluatorEngineTraits
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
};

// Single-patch evaluators.

template<>
struct EvaluatorEngineTraits<ScalarEngineTag>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<ConstantFunction>
{ 
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<class Functor>
struct EvaluatorEngineTraits<IndexFunction<Functor> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<int Dim2, class Functor>
struct EvaluatorEngineTraits<IndexFunctionView<Dim2, Functor> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<Brick>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<BrickView>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<CompressibleBrick>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<CompressibleBrickView>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<Dynamic>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<>
struct EvaluatorEngineTraits<DynamicView>
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

template<class A1,class A2>
struct EvaluatorEngineTraits<IndirectionTag<A1,A2> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}

  // This is actually going to be complicated.
  // This code is wrong in general.
  typedef SinglePatchEvaluatorTag Evaluator_t;
};

// Remote-single-patch evaluators.

template<class Tag>
struct EvaluatorEngineTraits<Remote<Tag> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef RemoteSinglePatchEvaluatorTag Evaluator_t;
};

// Multi-patch evaluators.

template<class LayoutTag, class PatchTag>
struct EvaluatorEngineTraits<MultiPatch<LayoutTag, PatchTag> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef MultiPatchEvaluatorTag Evaluator_t;
};

template<class LayoutTag, class PatchTag, int Dim2>
struct EvaluatorEngineTraits<MultiPatchView<LayoutTag, PatchTag, Dim2> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef MultiPatchEvaluatorTag Evaluator_t;
};

// Remote-multi-patch evaluators.

template<class LayoutTag, class Tag>
struct EvaluatorEngineTraits<MultiPatch<LayoutTag, Remote<Tag> > >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef RemoteMultiPatchEvaluatorTag Evaluator_t;
};

template<class LayoutTag, class Tag, int Dim2>
struct EvaluatorEngineTraits<MultiPatchView<LayoutTag, Remote<Tag>, Dim2> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef RemoteMultiPatchEvaluatorTag Evaluator_t;
};

// Must do some indirection to handle forwarding engines.

template<class Eng, class Components>
struct EvaluatorEngineTraits<CompFwd<Eng, Components> >
{
  EvaluatorEngineTraits() {}
  ~EvaluatorEngineTraits() {}
  typedef EvaluatorEngineTraits<typename Eng::Tag_t> ET;
  typedef typename ET::Evaluator_t Evaluator_t;
};

// Must traverse the parse tree to figure out what to do with expression
// engines.

template<class Expr>
struct EvaluatorEngineTraits<ExpressionTag<Expr> >
{
 EvaluatorEngineTraits() {}
 ~EvaluatorEngineTraits() {}
 typedef typename ForEach<Expr, EvaluatorTypeTag, 
   EvaluatorCombineTag>::Type_t Evaluator_t;
};


//-----------------------------------------------------------------------------
// DistributionTraits
//
// DistributionTraits contains information about the way an engine is
// distributed in multiple contexts.
//-----------------------------------------------------------------------------

template<class ETag>
struct DistributionTraits
{
  enum { remote = false };
  typedef ReplicatedTag LayoutTag_t;
};

template<class ETag>
struct DistributionTraits<Remote<ETag> >
{
  enum { remote = true };
  typedef DistributedTag LayoutTag_t;
};


#endif     // POOMA_EVALUATOR_ENGINETRAITS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: EngineTraits.h,v $   $Author: richard $
// $Revision: 1.33 $   $Date: 2004/11/01 18:16:40 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
