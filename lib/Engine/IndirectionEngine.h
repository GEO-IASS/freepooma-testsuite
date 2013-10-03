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
// Classes:
// Engine<Dim,T,IndirectionTag<D1,E1,T2,E2> >
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_INDIRECTIONENGINE_H
#define POOMA_ENGINE_INDIRECTIONENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 * Work in progress!
 *
 * Indirection has been implemented in serial and basically works, but there is
 * a significant amount of work that needs to be done to get it to work in
 * parallel.
 */

//-----------------------------------------------------------------------------
// Engine<Dim,T,IndirectionTag<E1,E2> > (aka Indirection-Engine) is the
// engine that contains two engines and performs indirection.
//-----------------------------------------------------------------------------

#include "Engine/Engine.h"
#include "Engine/EngineFunctor.h"
#include "Engine/DataObject.h"
#include "Evaluator/RequestLocks.h"
#include "Pooma/View.h"

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template<int Dim> class DomainLayout;

//-----------------------------------------------------------------------------
// Tag class used to encode the type of an Expression for Engine.
// The exact form of this tag is work in progress.
//-----------------------------------------------------------------------------

template<class A1,class A2>
struct IndirectionTag
{ };


/**
 * The indirection engine combines two arrays to represent ind = a1(a2),
 * where ind(i,j) = a1(a2(i,j)).
 */

template<int Dim,class T,class A1,class A2>
class Engine<Dim,T,IndirectionTag<A1,A2> >
{
public:

  //---------------------------------------------------------------------------
  // Required typedefs:
  //
  // Currently the definition of ElementRef_t is potentially bogus, since it
  // should be based on A1::ElementRef_t.  A1 can be an Array, however,
  // which does not define ElementRef_t.

  typedef IndirectionTag<A1,A2>              Tag_t;
  typedef Engine<Dim,T,Tag_t>                Engine_t;
  typedef typename A1::Element_t             Element_t;
  typedef typename A1::ElementRef_t          ElementRef_t;
  typedef typename A2::Domain_t              Domain_t;
  typedef DomainLayout<Dim>                  Layout_t;
  typedef typename A1::Engine_t              Engine1_t;
  typedef typename A2::Engine_t              Engine2_t;

  enum { dimensions = Dim };
  enum { hasDataObject = Engine1_t::hasDataObject ||
	 Engine2_t::hasDataObject };
  enum { multiPatch = Engine1_t::multiPatch ||
	 Engine2_t::multiPatch };
  enum { dynamic = false };
  enum { zeroBased = Engine2_t::zeroBased };

  //---------------------------------------------------------------------------
  // 

  inline
  Engine(const A1 &array1,const A2 &array2)
    : array1_m(array1),array2_m(array2)
  {
    // array2 takes the inputs, so this engine needs to have array2's Dim.
    CTAssert(A2::dimensions == Dim);
  }

  //---------------------------------------------------------------------------
  // Copy constructor.

  inline
  Engine(const Engine_t &engine)
    : array1_m(engine.array1()),array2_m(engine.array2())
  { }

  //---------------------------------------------------------------------------
  // Subsetting Constructor.
  // 

  template<int OtherDim,class OtherA2, class Domain>
  inline
  Engine(const Engine<OtherDim,T,IndirectionTag<A1,OtherA2> > &e,
	 const Domain &d)
    : array1_m(e.array1()),array2_m(e.array2(),d)
  {
    // array2 takes the inputs, so this engine needs to have array2's Dim.
    CTAssert(A2::dimensions == Dim);
  }

  //---------------------------------------------------------------------------
  // Accessor functions.

  inline const A1 &array1() const { return array1_m; }
  inline A1 &array1() { return array1_m; }
  inline const A2 &array2() const { return array2_m; }
  inline A2 &array2() { return array2_m; }

  //---------------------------------------------------------------------------
  // Get a private copy...

  //  Engine_t &makeOwnCopy();

  //---------------------------------------------------------------------------
  // Accessor functions for a specific element. We recursively go through the
  // expression tree asking each node and leaf to evaluate itself and combine
  // results based on the OpCombine. 

  inline Element_t read(int i0) const 
  {
    return array1_m.read(array2_m.read(i0));
  }

  inline Element_t read(int i0, int i1) const 
  {
    return array1_m.read(array2_m.read(i0,i1));
  }

  inline Element_t read(int i0, int i1,int i2) const 
  {
    return array1_m.read(array2_m.read(i0,i1,i2));
  }

  inline Element_t read(int i0, int i1,int i2,int i3) const 
  {
    return array1_m.read(array2_m.read(i0,i1,i2,i3));
  }

  inline Element_t read(int i0, int i1, int i2, int i3, int i4) const 
  {
    return array1_m.read(array2_m.read(i0,i1,i2,i3,i4));
  }

  inline Element_t read(int i0, int i1, int i2, int i3, int i4, int i5) const 
  {
    return array1_m.read(array2_m.read(i0,i1,i2,i3,i4,i5));
  }

  inline Element_t read(int i0, int i1, int i2, int i3, int i4, int i5,
			int i6) const 
  {
    return array1_m.read(array2_m.read(i0,i1,i2,i3,i4,i5,i6));
  }

  template<class Domain>
  inline Element_t read(const Domain &loc) const 
  {
    return array1_m.read(array2_m.read(loc));
  }

  //---------------------------------------------------------------------------
  // Accessor functions for a specific element. Normally these have read-write
  // semantics. However, the notion of writing to an expression makes no sense
  // so we simply use the read functions.

  inline ElementRef_t operator()(int i0) const 
  {
    return array1_m(array2_m.read(i0));
  }

  inline ElementRef_t operator()(int i0, int i1) const 
  {
    return array1_m(array2_m.read(i0,i1));
  }

  inline ElementRef_t operator()(int i0, int i1,int i2) const 
  {
    return array1_m(array2_m.read(i0,i1,i2));
  }

  inline ElementRef_t operator()(int i0, int i1,int i2,int i3) const 
  {
    return array1_m(array2_m.read(i0,i1,i2,i3));
  }

  inline ElementRef_t operator()(int i0, int i1, int i2, int i3, int i4) const
  {
    return array1_m(array2_m.read(i0,i1,i2,i3,i4));
  }

  inline ElementRef_t operator()(int i0, int i1, int i2, int i3, int i4,
				 int i5) const 
  {
    return array1_m(array2_m.read(i0,i1,i2,i3,i4,i5));
  }

  inline ElementRef_t operator()(int i0, int i1, int i2, int i3, int i4,
				 int i5, int i6) const 
  {
    return array1_m(array2_m.read(i0,i1,i2,i3,i4,i5,i6));
  }

  template<class Domain>
  inline ElementRef_t operator()(const Domain &loc) const 
  {
    return array1_m(array2_m.read(loc));
  }

  //---------------------------------------------------------------------------

  inline const Domain_t& domain() const 
  { 
    return array2_m.domain();
  }

  //---------------------------------------------------------------------------
  // Return the first index value for the specified direction.
  
  inline int first(int i) const
  {
    return array2_m.first(i);
  }
  
private:

  // The contained arrays are stored here.

  A1 array1_m;
  A2 array2_m;

};

//-----------------------------------------------------------------------------
// NewEngine
//
// Here we define the standard subsetting operation for indirection engine.  We
// use View to subset the indirector array.  Note that we don't use
// View1<A1,NewA2_t>, since that could involve the array that contains this
// engine and introduce a circular dependency.  (We're allow to take a view of
// the contained array, since it contains a different engine.)
//-----------------------------------------------------------------------------

template<int Dim,class T,class A1,class A2,class Domain>
struct NewEngine<Engine<Dim,T,IndirectionTag<A1,A2> >,Domain>
{
  typedef typename View1<A2,Domain>::Type_t NewA2_t;
  enum { newDim = NewA2_t::dimensions };
  typedef Engine<newDim,T,IndirectionTag<A1,NewA2_t> > Type_t;
};

//---------------------------------------------------------------------------
// Specialization of  DataObjectRequest engineFunctor to pass the request to
// the contained engine.
//---------------------------------------------------------------------------

template<int Dim, class T, class A1, class A2, class RequestType>
struct EngineFunctor<Engine<Dim, T, IndirectionTag<A1, A2> >,
  DataObjectRequest<RequestType> >
{
  typedef typename DataObjectRequest<RequestType>::Type_t Type_t;
  typedef typename DataObjectRequest<RequestType>::Combine_t Combine_t;

  static Type_t
  apply(const Engine<Dim, T, IndirectionTag<A1, A2> > &engine,
	const DataObjectRequest<RequestType> &tag)
  {
    return Combine2<Type_t,Type_t,OpAdd,
      Combine_t>::combine(
			  engineFunctor(engine.array1().engine(), tag),
			  engineFunctor(engine.array2().engine(), tag),
			  Combine_t()
			  ); 

  }
};

template<int Dim, class T, class A1, class A2>
struct EngineFunctor<Engine<Dim, T, IndirectionTag<A1, A2> >,
  DataObjectRequest<WriteRequest> >
{
  typedef typename DataObjectRequest<WriteRequest>::Type_t Type_t;

  static Type_t
  apply(const Engine<Dim, T, IndirectionTag<A1, A2> > &engine,
	const DataObjectRequest<WriteRequest> &tag)
  {
    engineFunctor(engine.array1().engine(), tag);
    return engineFunctor(engine.array2().engine(),
			 DataObjectRequest<ReadRequest>(tag));
  }
};

template<int Dim, class T, class A1, class A2>
struct EngineFunctor<Engine<Dim, T, IndirectionTag<A1, A2> >,
  DataObjectRequest<WriteRelease> >
{
  typedef typename DataObjectRequest<WriteRelease>::Type_t Type_t;

  static Type_t
  apply(const Engine<Dim, T, IndirectionTag<A1, A2> > &engine,
	const DataObjectRequest<WriteRelease> &tag)
  {
    engineFunctor(engine.array1().engine(), tag);
    return engineFunctor(engine.array2().engine(),
			 DataObjectRequest<ReadRelease>(tag));
  }
};

#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IndirectionEngine.h,v $   $Author: richard $
// $Revision: 1.27 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
