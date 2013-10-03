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
// Brick-Engine non-inline template definitions.
//-----------------------------------------------------------------------------

#include "Engine/BrickEngine.h"

///////////////////////////////////////////////////////////////////////////////
//
// Brick-Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,Brick> constructors:
//
//   Engine(const Domain_t &domain)
//   Engine(const Domain_t &domain, const T &model)
//   Engine(const Node<Domain_t> &node)
//   Engine(const Layout_t &layout)
//
// Constructs a Brick-Engine holding type T elements with the
// multidimensional domain given by either a Domain_t object, a
// Node object, or a Layout_t object. These constructors allocate
// memory for the elements, and initialize the memory with
// ElementProperties::construct. 
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,Brick>::Engine(const Domain_t &dom)
  : Base_t(dom), dataBlock_m(dom.size()), data_m(dataBlock_m.currentPointer())
{ }

template <int Dim, class T>
Engine<Dim,T,Brick>::Engine(const Node<Domain_t> &node)
  : Base_t(node),
    dataBlock_m(node.allocated().size(), node.affinity(), 
           typename DataBlockPtr<T>::WithAffinity_t()),
    data_m(dataBlock_m.currentPointer())
{ }

template <int Dim, class T>
Engine<Dim,T,Brick>::Engine(const Layout_t &layout)
  : Base_t(layout), dataBlock_m(layout.domain().size()),
    data_m(dataBlock_m.currentPointer())
{ }

template <int Dim, class T>
Engine<Dim,T,Brick>::Engine(const Domain_t &dom, const T& model)
  : Base_t(dom), dataBlock_m(dom.size(), model),
    data_m(dataBlock_m.currentPointer())
{ }

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,Brick>(T * foreignData, const Interval<Dim> &domain)
//
// Constructs a Brick-Engine holding type T elements with the
// multidimensional domain given by Interval<Dim>.  The actual data is
// provided in an external buffer pointed to by foreignData.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,Brick>::Engine(T * foreignData, const Domain_t &dom)
  : Base_t(dom), dataBlock_m(foreignData, dom.size()),
    data_m(dataBlock_m.currentPointer())
{ }

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,Brick>(const Engine<Dim,T,Brick> &)
//
// Copy constructor for Brick-Engine.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,Brick>::Engine(const This_t &modelEngine)
  : Base_t(modelEngine), dataBlock_m(modelEngine.dataBlock_m),
    data_m(modelEngine.data_m)
{ }

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,Brick> & operator=(const Engine<Dim,T,Brick> &)
//
// Assignment operator for Brick-Engines.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,Brick> & Engine<Dim,T,Brick>::operator=(const This_t &modelEngine)
{
  // Can skip the rest if we're trying to assign to ourselves

  if (this == &modelEngine) return *this;

  // Copy the base and the data block
  
  Base_t::operator=(modelEngine);
  dataBlock_m = modelEngine.dataBlock_m;
  data_m = modelEngine.data_m;
  PAssert(dataBlock_m.isAtBeginning());
  return *this;
}

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,Brick> & makeOwnCopy()
//
// Causes the Brick-Engine to obtain a private copy of the data
// that it refers to.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,Brick> &Engine<Dim,T,Brick>::makeOwnCopy()
{
  if (dataBlock_m.isValid() && dataBlock_m.count() > 1) 
    {
      PAssert(dataBlock_m.isAtBeginning());
      dataBlock_m.makeOwnCopy();
      data_m = dataBlock_m.currentPointer();
    }

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// BrickView Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Default constructor is required for containers.
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,BrickView>::
Engine()
  : Base_t(), dataBlock_m(), data_m(0)
{ }

//-----------------------------------------------------------------------------
//
// Engine(const Engine<Dim,T,BrickView<Dim2> > &);
//
// Copy constructor.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,BrickView>::
Engine(const This_t &modelEngine)
  : Base_t(modelEngine), dataBlock_m(modelEngine.dataBlock_m), data_m(dataBlock_m.currentPointer())
{ }

// What is this for again???

template <int Dim, class T>
Engine<Dim,T,BrickView>::
Engine(const This_t &modelEngine, const EngineConstructTag &)
  : Base_t(modelEngine), dataBlock_m(modelEngine.dataBlock_m), data_m(dataBlock_m.currentPointer())
{ }

//-----------------------------------------------------------------------------
//
// Engine<Dim,T,BrickView> & operator=(const Engine<Dim,T,BrickView> &)
//
// Assignment operator.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,BrickView> &
Engine<Dim,T,BrickView>::operator=(const This_t &modelEngine)
{
  if (this == &modelEngine) return *this;
  Base_t::operator=(modelEngine);
  dataBlock_m = modelEngine.dataBlock_m;
  data_m = modelEngine.data_m;
  return *this;
}

//-----------------------------------------------------------------------------
//
// Engine(const Engine<Dim,T,BrickView<Dim2> > &);
//
// Construct from compressible brick.
//
//-----------------------------------------------------------------------------

// Change this when CompressibleBrick is modified to use BaseDomain!

template <int Dim, class T>
Engine<Dim,T,BrickView>::
Engine(const Engine<Dim,T,CompressibleBrick> &model)
  : Base_t(model, false)
{
  dataBlock_m = DataBlockPtr<T>(model.dataBlock(),this->baseOffset());
  data_m = dataBlock_m.currentPointer();
}

//-----------------------------------------------------------------------------
//
// Engine(const Engine<Dim,T,BrickView<Dim2> > &);
//
// Construct from compressible brick view.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,BrickView>::
Engine(const Engine<Dim,T,CompressibleBrickView> &model)
  : Base_t(model, false)
{
  dataBlock_m = DataBlockPtr<T>(model.dataBlock(),this->baseOffset());
  data_m = dataBlock_m.currentPointer();
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BrickEngine.cpp,v $   $Author: richi $
// $Revision: 1.80 $   $Date: 2004/11/10 21:57:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
