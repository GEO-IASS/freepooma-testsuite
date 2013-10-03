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
// CompressibleBrick-Engine non-inline template definitions.
//-----------------------------------------------------------------------------

#include "Engine/CompressibleBrick.h"

///////////////////////////////////////////////////////////////////////////////
// namespace Pooma {


///////////////////////////////////////////////////////////////////////////////
//
// CompressibleBrick-Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//
// Engine<Dim,T,CompressibleBrick> constructors:
//
//   Engine(const Domain_t &domain)
//   Engine(const Domain_t &domain, const T &model)
//   Engine(const Node<Domain_t> &node)
//   Engine(const Layout_t &layout)
//
// Constructs a CompressibleBrick-Engine representing a Dim-dimensional
// brick (Fortran storage order) of elements of type T. The domain can
// be specified directly or by passing a Node or Layout_t object. If
// initializing with a Domain_t, one can optionally pass in a model 
// element, which will be used to initialize storage. CompressibleBricks
// are always born compressed, storing only a single value, which is 
// initialized with ElementProperties::construct.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick>::Engine(const Domain_t &domain)
  : Base_t(domain), cblock_m(domain.size())
{
  init();
}

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick>::Engine(const Node<Domain_t> &node)
  : Base_t(node.allocated()), 
    cblock_m(node.allocated().size(), node.affinity())
{
  init();
}

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick>::Engine(const Layout_t &layout)
  : Base_t(layout.domain()), cblock_m(layout.domain().size())
{
  init();
}

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick>::Engine(const Domain_t &domain, const T& model)
  : Base_t(domain), cblock_m(domain.size(),-1,model)
{
  init();
}

//-----------------------------------------------------------------------------
//
// Engine<D,T,CompressibleBrick> copy constructor
//
// We don't invoke the BrickBase copy constructor in the initializer
// list as we need to have the cblock locked before copying the 
// base class information. Instead, we invoke the BrickBase assignment
// operator after locking. 
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick>::
Engine(const Engine<Dim,T,CompressibleBrick> &modelEngine)
  : cblock_m(modelEngine.cblock_m)
{
  // Lock the controller so that compression can't occur
  // while we're copying.
  
  cblock_m.lock();
  
  data0_m = modelEngine.data0_m; // Must copy while locked.
  
  Base_t::operator=(modelEngine);
  
  if (cblock_m.isControllerValidUnlocked()) cblock_m.attach(this);
  
  cblock_m.unlock(); 
}


//-----------------------------------------------------------------------------
//
// Engine<D,T,CompressibleBrick> & Engine<D,T,CompressibleBrick>::
//       operator=(const Engine<D,T,CompressibleBrick> &)
//
// Assignment operator for CompressibleBrick-Engines.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick> & 
Engine<Dim,T,CompressibleBrick>::
operator=(const Engine<Dim,T,CompressibleBrick> &modelEngine)
{
  if (this != &modelEngine) 
    {
      // This only works if the RHS has a valid controller pointer.
      // (Perhaps this should be legal, but I don't want to figure that
      // out right now.)
      
      PAssert(modelEngine.cblock_m.isControllerPtrValid());
      
      // Lock the new cblock until we're done copying.
      
      modelEngine.cblock_m.lock();
      
      // Lock the old one and disable notification.
      // (Only do this if the RefCountedPtr actually points to something.)
      
      if (cblock_m.isControllerPtrValid())
        {
          cblock_m.lock();
          if (cblock_m.isControllerValidUnlocked())
            {
              cblock_m.detach(this);
            }
          cblock_m.unlock();
        }
      
      // This just copies the RCPtr<CBC> so it can be done while locked.
      
      cblock_m = modelEngine.cblock_m;
      
      // Lock our own mutex to ensure that no on else tries to copy
      // or use these strides/data while this update is occuring.
      // (It is important to lock the CBC first as that is the order 
      // when notify is called())
       
      lock();
      
      data0_m = modelEngine.data0_m;

      Base_t::operator=(modelEngine);
            
      unlock();
      
      if (cblock_m.isControllerValidUnlocked()) cblock_m.attach(this);
      
      // Unlock our cblock (which is also modelEngine's cblock):
        
      cblock_m.unlock();
    }
  return *this;
}

template <int Dim, class T> 
void Engine<Dim,T,CompressibleBrick>::init()
{
  // resetDataAndStrides gets compression dependent data from
  // the CBC. 
  // This is only called by CompressibleBrick constructors, so there
  // can't be any other viewers, and thus we don't need to lock
  // the CBC.
  
  resetDataAndStrides();
  PAssert(cblock_m.isControllerValidUnlocked());
  cblock_m.attach(this);
}

//-----------------------------------------------------------------------------
//
// Out-of-line destructor to shorten compile times.
// (It's ~2% effect on compile time, so re-inline it
// if the performance cost is big.)
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick>::~Engine()
{
  if (data0_m)
    {
      cblock_m.lock();
      if (cblock_m.isControllerValidUnlocked())
        {
          cblock_m.detach(this);
        }
      cblock_m.unlock();
    }
}


//-----------------------------------------------------------------------------
//
// Engine<Dim,T,CompressibleBrick> & makeOwnCopy()
//
// Causes the CompressibleBrick-Engine to obtain a private copy of the data
// that it refers to.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrick> &Engine<Dim,T,CompressibleBrick>::makeOwnCopy()
{
  // JIM: This is probably not thread safe??? 
  // There is a race from checking isShared to getting into cblock's
  // makeOwnCopy, which is thread safe. As a result, this should only
  // be called after a blockAndEvaluate() to ensure that nobody else
  // is messing with the underlying CBC while this is
  // occuring. (Logically, this is necessary anyway since you probably
  // want a copy of the data that results from all previous
  // computations having taken place.)  Also, as mentioned elsewhere,
  // the current implementation results in copying uncompressed data
  // in the parse thread, which will result in incorrect memory
  // affinity.
  
  if (cblock_m.isControllerValidUnlocked() && cblock_m.isShared()) 
    {
      cblock_m.detach(this);
      cblock_m.makeOwnCopy();
      cblock_m.attach(this);

      data0_m = cblock_m.data() + (cblock_m.compressed() ? 0 : this->baseOffset());
    }

  return *this;
}


//-----------------------------------------------------------------------------
//
// Notify function for the Observer<T*> base class:
// Compressible bricks observe the CompressibleBlock
// which notifies us when the data becomes compressed or uncompressed.
// The notification comes with a pointer to the new data.
//  
// Note: The CBC is locked when this is called, so we don't have to
// worry about contention for changing strides_m/data0_m, but we
// do need to make sure no one tries to make a copy of these data
// while they are being changed. Thus we lock our mutex.
//
// Also note: A corallary is that if you're going to lock both
// the engine and the CBC, always lock the CBC first. Otherwise
// there is a potential deadlock!
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
void Engine<Dim,T,CompressibleBrick>::
notify(T* &data, const ObserverEvent &event)
{
  switch (event.event())
    {
    default:
    case CompressibleBlock<T>::notifyDestruct:
      // cblock has destructed. this should never happen
      // if the engine still exists.
      PAssert(false);
      break;
    case CompressibleBlock<T>::notifyUncompress: 
      lock();
      this->restoreStrides();
      data0_m = data + this->baseOffset();
      unlock();
      break;
    case CompressibleBlock<T>::notifyCompress:
      lock();
      this->zeroStrides();
      data0_m = data;
      unlock();
      break;
    }
}

//-----------------------------------------------------------------------------
//
// resetDataAndStrides is a utility function used by various
// constructors.  It sets the strides based on the compression
// status, and sets the pointer to the data.  Once the Compressible
// Brick has been created, the strides and data pointer are updated
// by the notify function.
//
// NOTE: the cblock must be locked before this function is called.
//-----------------------------------------------------------------------------

template <int Dim, class T>
void Engine<Dim,T,CompressibleBrick>::resetDataAndStrides()
{
  if (cblock_m.compressed())
    {
      this->zeroStrides();
      data0_m = cblock_m.data();
    }
  else
    {
      this->restoreStrides();
      data0_m = cblock_m.data() + this->baseOffset();
    }
}

//-----------------------------------------------------------------------------
//
// Engine<Dim, T, CompressibleBrick>::
// elementsCompressed() const
//
// Return the number of compressed elements.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
long Engine<Dim,T,CompressibleBrick>::
elementsCompressed() const
{
  // If we are compressed, the number of compressed elements is the size
  // of our domain; otherwise, it is zero.
   
  if (compressed())
    return domain().size();
  else
    return 0L;
}

///////////////////////////////////////////////////////////////////////////////
//
// CompressibleBrickView-Engine Member Functions
//
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//
// Out-of-line destructor to shorten compile times.
// (It's ~2% effect on compile time, so re-inline it
// if the performance cost is big.)
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrickView>::
~Engine()
{
  cblock_m.lock();
  if (cblock_m.isControllerValidUnlocked())
    {
      cblock_m.detach(this);
    }
  cblock_m.unlock();
}

//-----------------------------------------------------------------------------
//
// Engine<D,T,CompressibleBrickView> & Engine<D,T,CompressibleBrickView>::
//       operator=(const Engine<D,T,CompressibleBrickView> &)
//
// Assignment operator for CompressibleBrickView-Engines.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrickView> & 
Engine<Dim,T,CompressibleBrickView>::
operator=(const Engine<Dim,T,CompressibleBrickView> &modelEngine)
{
  if (this != &modelEngine) 
    {
      // This only works if the RHS has a valid controller pointer.
      // (Perhaps this should be legal, but I don't want to figure that
      // out right now.)
      
      PAssert(modelEngine.cblock_m.isControllerPtrValid());
      
      // Lock the new cblock until we're done copying.
      
      modelEngine.cblock_m.lock();
      
      // Lock the old one and disable notification.
      // (Only do this if the RefCountedPtr actually points to something.)
      
      if (cblock_m.isControllerPtrValid())
        {
          cblock_m.lock();
          if (cblock_m.isControllerValidUnlocked())
            {
              cblock_m.detach(this);
            }
          cblock_m.unlock();
        }
      
      // This just copies the RCPtr<CBC> so it can be done while locked.
      
      cblock_m = modelEngine.cblock_m;
      entire_m = modelEngine.entire_m;
     
        
      // Lock our own mutex to ensure that no on else tries to copy
      // or use these strides/data while this update is occuring.
      // (It is important to lock the CBC first as that is the order 
      // when notify is called())
       
      lock();
      
      data0_m  = modelEngine.data0_m;
      Base_t::operator=(modelEngine);
      
      unlock();
      
      if (cblock_m.isControllerValidUnlocked()) cblock_m.attach(this);
      
      // Unlock our cblock (which is also modelEngine's cblock):
        
      cblock_m.unlock();
    }
  return *this;
}


//-----------------------------------------------------------------------------
//
// Engine<D,T,CompressibleBrickView>::
// Engine(const Engine<D,T,CompressibleBrickView> &)
//
// Copy constructor for CompressibleBrickView-Engine.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
Engine<Dim,T,CompressibleBrickView>::
Engine(const Engine<Dim,T,CompressibleBrickView> &modelEngine)
  : cblock_m(modelEngine.cblock_m),
    entire_m(modelEngine.entire_m)
{
  // Lock the controller so the RHS's compression state doesn't change.
  
  cblock_m.lock();
    
  // This being a constructor, no-one else can try to use our strides
  // and data0 until we're done, so locking our mutex is unnecessary.
  
  data0_m = modelEngine.data0_m;
  Base_t::operator=(modelEngine);
    
  if (cblock_m.isControllerValidUnlocked()) cblock_m.attach(this);
  
  cblock_m.unlock();  
}

template <int Dim, class T>
Engine<Dim,T,CompressibleBrickView>::
Engine(const Engine<Dim,T,CompressibleBrickView> &modelEngine,
       const EngineConstructTag &)
  : cblock_m(modelEngine.cblock_m),
    entire_m(modelEngine.entire_m)
{
  // Lock the controller so the RHS's compression state doesn't change.
  
  cblock_m.lock();
    
  // This being a constructor, no-one else can try to use our strides
  // and data0 until we're done, so locking our mutex is unnecessary.
  
  data0_m = modelEngine.data0_m;
  Base_t::operator=(modelEngine);
    
  if (cblock_m.isControllerValidUnlocked()) cblock_m.attach(this);
  
  cblock_m.unlock();  
}

//-----------------------------------------------------------------------------
//
// Engine<Dim, T, CompressibleBrickView>::
// elementsCompressed() const
//
// Return the number of compressed elements.
//
//-----------------------------------------------------------------------------

template <int Dim, class T>
long Engine<Dim,T,CompressibleBrickView>::
elementsCompressed() const
{
  // If we are compressed, the number of compressed elements is the size
  // of our domain; otherwise, it is zero.
   
  if (compressed())
    return domain().size();
  else
    return 0L;
}

// } // namespace Pooma
///////////////////////////////////////////////////////////////////////////////

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CompressibleBrick.cpp,v $   $Author: richard $
// $Revision: 1.28 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
