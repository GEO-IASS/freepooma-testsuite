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
//   ConstantFunction      - Tag class for defining a constant-function-engine.
//   Engine                - Specialization for ConstantFunction
//   NewEngine             - Specializations for ConstantFunction
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_CONSTANTFUNCTIONENGINE_H
#define POOMA_ENGINE_CONSTANTFUNCTIONENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 * Constant-function-engine objects provide a way to make a scalar behave like
 * an array.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/SliceDomain.h"
#include "Engine/Engine.h"
#include "Layout/INode.h"
#include "Layout/DomainLayout.h"
#include "PETE/ErrorType.h"


/**
 * ConstantFunction is just a tag class for the constant-function engine,
 * which makes a scalar look like an Array.
 */

struct ConstantFunction 
{ };

/**
 * Engine<Dim, T, ConstantFunction> is a specialization of Engine for 
 * ConstantFunction.
 *
 * This does all of the usual Engine things:
 *  - Typedefs for the tag, element types, domain and dimensions.
 *  - Operator() with integers to evaluate elements quickly.
 *  - Operator() with a domain to subset.
 *  - Accessor for the domain.
 */

template<int Dim, class T>
class Engine<Dim, T, ConstantFunction>
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants

  typedef ConstantFunction                   Tag_t;
  typedef Engine<Dim, T, ConstantFunction>   This_t;
  typedef This_t                             Engine_t;
  typedef Interval<Dim>                      Domain_t;
  typedef DomainLayout<Dim>                  Layout_t;
  typedef T                                  Element_t;
  typedef ErrorType                          ElementRef_t;

  enum { dimensions = Dim };
  enum { hasDataObject = false };
  enum { dynamic = false };
  enum { zeroBased = false };
  enum { multiPatch = false };

  //---------------------------------------------------------------------------
  // Default constructor.

  Engine() { }
  
  //---------------------------------------------------------------------------
  // Construct from a domain object.

  explicit Engine(const Domain_t &domain, T val = T())
  : val_m(val), domain_m(domain)
  { 
    for (int d = 0; d < Dim; ++d)
      firsts_m[d] = domain[d].first();
  }

  template<class Layout>
  explicit Engine(const Layout &layout, T val = T())
  : val_m(val), domain_m(layout.domain())
  { 
    for (int d = 0; d < Dim; ++d)
      firsts_m[d] = domain_m[d].first();
  }

  //---------------------------------------------------------------------------
  // Copy constructor.

  Engine(const Engine<Dim, T, ConstantFunction> &model)
  : val_m(model.constant()), domain_m(model.domain())
  {
    for (int d = 0; d < Dim; ++d)
      {
        firsts_m[d] = model.firsts_m[d];
      }
  }

  //---------------------------------------------------------------------------
  // Construct from various sorts of domains (e.g., take a view).

  template<class DT>
  Engine(const Engine<Dim, T, ConstantFunction> &e, const Domain<Dim, DT> &dom)
  : val_m(e.constant()), domain_m(Pooma::NoInit())
  {
    const typename DT::Domain_t &domain = dom.unwrap();
    for (int d = 0; d < Dim; ++d)
      {
        domain_m[d]  = Interval<1>(domain[d].length());
        firsts_m[d] = 0;
      }
  }

  template<int Dim2, class DT>
  Engine(const Engine<Dim2, T, ConstantFunction> &e, 
    const SliceDomain<DT> &dom)
  : val_m(e.constant()), domain_m(Pooma::NoInit())
  {
    // The domain's dimension should match ours.
    
    CTAssert(DT::sliceDimensions == Dim);
    CTAssert(DT::dimensions == Dim2);

    const typename DT::SliceDomain_t &domain = dom.sliceDomain();
    for (int d = 0; d < Dim; ++d)
      {
        domain_m[d]  = Interval<1>(domain[d].length());
        firsts_m[d] = 0;
      }
  }

  template<class Domain>
  Engine(const Engine<Dim, T, ConstantFunction> &e, const Node<Domain> &node)
  : val_m(e.constant()), domain_m(Pooma::NoInit())
  {
    // The nodes's dimension should match ours.
    
    CTAssert(Domain::dimensions == Dim);

    const Domain &domain = node.domain();
    for (int d = 0; d < Dim; ++d)
      {
        domain_m[d]  = Interval<1>(domain[d].length());
        firsts_m[d] = 0;
      }
  }

  Engine(const Engine<Dim, T, ConstantFunction> &e, const INode<Dim> &inode)
  : val_m(e.constant()), domain_m(Pooma::NoInit())
  {
    const typename INode<Dim>::Domain_t &domain = inode.domain();
    for (int d = 0; d < Dim; ++d)
      {
        domain_m[d]  = Interval<1>(domain[d].length());
        firsts_m[d] = 0;
      }
  }
    
  //---------------------------------------------------------------------------
  // Element access via ints for speed.
  // We only need read() functions since this engine should only be used in
  // a read-only array.

  inline Element_t read(int) const 
    {
      return val_m;
    }
  inline Element_t read(int, int) const 
    {
      return val_m;
    }
  inline Element_t read(int, int, int) const 
    {
      return val_m;
    }
  inline Element_t read(int, int, int, int) const 
    {
      return val_m;
    }
  inline Element_t read(int, int, int, int, int) const 
    {
      return val_m;
    }
  inline Element_t read(int, int, int, int, int, int) const 
    {
      return val_m;
    }
  inline Element_t read(int, int, int, int, int, int, int) const 
    {
      return val_m;
    }
  inline Element_t read(const Loc<Dim> &) const
    {
      return val_m;
    }

  //---------------------------------------------------------------------------
  // Return the domain.

  const Domain_t &domain() const { return domain_m; }

  //---------------------------------------------------------------------------
  // Return a layout.

  inline Layout_t layout() const { return Layout_t(domain_m); }

  // Return the first value for the specified direction.
  
  inline int first(int i) const
  {
    PAssert(i >= 0 && i < Dim);
    return firsts_m[i];
  }

  //---------------------------------------------------------------------------
  // Accessors/modifiers.

  T constant() const { return val_m; }
  void setConstant(T val) { val_m = val; }

private:

  T val_m;
  Domain_t domain_m;
  int firsts_m[Dim];
};

/**
 * NewEngine<Engine,SubDomain>
 *
 * Specializations of NewEngine for subsetting a constant-function-engine with
 * an arbitrary domain. 
 */

template <int Dim, class T>
struct NewEngine<Engine<Dim, T, ConstantFunction>, Interval<Dim> >
{
  typedef Engine<Dim, T, ConstantFunction> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim, T, ConstantFunction>, Range<Dim> >
{
  typedef Engine<Dim, T, ConstantFunction> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,ConstantFunction>, SliceInterval<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,ConstantFunction> Type_t;
};

template <int Dim, class T, int SliceDim>
struct NewEngine<Engine<Dim,T,ConstantFunction>, SliceRange<Dim,SliceDim> >
{
  typedef Engine<SliceDim,T,ConstantFunction> Type_t;
};

template <int Dim, class T, class Domain>
struct NewEngine<Engine<Dim, T, ConstantFunction>, Node<Domain> >
{
  typedef Engine<Dim, T, ConstantFunction> Type_t;
};

template <int Dim, class T>
struct NewEngine<Engine<Dim, T, ConstantFunction>, INode<Dim> >
{
  typedef Engine<Dim, T, ConstantFunction> Type_t;
};

#endif // POOMA_ENGINE_CONSTANTFUNCTIONENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ConstantFunctionEngine.h,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
