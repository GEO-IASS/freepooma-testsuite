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
//   IndexFunction         - Tag class for defining an index-function-engine.
//   IndexFunctionView     - Tag class for defining an view of an
//                           index-function-engine.
//   Engine                - Specialization for IndexFunction
//   NewEngine             - Specializations for IndexFunction
//-----------------------------------------------------------------------------

#ifndef POOMA_ENGINE_INDEXFUNCTIONENGINE_H
#define POOMA_ENGINE_INDEXFUNCTIONENGINE_H

/** @file
 * @ingroup Engine
 * @brief
 * Index-function-engine objects provide a way to make a function of indices
 * work like an array.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/SliceDomain.h"
#include "Engine/Engine.h"
#include "Engine/ViewEngine.h"
#include "Layout/INode.h"
#include "Layout/DomainLayout.h"
#include "PETE/PETE.h"
#include "Utilities/ViewIndexer.h"


/**
 * IndexFunction is just a tag class for the index-function-engine engine,
 * which makes a function of indices look like an array. It takes a Functor
 * class type as a template argument. This functor is what turns indices
 * into function values.
 */

template<class Functor>
struct IndexFunction 
{ };

/**
 * IndexFunctionView is the view analog of IndexFunction.
 * In addition to the function, this class includes the original dimension.
 */

template<int Dim, class Functor>
struct IndexFunctionView
{ };

/**
 * Engine<Dim, T, IndexFunction<Functor> > is a specialization of Engine for 
 * IndexFunction.
 *
 * This does all of the usual Engine things: 
 * - Typedefs for the tag, element types, domain and dimensions.
 * - Operator() with integers to evaluate elements quickly.
 * - Operator() with a domain to subset.
 * - Accessor for the domain.
 */

template<int Dim, class T, class Functor>
class Engine<Dim, T, IndexFunction<Functor> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants

  typedef IndexFunction<Functor>                   Tag_t;
  typedef Engine<Dim, T, Tag_t>                    This_t;
  typedef This_t                                   Engine_t;
  typedef Interval<Dim>                            Domain_t;
  typedef DomainLayout<Dim>                        Layout_t;
  typedef T                                        Element_t;
  typedef ErrorType                                ElementRef_t;

  enum { dimensions = Dim };
  enum { hasDataObject = false };
  enum { dynamic = false };
  enum { zeroBased = false };
  enum { multiPatch = false };

  //---------------------------------------------------------------------------
  /// Default constructor (allows subsequent initialization of domain/functor).
  
  Engine() { }

  //---------------------------------------------------------------------------
  /// Construct from a domain/layout object and an optional Functor object.

  explicit Engine(const Domain_t &domain, const Functor &f = Functor())
  : funct_m(f), domain_m(domain)
  { 
  }

  template<class Layout>
  explicit Engine(const Layout &layout, const Functor &f = Functor())
  : funct_m(f), domain_m(layout.domain())
  { 
  }

  //---------------------------------------------------------------------------
  /// Construct from another index-function-engine.

  Engine(const This_t &model)
  : funct_m(model.functor()), domain_m(model.domain())
  { 
  }

  //---------------------------------------------------------------------------
  /// Assign one index-function-engine to another.

  This_t &operator=(const This_t &rhs)
  {
    domain_m = rhs.domain();
    funct_m = rhs.functor();
    
    return *this;
  }
      
  //---------------------------------------------------------------------------
  /// @name Element access via ints for speed.
  //@{

  inline Element_t read(int i0) const 
    {
      return funct_m(i0);
    }
  inline Element_t read(int i0, int i1) const 
    {
      return funct_m(i0, i1);
    }
  inline Element_t read(int i0, int i1, int i2) const 
    {
      return funct_m(i0, i1, i2);
    }
  inline Element_t read(int i0, int i1, int i2, int i3) const 
    {
      return funct_m(i0, i1, i2, i3);
    }
  inline Element_t read(int i0, int i1, int i2, int i3, int i4) const 
    {
      return funct_m(i0, i1, i2, i3, i4);
    }
  inline Element_t read(int i0, int i1, int i2, int i3, int i4, 
    int i5) const 
    {
      return funct_m(i0, i1, i2, i3, i4, i5);
    }
  inline Element_t read(int i0, int i1, int i2, int i3, int i4, 
    int i5, int i6) const 
    {
      return funct_m(i0, i1, i2, i3, i4, i5, i6);
    }
  inline Element_t read(const Loc<1> &loc) const
    {
      return funct_m(loc[0].first());
    }
  inline Element_t read(const Loc<2> &loc) const
    {
      return funct_m(loc[0].first(), loc[1].first());
    }
  inline Element_t read(const Loc<3> &loc) const
    {
      return funct_m(loc[0].first(), loc[1].first(), loc[2].first());
    }
  inline Element_t read(const Loc<4> &loc) const
    {
      return funct_m(loc[0].first(), loc[1].first(), loc[2].first(), 
        loc[3].first());
    }
  inline Element_t read(const Loc<5> &loc) const
    {
      return funct_m(loc[0].first(), loc[1].first(), loc[2].first(), 
        loc[3].first(), loc[4].first());
    }
  inline Element_t read(const Loc<6> &loc) const
    {
      return funct_m(loc[0].first(), loc[1].first(), loc[2].first(), 
        loc[3].first(), loc[4].first(), loc[5].first());
    }
  inline Element_t read(const Loc<7> &loc) const
    {
      return funct_m(loc[0].first(), loc[1].first(), loc[2].first(), 
        loc[3].first(), loc[4].first(), loc[5].first(), loc[6].first());
    }

  //@}

  //---------------------------------------------------------------------------
  /// Return/set the domain. Also, return the base domain.

  inline const Domain_t &domain() const { return domain_m; }
  void setDomain(const Domain_t &dom) { domain_m = dom; }

  //---------------------------------------------------------------------------
  /// Return the first index value for the specified direction.
  
  inline int first(int i) const
  {
    PAssert(i >= 0 && i < Dim);
    return domain_m[i].first();
  }

  //---------------------------------------------------------------------------
  /// Returns the layout, which is constructed as a DomainLayout.

  Layout_t layout() const
  {
    return Layout_t(domain_m);
  }

  //---------------------------------------------------------------------------
  /// Accessor/modifier.

  const Functor &functor() const { return funct_m; }
  void setFunctor(const Functor &f) { funct_m = f; }

private:

  Functor funct_m;
  Domain_t domain_m;
};


/**
 * NewEngine<Engine,SubDomain>
 *
 * Specializations of NewEngine for subsetting a index-function-engines with
 * an arbitrary domain. 
 */

template <int Dim, class T, class Functor>
struct NewEngine<Engine<Dim, T, IndexFunction<Functor> >, Interval<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<Dim, IndexFunction<Functor> > > Type_t;
};

template <int Dim, class T, class Functor>
struct NewEngine<Engine<Dim, T, IndexFunction<Functor> >, Range<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<Dim, IndexFunction<Functor> > > Type_t;
};

template <int Dim, class T, class Functor, int SliceDim>
struct NewEngine<Engine<Dim, T, IndexFunction<Functor> >, 
  SliceInterval<Dim, SliceDim> >
{
  typedef Engine<SliceDim, T, ViewEngine<Dim, IndexFunction<Functor> > > Type_t;
};

template <int Dim, class T, class Functor, int SliceDim>
struct NewEngine<Engine<Dim, T, IndexFunction<Functor> >, 
  SliceRange<Dim, SliceDim> >
{
  typedef Engine<SliceDim, T, ViewEngine<Dim, IndexFunction<Functor> > > Type_t;
};

template <int Dim, class T, class Functor, class Domain>
struct NewEngine<Engine<Dim, T, IndexFunction<Functor> >, Node<Domain> >
{
  typedef Engine<Dim, T, ViewEngine<Dim, IndexFunction<Functor> > > Type_t;
};

template <int Dim, class T, class Functor>
struct NewEngine<Engine<Dim, T, IndexFunction<Functor> >, INode<Dim> >
{
  typedef Engine<Dim, T, ViewEngine<Dim, IndexFunction<Functor> > > Type_t;
};

#endif // POOMA_ENGINE_INDEXFUNCTIONENGINE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: IndexFunctionEngine.h,v $   $Author: richard $
// $Revision: 1.28 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
