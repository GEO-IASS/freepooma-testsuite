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
//   DynamicArray<T, EngineTag>
//   ElementProperties< DynamicArray<T, EngineTag> > 
//   CreateLeaf< DynamicArray<T, EngineTag> >
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Array
 * @brief
 * Dynamic arrays.
 */

#ifndef POOMA_DYNAMIC_ARRAY_DYNAMIC_ARRAY_H
#define POOMA_DYNAMIC_ARRAY_DYNAMIC_ARRAY_H

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

template<class T, class EngineTag> class DynamicArray;

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

// The following hack was put in the POOMA 2.3 to avoid a template
// instantiation error in KCC 3.4d.  The problem has been reported to KAI,
// but until it is fixed, the following odd ordering of includes fixes the
// problem.  The correct set of includes should be:
//
// include "Array/Array.h"
// include "DynamicArray/DynamicArrayOperators.h"
// include "DynamicArray/PoomaDynamicArrayOperators.h"
// include "DynamicArray/VectorDynamicArrayOperators.h"
// include "Layout/DynamicEvents.h"
// include "Engine/EnginePatch.h"

template<int Dim, class T, class EngineTag> class Array;

#include "PETE/PETE.h"
#include "Pooma/PETEExtras.h"

#include "DynamicArray/DynamicArrayOperators.h"
#include "DynamicArray/PoomaDynamicArrayOperators.h"
#include "DynamicArray/VectorDynamicArrayOperators.h"

#include "Array/Array.h"

#include "Layout/DynamicEvents.h"
#include "Engine/EnginePatch.h"

#include "Domain/IteratorPairDomain.h"

//-----------------------------------------------------------------------------
// Prototype for the assign function used to assign a DynamicArray to an Array.
//
// Prototypes defined here:
//   Array = DynamicArray
//-----------------------------------------------------------------------------

template<int Dim, class T, class EngineTag,
         class OtherT, class OtherEngineTag, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs,
       const DynamicArray<OtherT, OtherEngineTag> &rhs,
       const Op &op);


//-----------------------------------------------------------------------------
// View specializations for DynamicArray.
//-----------------------------------------------------------------------------

/// General versions.  Just defer to Array<1>.

template<class T, class EngineTag, class Sub1>
struct View1<DynamicArray<T, EngineTag>, Sub1>
{
  // Convenience typedef for the thing we're taking a view of.
  
  typedef DynamicArray<T, EngineTag> Subject_t;
  typedef Array<1, T, EngineTag> Base_t;

  // Deduce domains for the output type.
  
  typedef typename Subject_t::Domain_t Domain_t;
  typedef TemporaryNewDomain1<Domain_t, Sub1> NewDomain_t;
  typedef typename NewDomain_t::SliceType_t SDomain_t;
  
  // Deduce appropriate version of implementation to dispatch to.
  
  enum { sv = DomainTraits<SDomain_t>::singleValued };
  typedef View1Implementation<Base_t, SDomain_t, sv> Dispatch_t;

  // The optimized domain combiner:
  
  typedef CombineDomainOpt<NewDomain_t, sv> Combine_t;
  
  // The return type.
  
  typedef typename Dispatch_t::Type_t Type_t;

  // A function that creates the view.
  
  inline static
  Type_t make(const Subject_t &a, const Sub1 &s1)
  {
    SDomain_t s(Combine_t::make(a, s1));
    return Dispatch_t::make(a, s);
  }
};

template<class T1, class E1, class T2, class E2>
struct View1<DynamicArray<T1, E1>, Array<1, T2, E2> >
{
  typedef DynamicArray<T1, E1> Array1_t;
  typedef Array<1, T1, E1> BaseArray1_t;
  typedef Array<1, T2, E2> Array2_t;

  typedef typename View1<BaseArray1_t, Array2_t>::Type_t Type_t;

  inline static
  Type_t make(const Array1_t &a, const Array2_t &s)
  {
    return View1<BaseArray1_t, Array2_t>::make(a, s);
  }
};

template<class T1, class E1, class T2, class E2>
struct View1<DynamicArray<T1, E1>, DynamicArray<T2, E2> >
{
  typedef DynamicArray<T1, E1> Array1_t;
  typedef Array<1, T1, E1> BaseArray1_t;
  typedef DynamicArray<T2, E2> Array2_t;
  typedef Array<1, T2, E2> BaseArray2_t;

  typedef typename View1<BaseArray1_t, BaseArray2_t>::Type_t Type_t;

  inline static
  Type_t make(const Array1_t &a, const Array2_t &s)
  {
    return View1<BaseArray1_t, BaseArray2_t>::make(a, s);
  }
};

template<int D, class T1, class E1, class T2, class E2>
struct View1<Array<D, T1, E1>, DynamicArray<T2, E2> >
{
  typedef Array<D, T1, E1> Array1_t;
  typedef DynamicArray<T2, E2> Array2_t;
  typedef Array<1, T2, E2> BaseArray2_t;

  typedef typename View1<Array1_t, BaseArray2_t>::Type_t Type_t;

  inline static
  Type_t make(const Array1_t &a, const Array2_t &s)
  {
    return View1<Array1_t, BaseArray2_t>::make(a, s);
  }
};

//-----------------------------------------------------------------------------
// Patch specialization for DynamicArray.
//-----------------------------------------------------------------------------

template<class T, class EngineTag>
struct Patch<DynamicArray<T, EngineTag> >
{
  typedef DynamicArray<T, EngineTag> Subject_t;
  typedef Array<1, T, EngineTag> Base_t;
  typedef typename Patch<Base_t>::Type_t Type_t;

  inline static
  Type_t make(const Subject_t &subject, int i)
  {
    return Patch<Base_t>::make(subject, i);
  }
};

//-----------------------------------------------------------------------------
// ComponentView specialization for DynamicArray.
//-----------------------------------------------------------------------------

template <class Components, class T, class EngineTag>
struct ComponentView< Components, DynamicArray<T, EngineTag> >
{
  // Convenience typedef for the thing we're taking a component view of.
  
  typedef DynamicArray<T, EngineTag> Subject_t;
  typedef Array<1, T, EngineTag> Base_t;

  // The output type.
  
  typedef typename ComponentView<Components, Base_t>::Type_t Type_t;

  // A function that creates the view.
  
  inline static
  Type_t make(const Subject_t &a, const Components &loc)
  {
    return ComponentView<Components, Base_t>::make(a, loc);
  }
};

//-----------------------------------------------------------------------------
// DynamicArray
//-----------------------------------------------------------------------------

/**
 * A DynamicArray is a read-write array with extra create/destroy methods.
 * It can act just like a regular Array, but can have a dynamically-changing
 * domain.  Create and destroy methods will preserve the values of elements
 * that remain after these operations.  It is by definition 1-Dimensional,
 * and so it does not have a Dim parameter.  It provides the following extra
 * interface, beyond that of the standard Array class:
 *   -  void create(int num)
 *   -  void create(int num, PatchID_t patch)
 *   -  void destroy(const Domain &killlist, const DeleteMethod &method)
 *   -  void destroy(const Domain &killlist, PatchID_t patch,
 *                   const DeleteMethod &method)
 *   -  void copy(const Domain &dom, PatchID_t topatch)
 *   -  void copy(const Domain &dom, PatchID_t frompatch, PatchID_t topatch)
 *   -  void sync()
 *
 * BackFill and ShiftUp are tag classes used to indicate how elements should
 * be deleted - either by back-filling (moving elements from the bottom up)
 * or shift-up (just like the "erase" method in the STL vector class).
 *
 * shrink() will make sure the engine within a DynamicArray is not using any
 * more memory than necessary.
 *
 * sync() is something that a user should call if they perform some dynamic
 * create/destroy operations, and then want to use the DynamicArray in
 * expressions that will require knowledge of the global domain of the system.
 * Normally, create/destroy ops only modify the domain information for the
 * patches within the engine that are local to a context.  sync() will
 * syncronize with other contexts to make sure all contexts have up-to-date
 * domain information for the whole domain and all the patches, even those
 * that are on another context.
 */

template<class T = POOMA_DEFAULT_ELEMENT_TYPE,
  class EngineTag = POOMA_DEFAULT_DYNAMIC_ENGINE_TYPE>
class DynamicArray : public Array<1, T, EngineTag>
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  //---------------------------------------------------------------------------
  // Base_t is a convenience typedef referring to this class's base class.
  // This_t is a convenience typedef referring to this class.

  typedef Array<1, T, EngineTag>           Base_t;
  typedef DynamicArray<T, EngineTag>       This_t;
  typedef typename Base_t::Engine_t        Engine_t;
  typedef typename Base_t::Element_t       Element_t;
  typedef typename Base_t::ElementRef_t    ElementRef_t;
  typedef typename Base_t::Domain_t        Domain_t;
  typedef typename Engine_t::Layout_t      Layout_t;
  typedef typename Layout_t::PatchID_t     PatchID_t;
  typedef typename Layout_t::CreateSize_t  CreateSize_t;
  typedef EngineTag                        EngineTag_t;

  enum { dynamic = Engine_t::dynamic };


  //===========================================================================
  // Constructors
  //===========================================================================

  //---------------------------------------------------------------------------
  /// Default constructor for DynamicArray. Exists so this can be resized
  /// or given another layout.

  DynamicArray()
  {
    CTAssert(dynamic == true);
  }

  //---------------------------------------------------------------------------
  // Engine DynamicArray constructors. 

  explicit DynamicArray(const Engine_t &modelEngine)
    : Array<1, T, EngineTag>(modelEngine)
  {
    CTAssert(dynamic == true);
  }

  template<class T2, class EngineTag2>
  explicit DynamicArray(const Engine<1, T2, EngineTag2> &engine)
    : Array<1, T, EngineTag>(engine) 
  { 
    CTAssert(dynamic == true);
  }

  template<class T2, class EngineTag2, class Initializer>
  DynamicArray(const Engine<1, T2, EngineTag2> &engine, const Initializer &init)
    : Array<1, T, EngineTag>(engine, init) 
  { 
    CTAssert(dynamic == true);
  }

  /// One way for a user to construct a DynamicArray is to use another as a 
  /// model. The non-templated version is the copy constructor.
  /// For the templated versions to work, Engine_t must possess a constructor
  /// that takes an OtherEngine and, perhaps, an OtherDomain.

  DynamicArray(const This_t &model)
    : Array<1, T, EngineTag>(model.engine())
  {
    CTAssert(dynamic == true);
  }

  template<class OtherT, class OtherEngineTag, class OtherDomain>
  DynamicArray(const DynamicArray<OtherT, OtherEngineTag> &model, 
	       const OtherDomain &domain)
    : Array<1, T, EngineTag>(model.engine(), domain)
  {
    CTAssert(dynamic == true);
  }

  /// All of these constructors pass domain information to the engine. 
  /// Since DynamicArray is by definition 1-D, we only need the one-argument
  /// version.  These constructors call the default constructor
  /// for Element_t. 

  template<class Sub1>
  explicit DynamicArray(const Sub1 &s1)
    : Array<1, T, EngineTag>(s1)
  {
    CTAssert(dynamic == true);
  }

  /// All of these constructors pass domain information to the engine along
  /// with a model element, which is used to initialize all array elements.
  /// Therefore, it is assumed that Element_t either has deep-copy semantics
  /// by default, or can be made to have it using a wrapper class.
  /// Since DynamicArray is by definition 1-D, we only need the one-argument
  /// version.

  template<class Sub1>
  DynamicArray(const Sub1 &s1, const ModelElement<Element_t> &model)
    : Array<1, T, EngineTag>(s1, model)
  {
    CTAssert(dynamic == true);
  }


  //===========================================================================
  // Destructor
  //===========================================================================
  
  //---------------------------------------------------------------------------
  // The destructor is trivial since the base class is the only thing that
  // owns data.

  ~DynamicArray()
  {
  }


  //===========================================================================
  // Accessors and mutators
  //===========================================================================

  /// Get a reference to this same object, but cast to the Array base class
   
  inline Base_t &array()
  {
    return *this;
  }

  inline const Base_t &array() const
  {
    return *this;
  }

  inline Base_t &arrayAll()
  {
    return *this;
  }

  inline const Base_t &arrayAll() const
  {
    return *this;
  }

  /// Return a reference to the layout for this array.

  inline Layout_t &layout()
  {
    return this->engine().layout();
  }

  inline const Layout_t &layout() const
  {
    return this->engine().layout();
  }


  //===========================================================================
  // Dynamic interface methods.
  //===========================================================================

  ///@name Dynamic interface methods
  /// Dynamic interface methods are used to create or
  /// destroy elements in the array (this is what makes this class dynamic).
  /// These operations are passed on to the engine, so this class must
  /// be used with engines that support dynamic operations.  These
  /// operations will modify the domain of the array, and will rearrange
  /// values, but will not erase the values of existing elements that are
  /// meant to remain after the create/destroy operations.
  ///
  /// NOTE: Create/destroy operations only properly update the domain
  /// information for patches that are on the local context.  If the user
  /// needs to have all contexts have correct layout information about the
  /// dynamic array following create/destroy operations, they must call
  /// the sync() method, in an SPMD operations (that is, all contexts must
  /// call sync() at the same general time).

  //@{

  /// Create new elements for a DynamicArray, by extending the current domain
  /// on the local context by the requested number of elements, or the
  /// specified local patch.
  /// 'local' means on this same context.  The patch is refered to
  /// by local index, from 0 ... # local patches - 1.

  inline void create(CreateSize_t num)
  {
    this->engine().create(num);
  }

  inline void create(CreateSize_t num, PatchID_t patch)
  {
    this->engine().create(num, patch);
  }

  /// Delete the elements within our domain specified by the given
  /// Range, using either a backfill mechanism or a shift-up mechanism
  /// to delete the data.  The second argument should be one of the
  /// following types:
  ///  - BackFill() will move elements from the bottom up to fill the holes.
  ///  - ShiftUp() will shift elements up to fill in holes.
  /// The domain must be within the total domain of the DynamicArray.

  template <class Dom>
  inline void destroy(const Dom &killlist, BackFill method)
  {
    this->engine().destroy(killlist, method);
  }

  template <class Dom>
  inline void destroy(const Dom &killlist, ShiftUp method)
  {
    this->engine().destroy(killlist, method);
  }

  /// Use the default destroy method, BackFill.

  template <class Dom>
  inline void destroy(const Dom &killlist)
  {
    this->engine().destroy(killlist);
  }

  /// Versions that take a pair of random-access iterators. 

  template <class Iter>
  inline void destroy(Iter begin, Iter end, BackFill method)
  {
    Pooma::IteratorPairDomain<Iter> dom(begin, end);
    this->engine().destroy(dom, method);
  }

  template <class Iter>
  inline void destroy(Iter begin, Iter end, ShiftUp method)
  {
    Pooma::IteratorPairDomain<Iter> dom(begin, end);
    this->engine().destroy(dom, method);
  }

  template <class Iter>
  inline void destroy(Iter begin, Iter end)
  {
    Pooma::IteratorPairDomain<Iter> dom(begin, end);
    this->engine().destroy(dom);
  }

  /// Delete the elements within the specific local domain for the given
  /// patch.  The domain values in this case should all be zero-based,
  /// so that they are relative to the first element of the specified
  /// local patch.  The domain values should all be contained within the
  /// specified local patch as well.  To perform cross-patch destroy's,
  /// use the form where you do not specify a local patch number.

  template <class Dom>
  inline void destroy(const Dom &killlist, PatchID_t frompatch,
		      BackFill method)
  {
    this->engine().destroy(killlist, frompatch, method);
  }

  template <class Dom>
  inline void destroy(const Dom &killlist, PatchID_t frompatch,
		      ShiftUp method)
  {
    this->engine().destroy(killlist, frompatch, method);
  }

  template <class Dom>
  inline void destroy(const Dom &killlist, PatchID_t frompatch)
  {
    this->engine().destroy(killlist, frompatch, BackFill());
  }

  /// For the patch-specific destroy operations, the iterators must be
  /// convertible to "const int *" iterators as that is the type for
  /// which the underlying DynamicEventDomain objects are specialized.

  template <class Iter>
  inline void destroy(Iter begin, Iter end, PatchID_t frompatch,
                      BackFill method)
  {
    Pooma::IteratorPairDomain<const int *> dom(begin, end);
    this->engine().destroy(dom, frompatch, method);
  }

  template <class Iter>
  inline void destroy(Iter begin, Iter end, PatchID_t frompatch,
                      ShiftUp method)
  {
    Pooma::IteratorPairDomain<const int *> dom(begin, end);
    this->engine().destroy(dom, frompatch, method);
  }

  template <class Iter>
  inline void destroy(Iter begin, Iter end, PatchID_t frompatch)
  {
    Pooma::IteratorPairDomain<const int *> dom(begin, end);
    this->engine().destroy(dom, frompatch, BackFill());
  }

  /// Copy all elements of domain n to the end of the last patch or
  /// to the end of the specified patch.

  template<class Dom>
  inline void copy(const Dom &copylist)
  {
    this->engine().copy(copylist);
  }

  template<class Dom>
  inline void copy(const Dom &copylist, PatchID_t patch)
  {
    this->engine().copy(copylist, patch);
  }

  /// Copy all elements from the specified local patch to the end of
  /// the local to-patch.  The domain values in this case should all be zero-
  /// based,/ so that they are relative to the first element of the specified
  /// local patch.  The domain values should all be contained within the
  /// specified local patch as well.  To perform cross-patch copies,
  /// use the form where you do not specify a local from-patch number.

  template<class Dom>
  inline void copy(const Dom &copylist, PatchID_t frompatch, PatchID_t topatch)
  {
    this->engine().copy(copylist, frompatch, topatch);
  }

  /// Synchronize all the contexts to update their domain information.
  /// This should be used after create/destroy operations have modified
  /// the domain of local context's data, and all contexts must be told
  /// of the new situation.  This should be an SPMD call.

  void sync()
  {
    this->engine().sync();
  }

  //@}

  //---------------------------------------------------------------------------
  // Copy assignment operators.  Just use the versions from the base class.
  //
  // Note: PARTIAL ORDERING is required for the scalar assignment operator 
  // to be distinguishable from the others.

  inline This_t &operator=(const DynamicArray<T, EngineTag> &rhs)
  {
    Base_t::operator=(static_cast<const Base_t &>(rhs));
    return *this;
  }

  inline const This_t &operator=(const DynamicArray<T, EngineTag> &rhs) const
  {
    Base_t::operator=(static_cast<const Base_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator=(const T1 &rhs) const
  {
    Base_t::operator=(rhs);
    return *this;
  }


  //---------------------------------------------------------------------------
  // Op-assignment operators. 
  //
  // Note: PARTIAL ORDERING is required for the scalar assignment operator 
  // to be distinguishable from the others.

  // Addition.

  template<class OtherT, class OtherETag>
  const This_t &operator+=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator+=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator+=(const T1 &rhs) const
  {
    Base_t::operator+=(rhs);
    return *this;
  }

  // Subtraction.

  template<class OtherT, class OtherETag>
  const This_t &operator-=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator-=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator-=(const T1 &rhs) const
  {
    Base_t::operator-=(rhs);
    return *this;
  }

  // Multiplication.

  template<class OtherT, class OtherETag>
  const This_t &operator*=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator*=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator*=(const T1 &rhs) const
  {
    Base_t::operator*=(rhs);
    return *this;
  }

  // Division.

  template<class OtherT, class OtherETag>
  const This_t &operator/=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator/=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator/=(const T1 &rhs) const
  {
    Base_t::operator/=(rhs);
    return *this;
  }

  // Modulus.

  template<class OtherT, class OtherETag>
  const This_t &operator%=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator%=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator%=(const T1 &rhs) const
  {
    Base_t::operator%=(rhs);
    return *this;
  }

  // Bitwise-Or.

  template<class OtherT, class OtherETag>
  const This_t &operator|=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator|=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator|=(const T1 &rhs) const
  {
    Base_t::operator|=(rhs);
    return *this;
  }

  // Bitwise-And.

  template<class OtherT, class OtherETag>
  const This_t &operator&=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator&=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator&=(const T1 &rhs) const
  {
    Base_t::operator&=(rhs);
    return *this;
  }

  // Bitwise-Xor.

  template<class OtherT, class OtherETag>
  const This_t &operator^=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator^=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator^=(const T1 &rhs) const
  {
    Base_t::operator^=(rhs);
    return *this;
  }

  // Left shift.

  template<class OtherT, class OtherETag>
  const This_t &operator<<=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator<<=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator<<=(const T1 &rhs) const
  {
    Base_t::operator<<=(rhs);
    return *this;
  }

  // Right shift.

  template<class OtherT, class OtherETag>
  const This_t &operator>>=(const DynamicArray<OtherT, OtherETag> &rhs) const
  {
    typedef Array<1,OtherT,OtherETag> Array_t;
    Base_t::operator>>=(static_cast<const Array_t &>(rhs));
    return *this;
  }

  template<class T1>
  const This_t &operator>>=(const T1 &rhs) const
  {
    Base_t::operator>>=(rhs);
    return *this;
  }
};


//---------------------------------------------------------------------------
/// Implementation for assignment operators.  We define here the
/// extra versions of assign that tell how to assign
///    Array = DynamicArray
/// the other versions of assign are handled by routines in the base class.
///
/// Note: PARTIAL ORDERING is required for the scalar assignment operator 
/// to be distinguishable from the others.
//---------------------------------------------------------------------------

template<int Dim, class T, class EngineTag,
         class OtherT, class OtherEngineTag, class Op>
inline const Array<Dim, T, EngineTag> &
assign(const Array<Dim, T, EngineTag> &lhs,
       const DynamicArray<OtherT, OtherEngineTag> &rhs,
       const Op &op)
{
  // Cast the rhs to be an Array, and then invoke the version of
  // assign for Array = Array.

  typedef Array<1, OtherT, OtherEngineTag> Array_t;

  const Array_t &arhs = rhs;

  return assign(lhs, arhs, op);
}


//-----------------------------------------------------------------------------
/// Traits class telling RefCountedBlockPointer that this class has
/// shallow semantics and a makeOwnCopy method.
//-----------------------------------------------------------------------------

template <class T, class EngineTag>
struct ElementProperties< DynamicArray<T, EngineTag> > 
  : public MakeOwnCopyProperties< DynamicArray<T, EngineTag> >
{
};


//-----------------------------------------------------------------------------
/// A traits class that tells PETE how to make a DynamicArray into
/// an expression element.  DynamicArray just returns a reference to
/// itself, but cast as an Array.
//-----------------------------------------------------------------------------

template<class T, class EngineTag>
struct CreateLeaf< DynamicArray<T, EngineTag> >
{
  typedef DynamicArray<T, EngineTag>         Input_t;
  typedef Reference<Array<1, T, EngineTag> > Leaf_t;
  typedef Reference<Array<1, T, EngineTag> > Return_t;
  inline static
  Return_t make(const Input_t &a)
  {
    return Return_t(a);
  }
};

//-----------------------------------------------------------------------------
/// Generalized Engine Functors.
//-----------------------------------------------------------------------------

template<class T, class E, class Tag>
struct LeafFunctor<DynamicArray<T, E>, EngineFunctorTag<Tag> >
{
  typedef typename DynamicArray<T,E>::Engine_t Engine_t;
  typedef typename EngineFunctor<Engine_t,Tag>::Type_t Type_t;
  inline static
  Type_t apply(const DynamicArray<T, E> &array,
	       const EngineFunctorTag<Tag> &tag)
  {
    return EngineFunctor<Engine_t,Tag>::apply(array.engine(), tag.tag());
  }
};


#endif

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DynamicArray.h,v $   $Author: richard $
// $Revision: 1.35 $   $Date: 2004/11/01 18:16:34 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
