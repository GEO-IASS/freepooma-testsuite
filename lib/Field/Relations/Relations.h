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
// ----------------------------------------------------------------------
// Relation: intermediate templated base class for all relations.
// RelationCategory: base class for relation categories.
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Relations
 * @brief
 * Relation basics.
 */

#ifndef POOMA_FIELD_RELATIONS_RELATIONS_H
#define POOMA_FIELD_RELATIONS_RELATIONS_H

#include "Field/Relations/RelationBases.h"  // Base classes
#include "Field/Relations/RelationList.h"   // Class member

/**
 * InfluenceRelation supports the relation package by allowing fields that
 * have been modified---resulting in their dirty flag being set---to notify
 * dependent fields and set their relations' dirty flags as well.
 *
 * An InfluenceRelation should be added to the relation list of all fields
 * that influence other fields. This typically occurs when the field appears
 * on the RHS of a relation. The field being influenced, typically on the
 * LHS of a relation, should be the argument for the constructor.
 */

class InfluenceRelation : public RelationListItem
{
public:

  //---------------------------------------------------------------------------
  // Constructors.
  
  // Set the priority and store a weak reference to the target's relation list.

  template<class Target>
  InfluenceRelation(const Target &t)
  : RelationListItem(), 
    list_m(const_cast<RelationList*>(&t.fieldEngine().relations()))
    { 
      setPriority(100u);
    }

  // Copy constructor.
  
  InfluenceRelation(const InfluenceRelation &model)
  : RelationListItem(model), list_m(model.list_m)
    { }
    
  //---------------------------------------------------------------------------
  // Destructor. We maintain a weak reference to the relation list, so we have
  // nothing to do in the destructor.

  InfluenceRelation() { }

  //---------------------------------------------------------------------------
  // Methods.

  // All of the action is in setDirty(). If we are already dirty, there's
  // nothing to do. If clean, we set our dirty flag and then those in the
  // target's relation list.
    
  void setDirty()
    {
      if (!dirty())
        {
          RelationListItem::setDirty();
          list_m->setDirty();
        }
    }
      
  // We do nothing in response to an apply.

  void apply() { } 

private:

  RelationList *list_m;
};


/** Relation0 is a template used to construct relations, such as boundary
 * conditions, that do not depend on additional fields. The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L>
 *   RelationFunctor(const L &, const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L>
 *   void operator()(const L &) const
 *
 * where L is a Field. The constructor should use the arguments to
 * initialize itself and the function should apply the relation.
 */

template<class Target, class RelationFunctor>
class Relation0: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation0(const Target &t, const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation0() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    return new Relation0<Target, RelationFunctor>(target, this->functor_m);
  }

};


/** Relation1 is a template used to construct relations that depend on one
 * additional field (e.g., a = b). The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L, class R1, class R2>
 *   RelationFunctor(const L &, const R1 &, 
 *                   const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L, class R1>
 *   void operator()(const L &const R1 &) const
 *
 * where L and R1 are Fields. The constructor should use the arguments to
 * initialize itself and the function should apply the relation.
 */

template<class Target, class R1, class RelationFunctor>
class Relation1: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation1(const Target &t, const R1 &r, 
    const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f), r1_m(r)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation1() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m, r1_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    r1_m.addRelation(new InfluenceRelation(target));
    
    return new Relation1<Target, R1, RelationFunctor>
                 (target, r1_m, this->functor_m);
  }
  
protected:

  R1 r1_m;
};


/** Relation2 is a template used to construct relations that depend on two
 * additional fields (e.g., a = b + c). The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L, class R1, class R2>
 *   RelationFunctor(const L &, const R1 &, const R2 &, 
 *                   const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L, class R1, class R2>
 *   void operator()(const L &const R1 &, const R2 &) const
 *
 * where L, R1, and R2 are Fields. The constructor should use the arguments to
 * initialize itself and the function should apply the relation.
 */

template<class Target, class R1, class R2, class RelationFunctor>
class Relation2: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation2(const Target &t, const R1 &r1, const R2 &r2, 
    const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f), r1_m(r1), r2_m(r2)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation2() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m, r1_m, r2_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    r1_m.addRelation(new InfluenceRelation(target));
    r2_m.addRelation(new InfluenceRelation(target));
    
    return new Relation2<Target, R1, R2, RelationFunctor>
                 (target, r1_m, r2_m, this->functor_m);
  }
  
protected:

  R1 r1_m;
  R2 r2_m;
};


/** Relation3 is a template used to construct relations that depend on three
 * additional fields (e.g., a = b + c + d). The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L, class R1, class R2, class R3>
 *   RelationFunctor(const L &, const R1 &, const R2 &, const R3 &,
 *                   const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L, class R1, class R2, class R3>
 *   void operator()(const L &const R1 &, const R2 &, const R3 &) const
 *
 * where L, R1, R2, and R3 are Fields. The constructor should use 
 * the arguments to initialize itself and the function should apply the 
 * relation.
 */

template<class Target, class R1, class R2, class R3, class RelationFunctor>
class Relation3: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation3(const Target &t, const R1 &r1, const R2 &r2, const R3 &r3,
            const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f), 
    r1_m(r1), 
    r2_m(r2),
    r3_m(r3)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation3() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m, r1_m, r2_m, r3_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    r1_m.addRelation(new InfluenceRelation(target));
    r2_m.addRelation(new InfluenceRelation(target));
    r3_m.addRelation(new InfluenceRelation(target));
    
    return new Relation3<Target, R1, R2, R3, RelationFunctor>
      (target, r1_m, r2_m, r3_m, this->functor_m);
  }
  
protected:

  R1 r1_m;
  R2 r2_m;
  R3 r3_m;
};


/** Relation4 is a template used to construct relations that depend on four
 * additional fields (e.g., a = b + c + d + e). The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L, class R1, class R2, class R3, class R4>
 *   RelationFunctor(const L &, const R1 &, const R2 &, const R3 &,
 *                   const R4 &, const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L, class R1, class R2, class R3, class R4>
 *   void operator()(const L &const R1 &, const R2 &, const R3 &,
 *                   const R4 &) const
 *
 * where L, R1, R2, R3, and R4 are Fields. The constructor should use 
 * the arguments to initialize itself and the function should apply the 
 * relation.
 */

template<class Target, class R1, class R2, class R3, class R4,
         class RelationFunctor>
class Relation4: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation4(const Target &t, const R1 &r1, const R2 &r2, const R3 &r3,
            const R4 &r4, const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f), 
    r1_m(r1), 
    r2_m(r2),
    r3_m(r3),
    r4_m(r4)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation4() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m, r1_m, r2_m, r3_m, r4_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    r1_m.addRelation(new InfluenceRelation(target));
    r2_m.addRelation(new InfluenceRelation(target));
    r3_m.addRelation(new InfluenceRelation(target));
    r4_m.addRelation(new InfluenceRelation(target));
    
    return new Relation4<Target, R1, R2, R3, R4, RelationFunctor>
      (target, r1_m, r2_m, r3_m, r4_m, this->functor_m);
  }
  
protected:

  R1 r1_m;
  R2 r2_m;
  R3 r3_m;
  R4 r4_m;
};


/** Relation5 is a template used to construct relations that depend on five
 * additional fields (e.g., a = b + c + d + e + f). The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L, class R1, class R2, class R3, class R4, class R5>
 *   RelationFunctor(const L &, const R1 &, const R2 &, const R3 &,
 *                   const R4 &, const R5 &, const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L, class R1, class R2, class R3, class R4, class R5>
 *   void operator()(const L &const R1 &, const R2 &, const R3 &,
 *                   const R4 &, const R5 &) const
 *
 * where L, R1, R2, R3, R4, and R5 are Fields. The constructor should use 
 * the arguments to initialize itself and the function should apply the 
 * relation.
 */

template<class Target, class R1, class R2, class R3, class R4, class R5,
         class RelationFunctor>
class Relation5: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation5(const Target &t, const R1 &r1, const R2 &r2, const R3 &r3,
            const R4 &r4, const R5 &r5, const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f), 
    r1_m(r1), 
    r2_m(r2),
    r3_m(r3),
    r4_m(r4),
    r5_m(r5)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation5() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m, r1_m, r2_m, r3_m, r4_m, r5_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    r1_m.addRelation(new InfluenceRelation(target));
    r2_m.addRelation(new InfluenceRelation(target));
    r3_m.addRelation(new InfluenceRelation(target));
    r4_m.addRelation(new InfluenceRelation(target));
    r5_m.addRelation(new InfluenceRelation(target));
    
    return new Relation5<Target, R1, R2, R3, R4, R5, RelationFunctor>
      (target, r1_m, r2_m, r3_m, r4_m, r5_m, this->functor_m);
  }
  
protected:

  R1 r1_m;
  R2 r2_m;
  R3 r3_m;
  R4 r4_m;
  R5 r5_m;
};


/** Relation6 is a template used to construct relations that depend on six
 * additional fields (e.g., a = b + c + d + e + f + g). The Target must be
 * a Field. The RelationFunctor must be Default Constructable and Assignable.
 * In addition, it must provide the constructor
 *
 *   template<class L, class R1, class R2, class R3, class R4, class R5,
 *            class R6>
 *   RelationFunctor(const L &, const R1 &, const R2 &, const R3 &,
 *                   const R4 &, const R5 &, const R6 &, 
 *                   const RelationFunctor &)
 *
 * and the member function
 *
 *   template<class L, class R1, class R2, class R3, class R4, class R5,
 *            class R6>
 *   void operator()(const L &const R1 &, const R2 &, const R3 &,
 *                   const R4 &, const R5 &, const R6 &) const
 *
 * where L, R1, R2, R3, R4, R5, and R6 are Fields. The constructor should use 
 * the arguments to initialize itself and the function should apply the 
 * relation.
 */

template<class Target, class R1, class R2, class R3, class R4, class R5,
         class R6, class RelationFunctor>
class Relation6: public RelationBase<Target, RelationFunctor>
{
public:

  //---------------------------------------------------------------------------
  // Constructors.

  Relation6(const Target &t, const R1 &r1, const R2 &r2, const R3 &r3,
            const R4 &r4, const R5 &r5, const R6 &r6, const RelationFunctor &f)
  : RelationBase<Target, RelationFunctor>(t, f), 
    r1_m(r1), 
    r2_m(r2),
    r3_m(r3),
    r4_m(r4),
    r5_m(r5),
    r6_m(r6)
  { }

  //---------------------------------------------------------------------------
  // Destructor.

  ~Relation6() 
  { }

  //---------------------------------------------------------------------------
  // Methods.

  // Apply the supplied functor.
  
  void apply()
  {
    this->functor_m(this->target_m, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m);
  }

  virtual RelationListItem *retarget(const Target &target) const
  {
    r1_m.addRelation(new InfluenceRelation(target));
    r2_m.addRelation(new InfluenceRelation(target));
    r3_m.addRelation(new InfluenceRelation(target));
    r4_m.addRelation(new InfluenceRelation(target));
    r5_m.addRelation(new InfluenceRelation(target));
    r6_m.addRelation(new InfluenceRelation(target));
    
    return new Relation6<Target, R1, R2, R3, R4, R5, R6, RelationFunctor>
      (target, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m, this->functor_m);
  }
  
protected:

  R1 r1_m;
  R2 r2_m;
  R3 r3_m;
  R4 r4_m;
  R5 r5_m;
  R6 r6_m;
};


/**
 * Relation functors supporting the use of function pointers.
 */

template<class L>
class RelationFunctionPtr0 {
public:

  RelationFunctionPtr0(void (*f)(const L &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr0(const RelationFunctionPtr0<L> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l)
  {
    f_m(l);
  }

private:

  void (*f_m)(const L &);
};

template<class L, class R1>
class RelationFunctionPtr1 {
public:

  RelationFunctionPtr1(void (*f)(const L &, const R1 &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr1(const RelationFunctionPtr1<L, R1> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1)
  {
    f_m(l, r1);
  }

private:

  void (*f_m)(const L &, const R1 &);
};

template<class L, class R1, class R2>
class RelationFunctionPtr2 {
public:

  RelationFunctionPtr2(void (*f)(const L &, const R1 &, const R2 &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr2(const RelationFunctionPtr2<L, R1, R2> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2)
  {
    f_m(l, r1, r2);
  }

private:

  void (*f_m)(const L &, const R1 &, const R2 &);
};

template<class L, class R1, class R2, class R3>
class RelationFunctionPtr3 {
public:

  RelationFunctionPtr3(void (*f)(const L &, 
    const R1 &, const R2 &, const R3 &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr3(
    const RelationFunctionPtr3<L, R1, R2, R3> &model)
  : f_m(model.f_m)
  { }
  RelationFunctionPtr3(
    const RelationFunctionPtr3<L, R1, R2, R3> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3)
  {
    f_m(l, r1, r2, r3);
  }

private:

  void (*f_m)(const L &, const R1 &, const R2 &, const R3 &);
};

template<class L, class R1, class R2, class R3, class R4>
class RelationFunctionPtr4 {
public:

  RelationFunctionPtr4(void (*f)(const L &, 
    const R1 &, const R2 &, const R3 &, const R4 &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr4(
    const RelationFunctionPtr4<L, R1, R2, R3, R4> &model)
  : f_m(model.f_m)
  { }
  RelationFunctionPtr4(
    const RelationFunctionPtr4<L, R1, R2, R3, R4> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3, const R4 &r4)
  {
    f_m(l, r1, r2, r3, r4);
  }

private:

  void (*f_m)(const L &, const R1 &, const R2 &, 
    const R3 &, const R4 &);
};

template<class L, 
  class R1, class R2, class R3, class R4, class R5>
class RelationFunctionPtr5 {
public:

  RelationFunctionPtr5(void (*f)(const L &, 
    const R1 &, const R2 &, const R3 &, const R4 &, const R5 &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr5(
    const RelationFunctionPtr5<L, R1, R2, R3, R4, R5> &model)
  : f_m(model.f_m)
  { }
  RelationFunctionPtr5(
    const RelationFunctionPtr5<L, R1, R2, R3, R4, R5> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3, const R4 &r4, const R5 &r5)
  {
    f_m(l, r1, r2, r3, r4, r5);
  }

private:

  void (*f_m)(const L &, const R1 &, const R2 &, 
    const R3 &, const R4 &, const R5 &);
};

template<class L, 
  class R1, class R2, class R3, class R4, class R5, class R6>
class RelationFunctionPtr6 {
public:

  RelationFunctionPtr6(void (*f)(const L &, 
    const R1 &, const R2 &, const R3 &, const R4 &, const R5 &, const R6 &)) 
  : f_m(f) 
  { }
  RelationFunctionPtr6(
    const RelationFunctionPtr6<L, R1, R2, R3, R4, R5, R6> &model)
  : f_m(model.f_m)
  { }
  RelationFunctionPtr6(
    const RelationFunctionPtr6<L, R1, R2, R3, R4, R5, R6> &init, const L &)
  : f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3, const R4 &r4, const R5 &r5, const R6 &r6)
  {
    f_m(l, r1, r2, r3, r4, r5, r6);
  }

private:

  void (*f_m)(const L &, const R1 &, const R2 &, 
    const R3 &, const R4 &, const R5 &, const R6 &);
};


/**
 * Relation functors supporting the use of member function pointers.
 */

template<class C, class L>
class RelationMemberPtr0 {
public:

  RelationMemberPtr0(const C &obj, void (C::*f)(const L &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr0(const RelationMemberPtr0<C, L> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr0(const RelationMemberPtr0<C, L> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l)
  {
    (obj_m.*f_m)(l);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &);
};

template<class C, class L, class R1>
class RelationMemberPtr1 {
public:

  RelationMemberPtr1(const C &obj, void (C::*f)(const L &, 
    const R1 &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr1(const RelationMemberPtr1<C, L, R1> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr1(const RelationMemberPtr1<C, L, R1> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1)
  {
    (obj_m.*f_m)(l, r1);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &, const R1 &);
};

template<class C, class L, class R1, class R2>
class RelationMemberPtr2 {
public:

  RelationMemberPtr2(const C &obj, void (C::*f)(const L &, 
    const R1 &, const R2 &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr2(const RelationMemberPtr2<C, L, R1, R2> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr2(const RelationMemberPtr2<C, L, R1, R2> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2)
  {
    (obj_m.*f_m)(l, r1, r2);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &, const R1 &, const R2 &);
};

template<class C, class L, class R1, class R2, class R3>
class RelationMemberPtr3 {
public:

  RelationMemberPtr3(const C &obj, void (C::*f)(const L &, 
    const R1 &, const R2 &, const R3 &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr3(
    const RelationMemberPtr3<C, L, R1, R2, R3> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr3(
    const RelationMemberPtr3<C, L, R1, R2, R3> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3)
  {
    (obj_m.*f_m)(l, r1, r2, r3);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &, const R1 &, const R2 &, const R3 &);
};

template<class C, class L, class R1, class R2, class R3, class R4>
class RelationMemberPtr4 {
public:

  RelationMemberPtr4(const C &obj, void (C::*f)(const L &, 
    const R1 &, const R2 &, const R3 &, const R4 &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr4(
    const RelationMemberPtr4<C, L, R1, R2, R3, R4> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr4(
    const RelationMemberPtr4<C, L, R1, R2, R3, R4> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3, const R4 &r4)
  {
    (obj_m.*f_m)(l, r1, r2, r3, r4);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &, const R1 &, const R2 &, 
    const R3 &, const R4 &);
};

template<class C, class L, 
  class R1, class R2, class R3, class R4, class R5>
class RelationMemberPtr5 {
public:

  RelationMemberPtr5(const C &obj, void (C::*f)(const L &, 
    const R1 &, const R2 &, const R3 &, const R4 &, const R5 &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr5(
    const RelationMemberPtr5<C, L, R1, R2, R3, R4, R5> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr5(
    const RelationMemberPtr5<C, L, R1, R2, R3, R4, R5> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3, const R4 &r4, const R5 &r5)
  {
    (obj_m.*f_m)(l, r1, r2, r3, r4, r5);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &, const R1 &, const R2 &, 
    const R3 &, const R4 &, const R5 &);
};

template<class C, class L, 
  class R1, class R2, class R3, class R4, class R5, class R6>
class RelationMemberPtr6 {
public:

  RelationMemberPtr6(const C &obj, void (C::*f)(const L &, 
    const R1 &, const R2 &, const R3 &, const R4 &, const R5 &, const R6 &)) 
  : obj_m(obj), f_m(f) 
  { }
  RelationMemberPtr6(
    const RelationMemberPtr6<C, L, R1, R2, R3, R4, R5, R6> &model)
  : obj_m(model.obj_m), f_m(model.f_m)
  { }
  RelationMemberPtr6(
    const RelationMemberPtr6<C, L, R1, R2, R3, R4, R5, R6> &init, const L &)
  : obj_m(init.obj_m), f_m(init.f_m) 
  { }
    
  inline void operator()(const L &l, const R1 &r1, const R2 &r2, 
    const R3 &r3, const R4 &r4, const R5 &r5, const R6 &r6)
  {
    (obj_m.*f_m)(l, r1, r2, r3, r4, r5, r6);
  }

private:

  C obj_m;
  void (C::*f_m)(const L &, const R1 &, const R2 &, 
    const R3 &, const R4 &, const R5 &, const R6 &);
};
      

/**
 * RelationFunctorTraits is used to specify characteristics of a Relation
 * functor. Currently, the only one defined is defaultPriority.
 */

template<class RelationFunctor>
struct RelationFunctorTraits {
  
  enum { defaultPriority = 0 }; 

};


namespace Pooma {

  ///@name Standalone functions for creating relations.
  //@{

  //---------------------------------------------------------------------------
  // Functor versions
  
  template<class RelationFunctor, class L>
  void newRelation(const RelationFunctor &f, const L &l)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            RelationListItem *r = new Relation0<L, RelationFunctor>(lsub, f);
            r->setPriority(RelationFunctorTraits<RelationFunctor>::defaultPriority);
            lsub.addRelation(r);
          }
      }
  }
  
  template<class RelationFunctor, class L, class R1>
  void newRelation(const RelationFunctor &f, const L &l, 
                   const R1 &r1)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            const R1 &r1sub = r1.subField(m, c);
            
            r1sub.addRelation(new InfluenceRelation(lsub));
  
            RelationListItem *r = new 
              Relation1<L, R1, RelationFunctor>
                (lsub, r1sub, f);
            lsub.addRelation(r);
          }
      }
  }
  
  template<class RelationFunctor, class L, class R1, class R2>
  void newRelation(const RelationFunctor &f, const L &l, 
                   const R1 &r1, const R2 &r2)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            const R1 &r1sub = r1.subField(m, c);
            const R2 &r2sub = r2.subField(m, c);
            
            r1sub.addRelation(new InfluenceRelation(lsub));
            r2sub.addRelation(new InfluenceRelation(lsub));
  
            RelationListItem *r = 
              new Relation2<L, R1, R2, RelationFunctor>(lsub, r1sub, r2sub, f);
            lsub.addRelation(r);
          }
      }
  }
  
  template<class RelationFunctor, class L, class R1, class R2, class R3>
  void newRelation(const RelationFunctor &f, const L &l, 
                   const R1 &r1, const R2 &r2, const R3 &r3)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            const R1 &r1sub = r1.subField(m, c);
            const R2 &r2sub = r2.subField(m, c);
            const R3 &r3sub = r3.subField(m, c);
            
            r1sub.addRelation(new InfluenceRelation(lsub));
            r2sub.addRelation(new InfluenceRelation(lsub));
            r3sub.addRelation(new InfluenceRelation(lsub));
  
            RelationListItem *r = new 
              Relation3<L, R1, R2, R3, RelationFunctor>
                (lsub, r1sub, r2sub, r3sub, f);
            lsub.addRelation(r);
          }
      }
  }
  
  template<class RelationFunctor, class L, class R1, class R2, class R3,
    class R4>
  void newRelation(const RelationFunctor &f, const L &l, 
                   const R1 &r1, const R2 &r2, 
                   const R3 &r3, const R4 &r4)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            const R1 &r1sub = r1.subField(m, c);
            const R2 &r2sub = r2.subField(m, c);
            const R3 &r3sub = r3.subField(m, c);
            const R4 &r4sub = r4.subField(m, c);
            
            r1sub.addRelation(new InfluenceRelation(lsub));
            r2sub.addRelation(new InfluenceRelation(lsub));
            r3sub.addRelation(new InfluenceRelation(lsub));
            r4sub.addRelation(new InfluenceRelation(lsub));
  
            RelationListItem *r = new 
              Relation4<L, R1, R2, R3, R4, RelationFunctor>
                (lsub, r1sub, r2sub, r3sub, r4sub, f);
            lsub.addRelation(r);
          }
      }
  }
  
  template<class RelationFunctor, class L, class R1, class R2, class R3,
    class R4, class R5>
  void newRelation(const RelationFunctor &f, const L &l, 
                   const R1 &r1, const R2 &r2, 
                   const R3 &r3, const R4 &r4, const R5 &r5)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            const R1 &r1sub = r1.subField(m, c);
            const R2 &r2sub = r2.subField(m, c);
            const R3 &r3sub = r3.subField(m, c);
            const R4 &r4sub = r4.subField(m, c);
            const R5 &r5sub = r5.subField(m, c);
            
            r1sub.addRelation(new InfluenceRelation(lsub));
            r2sub.addRelation(new InfluenceRelation(lsub));
            r3sub.addRelation(new InfluenceRelation(lsub));
            r4sub.addRelation(new InfluenceRelation(lsub));
            r5sub.addRelation(new InfluenceRelation(lsub));
  
            RelationListItem *r = new 
              Relation5<L, R1, R2, R3, R4, R5, RelationFunctor>
                (lsub, r1sub, r2sub, r3sub, r4sub, r5sub, f);
            lsub.addRelation(r);
          }
      }
  }
  
  template<class RelationFunctor, class L, class R1, class R2, class R3,
    class R4, class R5, class R6>
  void newRelation(const RelationFunctor &f, const L &l, 
                   const R1 &r1, const R2 &r2, 
                   const R3 &r3, const R4 &r4, const R5 &r5, const R6 &r6)
  {
    for (int m = 0; m < l.numMaterials(); ++m)
      {
        for (int c = 0; c < l.centeringSize(); ++c)
          {
            const L &lsub = l.subField(m, c);
            const R1 &r1sub = r1.subField(m, c);
            const R2 &r2sub = r2.subField(m, c);
            const R3 &r3sub = r3.subField(m, c);
            const R4 &r4sub = r4.subField(m, c);
            const R5 &r5sub = r5.subField(m, c);
            const R6 &r6sub = r6.subField(m, c);
            
            r1sub.addRelation(new InfluenceRelation(lsub));
            r2sub.addRelation(new InfluenceRelation(lsub));
            r3sub.addRelation(new InfluenceRelation(lsub));
            r4sub.addRelation(new InfluenceRelation(lsub));
            r5sub.addRelation(new InfluenceRelation(lsub));
            r6sub.addRelation(new InfluenceRelation(lsub));
  
            RelationListItem *r = new 
              Relation6<L, R1, R2, R3, R4, R5, R6, RelationFunctor>
                (lsub, r1sub, r2sub, r3sub, r4sub, r5sub, r6sub, f);
            lsub.addRelation(r);
          }
      }
  }

  //---------------------------------------------------------------------------
  // Function pointer versions
  
  template<class L>
  RelationFunctionPtr0<L> 
  functionPtr(void (*f)(const L &))
  {
    return RelationFunctionPtr0<L>(f);
  }
  
  template<class L, class R1>
  RelationFunctionPtr1<L, R1> 
  functionPtr(void (*f)(const L &, const R1 &))
  {
    return RelationFunctionPtr1<L, R1>(f);
  }
  
  template<class L, class R1, class R2>
  RelationFunctionPtr2<L, R1, R2>
  functionPtr(void (*f)(const L &, const R1 &, const R2 &))
  {
    return RelationFunctionPtr2<L, R1, R2>(f);
  }
  
  template<class L, class R1, class R2, class R3>
  RelationFunctionPtr3<L, R1, R2, R3>
  functionPtr(void (*f)(const L &, const R1 &, const R2 &, const R3 &))
  {
    return RelationFunctionPtr3<L, R1, R2, R3>(f);
  }
   
  template<class L, class R1, class R2, class R3, class R4>
  RelationFunctionPtr4<L, R1, R2, R3, R4>
  functionPtr(void (*f)(const L &, const R1 &, const R2 &, const R3 &, 
    const R4 &r4))
  {
    return RelationFunctionPtr4<L, R1, R2, R3, R4>(f);
  }
 
  template<class L, class R1, class R2, class R3,
    class R4, class R5>
  RelationFunctionPtr5<L, R1, R2, R3, R4, R5>
  functionPtr(void (*f)(const L &, const R1 &, const R2 &, const R3 &, 
    const R4 &r4, const R5 &r5))
  {
    return RelationFunctionPtr5<L, R1, R2, R3, R4, R5>(f);
  }
  
  template<class L, class R1, class R2, class R3,
    class R4, class R5, class R6>
  RelationFunctionPtr6<L, R1, R2, R3, R4, R5, R6>
  functionPtr(void (*f)(const L &, const R1 &, const R2 &, const R3 &, 
    const R4 &r4, const R5 &r5, const R6 &r6))
  {
    return RelationFunctionPtr6<L, R1, R2, R3, R4, R5, R6>(f);
  }

  //---------------------------------------------------------------------------
  // Member function pointer versions
  
  template<class C, class L>
  RelationMemberPtr0<C, L> 
  memberPtr(const C &obj, void (C::*f)(const L &))
  {
    return RelationMemberPtr0<C, L>(obj, f);
  }
  
  template<class C, class L, class R1>
  RelationMemberPtr1<C, L, R1> 
  memberPtr(const C &obj, void (C::*f)(const L &, const R1 &))
  {
    return RelationMemberPtr1<C, L, R1>(obj, f);
  }
  
  template<class C, class L, class R1, class R2>
  RelationMemberPtr2<C, L, R1, R2> 
  memberPtr(const C &obj,  void (C::*f)(const L &, const R1 &, const R2 &))
  {
    return RelationMemberPtr2<C, L, R1, R2>(obj, f);
  }
  
  template<class C, class L, class R1, class R2, class R3>
  RelationMemberPtr3<C, L, R1, R2, R3> 
  memberPtr(const C &obj, 
    void (C::*f)(const L &, const R1 &, const R2 &, const R3 &))
  {
    return RelationMemberPtr3<C, L, R1, R2, R3>(obj, f);
  }
  
  template<class C, class L, class R1, class R2, class R3,
    class R4>
  RelationMemberPtr4<C, L, R1, R2, R3, R4> 
  memberPtr(const C &obj, 
    void (C::*f)(const L &, const R1 &, const R2 &, const R3 &, 
    const R4 &r4))
  {
    return RelationMemberPtr4<C, L, R1, R2, R3, R4>(obj, f);
  }
  
  template<class C, class L, class R1, class R2, class R3,
    class R4, class R5>
  RelationMemberPtr5<C, L, R1, R2, R3, R4, R5> 
  memberPtr(const C &obj, 
    void (C::*f)(const L &, const R1 &, const R2 &, const R3 &, 
    const R4 &r4, const R5 &r5))
  {
    return RelationMemberPtr5<C, L, R1, R2, R3, R4, R5>(obj, f);
  }
  
  template<class C, class L, class R1, class R2, class R3,
    class R4, class R5, class R6>
  RelationMemberPtr6<C, L, R1, R2, R3, R4, R5, R6> 
  memberPtr(const C &obj, 
    void (C::*f)(const L &, const R1 &, const R2 &, const R3 &, 
    const R4 &r4, const R5 &r5, const R6 &r6))
  {
    return RelationMemberPtr6<C, L, R1, R2, R3, R4, R5, R6>(obj, f);
  }

  //@}
  
} // namespace Pooma

#endif // POOMA_FIELD_RELATIONS_RELATIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Relations.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
