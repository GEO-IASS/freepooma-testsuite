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
// View1<Object,Domain>
// View2<Object,Domain1,Domain2>
// ...
// View7<Object,Domain1,...,Domain7>
// ComponentView[2-7]<Object>
//-----------------------------------------------------------------------------

#ifndef POOMA_POOMA_VIEW_H
#define POOMA_POOMA_VIEW_H

/** @file
 * @ingroup Pooma
 * @brief
 * View1<Object,Domain>::Type_t is defined for several Pooma objects, and
 * tells you the type of a(b).  You should be able to write code like:
 *
 * <PRE>
 * A a;
 * B b;
 * typename View1<A,B>::Type_t c = a(b);
 * </PRE>
 *
 * ViewN will also give the type for the read() member function in arrays
 * and fields.
 *
 * To define the view properties for a new class, you should specialize
 * View for that class.  Reliance on the partial specialization ordering
 * rules for more than one argument is a very bad thing, so never define
 * View for a particular domain and general A.  (It will be common for
 * us to define View for general domains with a particular A, so
 * specializations with general A would also require complete specializations
 * in order to avoid ambiguity.)
 */


/** View0 enables you to write:
 *
 * <PRE>
 * A a;
 * typename View0<A>::Type_t c = a();
 * </PRE>
 */

template<class Thing>
struct View0
{
};


/** View1 enables you to write:
 *
 * <PRE>
 * A a;
 * B b;
 * typename View1<A,B>::Type_t c = a(b);
 * </PRE>
 */

template<class Thing, class Sub>
struct View1
{
};


/** View2 enables you to write:
 *
 * <PRE>
 * A a;
 * B b;
 * C c;
 * typename View2<A,B,C>::Type_t d = a(b,c);
 * </PRE>
 */

template<class Thing, class Sub1, class Sub2>
struct View2
{
};


/** View3 enables you to write:
 *
 * <PRE>
 * A a;
 * B b;
 * C c;
 * D d;
 * typename View2<A,B,C,D>::Type_t e = a(b,c,d);
 * </PRE>
 */

template<class Thing, class Sub1, class Sub2, class Sub3>
struct View3
{
};


/** View4 enables you to write:
 *
 * <PRE>
 * A a;
 * B b;
 * C c;
 * D d;
 * E e;
 * typename View2<A,B,C,D,E>::Type_t thing = a(b,c,d,e);
 * </PRE>
 */

template<class Thing, class Sub1, class Sub2, class Sub3, class Sub4>
struct View4
{
};


/** View5 enables you to write:
 *
 * <PRE>
 * A a;
 * B b;
 * C c;
 * D d;
 * E e;
 * F f;
 * typename View5<A,B,C,D,E,F>::Type_t thing = a(b,c,d,e,f);
 * </PRE>
 */

template<class Thing, class Sub1, class Sub2, class Sub3, class Sub4,
  class Sub5>
struct View5
{
};


/** View6 enables you to write:
 *
 * <PRE>
 * A a;
 * B b;
 * C c;
 * D d;
 * E e;
 * F f;
 * G g;
 * typename View6<A,B,C,D,E,F,G>::Type_t thing = a(b,c,d,e,f,g);
 * </PRE>
 */

template<class Thing, class Sub1, class Sub2, class Sub3, class Sub4,
  class Sub5, class Sub6>
struct View6
{
};


/** View7 enables you to write:
 * <PRE>
 * A a;
 * B b;
 * C c;
 * D d;
 * E e;
 * F f;
 * G g;
 * H h;
 * typename View7<A,B,C,D,E,F,G,H>::Type_t thing = a(b,c,d,e,f,g,h);
 * </PRE>
 */

template<class Thing, class Sub1, class Sub2, class Sub3, class Sub4,
  class Sub5, class Sub6, class Sub7>
struct View7
{
};


/**
 * ComponentView<Components, Object>::Type_t is defined for several POOMA objects.
 * It tells you the type of a.comp(b).  You should be able to write code like:
 *
 * <PRE>
 * A a;
 * typename ComponentView<Loc<1>, A>::Type_t c = a.comp(i);
 * typename ComponentView<Loc<2>, A>::Type_t c = a.comp(i, j);
 * </PRE>
 *
 * The number gives the number of components you are applying.
 *
 * To define the view properties for a new class, you should specialize
 * ComponentView for that class.
 */

template<class Components, class Object>
struct ComponentView;

#endif     // POOMA_POOMA_VIEW_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: View.h,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:04 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
