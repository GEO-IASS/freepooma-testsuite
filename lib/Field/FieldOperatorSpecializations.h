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
// Function:
// min(Field<D,T,E>, Field<D,T,E>)
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_FIELDOPERATORSPECIALIZATIONS_H
#define POOMA_FIELD_FIELDOPERATORSPECIALIZATIONS_H

/** @file
 * @ingroup Field
 * @brief
 * min()/max() specializations
 *
 * The STL defines min() and max() functions with signatures
 * template<class T> T min(T,T);
 * If you generate a scope that includes both the std:: namespace,
 * and the namespace containing the Pooma versions of min() and max(),
 * then min(Field<1>,Field<1>) is ambiguous.  These specializations
 * disambiguate this common case.
 */

template<class D1,class T1,class E1>
inline typename MakeReturn<BinaryNode<FnMin,
  typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t,
  typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t> >::Expression_t
min(const Field<D1,T1,E1> & l,const Field<D1,T1,E1> & r)
{
  typedef BinaryNode<FnMin,
    typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t,
    typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t> Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(
    CreateLeaf<Field<D1,T1,E1> >::make(l),
    CreateLeaf<Field<D1,T1,E1> >::make(r)));
}

template<class D1,class T1,class E1>
inline typename MakeReturn<BinaryNode<FnMax,
  typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t,
  typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t> >::Expression_t
max(const Field<D1,T1,E1> & l,const Field<D1,T1,E1> & r)
{
  typedef BinaryNode<FnMax,
    typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t,
    typename CreateLeaf<Field<D1,T1,E1> >::Leaf_t> Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(
    CreateLeaf<Field<D1,T1,E1> >::make(l),
    CreateLeaf<Field<D1,T1,E1> >::make(r)));
}

#endif     // POOMA_FIELD_FIELDOPERATORSPECIALIZATIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldOperatorSpecializations.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:43 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
