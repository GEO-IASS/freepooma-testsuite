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
// Class:
// FunctorResult
//-----------------------------------------------------------------------------

#ifndef POOMA_POOMA_FUNCTORRESULT_H
#define POOMA_POOMA_FUNCTORRESULT_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Pooma
 * @brief
 * FunctorResult<Functor,T>::Type_t is the return type of a pooma functor
 * object.
 */


/**
 * Where Pooma uses functor objects, STL functors can typically be used,
 * where operator() is used to define the function.
 *
 * STL functors typically either give no type information about the argument
 * and return, or specify both through argument_type and result_type.
 * For Pooma, we prefer to use templated functors, so result_type is a
 * function of argument_type, which is why we use the class FunctorResult<>
 * to externally provide the mapping from argument_type to result_type.
 *
 * The default behaviour simply uses the argument_type for the result_type
 * which is sufficient for most STL function objects and most pooma function
 * objects.  A couple of examples where you would want to specialize this
 * class for different behaviour are:
 *
 * <PRE>
 * template<class T>
 * struct FunctorResult<logical_not<T>,T>
 * {
 *   typedef bool Type_t;
 * };
 *
 * template<int D, class T>
 * struct FunctorResult<MyNorm,Vector<D,T> >
 * {
 *   typedef T Type_t;
 * };
 * </PRE>
 */

template<class Functor, class ArgumentType>
struct FunctorResult
{
  typedef ArgumentType Type_t;
};

#endif     // POOMA_POOMA_FUNCTORRESULT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FunctorResult.h,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
