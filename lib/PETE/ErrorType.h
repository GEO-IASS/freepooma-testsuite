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

#ifndef PETE_PETE_ERRORTYPE_H
#define PETE_PETE_ERRORTYPE_H

//-----------------------------------------------------------------------------
// ErrorType
//
// ErrorType is a special type that can be used in traits computations to
// signify an illegal or undefined result.  Currently ErrorType is only used
// in CreateLeaf, in other places (EngineTraits for example) types are just
// left undefined to generate an error at compile time.
// ErrorType is used here for some types that end up not being used.  We want
// to avoid some compile time type computations, but don't want to generate an
// error.  For example, we define several operators:
// operator+(Array<D,T,E>,Array<D2,T2,E2>)
// operator+(T,Array<D2,T2,E2>)
// operator+(Array<D,T,E>,T2)
// Suppose we add two expressions (a+b+c)+(d+e+f).  The compiler will compute
// the return types to each operator+ before it decides that Array+Array is the
// most specialized.  For the second and third versions the return type
// would perform some hairy template metaprograms on Scalar<Array<...>> which
// never actually get used.  (This problem caused the compiler to self-destruct
// for expressions with sums of 15 or so terms.)  To avoid this problem, there
// are CreateLeaf specializations that just return ErrorType.
//-----------------------------------------------------------------------------

struct ErrorType
{
};

#endif // PETE_PETE_ERRORTYPE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ErrorType.h,v $   $Author: richard $
// $Revision: 1.2 $   $Date: 2004/11/01 18:16:56 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
