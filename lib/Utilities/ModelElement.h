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
//   ConstField
//-----------------------------------------------------------------------------

#ifndef POOMA_UTILITIES_MODELELEMENT_H
#define POOMA_UTILITIES_MODELELEMENT_H

/** @file
 * @ingroup Utilities
 * @brief
 * A wrapper class used to differentiate overloaded functions.
 */

template<class T>
class ModelElement
{
public:
  explicit ModelElement(const T &e) : e_m(e) { }
  ModelElement(const ModelElement<T> &m) : e_m(m.e_m) { }
  const T &element() const { return e_m; }
private:
  const T &e_m;
};

template<class T>
inline ModelElement<T> modelElement(const T &elem) 
  {
    return ModelElement<T>(elem);
  }

#endif // POOMA_UTILITIES_MODELELEMENT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ModelElement.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
