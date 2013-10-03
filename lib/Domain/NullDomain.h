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
//   NullDomain
//   ErrorDomain
// Functions:
//   contains(NullDomain, D)
//-----------------------------------------------------------------------------

#ifndef POOMA_DOMAIN_NULLDOMAIN_H
#define POOMA_DOMAIN_NULLDOMAIN_H

/** @file
 * @ingroup Domain
 * @brief
 * NullDomain and ErrorDomain special domains.
 */

//@{

/**
 * NullDomain and ErrorDomain are special "domains". ErrorDomains result when
 * someone tries an incorrect domain calculation.
 */

struct ErrorDomain
{
  ErrorDomain() {}
  ErrorDomain(const ErrorDomain&) {}
};

/**
 * NullDomain and ErrorDomain are special "domains".
 * NullDomain means a domain with nothing in it.
 */

struct NullDomain
{
  NullDomain() {}
  NullDomain(const NullDomain&) {}
};

//@}

/// The null-domain is contained by every domain.

template<class D>
bool contains(const NullDomain &, const D &)
{
  return true;
}

#endif // POOMA_DOMAIN_NULLDOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NullDomain.h,v $   $Author: richard $
// $Revision: 1.3 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
