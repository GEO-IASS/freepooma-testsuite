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
//   CombineDomainOpt
//-----------------------------------------------------------------------------

#ifndef POOMA_DOMAIN_COMBINEDOMAINOPT_H
#define POOMA_DOMAIN_COMBINEDOMAINOPT_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * CombineDomainOpt is a class that can be used to optimize the operation
 * NewDomainN<>::combineSlice(domain, s1, s2, ...).
 *
 * CombineDomainOpt is a class that can be used to optimize the operation
 * NewDomainN<>::combineSlice(domain, s1, s2, ...).
 *
 * Typically NewDomain is used by arrays to construct a view domain that
 * could be a slice, so typically you would call
 *
 * NewDomainN<>::combineSlice(a.totalDomain(), s1, s2);
 *
 * If the result is single-valued, the domain of a is not used, but the
 * function call a.totalDomain() may be hard to optimize away.  To avoid this
 * function call you can now say:
 *
 * typedef NewDomainN<...> NewDomain_t;
 * typedef typename NewDomain_t::SliceType_t SliceDomain_t;
 *
 * SliceDomain_t s(
 *     CombineDomainOpt<NewDomain_t, SliceDomain_t::singleValued>::
 *     make(a, s1, s2, ...)
 *                );
 *
 * If s is single-valued, the array a is never used.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

template<class NewDomain, bool sv>
struct CombineDomainOpt;

/**
 * Single valued specialization.
 * Since Locs can construct themselves from other Locs and
 * integers, we just use the constructors.  Another option is
 * return Loc<dim>(s1, s2, ...).
 */

template<class NewDomain>
struct CombineDomainOpt<NewDomain, true>
{
  typedef typename NewDomain::SliceType_t Type_t;

  template<class Array, class Sub1>
  inline static
  Type_t make(const Array &, const Sub1 &s1)
  {
    return Type_t(s1);
  }

  template<class Array, class Sub1, class Sub2>
  inline static
  Type_t make(const Array &, const Sub1 &s1, const Sub2 &s2)
  {
    return Type_t(s1, s2);
  }

  template<class Array, class Sub1, class Sub2, class Sub3>
  inline static
  Type_t make(const Array &,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3)
  {
    return Type_t(s1, s2, s3);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4>
  inline static
  Type_t make(const Array &,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4)
  {
    return Type_t(s1, s2, s3, s4);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4, class Sub5>
  inline static
  Type_t make(const Array &,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4, const Sub5 &s5)
  {
    return Type_t(s1, s2, s3, s4, s5);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4, class Sub5, class Sub6>
  inline static
  Type_t make(const Array &,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4, const Sub5 &s5, const Sub6 &s6)
  {
    return Type_t(s1, s2, s3, s4, s5, s6);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4, class Sub5, class Sub6, class Sub7>
  inline static
  Type_t make(const Array &,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4, const Sub5 &s5, const Sub6 &s6,
	      const Sub7 &s7)
  {
    return Type_t(s1, s2, s3, s4, s5, s6, s7);
  }
};

/**
 * Multi-valued version.  This one calls combineSlice to create
 * the final domain.
 */

template<class NewDomain>
struct CombineDomainOpt<NewDomain, false>
{
  typedef typename NewDomain::SliceType_t Type_t;

  template<class Array, class Sub1>
  inline static
  Type_t make(const Array &a, const Sub1 &s1)
  {
    return NewDomain::combineSlice(a.totalDomain(), s1);
  }

  template<class Array, class Sub1, class Sub2>
  inline static
  Type_t make(const Array &a, const Sub1 &s1, const Sub2 &s2)
  {
    return NewDomain::combineSlice(a.totalDomain(), s1, s2);
  }

  template<class Array, class Sub1, class Sub2, class Sub3>
  inline static
  Type_t make(const Array &a,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3)
  {
    return NewDomain::combineSlice(a.totalDomain(), s1, s2, s3);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4>
  inline static
  Type_t make(const Array &a,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4)
  {
    return NewDomain::combineSlice(a.totalDomain(), s1, s2, s3, s4);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4, class Sub5>
  inline static
  Type_t make(const Array &a,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4, const Sub5 &s5)
  {
    return NewDomain::combineSlice(a.totalDomain(), s1, s2, s3, s4, s5);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4, class Sub5, class Sub6>
  inline static
  Type_t make(const Array &a,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4, const Sub5 &s5, const Sub6 &s6)
  {
    return NewDomain::combineSlice(a.totalDomain(), 
      s1, s2, s3, s4, s5, s6);
  }

  template<class Array, class Sub1, class Sub2, class Sub3,
    class Sub4, class Sub5, class Sub6, class Sub7>
  inline static
  Type_t make(const Array &a,
	      const Sub1 &s1, const Sub2 &s2, const Sub3 &s3,
	      const Sub4 &s4, const Sub5 &s5, const Sub6 &s6,
	      const Sub7 &s7)
  {
    return NewDomain::combineSlice(a.totalDomain(), 
      s1, s2, s3, s4, s5, s6, s7);
  }
};

#endif     // POOMA_DOMAIN_COMBINEDOMAINOPT_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: CombineDomainOpt.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
