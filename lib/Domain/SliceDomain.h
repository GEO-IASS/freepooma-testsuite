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

#ifndef POOMA_DOMAIN_SLICE_DOMAIN_H
#define POOMA_DOMAIN_SLICE_DOMAIN_H

//-----------------------------------------------------------------------------
// Class:
// SliceDomain<DT>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * SliceDomain is a base class for all sliced domain objects.
 *
 * A sliced domain stores two pieces of information:
 *   -# A "total domain" of dimension TotalDim
 *   -# A "slice domain" of dimension SliceDim, with SliceDim < TotalDim.
 *
 * SliceDomain stores both domains, and provides accessors to get references
 * to them.  It does not have the full interface as regular domains, you
 * must get a reference to the relevant domain (total or slice) and then
 * use that as normal.  Unlike the regular Domain class, SliceDomain does
 * not have or need any 1D specializations, or any base class.
 * SliceDomain is templated on the DomainTraits<> class providing the
 * traits for the SliceDomain derived class.  It does not have a Dim
 * template parameter as Domain does.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainBase.h"
#include "Domain/DomainTraits.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"
#include <iosfwd>


//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------


/**
 * SliceDomain<DT> provides the bulk of the public interface for all 
 * sliced domain objects.  The template parameter should be DomainTraits<DT>,
 * where DT is whatever SliceDomain subclass is being set up (e.g.,
 * SliceInterval<D,SD>, etc).  The DomainTraits for SliceDomain subclasses
 * has a smaller set of requirements than for Domain subclasses, since
 * SliceDomain has a much more restricted set of capabilities and requirements.
 *
 * A SliceDomain is a special domain class that stores information on a
 * "Total" domain of dimension 'dimensions', and a "Slice" domain of dimension
 * 'sliceDimensions'.  The slice domain is a subset of the total domain, with
 * one or more dimensions of the total domain removed (sliced out).  This
 * arises when users want to select a subset of an Array or something, by
 * specifying a set of domains of different types - the "singleValued" domains
 * such as Loc or int should refer to sliced dimensions, while the other
 * domains should refer to regions of the domain space as normal.
 * SliceDomain stores the information on what the "total" domain is (the
 * simple combination of all the domain objects as normally occurs when
 * constructing new Domain's) and the "slice" domain (the total domain, with
 * the sliced dimensions taken out).
 *
 * SliceDomain is only meant for use in constructing the internal objects
 * which refer to a subset of another object, such as an Array which is a
 * slice of a larger-dimensional Array.  Thus, the user should not use
 * SliceDomain directly, and it has a much smaller interface.  SliceDomain<DT>
 * defines the following interface:
 * - typedef DT::Domain_t      Domain_t ........ the subclass of SliceDomain
 * - typedef DT::SliceDomain_t SliceDomain_t ... the slice domain
 * - typedef DT::TotalDomain_t TotalDomain_t ... the total domain
 * - static const int dimensions ............... the total # of dimensions
 * - static const int sliceDimensions .......... the # of dims in SliceDomain_t
 * - const SliceDomain_t &sliceDomain() const;
 * -       SliceDomain_t &sliceDomain();
 * - const TotalDomain_t &totalDomain() const;
 * -       TotalDomain_t &totalDomain();
 */

template<class DT>
class SliceDomain
{
public:
  //
  // Typedefs obtained from the traits template parameter.
  //

  typedef typename DT::Domain_t      Domain_t;
  typedef typename DT::SliceDomain_t SliceDomain_t;
  typedef typename DT::TotalDomain_t TotalDomain_t;
  
  //
  // Constructors.  SliceDomain has a default constructor, and the
  // domain objects stored here will be uninitialized.  Since slice domains
  // are only intended to be constructed right before being filled, this
  // should be OK.
  //

  SliceDomain()
    : slice_m(Pooma::NoInit()),
      domain_m(Pooma::NoInit()) {
    CTAssert(DT::sliceDimensions <= DT::dimensions);
    for (int d = 0; d < DT::dimensions; ++d)
      ignore_m[d] = true;
  }

  // skip initialization constructor ... if a Pooma::NoInit object is provided,
  // skip initialization.  This is the same as the default constructor, but is
  // provided here to keep the same interface as other domains.
  SliceDomain(const Pooma::NoInit &e)
    : slice_m(e), domain_m(e) {
    CTAssert(DT::sliceDimensions <= DT::dimensions);
    for (int d = 0; d < DT::dimensions; ++d)
      ignore_m[d] = true;
  }

  // copy constructor
  SliceDomain(const SliceDomain<DT> &sd)
    : slice_m(sd.slice_m),
      domain_m(sd.domain_m) {
    CTAssert(DT::sliceDimensions <= DT::dimensions);
    for (int d = 0; d < DT::dimensions; ++d)
      ignore_m[d] = sd.ignore_m[d];
  }

  // copy constructor from other slice domain object
  // (for SliceInterval -> SliceRange conversion.)
  template <class DTO>
  SliceDomain(const SliceDomain<DTO> &sd)
    : slice_m(sd.sliceDomain()),
      domain_m(sd.totalDomain()) {
    CTAssert(DT::sliceDimensions <= DT::dimensions);
    for (int d = 0; d < DT::dimensions; ++d)
      ignore_m[d] = sd.ignorable(d);
  }

  // unwrap this object back to its derived-class type
  Domain_t &unwrap() { return *static_cast<Domain_t *>(this); }

  // a const version of unwrap
  const Domain_t &unwrap() const {
    return *static_cast<const Domain_t *>(this);
  }

  //
  // Destructor.  Here, nothing to do
  //

  ~SliceDomain() { }

  //
  // SliceDomain accessors/modifiers.
  //

  // return a reference to the slice domain, which will be a subset
  // of the total domain's set of 1D domain objects.
  const SliceDomain_t &sliceDomain() const { return slice_m; }
  SliceDomain_t &sliceDomain() { return slice_m; }

  // return a reference to the 'full' domain
  const TotalDomain_t &totalDomain() const { return domain_m; }
  TotalDomain_t &totalDomain() { return domain_m; }
  
  // indicate that the given dimension in the total domain is not
  // ignorable (i.e., has not been sliced out)
  void cantIgnoreDomain(int d)
  {
    PAssert(d >= 0 && d < DT::dimensions);
    ignore_m[d] = false;
  }
  
  // handle to ignore_m so it can be set
  bool &ignorable(int d)
  {
    PAssert(d >= 0 && d < DT::dimensions);
    return ignore_m[d];
  }
  
  // return true if the given dimension in the total domain is ignorable
  // (i.e., has been sliced out)
  bool ignorable(int d) const
  {
    PAssert(d >= 0 && d < DT::dimensions);
    return ignore_m[d];
  }

  //
  // operator= ... just do memberwise copy
  //

  SliceDomain<DT> &operator=(const SliceDomain<DT> &sd) {
    slice_m  = sd.slice_m;
    domain_m = sd.domain_m;
    for (int d = 0; d < DT::dimensions; ++d)
      ignore_m[d] = sd.ignore_m[d];
    return *this;
  }
  
  // set the slice domain based on the state of the total domain.
  void setSliceFromTotal() {
    for (int d = 0, dt = 0; d < DT::dimensions; ++d)
      if (!ignore_m[d])
        slice_m[dt++] = domain_m[d];
  }

  //
  // I/O
  //

  // print a SliceDomain to a stream, in the format
  //   "[" first:last:stride, first:last:stride, ... first:last:stride "]"
  //   "==>" (same for slice dimensions)

  template<class Out>
  void print(Out &o) const {
    o << totalDomain() << "==>" << sliceDomain();
  }

private:
  // the slice domain
  SliceDomain_t slice_m;

  // the full domain
  TotalDomain_t domain_m;
  
  // bool array telling whether or not a particular slice in the
  // full domain is ignorable.
  bool ignore_m[DT::dimensions];
};

#define POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP,RETTYPE,LHSTYPE)	\
template <class T>							\
inline typename T::RETTYPE						\
OPFUNCTION(const SliceDomain<T> &d1, LHSTYPE d2) {			\
  typename T::RETTYPE ret(d1.unwrap());					\
  ret.totalDomain() OP##= d2;                                           \
  ret.setSliceFromTotal();                                              \
  return ret;								\
}									\
template <class T>							\
inline typename T::RETTYPE						\
OPFUNCTION(LHSTYPE d1, const SliceDomain<T> &d2) {			\
  typename T::RETTYPE ret(d2.unwrap());					\
  ret.totalDomain() = d2 OP ret.totalDomain();                          \
  ret.setSliceFromTotal();                                              \
  return ret;								\
}

#define POOMA_SLICEDOMAIN_ARITH(OPFUNCTION, OP, RETTYPE)		\
template <class T1, class T2>						\
inline typename T1::RETTYPE						\
OPFUNCTION(const SliceDomain<T1> &d1, const DomainBase<T2> &d2) {	\
  typename T1::RETTYPE ret(d1.unwrap());				\
  ret.totalDomain() OP##= d2.unwrap();                                  \
  ret.setSliceFromTotal();                                              \
  return ret;								\
}									\
template <class T1, class T2>						\
inline typename T1::RETTYPE						\
OPFUNCTION(const SliceDomain<T1> &d1, const SliceDomain<T2> &d2) {	\
  typename T1::RETTYPE ret(d1.unwrap());				\
  ret.totalDomain() OP##= d2.unwrap().totalDomain();                    \
  ret.setSliceFromTotal();                                              \
  return ret;								\
}									\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, char)		\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, unsigned char)	\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, short)		\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, unsigned short)	\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, int)		\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, unsigned int)	\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, long)		\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, unsigned long)	\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, float)		\
POOMA_SLICEDOMAIN_ARITH_SINGLE(OPFUNCTION,OP, RETTYPE, double)

//
// define the arithmetic operators
//

POOMA_SLICEDOMAIN_ARITH(operator+, +, Domain_t)
POOMA_SLICEDOMAIN_ARITH(operator-, -, Domain_t)
POOMA_SLICEDOMAIN_ARITH(operator*, *, Domain_t)
POOMA_SLICEDOMAIN_ARITH(operator/, /, Domain_t)


/// print a SliceDomain to a stream, in the format
///   "[" first:last:stride, first:last:stride, ... first:last:stride "]"
///   "==>" (same for slice dimensions)

template<class DT>
std::ostream& operator<<(std::ostream &o, const SliceDomain<DT> &dbase) {
  dbase.print(o);
  return o;
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_SLICE_DOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: SliceDomain.h,v $   $Author: richard $
// $Revision: 1.20 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
