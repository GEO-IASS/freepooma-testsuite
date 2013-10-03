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
// Domain<DT,Dim>
//-----------------------------------------------------------------------------

#ifndef POOMA_DOMAIN_DOMAIN_H
#define POOMA_DOMAIN_DOMAIN_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * Domain is a base class for all domain objects, but one which can be
 * specialized for N-dimensional (N>1), and 1-dimensional domain objects.
 * The first template parameter is a dimension, which is used 
 * to specialize this class to 1-D objects.
 * Its second template parameter should be a traits class that describes
 * all the characteristics of the domain object, and the dimension of the
 * object.  
 * This base class provides the implementation for most of the public
 * interface (other than constructors) for the domain objects.  A few functions
 * which are common to all Domain objects, regardless of whether they are
 * 1-D or N-D domains, are collected into the DomainBase class which is a
 * base class for Domain.
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

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * Domain<Dim,DT> provides the bulk of the public interface for all domain
 * objects. The first template parameter is a dimension that
 * is used to specialize this class to 1-D objects.
 * The second template parameter should be DomainTraits<DT>, where
 * DT is whatever Domain object is being set up (Loc<N>, Interval<N>, etc).
 * DomainTraits is a traits class which must be specialized for the particular
 * type of domain this Domain base class is setting up.  DomainTraits
 * must include a number of typedefs, static data, and static methods to 
 * specialize Domain to do the right thing for the different Domain objects.
 *
 * When a Domain is created, it will initialize its storage if the default
 * constructor is used.  However, if you wish to avoid the work of
 * initialization, you can use the constructor which takes a Pooma::NoInit
 * object.  In that case, storage space for the domain will be maintained,
 * but it will not be initialized.  This is useful if you know you will be
 * changing the values later, and do not want to spend the extra time filling
 * in zeros or something into the storage.
 *
 * There is a general version of Domain, for an dimension N, and one partially
 * specialized to the case where the dimension == 1.  The 1D specialization
 * adds extra functionality for 1D domain objects, which is not available
 * for multidimensional objects.  For example, Range<2> objects do not have
 * first(), last(), length(), etc. methods, but Range<1> objects do.
 *
 * Domain inherits from DomainBase<DT>, regardless of the dimension;
 * DomainBase provides the definition of all methods which do not depend
 * on what the dimension is.
 *
 * Domain<N,DT> defines the following operations for domain objects:
 *   - operator[]
 *   - operator+=, -=, *=, /=
 *   - operator++, operator-- : just like using += stride(), -= stride()
 *     for each dimension
 *   - int size()
 *   - bool empty()
 *   - operator <, == (other comparisons defined in DomainBase.h)
 *
 * Note that if some operation is not defined for a particular type of
 * domain (e.g., *=, which is not legal for Loc or Interval objects), the
 * DomainTraits class will be missing the particular method needed to
 * implement that operation, and a compile-time error will occur if the
 * user tries to employ that disallowed operation.  The error will be
 * extremely cryptic, of course, but it will occur.
 *
 * Domain<1,DT> defines, in addition to the Domain<N,DT> methods, these
 * extra interface functions for just 1D domains:
 *   - int first(), last(), stride(), min(), max()
 *   - long length()
 *
 * The actual data for the domain (such as what the endpoints are) is
 * kept in DomainBase.  The DomainTraits class defines what the type
 * of the storage should be.  By putting that storage in the base class,
 * we avoid the 'empty base class' penalty of many compilers.  The storage
 * object is named domain_m.
 */

template<int Dim, class DT>
class Domain : public DomainBase<DT>
{
  /// convenience typedef
  typedef DomainBase<DT>                 Base_t;

public:
  //
  /// Typedefs obtained from the DomainBase and DomainTraits.
  //

  typedef typename DT::Size_t            Size_t;
  typedef typename Base_t::Domain_t      Domain_t;
  typedef typename DT::OneDomain_t       OneDomain_t;

  /// Iterator typedefs.  An N-dimensional iterator is a forward
  /// iterator, it works only with operator++ (it does not have operator--).

  typedef typename Base_t::const_iterator const_iterator;
  typedef typename Base_t::iterator       iterator;

  /// Block iterator typedefs.  All domains use the DomainBlockIterator
  /// class to iterate through blocks defined by the domain points.  A
  /// block iterator is a forward iterator, it works only with operator++.

  typedef typename Base_t::const_blockIterator const_blockIterator;
  typedef typename Base_t::blockIterator       blockIterator;


  /// Domain has a default constructor, which only
  /// makes sure for now that the Dim parameter is consistent with the
  /// DT parameter.  The storage object is kept and initialized in DomainBase.

  inline
  Domain() {
    CTAssert(DT::dimensions == Dim && Dim > 0);
  }

  /// If an Pooma::NoInit object is given in the constructor, we skip
  /// initialization of our array of 1D domains.

  inline Domain(const Pooma::NoInit &d) : Base_t(d) {
    CTAssert(DT::dimensions == Dim && Dim > 0);
  }

  //
  // Destructor.  Here, nothing to do.
  //

  inline
  ~Domain() { }

  //
  // Domain accessors.
  //

  /// return the Nth element of this domain, using []
  inline
  const OneDomain_t &operator[](int d) const { return this->domain_m[d]; }

  /// return the Nth element of this domain, using []
  inline
  OneDomain_t &operator[](int d) { return this->domain_m[d]; }

  /// return the total size of the domain, which is the product
  /// of all the lengths of the 1D domains
  inline
  Size_t size() const {
    Size_t sz = this->domain_m[0].size();
    for (int i = 1; i < Dim; i++)
      sz *= this->domain_m[i].size();
    return sz;
  }

  /// return if this domain is empty, which reports whether any of the
  /// N 1-dimensional domains are empty.
  inline
  bool empty() const {
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i].empty())
        return true;
    return false;
  }

  /// return whether this domain has been initialized.  This is the
  /// same as saying it is not empty.
  inline
  bool initialized() const { return (!empty()); }


  /// @name Comparison operators ==, !=, <, >, <=, >=
  //@{

  template<class T>
  bool operator==(const T &d2) const {
    CTAssert(Dim == DomainTraits<T>::dimensions);
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i] != DomainTraits<T>::getDomain(d2, i))
        return false;
    return true;
  }

  template<class T>
  bool operator<(const T &d2) const {
    CTAssert(Dim == DomainTraits<T>::dimensions);
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i] >= DomainTraits<T>::getDomain(d2, i))
        return false;
    return true;
  }

  template<class T>
  bool operator!=(const T &d2) const {
    CTAssert(Dim == DomainTraits<T>::dimensions);
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i] != DomainTraits<T>::getDomain(d2, i))
        return true;
    return false;
  }

  template<class T>
  bool operator>(const T &d2) const {
    CTAssert(Dim == DomainTraits<T>::dimensions);
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i] <= DomainTraits<T>::getDomain(d2, i))
        return false;
    return true;
  }

  template<class T>
  bool operator<=(const T &d2) const {
    CTAssert(Dim == DomainTraits<T>::dimensions);
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i] > DomainTraits<T>::getDomain(d2, i))
        return false;
    return true;
  }

  template<class T>
  bool operator>=(const T &d2) const {
    CTAssert(Dim == DomainTraits<T>::dimensions);
    for (int i = 0; i < Dim; i++)
      if (this->domain_m[i] < DomainTraits<T>::getDomain(d2, i))
        return false;
    return true;
  }

  //@}

  /// @name Arithmetic accumulation operators
  /// These are only allowed to
  /// occur with domain objects which are single-valued and have the
  /// right number of dimensions (basically, Loc's and scalar's).
  ///
  /// All return a reference to this object, but cast down to the
  /// derived type (e.g., Loc<N> instead of Domain<DomainTraits<Loc<N>>>

  //@{

  template<class T>
  Domain_t &operator+=(const T &d2) 
  {
    CTAssert(DomainTraits<T>::dimensions == Dim ||
	     DomainTraits<T>::dimensions == 1 ||
	     Dim == 1);
    int d = DomainTraits<T>::dimensions > Dim ? 
      DomainTraits<T>::dimensions : Dim;
    for (int i = 0; i < d; i++)
      this->domain_m[i] += DomainTraits<T>::getPointDomain(d2, i);
    return this->unwrap();
  }

  template<class T>
  Domain_t &operator-=(const T &d2) {
    CTAssert(DomainTraits<T>::singleValued);
    CTAssert(DomainTraits<T>::dimensions == Dim ||
	     DomainTraits<T>::dimensions == 1);
    for (int i = 0; i < Dim; i++)
      this->domain_m[i] -= DomainTraits<T>::getPointDomain(d2, i);
    return this->unwrap();
  }

  template<class T>
  Domain_t &operator*=(const T &d2) {
    CTAssert(DomainTraits<T>::singleValued);
    CTAssert(DomainTraits<T>::dimensions == Dim ||
	     DomainTraits<T>::dimensions == 1);
    for (int i = 0; i < Dim; i++)
      this->domain_m[i] *= DomainTraits<T>::getPointDomain(d2, i);
    return this->unwrap();
  }

  template<class T>
  Domain_t &operator/=(const T &d2) {
    CTAssert(DomainTraits<T>::singleValued);
    CTAssert(DomainTraits<T>::dimensions == Dim ||
	     DomainTraits<T>::dimensions == 1);
    for (int i = 0; i < Dim; i++)
      this->domain_m[i] /= DomainTraits<T>::getPointDomain(d2, i);
    return this->unwrap();
  }

  //@}

private:
  // make the copy constructor and operator= private and undefined
  // so that they will not be generated and so it will be an error if
  // the user tries to use them.  The classes derived from Domain should
  // provide all the constructors and operator='s needed.
  Domain(const Domain<Dim,DT> &);
  void operator=(const Domain<Dim,DT> &);
};


/**
 * SetDomainFunctor is a simple wrapper around the setDomain method in
 * the DomainTraits class.  It is templated on the DomainTraits type,
 * the domain storage type, the type of domain being copied into the
 * domain, and a boolean indicating if the domain is a wildcard type or not.
 * If it is a wildcard, the set operation is skipped, which can be used
 * to save time.  A specialization for wildcard == true is provided which
 * just does nothing, instead of calling DT::setDomain.
 *
 * When a wildcard is to be used to determine the final domain, a separate
 * setWildcardDomain method is available which takes an extra user-supplied
 * reference domain.  This reference domain is used by the wildcard to
 * calculate what the true domain should be.
 */

template<class DT, class ST, class T, class UT, bool wildcard>
struct SetDomainFunctor
{
  inline
  static void setDomain(ST &domain, const T &newdom) {
    DT::setDomain(domain, newdom);
  }
  inline
  static void setWildcardDomain(ST &domain, const UT &, const T &newdom) {
    DT::setDomain(domain, newdom);
  }
};

template<class DT, class ST, class T, class UT>
struct SetDomainFunctor<DT, ST, T, UT, true>
{
  inline
  static void setDomain(ST &, const T &) { }
  inline
  static void setWildcardDomain(ST &domain, const UT &u, const T &newdom) {
    DT::setWildcardDomain(domain, u, newdom);
  }
};


/**
 * The 1D-specialized version of Domain, which acts much like the ND version
 * but also provides a number of new or redefined interface functions:
 *  - int first(), last(), stride(), min(), max()
 *  - long length()
 */

template<class DT>
class Domain<1, DT> : public DomainBase<DT>
{
  // convenience typedef
  typedef DomainBase<DT>                 Base_t;

public:
  //
  // Typedefs obtained from the DomainBase and DomainTraits
  //

  typedef typename DT::Size_t            Size_t;
  typedef typename DT::Element_t         Element_t;
  typedef typename Base_t::Domain_t      Domain_t;
  typedef typename Base_t::Storage_t     Storage_t;

  // Iterator typedefs.  An N-dimensional iterator is a forward
  // iterator, it works only with operator++ (it does not have operator--).

  typedef typename Base_t::const_iterator const_iterator;
  typedef typename Base_t::iterator       iterator;

  // Block iterator typedefs.  All domains use the DomainBlockIterator
  // class to iterate through blocks defined by the domain points.  A
  // block iterator is a forward iterator, it works only with operator++.

  typedef typename Base_t::const_blockIterator const_blockIterator;
  typedef typename Base_t::blockIterator       blockIterator;


  //
  // Constructors.
  //

  inline
  Domain() { }

  inline
  Domain(const Pooma::NoInit &d) : Base_t(d) { }

  //
  // Destructor.  Here, nothing to do.
  //

  inline
  ~Domain() { }

  //
  // Domain accessors.
  //

  // return the Nth element of this domain, using [].  Since this is a
  // 1D object, this  just returns this same object back using unwrap()
  inline
  Domain_t &operator[](int) { return this->unwrap(); }
  inline
  const Domain_t &operator[](int) const { return this->unwrap(); }

  // return the ith value of the domain.
  inline
  Element_t elem(int n) const { return DT::elem(this->domain_m, n); }
  inline
  Element_t operator()(int n) const { return DT::elem(this->domain_m, n); }

  // return the first and last points in the domain (the endpoints)
  inline
  Element_t first() const { return DT::first(this->domain_m); }
  inline
  Element_t last() const { return DT::last(this->domain_m); }

  // return the stride of the domain
  inline
  Element_t stride() const { return DT::stride(this->domain_m); }
  
  // return the length of the domain, which is the number of points
  // (including the endpoints) for the domain
  inline
  Size_t length() const { return DT::length(this->domain_m); }

  // return the min and max values of the domain endpoints
  inline
  Element_t min() const { return DT::min(this->domain_m); }
  inline
  Element_t max() const { return DT::max(this->domain_m); }

  // return the total size of the domain, which is the product
  // of all the lengths of the 1D domains
  inline
  Size_t size() const { return length(); }

  // return if this domain is empty.
  inline
  bool empty() const { return DT::empty(this->domain_m); }

  // return whether this domain has been initialized.  This is the
  // same as saying it is not empty.
  inline
  bool initialized() const { return (!empty()); }

  // return which loop this domain corresponds to (not all domains have this
  // kind of information, but for those cases default values will be used)
  inline
  int loop() const { return DT::loop(this->domain_m); }

  //
  // Domain modifiers.
  //

  // setDomain: for a 1D domain, this actually tries to change the current
  // domain settings to those of the given 1D domain.  If the given object
  // is not 1D, or if it does not have information that we require, it is
  // a compile-time error.
  template<class T>
  inline
  void setDomain(const T &newdom) {
    SetDomainFunctor<DT,Storage_t,T,T,DomainTraits<T>::wildcard>::
      setDomain(this->domain_m, newdom);
  }

  // setWildcardDomain: the same as setDomain, except that if the new domain
  // object is a wildcard, this will use a user-supplied reference 1D domain
  // to calculate what the proper domain should be.
  template<class UT, class T>
  inline
  void setWildcardDomain(const UT &u, const T &newdom) {
    SetDomainFunctor<DT,Storage_t,T,UT,DomainTraits<T>::wildcard>::
      setWildcardDomain(this->domain_m, u, newdom);
  }

  // setLoop: change which loop variable this dimension should correspond
  // to.  Some domain object may just ignore this information.
  inline
  void setLoop(int newloop) { DT::setLoop(this->domain_m, newloop); }

  //
  // Main comparison operators == and <
  //

  // equality comparison for a 1D domain.
  template<class T>
  bool operator==(const T &d2) const { 
    return DT::isEqualTo(this->domain_m, d2); 
  }

  // less-than comparison for a 1D domain.
  template<class T>
  bool operator<(const T &d2) const { 
    return DT::isLessThan(this->domain_m, d2); 
  }

  //
  // other comparison operators, built using == and <
  //

  template<class T>
  bool operator!=(const T &d2) const { return !(*this == d2); }

  template<class T>
  bool operator>(const T &d2) const { return !(*this < d2 || *this == d2); }

  template<class T>
  bool operator<=(const T &d2) const { return (*this < d2 || *this == d2); }

  template<class T>
  bool operator>=(const T &d2) const { return !(*this < d2); }

  //
  // Arithmetic accumulation operators.  These are only allowed to
  // occur with domain objects which are single-valued and have the
  // right number of dimensions (basically, Loc's and int's).
  //
  // All return a reference to this object, but cast down to the
  // derived type (e.g., Loc<N> instead of Domain<DomainTraits<Loc<N>>>
  //

  template<class T>
  inline
  Domain_t &operator+=(const T &d2) { 
    DT::addAccum(this->domain_m, d2);
    return this->unwrap();
  }

  template<class T>
  inline
  Domain_t &operator-=(const T &d2) {
    DT::subtractAccum(this->domain_m, d2);
    return this->unwrap();
  }

  template<class T>
  Domain_t &operator*=(const T &d2) {
    DT::multiplyAccum(this->domain_m, d2);
    return this->unwrap();
  }

  template<class T>
  Domain_t &operator/=(const T &d2) {
    DT::divideAccum(this->domain_m, d2);
    return this->unwrap();
  }

private:
  // make the copy constructor and operator= private and undefined
  // so that they will not be generated and so it will be an error if
  // the user tries to use them.  The classes derived from Domain should
  // provide all the constructors and operator='s needed.
  Domain(const Domain<1,DT> &);
  void operator=(const Domain<1,DT> &);
};


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Domain.h,v $   $Author: richard $
// $Revision: 1.34 $   $Date: 2004/11/01 18:16:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
