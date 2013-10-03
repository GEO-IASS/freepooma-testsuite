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

#ifndef POOMA_DOMAIN_DOMAIN_BASE_H
#define POOMA_DOMAIN_DOMAIN_BASE_H

//-----------------------------------------------------------------------------
// Class:
// DomainBase<T>
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file 
 * @ingroup Domain
 * @brief
 * DomainBase is a common base class for all domain objects.
 *
 * The template
 * parameter T should be a traits class that describes all the
 * characteristics of the domain object, and the dimension of the object.
 * This base class provides a collection of all the functionality that is
 * common to all DomainBase-derived objects, regardless of whether they are
 * specialized to a specific number of dimensions or not.  For example, both
 * Domain<N, DomainTraits<Loc<N>>> and Domain<1, DomainTraits<Loc<N>>>
 * use DomainBase as a base class.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainArithOpsTraits.h"
#include "Domain/DomainIterator.h"
#include "Domain/DomainBlockIterator.h"
#include "Domain/DomainTraits.h"
#include "Utilities/NoInit.h"
#include "Utilities/PAssert.h"
#include <iosfwd>

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

/**
 * DomainBase<DT> is the common base class for all Domain<Dim, DT> objects,
 * regardless of the value of Dim.  Since Domain<Dim,DT> is defined for
 * general Dim and also specialized to Dim==1, then it is useful to
 * collect common code for both cases into this base class.  It also
 * provides the standard typedefs and enumerations which all Domain's should
 * have, with values determined by looking at the DomainTraits class.
 * DT should be DomainTraits<Dom>, where Dom is a domain class like Loc<N>,
 * Range<N>, etc.  DomainTraits should be specialized for the different
 * domain types which are to be set up.
 *
 * The following typedefs and enums are defined here, basically just as
 * wrappers around the values from the DomainTraits class:
 *   - Domain_t : the domain type itself, on which the DomainTraits
 *     is templated, such as Loc<N> or Range<N>
 *   - OneDomain_t : a 1D version of the same domain type; if Domain_t is
 *     Range<N>, then OneDomain_t would be Range<1>
 *   - Element_t : the type of data used for the sequence values, mostly
 *     this is int
 *   - Storage_t : the type for actually storing the domain data itself.  A
 *     variable of type Storage_t is a member of DomainBase.
 *   - dimensions : a static const int, set to the dimensionality of the domain.
 *   - loopAware : if true, this domain stores info about which loop variable
 *     it refers to.  If false, this domain has no loop variable knowledge.
 *   - singleValued : if true, this domain is known at compile-time to refer to
 *     just one point, instead of a possible sequence.
 *   - unitStride : if true, this domain is known at compile-time to refer to a
 *     unit-strided sequence.
 *
 * This base class defines the following common interface methods for all
 * Domain objects regardless of their dimensionality:
 *   - unwrap() - cast this particular object down to the actual domain type;
 *     for example, DomainBase<DomainTraits<Loc<N>>> would be cast down
 *     to Loc<N>, and unwrap returns a reference to this cast to that type.
 *   - operator- : just like returning (*this * -1)
 *
 * This file also defines binary arithmetic operations +, -, *, /, and general
 * comparison operations <=, >=, !=, > (these use the < and == operators
 * defined in Domain<DT,Dim>, just like the +,-,*,/ operators used the +=,
 * -=, etc. operators in Domain<DT,Dim>).
 *
 * When a DomainBase is created, it will initialize its storage if the default
 * constructor is used.  However, if you wish to avoid the work of
 * initialization, you can use the constructor which takes a Pooma::NoInit
 * object.  In that case, storage space for the domain will be maintained,
 * but it will not be initialized.  This is useful if you know you will be
 * changing the values later, and do not want to spend the extra time filling
 * in zeros or something into the storage.
 *
 * Finally, this file defines the operations to print a Domain to an ostream;
 * the format for printing a domain is  "[" followed by first():last():stride()
 * for each dimension, followed by "]".  For example, a 2D Range with the
 * same sequence 1 ... 9 step 2 would be  "[1:9:2,1:9:2]"
 */

template<class DT>
class DomainBase
{
public:
  //
  // Typedefs obtained from the traits template parameter.
  //

  typedef typename DT::Domain_t      Domain_t;
  typedef typename DT::AskDomain_t   AskDomain_t;
  typedef typename DT::MultResult_t  MultResult_t;
  typedef typename DT::Storage_t     Storage_t;

  // Iterator typedefs.  An N-dimensional iterator is a forward
  // iterator, it works only with operator++ (it does not have operator--).

  typedef DomainIterator<Domain_t> const_iterator;
  typedef DomainIterator<Domain_t> iterator;

  // Block iterator typedefs.  All domains use the DomainBlockIterator
  // class to iterate through blocks defined by the domain points.  A
  // block iterator is a forward iterator, it works only with operator++.

  typedef DomainBlockIterator<Domain_t> const_blockIterator;
  typedef DomainBlockIterator<Domain_t> blockIterator;


  //
  // Constructors.  DomainBase has a default constructor, which
  // makes sure for now that the Dim parameter is consistent with the
  // DT parameter.  The domain traits class knows specifically how the
  // storage should be initialized.
  //

  inline
  DomainBase() {
    CTAssert(DT::dimensions > 0);
    DT::initializeStorage(domain_m);
  }

  //
  // If a Pooma::NoInit object is given in the constructor, we skip
  // initialization of our array of 1D domains.
  //

  inline
  DomainBase(const Pooma::NoInit &) {
    CTAssert(DT::dimensions > 0);
  }

  //
  // Destructor.  Here, nothing to do
  //

  inline
  ~DomainBase() { }

  //
  // DomainBase accessors.
  //

  // unwrap this object back to its derived-class type
  inline
  Domain_t &unwrap() { return *static_cast<Domain_t *>(this); }

  // a const version of unwrap
  const Domain_t &unwrap() const {
    return *static_cast<Domain_t *>(const_cast<DomainBase<DT> *>(this));
  }

  // return the first elements of the domain in another domain object
  inline
  AskDomain_t firsts() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].first();
    return retval;
  }

  // return the last elements of the domain in another domain object
  inline
  AskDomain_t lasts() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].last();
    return retval;
  }

  // return the stride of the domain in another domain object
  inline
  AskDomain_t strides() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].stride();
    return retval;
  }

  // return the lengths of the domain in another domain object
  inline
  AskDomain_t lengths() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].length();
    return retval;
  }

  // return the sizes of the 1D domains in another domain object
  inline
  AskDomain_t sizes() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].size();
    return retval;
  }

  // return the min values of the 1D domains in another domain object
  inline
  AskDomain_t mins() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].min();
    return retval;
  }

  // return the max values of the 1D domains in another domain object
  inline
  AskDomain_t maxes() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].max();
    return retval;
  }

  // return the loop values of the 1D domains in another domain object
  inline
  AskDomain_t loops() const {
    AskDomain_t retval = Pooma::NoInit();
    for (int i=0; i < DT::dimensions; ++i)
      retval[i] = unwrap()[i].loop();
    return retval;
  }

  //
  // Negation operator
  //

  MultResult_t operator-() const { 
    return (MultResult_t(unwrap()) * (-1)); 
  }

  //
  // Increment/decrement operators
  //

  Domain_t &operator++() {
    for (int i = 0; i < DT::dimensions; i++)
      domain_m[i] += domain_m[i].stride();
    return unwrap();
  }

  Domain_t &operator--() {
    for (int i = 0; i < DT::dimensions; i++)
      domain_m[i] -= domain_m[i].stride();
    return unwrap();
  }

  //
  // Iterator accessor functions for 1-D Domain.
  //

  // return begin and end iterators.  const_iterator and iterator
  // are the same here, so we only need one version.

  const_iterator begin() const { return const_iterator(unwrap()); }
  const_iterator end() const { return const_iterator(unwrap(),
						     unwrap().size()); }

  // return begin and end block iterators.  const_blockiterator and
  // blockIterator are the same here, so we only need one version.

  const_blockIterator beginBlock() const{return const_blockIterator(unwrap());}
  const_blockIterator endBlock() const { return const_blockIterator(); }

  //
  // I/O
  //

  // print a domain to a stream, in the format
  //   "[" first:last:stride, first:last:stride, ... first:last:stride "]"

  template<class Out>
  void print(Out &o) const {
    const Domain_t &d = unwrap();
    o << "[";
    for (int i=0; i < DT::dimensions; ++i) {
      o << d[i].first() << ":" << d[i].last() << ":" << d[i].stride();
      if (i < (DT::dimensions-1))
        o << ",";
    }
    o << "]";
  }

protected:
  // The storage for the domain data.  We put it here in the base
  // class so that the base class is not empty, and so that the
  // other methods implemented in the base class which manipulate the
  // data can actually see the storage.
  Storage_t domain_m;

private:
  // make the copy constructor and operator= private and undefined
  // so that they will not be generated and so it will be an error if
  // the user tries to use them.  The classes derived from DomainBase should
  // provide all the constructors and operator='s needed.
  DomainBase(const DomainBase<DT> &);
  void operator=(const DomainBase<DT> &);
};


//////////////////////////////////////////////////////////////////////
//
// Inline implementations of arithmetic operators for DomainBase objects.
// These build on the accumulation operators +=, -=, *=, etc. that each
// Domain must defined in a derived class, and the comparison operators
// ==, !=, <, etc..  They work with two domain objects, or a domain and
// a scalar.
//
// Note that for many operators, the Domain class defines templated
// functions which allow for the domain object on the LHS, and an
// arbitrary type on the RHS.  But this does not cover the case of having
// an arbitrary type on the LHS, and a domain object on the RHS.  The
// following global functions try to cover that possibility.  Unfortunately,
// due to template resolution ambiguities, we cannot have ALL possible types
// on the LHS, only  1) other domain types, and 2) basic scalars.
//
//////////////////////////////////////////////////////////////////////

//
// first define some helper macros, to get the functions defined
// for different scalar types on the LHS.  It would be great if
// we could just have a templated type for the LHS, but this leads
// to errors in 'multiple functions defined for a single operation'.
// So, we're required to enumerate the possibilities for the LHS of these
// operations.  Note that for a domain on the LHS of the operation,
// the RHS type IS templated, so you can have any type on the RHS.
// But if the scalar appears on the LHS, you can only have the types
// listed here
//

#define POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, LHSTYPE)	\
template<class T>						\
inline bool							\
OPFUNCTION(LHSTYPE d1, const DomainBase<T> &d2) {		\
  return (d2.unwrap() OP d1);					\
}


#define POOMA_DOMAIN_OPERATOR(OPFUNCTION, OP)			\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, char)		\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, unsigned char)	\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, short)		\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, unsigned short)	\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, int)		\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, unsigned int)	\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, long)		\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, unsigned long)	\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, float)		\
POOMA_DOMAIN_OPERATOR_SINGLE(OPFUNCTION, OP, double)


#define POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, LHSTYPE)	\
template<class T>							\
inline typename T::RETTYPE						\
OPFUNCTION(const DomainBase<T> &d1, LHSTYPE d2) {			\
  typename T::RETTYPE retval(d1.unwrap());				\
  return (retval OP d2);						\
}									\
template<class T>							\
inline typename T::RETTYPE						\
OPFUNCTION(LHSTYPE d1, const DomainBase<T> &d2) {			\
  typename T::RETTYPE retval(d2.unwrap());				\
  return (retval OP d1);						\
}

#define POOMA_DOMAIN_ARITH(OPFUNCTION, OP, RETTYPE)                     \
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, char)		\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, unsigned char)	\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, short)		\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, unsigned short)	\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, int)		        \
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, unsigned int)	\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, long)		\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, unsigned long)	\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, float)		\
POOMA_DOMAIN_ARITH_SINGLE(OPFUNCTION, OP, RETTYPE, double)


template <class T1, class T2, bool SV>
struct DomPair {

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::AddResult_t
  add(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {									
    typename DomainArithOpsTraits<typename T1::Domain_t,
                                  typename T2::Domain_t>::AddResult_t   
      retval(d1.unwrap());   
    return (retval += d2.unwrap());					
  }

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::SubResult_t
  sub(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {
    typename DomainArithOpsTraits<typename T1::Domain_t,            	
                                  typename T2::Domain_t>::SubResult_t   
      retval(d1.unwrap());   
    return (retval -= d2.unwrap());					
  }

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::MultResult_t
  mult(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {									
    typename DomainArithOpsTraits<typename T1::Domain_t,            	
                                  typename T2::Domain_t>::MultResult_t   
      retval(d1.unwrap());   
    return (retval *= d2.unwrap());					
  }

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::MultResult_t
  div(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {									
    typename DomainArithOpsTraits<typename T1::Domain_t,            	
                                  typename T2::Domain_t>::MultResult_t   
      retval(d1.unwrap());   
    return (retval /= d2.unwrap());					
  }
};


template <class T1, class T2>
struct DomPair<T1,T2,false> {

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::AddResult_t
  add(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {                                                                	
    typename DomainArithOpsTraits<typename T1::Domain_t,           
                                  typename T2::Domain_t>::AddResult_t
      retval(d2.unwrap());     
    return (retval += d1.unwrap());					
  }

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::SubResult_t
  sub(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {                                                                	
    typename DomainArithOpsTraits<typename T1::Domain_t,           
                                  typename T2::Domain_t>::SubResult_t
      retval(d2.unwrap());
    retval = -retval;

    return (retval += d1.unwrap());					
  }

  static typename DomainArithOpsTraits<typename T1::Domain_t,
                                       typename T2::Domain_t>::MultResult_t
  mult(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
  {                                                                	
    typename DomainArithOpsTraits<typename T1::Domain_t,           
                                  typename T2::Domain_t>::MultResult_t
      retval(d2.unwrap());     
    return (retval *= d1.unwrap());					
  }

// It makes no sense to divide by a non singlevalued domain. 
// static typename DomainArithOpsTraits<typename T1::Domain_t,          
//                                      typename T2::Domain_t>::MultResult_t
// div(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
// {                                                                	
//   typename DomainArithOpsTraits<typename T1::Domain_t,           
//                                 typename T2::Domain_t>::MultResult_t
//     retval(1/d2.unwrap());  
//   What does it mean to say 1/Range<2> ?????!?!?! 
//   return (retval *= d1.unwrap());				       
// }
};





template<class T1, class T2>						
inline typename DomainArithOpsTraits<typename T1::Domain_t,          
                                     typename T2::Domain_t>::AddResult_t
operator+(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
{
  return DomPair<T1,T2,T2::singleValued>::add(d1,d2);
}

template<class T1, class T2>						
inline typename DomainArithOpsTraits<typename T1::Domain_t,          
                                     typename T2::Domain_t>::SubResult_t
operator-(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
{
  return DomPair<T1,T2,T2::singleValued>::sub(d1,d2);
}

template<class T1, class T2>						
inline typename DomainArithOpsTraits<typename T1::Domain_t,          
                                     typename T2::Domain_t>::MultResult_t
operator*(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
{
  return DomPair<T1,T2,T2::singleValued>::mult(d1,d2);
}

template<class T1, class T2>						
inline typename DomainArithOpsTraits<typename T1::Domain_t,          
                                     typename T2::Domain_t>::MultResult_t
operator/(const DomainBase<T1> &d1, const DomainBase<T2> &d2)
{
  return DomPair<T1,T2,T2::singleValued>::div(d1,d2);
}



//
// define the comparison operators
//

POOMA_DOMAIN_OPERATOR(operator==, ==)
POOMA_DOMAIN_OPERATOR(operator!=, !=)
POOMA_DOMAIN_OPERATOR(operator<,  <)
POOMA_DOMAIN_OPERATOR(operator>,  >)
POOMA_DOMAIN_OPERATOR(operator<=, <=)
POOMA_DOMAIN_OPERATOR(operator>=, >=)

//
// define the scalar arithmetic operators
//

POOMA_DOMAIN_ARITH(operator+, +=, AddResult_t)
POOMA_DOMAIN_ARITH(operator-, -=, AddResult_t)
POOMA_DOMAIN_ARITH(operator*, *=, MultResult_t)
POOMA_DOMAIN_ARITH(operator/, /=, MultResult_t)


//////////////////////////////////////////////////////////////////////
//
// Routines to convert between a Domain-like object and a
// Vector-like object.  These convert from the first argument into
// the second.
//
//////////////////////////////////////////////////////////////////////

template<class D, class V>
inline void DomainToVector(const D &dom, V &vec)
{
  // Loop over the number of elements in the type D, and assign
  // them to elements in vec.

  for (int i=0; i < DomainTraits<D>::dimensions; ++i)
    vec(i) = dom[i].first();
}


template<class V, class D>
inline void VectorToDomain(const V &vec, D &dom)
{
  // Loop over the number of elements in the type D, and assign to
  // dom from the first argument.

  for (int i=0; i < DomainTraits<D>::dimensions; ++i)
    dom[i] = DomainTraits<D>::OneDomain_t(vec(i), vec(i));
}


//////////////////////////////////////////////////////////////////////
//
// print a domain to a stream, in the format
//   "[" first:last:stride, first:last:stride, ... first:last:stride "]"
//
//////////////////////////////////////////////////////////////////////

template<class DT>
std::ostream& operator<<(std::ostream &o, const DomainBase<DT> &dbase)
{
  dbase.print(o);
  return o;
}


//////////////////////////////////////////////////////////////////////

#endif     // POOMA_DOMAIN_DOMAIN_BASE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DomainBase.h,v $   $Author: richard $
// $Revision: 1.32 $   $Date: 2004/11/01 18:16:31 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
