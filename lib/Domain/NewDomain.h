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

#ifndef POOMA_DOMAIN_NEWDOMAIN_H
#define POOMA_DOMAIN_NEWDOMAIN_H

//-----------------------------------------------------------------------------
// Class:
// NewDomain structs
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Domain
 * @brief
 * A set of simple structs which tell how to combine different Domain
 * objects together.
 *
 * They are named NewDomain1 ... NewDomain7, are
 * templated on from 1 ... 7 different domain types, and provide the
 * following interface:
 *  - NewDomain2<T1,T2>::Type_t newdom; // resulting type when Domains combined
 *  - NewDomain2<T1,T2>::SliceType_t slicedom; // type for 'sliced' Dom's
 *  - newdom = NewDomain2<T1,T2>::combine(a,b); // combine a & b, return combo
 *  - NewDomain2<T1,T2>::fill(newdom, a, b); // combine a & b into newdom
 *  - slicedom = NewDomain2<T1,T2>::combineSlice(a,b); // 'slice' a & b
 *  - NewDomain2<T1,T2>::fillSlice(slicedom, a, b); // 'slice' into slicedom
 *
 * similarly for NewDomain1, and NewDomain3 ... NewDomain7
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/DomainTraits.h"
#include "Utilities/PAssert.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

template<int Dim>               class Loc;
template<int Dim>               class Interval;
template<int Dim>               class Range;
template<int Dim>               class Grid;
template<int Dim>               class AllDomain;
template<int Dim>               class LeftDomain;
template<int Dim>               class RightDomain;
template<int Dim, int SliceDim> class SliceInterval;
template<int Dim, int SliceDim> class SliceRange;
template<int Dim, class T>      class Region;
template<class T>               class IndirectionList;


/**
 * CombineDomain is a utility class used only by the 'fill*' methods
 * of NewDomain2 ... NewDomain7.  It is templated on the following params:
 *  - RT : the return domain type, which is the domain we are constructing
 *  - CT : the domain type being combined with RT
 *  - DS : the domain index at which we should start adding CT into RT
 * CombineDomain takes a domain object of type CT, and combines it into
 * an existing domain object RT.  The domain RT is then changed.  DS provides
 * the 'beginning' index at which CT should be added, for example if RT
 * is a dim-5 domain, CT a dim-2 domain, and DS is 1, then RT[1] is
 * set to be CT[0], and RT[2] is set to be CT[1].
 *
 * For "slice" rules (used by fillSlice methods), a separate CombineSliceDomain
 * class is available.  The complication here is that SliceDomain's store two
 * pieces of information (the "total" domain, and the "slice" domain), instead
 * of just one (the domain).  To store these values, CombineSliceDomain takes
 * a slightly different set of parameters:
 *  - RT: the return domain type, which is the domain we are constructing
 *  - CT: the domain type being combined with RT
 *  - DS: the "total" domain index at which we should start adding CT into RT
 *  - SliceDS: the "slice" domain index at which we should start adding CT into RT
 *  - incl: if this is false, then either RT is not a SliceDomain, or the
 *          type CT is a "singleValued" domain which should not be included in
 *          the slice domain.
 *          If this is true, then RT is a SliceDomain, and CT is a
 *          domain such as Interval or Range that should be included in both
 *          the total domain and the slice domain.  To do this, the class
 *          is partially specialized to the case where incl=true.
 *          Those specializations fill in CT into both the total
 *          domain of RT and the slice domain of RT.
 *  - wildcard: if this is true, the domain being combined, of type CT, is a
 *          wildcard domain, so the user must provide a "reference" domain
 *          with which to combine.  The wildcard takes the reference domain
 *          and copies or modifies it based on the type of wildcard.
 *          if this is false, the CT domain is not a wildcard, and the user-
 *          supplied reference domain is ignored.
 */

//
// the general version of CombineDomain.   By using DomainTraits::getDomain
// on rt, this will work for even if rt is a 1D domain.
//

template<class RT, class CT, int DS>
struct CombineDomain {
  enum { DRT = DomainTraits<RT>::dimensions };
  enum { DCT = DomainTraits<CT>::dimensions };

  static void combine(RT &rt, const CT& ct) {
    CTAssert(DS >= 0);
    CTAssert(DRT > (DS + DCT - 1));
    for (int i=0; i < DCT; ++i)
      DomainTraits<RT>::getDomain(rt, DS + i).setDomain(
	                        DomainTraits<CT>::getDomain(ct, i));
  }
};


/**
 * the general version of CombineSliceDomainWC ... this class is used
 * by CombineSliceDomain to set up slice domain objects.  It is similar
 * to CombineDomain except that it will on occasion fill in a separate
 * 'slice domain'.  It can also use wildcard domains, if it is given a
 * reference domain which is used by the wildcard to determine the correct
 * domain.  If no slicing is being done, and no wildcard is used, this
 * just fills in the reduced 'slice domain'.
 */

template<class RT, class UT, class CT, int DS, int SliceDS,
         bool incl, bool wc>
struct CombineSliceDomainWC {
  enum { DRT = DomainTraits<RT>::dimensions };
  enum { DCT = DomainTraits<CT>::dimensions };

  static void combine(RT &rt, const UT &, const CT& ct) {
    CTAssert(DS >= 0);
    CTAssert(DRT > (DS + DCT - 1));
    for (int i=0; i < DCT; ++i)
      DomainTraits<RT>::getDomain(rt, DS + i).setDomain(
	                        DomainTraits<CT>::getPointDomain(ct, i));
  }
};

/**
 * specialization of CombineSliceDomainWC in which we fill in slice
 * domain values and full domain values, but without using wildcards
 */

template<class RT, class UT, class CT, int DS, int SliceDS>
struct CombineSliceDomainWC<RT,UT,CT,DS,SliceDS,true,false> {
  enum { DRT = DomainTraits<RT>::dimensions };
  enum { DCT = DomainTraits<CT>::dimensions };

  static void combine(RT &rt, const UT &, const CT& ct) {
    CTAssert(DS >= 0 && SliceDS >= 0);
    CTAssert(DRT > (DS + DCT - 1));
    for (int i=0; i < DCT; ++i) {
      DomainTraits<RT>::getDomain(rt, DS + i).setDomain(
	DomainTraits<CT>::getPointDomain(ct, i));
      DomainTraits<RT>::setIgnorable(rt, DS + i,
        DomainTraits<CT>::getIgnorable(ct, i));
    }
    rt.setSliceFromTotal();
  }
};

/**
 * specialization of CombineSliceDomainWC in which we only fill in total
 * domain values, using wildcards
 */

template<class RT, class UT, class CT, int DS, int SliceDS>
struct CombineSliceDomainWC<RT,UT,CT,DS,SliceDS,false,true> {
  enum { DRT = DomainTraits<RT>::dimensions };
  enum { DUT = DomainTraits<UT>::dimensions };
  enum { DCT = DomainTraits<CT>::dimensions };

  static void combine(RT &rt, const UT &u, const CT& ct) {
    CTAssert(DS >= 0);
    CTAssert(DRT > (DS + DCT - 1));
    CTAssert(DUT == DRT);
    for (int i=0; i < DCT; ++i)
      DomainTraits<RT>::getDomain(rt, DS + i).setWildcardDomain(
	DomainTraits<UT>::getPointDomain(u, DS + i),
	DomainTraits<CT>::getPointDomain(ct, i));
  }
};

/**
 * specialization of CombineSliceDomainWC in which we fill in slice
 * domain values and full domain values, using wildcards
 */

template<class RT, class UT, class CT, int DS, int SliceDS>
struct CombineSliceDomainWC<RT,UT,CT,DS,SliceDS,true,true> {
  enum { DRT = DomainTraits<RT>::dimensions };
  enum { DUT = DomainTraits<UT>::dimensions };
  enum { DCT = DomainTraits<CT>::dimensions };

  static void combine(RT &rt, const UT &u, const CT& ct) {
    CTAssert(DS >= 0 && SliceDS >= 0);
    CTAssert(DRT > (DS + DCT - 1));
    CTAssert((int)DUT == DRT);
    for (int i=0; i < DCT; ++i) {
      DomainTraits<RT>::getDomain(rt, DS + i).setWildcardDomain(
	DomainTraits<UT>::getPointDomain(u, DS + i),
	DomainTraits<CT>::getPointDomain(ct, i));
      DomainTraits<RT>::getSliceDomain(rt, SliceDS + i).setWildcardDomain(
	DomainTraits<UT>::getPointDomain(u, DS + i),
	DomainTraits<CT>::getPointDomain(ct, i));
      DomainTraits<RT>::cantIgnoreDomain(rt, DS + i);
    }
  }
};

/**
 * the general version of CombineSliceDomain ... by default, it just
 * does the same thing as CombineDomain, except, for domains which store
 * a slice, it will fill in a second 'total' domain with the extra info
 * about the domains that are sliced out.  If the boolean type 'incl'
 * is true, there is a specialization here to also fill in
 * these slice dimensions.  If the combining domain is a wildcard, then
 * we use a separate method to fill in the domain using the user-supplied
 * reference domain.  To get all this done, CombineSliceDomain defers to
 * a separate CombineSliceDomainWC struct which has an extra boolean
 * template param 'wildcard' indicating whether to use wildcard set routines.
 */

template<class RT, class UT, class CT, int DS, int SliceDS, bool incl>
struct CombineSliceDomain {
  static void combine(RT &rt, const UT &u, const CT& ct) {
    CombineSliceDomainWC<RT,UT,CT,DS,SliceDS,incl,
      DomainTraits<CT>::wildcard>::combine(rt, u, ct);
  }
};


/** @defgroup NewDomain
 * @ingroup Domain
 * NewDomain2 ... NewDomain7 are simple helper structs which are used to
 * combine Domain objects of different types together (when possible).
 * They are basically traits classes used to tell what the type is when
 * you combine several domains together, and what the combined object
 * looks like.
 *
 * If you have N Domain objects to combine (where N is from 1 to 7),
 * the type of the Domain object which results when you combine the
 * Domain objects together is (for, say, combining 4 Domains of types
 * T1, T2, T3, and T4 together)
 *    NewDomain4<T1,T2,T3,T4>::Type_t
 *
 * The new Domain will be one of the Domain types listed in the include and
 * forward ref sections above.  Generally, combining a more specific
 * type of domain with a more general type of domain (say, a Loc with a
 * Range) results in return domain type of the more general category.  E.g.,
 * in the case of combining a Loc with Range, the result would be a Range.
 *
 * The dimension of the new Domain will be the sum of the dimensions of
 * the input objects.  The dimensions will be filled in from left to right
 * from the arguments, going from 1 ... Dim for the dimensions.  For example,
 * if you combine three domains, the first of Dim 2, the next of Dim 1,
 * the third of Dim 2, the result will be of Dim=5, and will have Dim 1&2
 * taken from the first argument, Dim 3 from the second, and Dim 4&5 from
 * the third.
 *
 * To combine domains together, there are two possibilities:
 * -# combine them into a separate object, returned by the 'combine' method
 *    This takes N domain objects, and puts them into a new domain object
 *    of type NewDomainN<T1,T2,...,TN>::Type_t.  The method
 *        Type_t NewDomainN<...>::combine(a,b,c,d,...)
 *    will combine the elements a, b, c, ... and return a new domain object
 * -# fill them into a provided domain object, using the 'fill' method
 *    This does a similar thing as 'combine', but is more efficient since
 *    the user provides (in the first argument to 'fill') the domain object
 *    into which the elements should be combined.  The syntax is
 *        void NewDomainN<...>::fill(retval,a,b,c,d,...)
 *    where a,b,c,d,... are the domain object to combine together, and
 *    retval is the domain object which should hold the result.  Note that
 *    retval is a templated item, it does not have to be the same type as
 *    NewDomainN<...>::Type_t.  fill returns a reference to the retval
 *    item.
 *
 * Both of these methods will perform 1D assigments of pieces of the domains
 * to the resulting domain object, using the setDomain method of a 1D
 * domain object.
 *
 * Determining the type of the combined domain is done pair-wise, starting
 * with the first two types, then combining the result of that with the
 * third type, etc.  So, most of the work here is in NewDomain2, which
 * must be specialized for the different possible pairwise combinations.
 * The other NewDomainN's then use NewDomain2 and NewDomain(N-1)'s.
 *
 * There is actually more than one way to combine Domain objects together.
 * Each method has an associated set of typedef and combine/fill methods.
 * The above discussion was just for one of these methods; the different
 * combining methods are:
 *   -# "Domain/Array constructor" rules:
 *      this combining method is applicable to
 *      the situation that exists when combining domain objects together
 *      in the constructors/operator='s of domain objects.  In this case,
 *      the total dimension of the resulting object is the sum of the
 *      dimensions of the arguments.  Loc<1> objects are treated like
 *      single points, but int's are treated like Interval<1> objects of
 *      size 0 ... int value - 1.  For this case,
 *      the typedefs and combine/fill methods are:
 *       -  NewDomainN::Type_t;                               // typedef
 *       -  NewDomainN::Type_t NewDomainN::combine(a,b,...);  // combine
 *       -  void NewDomainN::fill(retval, a, b, ...);         // fill
 *   -# "Slice" rules:
 *      this combining method is applicable to the situation
 *      that exists when combining domain objects together in the operator()'s
 *      of Array objects.  In this case, the resulting domain may represent
 *      a "slice" of the full domain which would have been created if
 *      "domain constructor" rules had been applied.  The slice is actually
 *      a lower-dimensional subdomain of the complete domain.  Specifying
 *      a Loc or int for one or more of the dimensions in the operator() will
 *      result in that dimension being "sliced out", so that the total dim
 *      of the slice = sum of dimensions of Interval, Range objects
 *      - sum dimensions of Loc and int objects.  If slice rules are used,
 *      and any Loc or int objects are being combined in, the resulting
 *      domain object will be a SliceDomain subclass (e.g., SliceInterval,
 *      SliceRange, whatever is the most general based on the other combined
 *      domain types).  This SliceDomain subclass is templated on the total
 *      dimension, and the slice dimension, where slicedim < totaldim.  If
 *      slice rules are used, but no singleValue'd domains are used, then
 *      regular "domain constructor" rules apply.  The pairwise combination
 *      mechanism only allows you to combine a SliceDomain subclass with
 *      a non-SliceDomain, which means that you cannot combine, say,
 *      a SliceInterval with a SliceRange.  Since slice domain objects are
 *      only meant for internal use when constructing a subdimensional
 *      Array from a larger dimensional Array, this is all that is needed.
 *      For this case, the typedefs and methods are:
 *       - NewDomainN::SliceType_t;                                  // typedef
 *       - NewDomainN::SliceType_t NewDomainN::combineSlice(a,b,..); // combine
 *       - void NewDomainN::fillSlice(retval, a, b, ...);            // fill
 *
 * So, to use NewDomainN objects, first determine what kind of combining
 * rules you need in your particular context, and then use the relevant
 * set of typedef and combine/fill methods from the list above.
 */



/** @addtogroup NewDomain
 * A simple base class for all specialized NewDomain2 objects, defining
 * the functions which all have in common.  This is most easily done through
 * a base class, since we need to define several different specializations
 * of NewDomain2.
 */

template<class T1, class T2, class TCombine, class TSCombine>
struct NewDomain2Base
{
  // the Type_t typedef is provided by the derived class as the TCombine param,
  // the SliceType_t typedef is provided via the TSCombine param.
  typedef TCombine  Type_t;
  typedef TSCombine SliceType_t;

  // static data 
  enum { S2 = DomainTraits<T1>::dimensions };
  enum { DX1 = DomainTraits<T1>::sliceDimensions };
  enum { DX2 = DomainTraits<T2>::sliceDimensions };

  // combine is used to combine two items together
  inline static Type_t combine(const T1 &a, const T2 &b)
    {
      Type_t retval = Pooma::NoInit();
      return fill(retval, a, b);
    }

  // fill does the same as combine, but fills info into an arbitrary type
  // RT instead of into a newly-created instance of type Type_t
  template<class RT>
  inline static RT &fill(RT &retval, const T1 &a, const T2 &b)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      CombineDomain<RT,T2,S2>::combine(retval,b);
      return retval;
    }

  // combineSlice is used to combine items together using the 'slice' domain.
  // It uses a given 'user' domain which, if wildcards are being combined
  // here, is used to fill in the resulting domain.  It should have the
  // same number of dimensions as the SliceType_t domain.
  template<class UT>
  inline static SliceType_t combineSlice(const UT &u,
					 const T1 &a, const T2 &b)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a, b);
    }

  // fillSlice does the same as combineSlice, but fills info into an arbitrary
  // type RT instead of into a newly-created instance of type SliceType_t.
  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a, const T2 &b)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::
	combine(retval,u,a);
      CombineSliceDomain<RT,UT,T2,S2,DX1,(DX2>0 && RDX)>::
	combine(retval,u,b);
      return retval;
    }
};


/** @class NewDomain2
 * @addtogroup NewDomain
 * The general version of NewDomain2.  The allowed versions of NewDomain2
 * are given as partial specializations of this general case.  The general
 * case assumes T1 and T2 are single-valued domains for which DomainTraits
 * exist, that combine together to form Interval's (or Loc's for slice
 * combine rules).
 */

// first, create a simple struct used to add two dimensions together at
// compile time
template<class T1, class T2>
struct AddNewDomain2Dimensions
{
  enum { dimensions = 
    DomainTraits<T1>::dimensions + DomainTraits<T2>::dimensions };
};


template<class T1, class T2>
struct NewDomain2
  : public NewDomain2Base<T1, T2,
                          Interval<AddNewDomain2Dimensions<T1,T2>::dimensions>,
                          Loc<AddNewDomain2Dimensions<T1,T2>::dimensions> >
{ };


//-----------------------------------------------------------------------------
// Specific versions of NewDomain2, for all the allowed combinations.
//-----------------------------------------------------------------------------


// macros for use in defining NewDomain2.  The first sets up the
// combination of a domain with itself and with an int or Loc.  The second
// sets up the combination of a domain with another domain of a different
// type.  The third sets up the combination of a SliceDomain with another
// domain type.
//
// The rules for combining domains using "domain/array constructor" rules are:
//   1. int's, char's, shorts's, long's are treated as Interval<1>(int val),
//      e.g., 7 ==> Interval<1>(7)
//   2. The resulting domain is the most general which can store all the
//      combined objects.
//   3. The final number of dimensions = sum of dimension of combined domains
//
// The rules for combining slice domains deserve special attention.  They are:
//   1. Slice domains can only be combined with non-slice domains
//   2. Combining a slice domain with an int or Loc (or any single-valued
//      domain) increases the total dimension of the domain, but not the
//      slice dimension
//   3. Combining a slice domain with any other non-slice domain increases
//      both the total dimension and the slice dimension

#define POOMA_NEWDOMAIN_SAME_SCALAR(DOM,SLICEDOM,S)			   \
template <int D>							   \
struct NewDomain2< DOM<D>, S >						   \
  : public NewDomain2Base< DOM<D>, S, DOM<D+1>, SLICEDOM<D+1,D> > { };	   \
template <int D>							   \
struct NewDomain2< S, DOM<D> >						   \
  : public NewDomain2Base< S, DOM<D>, DOM<D+1>, SLICEDOM<D+1,D> > { };	   \
template <int D>							   \
struct NewDomain2< DOM<D>, unsigned S >					   \
  : public NewDomain2Base< DOM<D>, unsigned S, DOM<D+1>, SLICEDOM<D+1,D> > \
{ };                                                                       \
template <int D>							   \
struct NewDomain2< unsigned S, DOM<D> >					   \
  : public NewDomain2Base< unsigned S, DOM<D>, DOM<D+1>, SLICEDOM<D+1,D> > \
{ };


#define POOMA_NEWDOMAIN_SAME(DOM,SLICEDOM)				   \
template <int D1, int D2>						   \
struct NewDomain2< DOM<D1>, DOM<D2> >					   \
  : public NewDomain2Base< DOM<D1>, DOM<D2>, DOM<D1+D2>, DOM<D1+D2> > { }; \
template <int D1, int D2>						   \
struct NewDomain2< DOM<D1>, Loc<D2> >					   \
  : public NewDomain2Base< DOM<D1>, Loc<D2>, DOM<D1+D2>,                   \
                           SLICEDOM<D1+D2,D1> > { };                       \
template <int D1, int D2>						   \
struct NewDomain2< Loc<D2>, DOM<D1> >					   \
  : public NewDomain2Base< Loc<D2>, DOM<D1>, DOM<D1+D2>,                   \
                           SLICEDOM<D1+D2,D1> > { };                       \
POOMA_NEWDOMAIN_SAME_SCALAR(DOM,SLICEDOM,char)				   \
POOMA_NEWDOMAIN_SAME_SCALAR(DOM,SLICEDOM,short)			           \
POOMA_NEWDOMAIN_SAME_SCALAR(DOM,SLICEDOM,int)				   \
POOMA_NEWDOMAIN_SAME_SCALAR(DOM,SLICEDOM,long)


#define POOMA_NEWDOMAIN_OTHER(DOM1,DOM2)				   \
template <int D1, int D2>						   \
struct NewDomain2< DOM1<D1>, DOM2<D2> >					   \
  : public NewDomain2Base< DOM1<D1>, DOM2<D2>, DOM1<D1+D2>, DOM1<D1+D2> >  \
{ };                                                                       \
template <int D1, int D2>						   \
struct NewDomain2< DOM2<D1>, DOM1<D2> >					   \
  : public NewDomain2Base< DOM2<D1>, DOM1<D2>, DOM1<D1+D2>, DOM1<D1+D2> >  \
{ };



#define POOMA_NEWDOMAIN_SLICE_SAME_SCALAR(SLICEDOM,S)			   \
template <int D1, int DS1>						   \
struct NewDomain2< SLICEDOM<D1,DS1>, S >				   \
  : public NewDomain2Base< SLICEDOM<D1,DS1>, S ,SLICEDOM<D1+1,DS1>,	   \
                           SLICEDOM<D1+1,DS1> > { };			   \
template <int D1, int DS1>						   \
struct NewDomain2< SLICEDOM<D1,DS1>, unsigned S >			   \
  : public NewDomain2Base< SLICEDOM<D1,DS1>, unsigned S,                   \
                           SLICEDOM<D1+1,DS1>, SLICEDOM<D1+1,DS1> > { };   \
template <int D1, int DS1>						   \
struct NewDomain2< S, SLICEDOM<D1,DS1> >				   \
  : public NewDomain2Base< S, SLICEDOM<D1,DS1>, SLICEDOM<D1+1,DS1>,	   \
                           SLICEDOM<D1+1,DS1> > { };			   \
template <int D1, int DS1>						   \
struct NewDomain2< unsigned S, SLICEDOM<D1,DS1> >			   \
  : public NewDomain2Base< unsigned S, SLICEDOM<D1,DS1>,                   \
                           SLICEDOM<D1+1,DS1>, SLICEDOM<D1+1,DS1> > { };


#define POOMA_NEWDOMAIN_SLICE_SAME(SLICEDOM)				   \
template <int D1, int DS1, int D2>					   \
struct NewDomain2< SLICEDOM<D1,DS1>, Loc<D2> >				   \
  : public NewDomain2Base< SLICEDOM<D1,DS1>, Loc<D2>, SLICEDOM<D1+D2,DS1>, \
                           SLICEDOM<D1+D2,DS1> > { };			   \
template <int D1, int DS1, int D2>					   \
struct NewDomain2< Loc<D2>, SLICEDOM<D1,DS1> >				   \
  : public NewDomain2Base< Loc<D2>, SLICEDOM<D1,DS1>, SLICEDOM<D1+D2,DS1>, \
                           SLICEDOM<D1+D2,DS1> > { };			   \
POOMA_NEWDOMAIN_SLICE_SAME_SCALAR(SLICEDOM,char)			   \
POOMA_NEWDOMAIN_SLICE_SAME_SCALAR(SLICEDOM,short)			   \
POOMA_NEWDOMAIN_SLICE_SAME_SCALAR(SLICEDOM,int)			           \
POOMA_NEWDOMAIN_SLICE_SAME_SCALAR(SLICEDOM,long)


#define POOMA_NEWDOMAIN_SLICE_OTHER(DOM1,DOM2,SLICEDOM)	                   \
template <int D1, int DS1, int D2>				           \
struct NewDomain2< DOM1<D1,DS1>, DOM2<D2> >			           \
  : public NewDomain2Base< DOM1<D1,DS1>, DOM2<D2>, SLICEDOM<D1+D2,DS1+D2>, \
			   SLICEDOM<D1+D2,DS1+D2> > { };		   \
template <int D1, int DS1, int D2>				           \
struct NewDomain2< DOM2<D2>, DOM1<D1,DS1> >			           \
  : public NewDomain2Base< DOM2<D2>, DOM1<D1,DS1>, SLICEDOM<D1+D2,DS1+D2>, \
			   SLICEDOM<D1+D2,DS1+D2> > { };


//
// Range with others
//

POOMA_NEWDOMAIN_SAME(Range,SliceRange)
POOMA_NEWDOMAIN_OTHER(Range,Interval)
POOMA_NEWDOMAIN_OTHER(Range,AllDomain)
POOMA_NEWDOMAIN_OTHER(Range,LeftDomain)
POOMA_NEWDOMAIN_OTHER(Range,RightDomain)
POOMA_NEWDOMAIN_SLICE_SAME(SliceRange)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceRange,Range,SliceRange)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceRange,Interval,SliceRange)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceRange,AllDomain,SliceRange)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceRange,LeftDomain,SliceRange)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceRange,RightDomain,SliceRange)

//
// Interval with others
//

POOMA_NEWDOMAIN_SAME(Interval,SliceInterval)
POOMA_NEWDOMAIN_OTHER(Interval,AllDomain)
POOMA_NEWDOMAIN_OTHER(Interval,LeftDomain)
POOMA_NEWDOMAIN_OTHER(Interval,RightDomain)
POOMA_NEWDOMAIN_SLICE_SAME(SliceInterval)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceInterval,Interval,SliceInterval)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceInterval,Range,SliceRange)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceInterval,AllDomain,SliceInterval)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceInterval,LeftDomain,SliceInterval)
POOMA_NEWDOMAIN_SLICE_OTHER(SliceInterval,RightDomain,SliceInterval)

//
// Wildcards with themselves
//

POOMA_NEWDOMAIN_SAME(AllDomain,SliceInterval)
POOMA_NEWDOMAIN_SAME(LeftDomain,SliceInterval)
POOMA_NEWDOMAIN_SAME(RightDomain,SliceInterval)
POOMA_NEWDOMAIN_OTHER(AllDomain,LeftDomain)
POOMA_NEWDOMAIN_OTHER(AllDomain,RightDomain)
POOMA_NEWDOMAIN_OTHER(LeftDomain,RightDomain)


//
// Grid with others
//

POOMA_NEWDOMAIN_SAME(Grid,SliceRange)
POOMA_NEWDOMAIN_OTHER(Grid,Range)
POOMA_NEWDOMAIN_OTHER(Grid,Interval)
POOMA_NEWDOMAIN_OTHER(Grid,AllDomain)
POOMA_NEWDOMAIN_OTHER(Grid,LeftDomain)
POOMA_NEWDOMAIN_OTHER(Grid,RightDomain)

//
// Grid with IndirectionList
//

template<int D>
struct NewDomain2< Grid<D>, IndirectionList<int> >
  : public NewDomain2Base<Grid<D>,IndirectionList<int>,Grid<D+1>,Grid<D+1> >{};

template<int D>
struct NewDomain2< IndirectionList<int>, Grid<D> >
  : public NewDomain2Base<IndirectionList<int>,Grid<D>,Grid<D+1>,Grid<D+1> >{};


//
// complete specializations for Loc with Loc.  Other combinations of Loc
// with a scalar are covered by the default case.
//

template<int D1, int D2>
struct NewDomain2< Loc<D1>, Loc<D2> >
  : public NewDomain2Base<Loc<D1>, Loc<D2>, Loc<D1+D2>, Loc<D1+D2> > { };


//
// macros for use in defining combinations of continuous domains
//

#define POOMA_NEWDOMAIN_CONTINUOUS_SAME(DOM)			           \
template <int D1, class T1, int D2, class T2>			           \
struct NewDomain2< DOM<D1,T1>, DOM<D2,T2> >			           \
  : public NewDomain2Base< DOM<D1,T1>, DOM<D2,T2>, DOM<D1+D2,T1>,	   \
                           DOM<D1+D2,T1> > { };


#define POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(DOM,S)			   \
template <int D1, class T1>						   \
struct NewDomain2< DOM<D1,T1>, S >					   \
  : public NewDomain2Base< DOM<D1,T1>, S, DOM<D1+1,T1>, DOM<D1+1,T1> >     \
{ };	                                                                   \
template <int D1, class T1>						   \
struct NewDomain2< S, DOM<D1,T1> >					   \
  : public NewDomain2Base< S, DOM<D1,T1>, DOM<D1+1,T1>, DOM<D1+1,T1> >     \
{ };


#define POOMA_NEWDOMAIN_CONTINUOUS_OTHER(DOM1,DOM2)		           \
template <int D1, class T1, int D2>				           \
struct NewDomain2< DOM1<D1,T1>, DOM2<D2> >			           \
  : public NewDomain2Base< DOM1<D1,T1>, DOM2<D2>, DOM1<D1+D2,T1>,	   \
                           DOM1<D1+D2,T1> > { };			   \
template <int D1, class T1, int D2>				           \
struct NewDomain2< DOM2<D2>, DOM1<D1,T1> >			           \
  : public NewDomain2Base< DOM2<D2>, DOM1<D1,T1>, DOM1<D1+D2,T1>,	   \
                           DOM1<D1+D2,T1> > { };


#define POOMA_NEWDOMAIN_JUST_SCALAR_SAME(S,DOM)		                   \
template <>							           \
struct NewDomain2<S, S>						           \
  : public NewDomain2Base< S, S, DOM<2,S>, DOM<2,S> > { };


#define POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(S1,S2,DOM,S3)	                   \
template <>							           \
struct NewDomain2<S1, S2>					           \
  : public NewDomain2Base< S1, S2, DOM<2,S3>, DOM<2,S3> > { };	           \
template <>							           \
struct NewDomain2<S2, S1>					           \
  : public NewDomain2Base< S2, S1, DOM<2,S3>, DOM<2,S3> > { };


#define POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(S,DOM1,DOM2)		           \
template <int D1>							   \
struct NewDomain2< S, DOM1<D1> >					   \
  : public NewDomain2Base< S, DOM1<D1>, DOM2<D1+1,S>, DOM2<D1+1,S> > { };  \
template <int D1>							   \
struct NewDomain2< DOM1<D1>, S >					   \
  : public NewDomain2Base< DOM1<D1>, S, DOM2<D1+1,S>, DOM2<D1+1,S> > { };

//
// combinations involving Region
//

POOMA_NEWDOMAIN_CONTINUOUS_SAME(Region)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,char)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,unsigned char)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,short)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,unsigned short)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,int)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,unsigned int)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,long)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,unsigned long)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,float)
POOMA_NEWDOMAIN_CONTINUOUS_SCALAR(Region,double)
POOMA_NEWDOMAIN_CONTINUOUS_OTHER(Region,Range)
POOMA_NEWDOMAIN_CONTINUOUS_OTHER(Region,Interval)
POOMA_NEWDOMAIN_CONTINUOUS_OTHER(Region,Loc)
POOMA_NEWDOMAIN_CONTINUOUS_OTHER(Region,AllDomain)
POOMA_NEWDOMAIN_CONTINUOUS_OTHER(Region,LeftDomain)
POOMA_NEWDOMAIN_CONTINUOUS_OTHER(Region,RightDomain)

//
// combinations involving scalars other than int
//

POOMA_NEWDOMAIN_JUST_SCALAR_SAME(double,Region)
POOMA_NEWDOMAIN_JUST_SCALAR_SAME(float,Region)

POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,char,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,unsigned char,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,short,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,unsigned short,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,int,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,unsigned int,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,long,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,unsigned long,Region,double)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(double,float,Region,double)

POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,char,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,unsigned char,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,short,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,unsigned short,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,int,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,unsigned int,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,long,Region,float)
POOMA_NEWDOMAIN_JUST_SCALAR_OTHER(float,unsigned long,Region,float)

POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(double,Loc,Region)
POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(double,Interval,Region)
POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(double,Range,Region)
POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(float,Loc,Region)
POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(float,Interval,Region)
POOMA_NEWDOMAIN_JUST_SCALAR_DOMAIN(float,Range,Region)


/** @addtogroup NewDomain
 * NewDomainNBase: a base class for all NewDomainN classes for N >= 3.
 * This just simplifies the work of specifying all the typedefs for the
 * combination results.  The implementations of the combine and fill
 * routines are left to the NewDomainN classes which are derived from
 * NewDomainNBase.
 *    T = final type in list of types for NewDomainN.  Each NewDomainN is
 * templated on N types T1 ... TN.  T should be the same as TN.
 *   ND = type of lower-dimensional NewDomainN class which the subclass will
 * build on.  For example, NewDomain4<T1,..T4> should be derived from
 * NewDomainNBase<NewDomain3<T1,T2,T3>, T4>.  Thus, in this case, we have
 * ND ==> NewDomain3<T1,T2,T3>, and T ==> T4.
 */

template<class ND, class T>
struct NewDomainNBase
{
  typedef typename ND::Type_t      PrevType_t;
  typedef typename ND::SliceType_t PrevSliceType_t;

  typedef typename NewDomain2<PrevType_t,T>::Type_t           Type_t;
  typedef typename NewDomain2<PrevSliceType_t,T>::SliceType_t SliceType_t;
};


/** @addtogroup NewDomain
 * NewDomain1 ... NewDomain7: these build on each other pairwise, eventually
 * returning back to NewDomain2 for the final determination.  NewDomain1
 * does no actual combining, just a conversion; and NewDomain2 is defined
 * already above.
 */

template<class T1>
struct NewDomain1
{
  // typedefs
  typedef typename DomainTraits<T1>::Domain_t          Type_t;
  typedef typename DomainTraits<T1>::NewDomain1_t      SliceType_t;

  // static data
  enum { DX1 = DomainTraits<T1>::sliceDimensions };

  inline static Type_t combine(const T1 &a)
    {
      Type_t retval = Pooma::NoInit();
      return fill(retval, a);
    }

  template<class RT>
  inline static RT &fill(RT &retval, const T1 &a)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      return retval;
    }

  template<class UT>
  inline static SliceType_t combineSlice(const UT &u, const T1 &a)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a);
    }

  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::combine(retval,u,a);
      return retval;
    }
};


template<class T1, class T2, class T3>
struct NewDomain3 : public NewDomainNBase<NewDomain2<T1,T2>, T3>
{
  // typedefs
  typedef typename NewDomainNBase<NewDomain2<T1,T2>, T3>::Type_t Type_t;
  typedef typename NewDomainNBase<NewDomain2<T1,T2>, T3>::SliceType_t SliceType_t;
  
  // static data
  enum { S2 = DomainTraits<T1>::dimensions };
  enum { S3 = S2 + DomainTraits<T2>::dimensions };
  enum { DX1 = DomainTraits<T1>::sliceDimensions };
  enum { DX2 = DomainTraits<T2>::sliceDimensions };
  enum { DX3 = DomainTraits<T3>::sliceDimensions };

  inline static Type_t combine(const T1 &a, const T2 &b, const T3 &c)
    {
      Type_t retval = Pooma::NoInit();
      return fill(retval, a, b, c);
    }  

  template<class RT>
  inline static RT &fill(RT &retval,
			 const T1 &a, const T2 &b, const T3 &c)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      CombineDomain<RT,T2,S2>::combine(retval,b);
      CombineDomain<RT,T3,S3>::combine(retval,c);
      return retval;
    }

  template<class UT>
  inline static SliceType_t combineSlice(const UT &u,
					 const T1 &a, const T2 &b, const T3 &c)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a, b, c);
    }

  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a, const T2 &b, const T3 &c)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::combine(retval,u,a);
      CombineSliceDomain<RT,UT,T2,S2,DX1,(DX2>0 && RDX)>::combine(retval,u,b);
      CombineSliceDomain<RT,UT,T3,S3,DX1+DX2,(DX3>0 && RDX)>::
        combine(retval,u,c);
      return retval;
    }
};


template<class T1, class T2, class T3, class T4>
struct NewDomain4 : public NewDomainNBase<NewDomain3<T1,T2,T3>, T4>
{
  // typedefs
  typedef typename NewDomainNBase<NewDomain3<T1,T2,T3>, T4>::Type_t Type_t;
  typedef typename NewDomainNBase<NewDomain3<T1,T2,T3>, T4>::SliceType_t 
    SliceType_t;
  
  // static data
  enum { S2 = DomainTraits<T1>::dimensions };
  enum { S3 = S2 + DomainTraits<T2>::dimensions };
  enum { S4 = S3 + DomainTraits<T3>::dimensions };
  enum { DX1 = DomainTraits<T1>::sliceDimensions };
  enum { DX2 = DomainTraits<T2>::sliceDimensions };
  enum { DX3 = DomainTraits<T3>::sliceDimensions };
  enum { DX4 = DomainTraits<T4>::sliceDimensions };

  inline static Type_t combine(const T1 &a, const T2 &b, const T3 &c,
			       const T4 &d)
    {
      Type_t retval = Pooma::NoInit();
      return fill(retval, a, b, c, d);
    }

  template<class RT>
  inline static RT &fill(RT &retval,
			 const T1 &a, const T2 &b, const T3 &c,
			 const T4 &d)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      CombineDomain<RT,T2,S2>::combine(retval,b);
      CombineDomain<RT,T3,S3>::combine(retval,c);
      CombineDomain<RT,T4,S4>::combine(retval,d);
      return retval;
    }

  template<class UT>
  inline static SliceType_t combineSlice(const UT &u,
					 const T1 &a, const T2 &b, const T3 &c,
					 const T4 &d)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a, b, c, d);
    }

  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a, const T2 &b, const T3 &c,
			      const T4 &d)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::combine(retval,u,a);
      CombineSliceDomain<RT,UT,T2,S2,DX1,(DX2>0 && RDX)>::combine(retval,u,b);
      CombineSliceDomain<RT,UT,T3,S3,DX1+DX2,(DX3>0 && RDX)>::
        combine(retval,u,c);
      CombineSliceDomain<RT,UT,T4,S4,DX1+DX2+DX3,(DX4>0 && RDX)>::
        combine(retval,u,d);
      return retval;
    }
};


template<class T1, class T2, class T3, class T4, class T5>
struct NewDomain5 : public NewDomainNBase<NewDomain4<T1,T2,T3,T4>, T5>
{
  // typedefs
  typedef typename NewDomainNBase<NewDomain4<T1,T2,T3,T4>, T5>::Type_t Type_t;
  typedef typename NewDomainNBase<NewDomain4<T1,T2,T3,T4>, T5>::SliceType_t 
    SliceType_t;
  
  // static data
  enum { S2 = DomainTraits<T1>::dimensions };
  enum { S3 = S2 + DomainTraits<T2>::dimensions };
  enum { S4 = S3 + DomainTraits<T3>::dimensions };
  enum { S5 = S4 + DomainTraits<T4>::dimensions };
  enum { DX1 = DomainTraits<T1>::sliceDimensions };
  enum { DX2 = DomainTraits<T2>::sliceDimensions };
  enum { DX3 = DomainTraits<T3>::sliceDimensions };
  enum { DX4 = DomainTraits<T4>::sliceDimensions };
  enum { DX5 = DomainTraits<T5>::sliceDimensions };

  inline static Type_t combine(const T1 &a, const T2 &b, const T3 &c,
			       const T4 &d, const T5 &e)
    {
      Type_t retval = Pooma::NoInit();
      return fill(retval, a, b, c, d, e);
    }

  template<class RT>
  inline static RT &fill(RT &retval,
			 const T1 &a, const T2 &b, const T3 &c,
			 const T4 &d, const T5 &e)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      CombineDomain<RT,T2,S2>::combine(retval,b);
      CombineDomain<RT,T3,S3>::combine(retval,c);
      CombineDomain<RT,T4,S4>::combine(retval,d);
      CombineDomain<RT,T5,S5>::combine(retval,e);
      return retval;
    }

  template<class UT>
  inline static SliceType_t combineSlice(const UT &u,
					 const T1 &a, const T2 &b, const T3 &c,
					 const T4 &d, const T5 &e)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a, b, c, d, e);
    }

  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a, const T2 &b, const T3 &c,
			      const T4 &d, const T5 &e)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::combine(retval,u,a);
      CombineSliceDomain<RT,UT,T2,S2,DX1,(DX2>0 && RDX)>::combine(retval,u,b);
      CombineSliceDomain<RT,UT,T3,S3,DX1+DX2,(DX3>0 && RDX)>::
        combine(retval,u,c);
      CombineSliceDomain<RT,UT,T4,S4,DX1+DX2+DX3,(DX4>0 && RDX)>::
        combine(retval,u,d);
      CombineSliceDomain<RT,UT,T5,S5,DX1+DX2+DX3+DX4,(DX5>0 && RDX)>::
        combine(retval,u,e);
      return retval;
    }
};


template<class T1, class T2, class T3, class T4, class T5, class T6>
struct NewDomain6 : public NewDomainNBase<NewDomain5<T1,T2,T3,T4,T5>, T6>
{
  // typedefs
  typedef typename 
    NewDomainNBase<NewDomain5<T1,T2,T3,T4,T5>, T6>::Type_t Type_t;
  typedef typename 
    NewDomainNBase<NewDomain5<T1,T2,T3,T4,T5>, T6>::SliceType_t SliceType_t;
  
  // static data
  enum { S2 = DomainTraits<T1>::dimensions };
  enum { S3 = S2 + DomainTraits<T2>::dimensions };
  enum { S4 = S3 + DomainTraits<T3>::dimensions };
  enum { S5 = S4 + DomainTraits<T4>::dimensions };
  enum { S6 = S5 + DomainTraits<T5>::dimensions };
  enum { DX1 = DomainTraits<T1>::sliceDimensions };
  enum { DX2 = DomainTraits<T2>::sliceDimensions };
  enum { DX3 = DomainTraits<T3>::sliceDimensions };
  enum { DX4 = DomainTraits<T4>::sliceDimensions };
  enum { DX5 = DomainTraits<T5>::sliceDimensions };
  enum { DX6 = DomainTraits<T6>::sliceDimensions };

  inline static Type_t combine(const T1 &a, const T2 &b,
			       const T3 &c, const T4 &d,
			       const T5 &e, const T6 &f)
    {
      Type_t retval = Pooma::NoInit();
      return fill(retval, a, b, c, d, e, f);
    }

  template<class RT>
  inline static RT &fill(RT &retval,
			 const T1 &a, const T2 &b,
			 const T3 &c, const T4 &d,
			 const T5 &e, const T6 &f)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      CombineDomain<RT,T2,S2>::combine(retval,b);
      CombineDomain<RT,T3,S3>::combine(retval,c);
      CombineDomain<RT,T4,S4>::combine(retval,d);
      CombineDomain<RT,T5,S5>::combine(retval,e);
      CombineDomain<RT,T6,S6>::combine(retval,f);
      return retval;
    }

  template<class UT>
  inline static SliceType_t combineSlice(const UT &u,
					 const T1 &a, const T2 &b,
					 const T3 &c, const T4 &d,
					 const T5 &e, const T6 &f)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a, b, c, d, e, f);
    }

  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a, const T2 &b, const T3 &c,
			      const T4 &d, const T5 &e, const T6 &f)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::combine(retval,u,a);
      CombineSliceDomain<RT,UT,T2,S2,DX1,(DX2>0 && RDX)>::combine(retval,u,b);
      CombineSliceDomain<RT,UT,T3,S3,DX1+DX2,(DX3>0 && RDX)>::
        combine(retval,u,c);
      CombineSliceDomain<RT,UT,T4,S4,DX1+DX2+DX3,(DX4>0 && RDX)>::
        combine(retval,u,d);
      CombineSliceDomain<RT,UT,T5,S5,DX1+DX2+DX3+DX4,(DX5>0 && RDX)>::
        combine(retval,u,e);
      CombineSliceDomain<RT,UT,T6,S6,DX1+DX2+DX3+DX4+DX5,(DX6>0 && RDX)>::
        combine(retval,u,f);
      return retval;
    }
};


template<class T1, class T2, class T3, class T4, class T5, class T6, class T7>
struct NewDomain7 : public NewDomainNBase<NewDomain6<T1,T2,T3,T4,T5,T6>, T7>
{
  // typedefs
  typedef typename 
    NewDomainNBase<NewDomain6<T1,T2,T3,T4,T5,T6>, T7>::Type_t Type_t;
  typedef typename 
    NewDomainNBase<NewDomain6<T1,T2,T3,T4,T5,T6>, T7>::SliceType_t SliceType_t;
  
  // static data
  enum { S2 = DomainTraits<T1>::dimensions };
  enum { S3 = S2 + DomainTraits<T2>::dimensions };
  enum { S4 = S3 + DomainTraits<T3>::dimensions };
  enum { S5 = S4 + DomainTraits<T4>::dimensions };
  enum { S6 = S5 + DomainTraits<T5>::dimensions };
  enum { S7 = S6 + DomainTraits<T6>::dimensions };
  enum { DX1 = DomainTraits<T1>::sliceDimensions };
  enum { DX2 = DomainTraits<T2>::sliceDimensions };
  enum { DX3 = DomainTraits<T3>::sliceDimensions };
  enum { DX4 = DomainTraits<T4>::sliceDimensions };
  enum { DX5 = DomainTraits<T5>::sliceDimensions };
  enum { DX6 = DomainTraits<T6>::sliceDimensions };
  enum { DX7 = DomainTraits<T7>::sliceDimensions };

  inline static Type_t combine(const T1 &a, const T2 &b,
			       const T3 &c, const T4 &d,
			       const T5 &e, const T6 &f,
			       const T7 &g)
    {
      Type_t retval = Pooma::NoInit();
      NewDomain7<T1,T2,T3,T4,T5,T6,T7>::fill(retval, a, b, c, d, e, f, g);
      return fill(retval, a, b, c, d, e, f,g);
    }

  template<class RT>
  inline static RT &fill(RT &retval,
			 const T1 &a, const T2 &b,
			 const T3 &c, const T4 &d,
			 const T5 &e, const T6 &f,
			 const T7 &g)
    {
      CombineDomain<RT,T1,0>::combine(retval,a);
      CombineDomain<RT,T2,S2>::combine(retval,b);
      CombineDomain<RT,T3,S3>::combine(retval,c);
      CombineDomain<RT,T4,S4>::combine(retval,d);
      CombineDomain<RT,T5,S5>::combine(retval,e);
      CombineDomain<RT,T6,S6>::combine(retval,f);
      CombineDomain<RT,T7,S7>::combine(retval,g);
      return retval;
    }

  template<class UT>
  inline static SliceType_t combineSlice(const UT &u,
					 const T1 &a, const T2 &b,
					 const T3 &c, const T4 &d,
					 const T5 &e, const T6 &f,
					 const T7 &g)
    {
      SliceType_t retval = Pooma::NoInit();
      return fillSlice(retval, u, a, b, c, d, e, f, g);
    }

  template<class RT, class UT>
  inline static RT &fillSlice(RT &retval, const UT &u,
			      const T1 &a, const T2 &b,
			      const T3 &c, const T4 &d,
			      const T5 &e, const T6 &f,
			      const T7 &g)
    {
      enum { RDX = 
        DomainTraits<RT>::dimensions > DomainTraits<RT>::sliceDimensions };
      CombineSliceDomain<RT,UT,T1,0,0,(DX1>0 && RDX)>::combine(retval,u,a);
      CombineSliceDomain<RT,UT,T2,S2,DX1,(DX2>0 && RDX)>::combine(retval,u,b);
      CombineSliceDomain<RT,UT,T3,S3,DX1+DX2,(DX3>0 && RDX)>::
        combine(retval,u,c);
      CombineSliceDomain<RT,UT,T4,S4,DX1+DX2+DX3,(DX4>0 && RDX)>::
        combine(retval,u,d);
      CombineSliceDomain<RT,UT,T5,S5,DX1+DX2+DX3+DX4,(DX5>0 && RDX)>::
        combine(retval,u,e);
      CombineSliceDomain<RT,UT,T6,S6,DX1+DX2+DX3+DX4+DX5,(DX6>0 && RDX)>::
        combine(retval,u,f);
      CombineSliceDomain<RT,UT,T7,S7,DX1+DX2+DX3+DX4+DX5+DX6,(DX7>0 && RDX)>::
        combine(retval,u,g);
      return retval;
    }
};

/** @addtogroup NewDomain
 * TemporaryNewDomain1 is a fix for a deficiency in NewDomain1, which
 * synthesizes SliceType_t from a single parameter, which is incorrect. We
 * need the array's domain type in order to get the right answer when the
 * user specifies a single AllDomain<Dim>.
 */

template<class Domain, class Sub>
struct TemporaryNewDomain1
{
  typedef typename NewDomain1<Sub>::SliceType_t SliceType_t;
  static inline
  SliceType_t combineSlice(const Domain &d, const Sub &s)
  {
    return NewDomain1<Sub>::combineSlice(d, s);
  }
};

template<class Domain, int N>
struct TemporaryNewDomain1<Domain, AllDomain<N> > 
{
  typedef Domain SliceType_t;
  static inline
  const SliceType_t &combineSlice(const Domain &d, const AllDomain<N> &)
  {
    return d;
  }
};

#endif     // POOMA_DOMAIN_NEWDOMAIN_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NewDomain.h,v $   $Author: richard $
// $Revision: 1.36 $   $Date: 2004/11/01 18:16:32 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
