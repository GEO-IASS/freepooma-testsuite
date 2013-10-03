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
// array_test19.cpp some tests to verify correctnes of stencil objects
//-----------------------------------------------------------------------------

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Partition/UniformGridPartition.h"
#include "Layout/UniformGridLayout.h"
#include "Layout/GuardLayers.h"
#include "Engine/BrickEngine.h"
#include "Engine/MultiPatchEngine.h"
#include "Engine/Stencil.h"
#include "Array/Array.h"
#include "Tiny/Vector.h"
#include "Pooma/FunctorResult.h"

#include <iostream>
#include <complex>
#include <cmath>

template<class Array>
bool isSmall(const Array &a)
{
  double epsilon = 0.000000001;
  return (sum(a*a) < epsilon);
}

template<class Array, class Tester>
void checkArray(const Array &a, Tester &tester, char *comment)
{
  bool ok = isSmall(a);
  tester.check(ok);
  if (!ok)
  {
    tester.out() << "Failure from:" << comment << std::endl;
    tester.out() << a << std::endl;
  }
}

template<class Tester>
void checkFlag(bool ok, Tester &tester, char *comment)
{
  tester.check(ok);
  if (!ok)
  {
    tester.out() << "Failure from:" << comment << std::endl;
  }
}


class TwoPoint
{
public:
  TwoPoint() {}

  template <class A>
  inline
  typename A::Element_t
  operator()(const A& x, int i) const
  {
    return ( x(i-1) + x(i) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 0; }

private:
};

class ThreePoint
{
public:
  ThreePoint() {}

  template <class A>
  inline
  typename A::Element_t
  operator()(const A& x, int i) const
  {
    return ( x(i-1) + x(i) + x(i+1) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 1; }

private:
};

template<class T>
struct NormResult
{ };

template<int D, class T>
struct NormResult<Vector<D,T> >
{
  typedef T Type_t;
};

class NormThing
{
public:
  NormThing() {}

  template <class A>
  inline
  typename NormResult<typename A::Element_t>::Type_t
  operator()(const A& x, int i) const
  {
    return ( 0.5*(dot(x(i-1), x(i)) + dot(x(i), x(i+1))) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 1; }

private:
};

// To apply stencils that return a different type than
// they input, we must tell pooma the return type using
// FunctorResult.  (It's analagous to the STL result_type
// in functors.)

template<class T>
struct FunctorResult<NormThing,T>
{
  typedef typename NormResult<T>::Type_t Type_t;
};

class AsymDoof
{
public:
  AsymDoof() {}

  template <class A>
  inline
  typename A::Element_t
  operator()(const A& x, int i, int j) const
  {
    return ( (1.0/15.0) *
             ( x(i+1,j+1) + 2*x(i+1,j  ) + 3*x(i+1,j-1) +
               3*x(i  ,j+1) + x(i  ,j  ) + 2*x(i  ,j-1) +
               4*x(i-1,j+1) + 3*x(i-1,j  ) + 5*x(i-1,j-1) ) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 1; }  

private:
};

class DoofNinePt
{
public:
  DoofNinePt() {}

  template <class A>
  inline
  typename A::Element_t
  operator()(const A& x, int i, int j) const
  {
    return ( (1.0/9.0) *
             ( x(i+1,j+1) + x(i+1,j  ) + x(i+1,j-1) +
               x(i  ,j+1) + x(i  ,j  ) + x(i  ,j-1) +
               x(i-1,j+1) + x(i-1,j  ) + x(i-1,j-1) ) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 1; }  

private:
};

class DoofLower
{
public:
  DoofLower(complex<double> alpha)
    : alpha_m(alpha), alphaConj_m(conj(alpha))
  {}

  template <class A>
  inline
  complex<double>
  operator()(const A& x, int i, int j) const
  {
    return ( (1.0/3.0) *
             ( x(i  ,j  ) + alpha_m * x(i  ,j-1) +
               alphaConj_m * x(i-1,j  ) + x(i-1,j-1) ) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 0; }  

private:
  complex<double> alpha_m;
  complex<double> alphaConj_m;
};

class DoofUpper
{
public:
  DoofUpper(complex<double> alpha)
    : alpha_m(alpha), alphaConj_m(conj(alpha))
  {}

  template <class A>
  inline
  complex<double>
  operator()(const A& x, int i, int j) const
  {
    return ( (1.0/3.0) *
             ( x.read(i  ,j  ) + alpha_m * x.read(i  ,j+1) +
               alphaConj_m * x.read(i+1,j  ) + x.read(i+1,j+1) ) );
  }

  inline int lowerExtent(int) const { return 0; }
  inline int upperExtent(int) const { return 1; }  

private:
  complex<double> alpha_m;
  complex<double> alphaConj_m;
};

template<class T>
struct FunctorResult<DoofLower,T>
{
  typedef complex<double> Type_t;
};

template<class T>
struct FunctorResult<DoofUpper,T>
{
  typedef complex<double> Type_t;
};

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << argv[0] << ": Tests of stencil objects on arrays."
	       << std::endl;
  tester.out() << "------------------------------------------------"
	       << std::endl;

  Interval<1> i1(100),i2(200,299),i3(51,151);
  Range<1> r1(10,80,10);

  Array<1> initial(i1);

  initial = 1.0;

  Pooma::blockAndEvaluate();

  initial(4) = 2.0;
  initial(23) = 3.0;
  initial(52) = 4.0;
  initial(r1) = 5.0;
  initial(r1+3) = 6.0;

  Stencil<TwoPoint>   twoPoint;
  Stencil<ThreePoint> threePoint;

  Interval<1> out = threePoint.insetDomain(i1);

  Array<1> a1(i1),a2(i1),b1(i2),b2(i1);

  a1 = initial;
  b1 = initial;
  a2 = 0.0;
  b2 = 0.0;
  a2(out) = threePoint(a1);
  b2(out) = threePoint(b1);

  tester.check( isSmall(a2 - b2) );
  tester.out() << a2 << std::endl;
  tester.out() << b2 << std::endl;
  tester.out() << "Test #1:"
	       << (isSmall(a2 - b2) ? "passed" : "failed") << std::endl;

  Interval<1> view(20,40),v2(219,241);
  a2 = 0.0;
  b2 = 0.0;
  a2(view) = threePoint(a1,view);
  b2(view) = threePoint(b1(v2));

  tester.out() << a2 << std::endl;
  tester.out() << b2 << std::endl;
  checkArray(a2 - b2, tester, "test #2");
  //  tester.check( isSmall(a2 - b2) );
  //  tester.out() << "Test #2:"
  //	       << (isSmall(a2 - b2) ? "passed" : "failed") << std::endl;

  a2 = 0.0;
  b2 = 0.0;
  a2(r1) = threePoint(a1,r1);
  b2(out) = threePoint(b1);

  tester.check( isSmall( (a2 - b2)(r1) ) );
  tester.out() << a2(r1) << std::endl;
  tester.out() << b2(r1) << std::endl;
  tester.out() << "Test #3:"
	       << (isSmall(a2(r1) - b2(r1)) ? "passed" : "failed")
	       << std::endl;

  Array<1,Vector<3,double> > v(i1);
  Stencil<NormThing> normThing;

  v  = 1.0;

  b2 = 0.0;
  b2(out) = normThing(v);

  Pooma::blockAndEvaluate();

  double check1 = 0.5*(dot(v(23), v(24)) + dot(v(24), v(25)));

  checkFlag(check1 == b2(24), tester, "stencil with different return type");

  complex<double> alpha(0.5,0.5*sqrt(3.0));

  // Test of stencil of stencil.
  // The stencils DoofLower and DoofUpper are 2x2 stencil that
  // can be composed to form the 3x3 DoofNinePt stencil.

  Stencil<DoofNinePt> doof;
  Stencil<DoofLower> doofL(alpha);
  Stencil<DoofUpper> doofU(alpha);

  Interval<1> inew(10);
  Interval<2> d2(inew,inew);
  Interval<2> inset = doof.insetDomain(d2);

  Array<2,double,Brick> init(d2),version1(d2),version2(d2);

  init = 0.0;
  version1 = 0.0;
  version2 = 0.0;

  Pooma::blockAndEvaluate();

  init(3,3) = 2.0;

  version1(inset) = doof(init);
  version2(inset) = real(doofU(doofL(init)));

  checkArray(version1 - version2, tester, "stencil of stencil");
  checkArray(imag(doofU(doofL(init))), tester, "imag");

  // Now some tests of views of 2D stencils.

  Stencil<AsymDoof> doofA;

  Interval<1> isub(2,5);
  Interval<2> d3(isub,isub);

  version1(inset) = doofA(init);

  Array<2,double,Brick> v3(d3),v4(d3);

  v3 = doofA(init)(d3);
  v4 = version1(inset)(d3);

  tester.out() << v3 << std::endl;
  tester.out() << v4 << std::endl;

  checkArray(v3 - v4, tester, "interval view");

  Range<1> rsub(1,7,2);
  Range<2> r3(rsub,rsub);

  v3 = doofA(init)(r3);
  v4 = version1(inset)(r3);

  tester.out() << v3 << std::endl;
  tester.out() << v4 << std::endl;

  checkArray(v3 - v4, tester, "range view");

  Array<1,double,Brick> v5(isub),v6(isub);

  v5 = doofA(init)(rsub,2);
  v6 = version1(inset)(rsub,2);

  tester.out() << v5 << std::endl;
  tester.out() << v6 << std::endl;

  checkArray(v5 - v6, tester, "slice view");

  Interval<1> i4(0,1),i5(4,5),i6(2,3);
  Range<1> r4(0,2,2);

  v5(i5) = doofA(init)(rsub,2)(i4);
  v6(i5) = version1(inset)(rsub,2)(i4);

  tester.out() << v5 << std::endl;
  tester.out() << v6 << std::endl;

  checkArray(v5 - v6, tester, "view of slice view");

  v5(i6) = doofA(init)(rsub,2)(r4);
  v6(i6) = version1(inset)(rsub,2)(r4);

  tester.out() << v5 << std::endl;
  tester.out() << v6 << std::endl;

  checkArray(v5 - v6, tester, "range view of slice view");

  UniformGridPartition<2> partition(Loc<2>(2,2), GuardLayers<2>(1));
  UniformGridLayout<2>    layout(d2, partition,ReplicatedTag());

  Array<2,double,MultiPatch<UniformTag,Brick> > v7(layout), initm(layout);

  initm = 0.0;
  initm(d2) = init;

  v7 = 0.0;
  v7(d2) = doofA(initm);

  checkArray(v7(d2) - version1, tester, "multipatch stencil");

  v3 = doofA(initm,inset)(r3);
  v4 = version1(inset)(r3);

  tester.out() << v3 << std::endl;
  tester.out() << v4 << std::endl;

  checkArray(v3 - v4, tester, "range view of multipatch stencil");

  v5 = doofA(initm,inset)(rsub,2);
  v6 = version1(inset)(rsub,2);

  tester.out() << v5 << std::endl;
  tester.out() << v6 << std::endl;

  checkArray(v5 - v6, tester, "slice view of multipatch stencil");

  tester.out() << "------------------------------------------------"
	       << std::endl;

  int retval = tester.results("array_test19");

  Pooma::finalize();

  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test19.cpp,v $   $Author: richi $
// $Revision: 1.22 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
