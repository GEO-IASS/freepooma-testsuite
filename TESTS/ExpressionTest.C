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
// Test some basic expressions with different engines.
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Engine/UserFunction.h"
#include "Engine/Stencil.h"
#include "Tiny/Vector.h"
#include "Pooma/FunctorResult.h"
#include "Pooma/Fields.h"

#include <iostream>
#include <cmath>

#if POOMA_MESSAGING
  typedef DistributedTag LayoutTag_t;
  typedef Remote<Brick> BrickTag_t;
  typedef Remote<CompressibleBrick> CompBrickTag_t;
#else
  typedef ReplicatedTag LayoutTag_t;
  typedef Brick BrickTag_t;
  typedef CompressibleBrick CompBrickTag_t;
#endif

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
// These functions test various data-parallel statements against corresponding
// serial versions.  They're templated on the array type, so they can be used
// with different kinds of engines.
//
// Typically, the functions take 5 arrays.  The last array is used as an
// initial condition, and the first four are used as the left and right sides
// of the serial and data-parallel expressions.
//
//-----------------------------------------------------------------------------

struct Norm
{
  Norm() { }
  Norm(const Norm &) { }

  template<class T>
  inline
  T operator()(const T &a) const
  {
    return a*a;
  }

  template<int Dim, class T>
  inline
  T operator()(const Vector<Dim,T> &a) const
  {
    return dot(a,a);
  }
};

template<int Dim, class T>
struct FunctorResult<Norm,Vector<Dim,T> >
{
  typedef T Type_t;
};

template<class Array>
bool isSmall(const Array &a)
{
  //  UserFunction<Norm> norm;
  double epsilon = 0.000000001;
  return (sum(abs(a)) < epsilon);
}

template<class A2, class A4>
bool checkTest(Pooma::Tester& tester, int test,
	       const A2 &a2, const A4 &a4)
{
  Pooma::blockAndEvaluate();

  bool passed    = isSmall(a2 - a4);

  if (passed)
  {
    tester.out() << "Test #" << test << " passed." << std::endl;
  }
  else
  {
    tester.out() << "Test #" << test << " failed." << std::endl;
    tester.out() << "loop version:" << std::endl
		 << a2 << std::endl;
    tester.out() << "data-parallel version:" << std::endl
		 << a4 << std::endl;
  }

  return passed;
}

template<class A1,class A2,class A3,class A4, class AInit>
void test1(Pooma::Tester& tester, int test,
	   const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	   const AInit &initial, const Interval<1> &I)
{
  //---------------------------------------------------------------------
  // test #1
  // Simple expression with data-parallel stencil.

  int from = I.first();
  int to = I.last();
  int i;

  a1 = initial;
  a2 = initial;
  a3 = initial;
  a4 = initial;

  Pooma::blockAndEvaluate();

  for (i = from; i <= to; ++i)
  {
    a2(i) = initial(i) + a1(i-1) + a1(i);
  }

  a4(I) = initial(I) + a3(I-1) + a3(I);

  Pooma::blockAndEvaluate();

  tester.check(checkTest(tester,test,a2,a4));
}

class CosTimes
{
public:
  CosTimes(double x) : x_m(x) {}
  double operator()(double y) const { return cos(x_m*y); }

  CosTimes() : x_m(0.0) { }
  CosTimes(const CosTimes & a) : x_m(a.x_m) { }

private:
  double x_m;
};

template<class A1,class A2,class A3,class A4, class AInit>
void test2(Pooma::Tester& tester, int test,
	   const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	   const AInit &initial, const Interval<1> &I)
{
  //---------------------------------------------------------------------
  // test #2
  // UserFunction engine

  int from = I.first();
  int to = I.last();
  int i;

  UserFunction<CosTimes> cosTimes(0.15);

  a1 = initial;
  a2 = initial;
  a3 = initial;
  a4 = initial;

  Pooma::blockAndEvaluate();

  for (i = from;i <= to; ++i)
  {
    a2(i) = initial(i) + cos(0.15*(a1(i-1) + a1(i)));
  }

  a4(I) = initial(I) + cosTimes(a3(I-1) + a3(I));

  Pooma::blockAndEvaluate();

  tester.check(checkTest(tester, test, a2, a4));
}

template <int Dim>
struct TwoPt
{
  typedef double OutputElement_t;
  TwoPt() { }
  TwoPt(const TwoPt &m) : inputCentering_m(m.inputCentering_m) { }
  template <class FE>
  TwoPt(const FE& fe)
  {
    inputCentering_m = fe.centering();
  }

  Centering<Dim> outputCentering() const
  {
    return inputCentering_m;
  }
  Centering<Dim> inputCentering() const
  {
    return inputCentering_m;
  }

  template <class A>
  inline
  typename A::Element_t
  operator()(const A& x, int i) const
  {
    return ( x.read(i-1) + x.read(i) );
  }

  inline int lowerExtent(int) const { return 1; }
  inline int upperExtent(int) const { return 0; }

  Centering<Dim> inputCentering_m;
};

template <class F, class Dom>
typename FieldStencilSimple<TwoPt<F::dimensions>, F>::Type_t
twoPt(const F& expr, const Dom &domain)
{
  typedef FieldStencilSimple<TwoPt<F::dimensions>, F> Ret_t;
  return Ret_t::make(TwoPt<F::dimensions>(expr), expr, domain);
}

template<class A1,class A2,class A3,class A4, class AInit>
void test3(Pooma::Tester& tester, int test,
	   const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	   const AInit &initial, const Interval<1> &I)
{
  //---------------------------------------------------------------------
  // test #3
  // stencil engine

  int from = I.first();
  int to = I.last();
  int i;

  a1 = initial;
  a2 = initial;
  a3 = initial;
  a4 = initial;

  Pooma::blockAndEvaluate();

  for (i = from;i <= to; ++i)
  {
    a2(i) = initial(i) + a1(i-1) + a1(i);
  }

  a4(I) = initial(I) + twoPt(a3, I);

  Pooma::blockAndEvaluate();

  tester.check(checkTest(tester, test, a2, a4));
}

template<class A1,class A2,class A3,class A4, class AInit>
void test4(Pooma::Tester& tester, int test,
	   const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	   const AInit &initial, const Interval<1> &I)
{
  //---------------------------------------------------------------------
  // test #4
  // expression inside stencil

  int from = I.first();
  int to = I.last();
  int i;

  a1 = initial;
  a2 = initial;
  a3 = initial;
  a4 = initial;

  Pooma::blockAndEvaluate();

  for (i = from; i <= to; ++i)
  {
    a2(i) = initial(i) + 1.0 + a1(i-1) + 1.0 + a1(i);
  }

  a4(I) = initial(I) + twoPt(1.0 + a3, I);

  Pooma::blockAndEvaluate();

  tester.check(checkTest(tester, test, a2, a4));
}

template<class A1,class A2,class A3,class A4, class AInit>
void test5(Pooma::Tester& tester, int test,
	   const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	   const AInit &initial, const Interval<1> &I)
{
  //---------------------------------------------------------------------
  // test #5
  // component forward + user function

  int from = I.first();
  int to = I.last();
  int i;

  UserFunction<CosTimes> cosTimes(0.15);

  a1 = initial;
  a2 = initial;
  a3 = initial;
  a4 = initial;

  Pooma::blockAndEvaluate();

  for (i = from; i <= to; ++i)
  {
    a2(i)(1) = initial(i)(1) + cos(0.15 * (a1(i - 1)(1)));
  }

  a4.comp(1)(I) = initial.comp(1)(I) + cosTimes(a3.comp(1)(I - 1));

  Pooma::blockAndEvaluate();

  tester.check(checkTest(tester, test, a2.comp(1), a4.comp(1)));
}

template<class A1,class A2,class A3,class A4, class AInit>
void test6(Pooma::Tester& tester, int test,
	   const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	   const AInit &initial, const Interval<1> &I)
{
  //---------------------------------------------------------------------
  // test #6
  // component forward

  int from = I.first();
  int to = I.last();
  int i;

  a1 = initial;
  a2 = initial;
  a3 = initial;
  a4 = initial;

  Pooma::blockAndEvaluate();

  for (i = from; i <= to; ++i)
  {
    a2(i)(1) = initial(i)(1) + a1(i - 1)(1);
  }

  a4.comp(1)(I) = initial.comp(1)(I) + a3.comp(1)(I - 1);

  Pooma::blockAndEvaluate();

  tester.check(checkTest(tester, test, a2.comp(1), a4.comp(1)));
}

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Interval<1> physicalVertexDomain(10);
  Interval<1> cellInterior(2,7);
  
  Loc<1> blocks(2);
  UniformGridPartition<1> partition(blocks, GuardLayers<1>(1));   
  UniformGridLayout<1> layout(physicalVertexDomain, partition,
			      LayoutTag_t());
  
  tester.out() << layout << std::endl;
  tester.out() << layout.domain() << std::endl;
  
  // Now, we can declare a field.

  Centering<1> cell = canonicalCentering<1>(CellType, Continuous);

  typedef 
    Field<UniformRectilinearMesh<1>, double, MultiPatch<UniformTag, BrickTag_t> >
    Field_t;

  Field_t initial(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  Field_t a1(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  Field_t a2(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  Field_t a3(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  Field_t a4(cell, layout, Vector<1>(0.0), Vector<1>(1.0));

  tester.out() << "physical domain:" << a1.physicalDomain() << std::endl;
  tester.out() << "interior domain:" << cellInterior << std::endl;

  initial = 0.0;

  Pooma::blockAndEvaluate();

  initial(4) = 2.0;
  initial(7) = 3.0;

  test1(tester, 1, a1, a2, a3, a4, initial, cellInterior);
  // No UserFunction for Field yet.
  //  test2(tester, 2, a1, a2, a3, a4, initial, cellInterior);

  // Need to replace the stencil code above with Field Stencil code.
  test3(tester, 3, a1, a2, a3, a4, initial, cellInterior);
  //test4(tester, 4, a1, a2, a3, a4, initial, cellInterior);

  typedef 
    Field<UniformRectilinearMesh<1>, double, MultiPatch<UniformTag,
                                                    CompBrickTag_t> >
    CompressibleField_t;

  CompressibleField_t cinit(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  CompressibleField_t ca1(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  CompressibleField_t ca2(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  CompressibleField_t ca3(cell, layout, Vector<1>(0.0), Vector<1>(1.0));
  CompressibleField_t ca4(cell, layout, Vector<1>(0.0), Vector<1>(1.0));

  cinit = 0.0;

  Pooma::blockAndEvaluate();

  cinit(4) = 2.0;
  cinit(7) = 3.0;

  test1(tester, 1, ca1, ca2, ca3, ca4, cinit, cellInterior);
  //  test2(tester, 2, ca1, ca2, ca3, ca4, cinit, cellInterior);
  test3(tester, 3, ca1, ca2, ca3, ca4, cinit, cellInterior);
  //test4(tester, 4, ca1, ca2, ca3, ca4, cinit, cellInterior);


  int ret = tester.results("ExpressionTest");
  Pooma::finalize();
  return ret;
}




// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ExpressionTest.cpp,v $   $Author: richi $
// $Revision: 1.6 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
