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
// Tiny operations test
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Tiny/TinyMatrix.h"
#include "Tiny/Vector.h"
#include "Tiny/VectorTinyMatrix.h"

const int D=3;
const int N=100;
Pooma::Tester *tester;

template<int D1, int D2, class T, class E>
void fill( TinyMatrix<D1,D2,T,E>& x )
{
  for (int j=0; j<D1; ++j)
    for (int k=0; k<D2; ++k)
      x(j,k) = rand() % 1000000;
}

template<int D1, class T, class E>
void fill( Vector<D1,T,E>& x)
{
  for (int j=0; j<D1; ++j)
    x(j) = rand() % 1000000;
}

void testTinyMatrixDot()
{
  TinyMatrix<D,D> x,y;
  int i,j,k,n;

  double a = 0;
  for (n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      TinyMatrix<D,D> z = dot(x,y);
      for (i=0; i<D; ++i)
        for (j=0; j<D; ++j)
          for (k=0; k<D; ++k)
            z(i,k) -= x(i,j)*y(j,k);
      z = z*z;
      for (i=0; i<D; ++i)
        for (j=0; j<D; ++j)
          a += z(i,j);
    }

  tester->check("TinyMatrix dot", a == 0);
}

void testVectorDot()
{
  Vector<D> x,y;

  double a = 0;
  for (int n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      a += dot(x,y);
      for (int i=0; i<D; ++i)
        a -= x(i)*y(i);
    }

  tester->check("Vector dot", a == 0);
}

void testTinyMatrixEquality()
{
  TinyMatrix<D,D> x,y;

  bool eq, ans;  
  for (int n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      eq = x == y;
      ans = true;
      for (int j=0; j<D; ++j)
        for (int i=0; i<D; ++i)
        ans = ans && (x(i,j) == y(i,j));
    }

  tester->check("TinyMatrix equality", ans == eq);
}

void testVectorEquality()
{
  Vector<D> x,y;

  bool eq, ans;  
  for (int n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      eq = x == y;
      ans = true;
      for (int i=0; i<D; ++i)
        ans = ans && (x(i) == y(i));
    }

  tester->check("Vector equality", ans == eq);
}

void testVectorAdd()
{
  Vector<D> x,y;
  double a=0;

  for (int n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      Vector<D> b = x + y;
      for (int i=0; i<D; ++i)
        a += b(i) - (x(i)+y(i));
    }

  tester->check("Vector add", a == 0);
}

void testVectorScalar()
{
  Vector<D> x;
  double a=0;
  int i;

  for (int n=0; n<N; ++n)
    {
      fill(x);
      Vector<D> b = x + 1.0;
      for (i=0; i<D; ++i)
        a += b(i) - (x(i)+1.0);
      b = 1.0 + x;
      for (i=0; i<D; ++i)
        a += b(i) - (1.0+x(i));
    }
 
  tester->check("Vector scalar", a == 0);
}

void testTinyMatrixAdd()
{
  TinyMatrix<D,D> x,y;
  double a=0;

  for (int n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      TinyMatrix<D,D> b = x + y;
      for (int i=0; i<D; ++i)
        for (int j=0; j<D; ++j)
          a += b(i,j) - (x(i,j)+y(i,j));
    }

  tester->check("TinyMatrix add", a == 0);
}

void testTinyMatrixScalar()
{
  TinyMatrix<D,D> x;
  double a=0;
  int i,j;

  for (int n=0; n<N; ++n)
    {
      fill(x);
      TinyMatrix<D,D> b = x + 1.0;
      for (i=0; i<D; ++i)
        for (j=0; j<D; ++j)
          a += b(i,j) - (x(i,j)+1.0);
      b = 1.0 + x;
      for (i=0; i<D; ++i)
        for (j=0; j<D; ++j)
          a += b(i,j) - (1.0+x(i,j));
    }

  tester->check("TinyMatrix scalar", a == 0);
}

void testVectorNegate()
{
  Vector<D> x;
  double a=0;

  for (int n=0; n<N; ++n)
    {
      fill(x);
      Vector<D> b = -x;
      for (int i=0; i<D; ++i)
        a += b(i) + x(i);
    }

  tester->check("Vector negate", a == 0);
}

void testTinyMatrixNegate()
{
  TinyMatrix<D,D> x;
  double a=0;

  for (int n=0; n<N; ++n)
    {
      fill(x);
      TinyMatrix<D,D> b = -x;
      for (int i=0; i<D; ++i)
        for (int j=0; j<D; ++j)
          a += b(i,j) + x(i,j);
    }

  tester->check("TinyMatrix negate", a == 0);
}

void testVectorDotTinyMatrix()
{
  Vector<D> x;
  TinyMatrix<D,D> y;
  int i,j,n;

  double a = 0;
  for (n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      Vector<D> z = dot(x,y);
      for (i=0; i<D; ++i)
        for (j=0; j<D; ++j)
          z(j) -= x(i)*y(i,j);
      for (i=0; i<D; ++i)
        a += z(i)*z(i);
    }

  tester->check("Vector dot TinyMatrix", a == 0);
}

void testTinyMatrixDotVector()
{
  Vector<D> x;
  TinyMatrix<D,D> y;
  int i,j,n;

  double a = 0;
  for (n=0; n<N; ++n)
    {
      fill(x);
      fill(y);
      Vector<D> z = dot(y,x);
      for (i=0; i<D; ++i)
        for (j=0; j<D; ++j)
          z(i) -= y(i,j)*x(j);
      for (i=0; i<D; ++i)
        a += z(i)*z(i);
    }

  tester->check("TinyMatrix dot Vector", a == 0);
}

void testTinyMatrixDot2()
{
  const int D1 = 3;
  const int D2 = 2;
  const int D3 = 4;
  TinyMatrix<D1,D2> t1;
  TinyMatrix<D2,D3> t2;
  TinyMatrix<D1,D3> t3;

  double a = 0;
  for (int i=0; i<N; ++i)
    {
      fill(t1);
      fill(t2);
      t3 = dot(t1,t2);
      for (int i1=0; i1<D1; ++i1)
        for (int i3=0; i3<D3; ++i3)
          {
            double x = t3(i1,i3);
            for (int i2=0; i2<D2; ++i2)
              x -= t1(i1,i2)*t2(i2,i3);
            a += x*x;
          }
    }

  tester->check("TinyMatrix<3,2> dot TinyMatrix<2,4>", a == 0);
}

void testTinyMatrixDotVector2()
{
  const int D1 = 3;
  const int D2 = 2;
  TinyMatrix<D1,D2> t1;
  Vector<D1> v1;
  Vector<D2> v2;

  double a = 0;
  for (int i=0; i<N; ++i)
    {
      fill(t1);
      fill(v2);
      v1 = dot(t1,v2);
      for (int i1=0; i1<D1; ++i1)
        {
          double x = v1(i1);
          for (int i2=0; i2<D2; ++i2)
            x -= t1(i1,i2)*v2(i2);
          a += x*x;
        }
    }

  tester->check("TinyMatrix<3,2> dot Vector<2>", a == 0);
}

void testVectorDotTinyMatrix2()
{
  const int D1 = 3;
  const int D2 = 2;
  TinyMatrix<D1,D2> t1;
  Vector<D1> v1;
  Vector<D2> v2;

  double a = 0;
  for (int i=0; i<N; ++i)
    {
      fill(t1);
      fill(v1);
      v2 = dot(v1,t1);
      for (int i2=0; i2<D2; ++i2)
        {
          double x = v2(i2);
          for (int i1=0; i1<D1; ++i1)
            x -= v1(i1)*t1(i1,i2);
          a += x*x;
        }
    }

  tester->check("Vector<3> dot TinyMatrix<3,2>", a == 0);
}

void testVectorAccum()
{
  Vector<D> v1,v2,v3;

  double a = 0;
  for (int i=0; i<N; ++i)
    {
      fill(v1);
      fill(v2);
      v3 = v1;
      v1 += v2;
      for (int j=0; j<D; ++j)
        {
          double x = v1(j) - ( v3(j) + v2(j) );
          a += x*x;
        }

      fill(v1);
      v2 = v1;
      v1 += 73;
      for (int k=0; k<D; ++k)
        {
          double x = v1(k) - ( v2(k) + 73 );
          a += x*x;
        }
    }

  tester->check("Vector accum", a == 0);
}

void testTinyMatrixAccum()
{
  TinyMatrix<D,D> v1,v2,v3;
  int i,j,k;

  double a = 0;
  for (i=0; i<N; ++i)
    {
      fill(v1);
      fill(v2);
      v3 = v1;
      v1 += v2;
      for (j=0; j<D; ++j)
        for (k=0; k<D; ++k)
          {
            double x = v1(j,k) - ( v3(j,k) + v2(j,k) );
            a += x*x;
          }

      fill(v1);
      v2 = v1;
      v1 += 73;
      for (j=0; j<D; ++j)
        for (k=0; k<D; ++k)
          {
            double x = v1(j,k) - ( v2(j,k) + 73 );
            a += x*x;
          }
    }

  tester->check("TinyMatrix accum", a == 0);
}

void testNorm()
{
  const int D=3;
  const int N=100;
  Vector<D> x[N];

  for (int i=0; i<N; ++i)
    for (int j=0; j<D; ++j)
      x[i](j) = rand();

  int i;
  for (i=0; i<N; ++i)
    {
      double n0 = norm2(x[i]);
      double n1 = norm(x[i]);
      double n2 = 0;
      for (int j=0; j<D; ++j)
	n2 += x[i](j)*x[i](j);
      if ( fabs(sqrt(n2) - n1)/n1 > 1e-12 )
	break;
      if ( fabs(n2 - n0)/n0 > 1e-12 )
	break;
    }
  tester->check("norm/norm2", i == N);
}

#if POOMA_EXCEPTIONS
#undef POOMA_BOUNDS_CHECK
#define POOMA_BOUNDS_CHECK POOMA_YES
void testBoundsChecking()
{
  int ecnt = 0;
  try 
    {
      Vector<3, double> v;
      v(0) = 0.0; v(1) = 2.0; v(2) = -4.0;
      v(6) = 1.3;
    }
  catch (Pooma::Assertion &a)
    {
      a.print(tester->out());
      tester->out() << std::endl;
      ecnt++;
    }
  try 
    {
      TinyMatrix<2, 3, double> t;
      t(0, 0) = 0.0; t(0, 1) = 2.0; t(0, 2) = -4.0;
      t(1, 0) = 0.2; t(1, 1) = 2.6; t(1, 2) = -0.4;
      t(-1,0) = 1.3;
    }
  catch (Pooma::Assertion &a)
    {
      a.print(tester->out());
      tester->out() << std::endl;
      ecnt++;
    }

  tester->check("bounds checking", ecnt == 2);
}
#endif

int main(int argc, char **argv)
{
  Pooma::initialize(argc, argv);
  tester = new Pooma::Tester(argc, argv);

  testTinyMatrixDot();
  testVectorDot();
  testTinyMatrixEquality();
  testVectorEquality();
  testTinyMatrixAdd();
  testVectorAdd();
  testTinyMatrixNegate();
  testVectorNegate();
  testVectorScalar();
  testTinyMatrixScalar();
  testVectorDotTinyMatrix();
  testTinyMatrixDotVector();
  testTinyMatrixDot2();
  testTinyMatrixDotVector2();
  testVectorDotTinyMatrix2();
  testVectorAccum();
  testTinyMatrixAccum();
  testNorm();
#if POOMA_EXCEPTIONS
  testBoundsChecking();
#endif

  int ret = tester->results("t1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: t1.cpp,v $   $Author: richard $
// $Revision: 1.12 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
