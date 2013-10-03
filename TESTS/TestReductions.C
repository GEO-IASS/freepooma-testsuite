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
// Tiny Reduction Tests.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Tiny.h"
#include "Utilities/Tester.h"

// Answers: sum, prod, min, max, all, any, bitOr, bitAnd.

static int vans1[] = { 2, 2, 2, 2, 1, 1, 2, 2 };
static int vans2[] = { 3, 2, 1, 2, 1, 1, 3, 0 };
static int vans3[] = { 3, 0, 0, 2, 0, 1, 3, 0 };

static int taans1[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
static int taans2[] = { 0, 0, -2, 2, 0, 1, (2 | -2), 0 };
static int taans3[] = { 0, 0, -2, 2, 0, 1, (2 | 1 | -1 | -2), 0 };

static int tdans1[] = { 2, 2, 2, 2, 1, 1, 2, 2 };
static int tdans2[] = { 3, 0, 0, 2, 0, 1, 3, 0 };
static int tdans3[] = { 3, 0, 0, 2, 0, 1, 3, 0 };

static int tfans1[] = { 2, 2, 2, 2, 1, 1, 2, 2 };
static int tfans2[] = { 6, 4, 1, 2, 1, 1, 3, 0 };
static int tfans3[] = { 9, 0, 0, 2, 0, 1, 3, 0 };

static int tsans1[] = { 2, 2, 2, 2, 1, 1, 2, 2 };
static int tsans2[] = { 6, 4, 1, 2, 1, 1, 3, 0 };
static int tsans3[] = { 10, 0, 0, 2, 0, 1, 3, 0 };

static int tmans1[] = { 3, 2, 1, 2, 1, 1, 3, 0 };
static int tmans2[] = { 6, 0, 0, 2, 0, 1, 3, 0 };
static int tmans3[] = { 9, 8, 1, 2, 1, 1, 3, 0 };

template<int Dim>
void initialize(Vector<Dim, int> &v, Pooma::Tester &tester)
{
  for (int i = 0, k = 2; i < Dim; i++, k--)
    v(i) = k;
  tester.out() << v << std::endl;
}

template<int Dim>
void initialize(Tensor<Dim, int, Antisymmetric> &t, Pooma::Tester &tester)
{
  for (int i = 0; i < TensorStorageSize<Dim, Antisymmetric>::Size; i++)
    t(i) = 2 - (i % Dim);
  tester.out() << t << std::endl;
}

template<int Dim>
void initialize(Tensor<Dim, int, Diagonal> &t, Pooma::Tester &tester)
{
  for (int i = 0, k = 2; i < Dim; i++, k--)
    t(i, i) = k;
  tester.out() << t << std::endl;
}

template<int Dim>
void initialize(Tensor<Dim, int, Full> &t, Pooma::Tester &tester)
{
  for (int i = 0; i < Dim; i++)
    for (int j = 0, k = 2; j < Dim; j++, k--)
      t(i, (j + i) % Dim) = k;
  tester.out() << t << std::endl;
}

template<int Dim>
void initialize(Tensor<Dim, int, Symmetric> &t, Pooma::Tester &tester)
{
  for (int i = 0; i < TensorStorageSize<Dim, Symmetric>::Size; i++)
    t(i) = 2 - (i % Dim);
  tester.out() << t << std::endl;
}

template<int Dim1, int Dim2>
void initialize(TinyMatrix<Dim1, Dim2, int> &m, Pooma::Tester &tester)
{
  for (int i = 0; i < Dim1 * Dim2; i++)
    m(i) = 2 - (i % Dim2);
  tester.out() << m << std::endl;
}

template<class TinyObject>
void test(const TinyObject &o, Pooma::Tester &tester, 
  const int *answers)
{
  tester.check("sum", sum(o), answers[0]);
  tester.check("prod", prod(o), answers[1]);
  tester.check("min", min(o), answers[2]);
  tester.check("max", max(o), answers[3]);
  tester.check("all", all(o), answers[4] == 1);
  tester.check("any", any(o), answers[5] == 1);
  tester.check("bitOr", bitOr(o), answers[6]);
  tester.check("bitAnd", bitAnd(o), answers[7]);
}

void testVectors(Pooma::Tester &tester)
{
  tester.out() << "Vector tests:" << std::endl;
  
  Vector<1, int> v1; 
  Vector<2, int> v2; 
  Vector<3, int> v3;

  tester.out() << "1D" << std::endl;
  initialize(v1, tester); test(v1, tester, vans1);

  tester.out() << "2D" << std::endl;
  initialize(v2, tester); test(v2, tester, vans2);

  tester.out() << "3D" << std::endl;
  initialize(v3, tester); test(v3, tester, vans3);
}

template<class EngineTag>
void testTensors(const char *tag, Pooma::Tester &tester, 
  const int *ans1, const int *ans2, const int *ans3)
{
  tester.out() << tag << " Tensor tests:" << std::endl;
  
  Tensor<1, int, EngineTag> t1; 
  Tensor<2, int, EngineTag> t2; 
  Tensor<3, int, EngineTag> t3;

  tester.out() << "1D" << std::endl;
  initialize(t1, tester); test(t1, tester, ans1);

  tester.out() << "2D" << std::endl;
  initialize(t2, tester); test(t2, tester, ans2);

  tester.out() << "3D" << std::endl;
  initialize(t3, tester); test(t3, tester, ans3);
}

void testTinyMatrices(Pooma::Tester &tester)
{
  tester.out() << "TinyMatrix tests:" << std::endl;
  
  TinyMatrix<1, 2, int> m1; 
  TinyMatrix<2, 3, int> m2; 
  TinyMatrix<3, 2, int> m3;

  tester.out() << "1 x 2" << std::endl;
  initialize(m1, tester); test(m1, tester, tmans1);

  tester.out() << "2 x 3" << std::endl;
  initialize(m2, tester); test(m2, tester, tmans2);

  tester.out() << "3 x 2" << std::endl;
  initialize(m3, tester); test(m3, tester, tmans3);
}

int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  testVectors(tester);

  testTensors<Antisymmetric>("Antisymmetric", tester, taans1, taans2, taans3);
  testTensors<Diagonal>("Diagonal", tester, tdans1, tdans2, tdans3);
  testTensors<Full>("Full", tester, tfans1, tfans2, tfans3);
  testTensors<Symmetric>("Symmetric", tester, tsans1, tsans2, tsans3);

  testTinyMatrices(tester);

  int retval = tester.results("TestReductions");

  Pooma::finalize();

  return retval;  
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestReductions.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
