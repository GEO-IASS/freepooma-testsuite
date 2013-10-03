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

// ----------------------------------------------------------------------------
// Various tests of Tensor<D,double,Diagonal>
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"

//-----------------------------------------------------------------------------
// Main program:
//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc, argv);

  // --------------------------------------------------------------------------
  // 3D 
  // --------------------------------------------------------------------------
  Tensor<3,double,Full> t3f1(0.0, 3.0, 6.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0);
  tester.out() << "t3f1: " << t3f1 << std::endl;
  Tensor<3,double,Full> t3f2 = -t3f1;
  tester.out() << "t3f2: " << t3f2 << std::endl;

  Tensor<3,double,Diagonal> t3d1(1.0, 2.0, 3.0);
  tester.out() << "t3d1: " << t3d1 << std::endl;
  Tensor<3,double,Diagonal> t3d2(-1.0, -2.0, -3.0);
  tester.out() << "t3d2: " << t3d2 << std::endl;

  Tensor<3,double,Full> t3d1AsFull(1.0,0.0,0.0, 0.0,2.0,0.0, 0.0,0.0,3.0);
  tester.out() << "t3d1AsFull: " << t3d1AsFull << std::endl;
  Tensor<3,double,Full> t3d2AsFull = -t3d1AsFull;
  tester.out() << "t3d2AsFull: " << t3d2AsFull << std::endl;

  Tensor<3,double,Diagonal> t3d3(9.0, 9.0, 9.0), t3d4(9.0, 9.0, 9.0);

  t3d3 = t3d1 + t3d2;
  tester.out() << "t3d3 = t3d1 + t3d2: " << t3d3 << std::endl;
  tester.check("t3d3", t3d3, Tensor<3,double,Diagonal>(0.0));
  tester.check("t3d3 against Full", 
               (t3d3 == Tensor<3,double,Diagonal>(0.0)));

  Tensor<3,double,Full> t3f3(99.9), t3f4(99.9), t3f5(99.9), t3f6(99.9);

  t3f3 = t3f1 + t3f2; // No need to check results here; done in TestTensors

  t3f4 = t3d1 + t3d2;
  tester.out() << "t3f4 = t3d1 + t3d2: " << t3f4 << std::endl;
  tester.check("t3f4", (t3f4 == t3d3));

  t3f5 = t3f1 + t3d2;
  tester.out() << "t3f5 = t3f1 + t3d2: " << t3f5 << std::endl;
  tester.check("t3f5", t3f5, t3f1 + t3d2AsFull);

  t3f6 = t3d2 + t3f1;
  tester.out() << "t3f6 = t3d2 + t3f1: " << t3f6 << std::endl;
  tester.check("t3f6", t3f6, t3f1 + t3d2AsFull);

  t3f6 -= t3f1;
  tester.out() << "t3f6 -= t3f1: " << t3f6 << std::endl;
  tester.check("t3f6", t3f6, t3d2AsFull);

  t3d4 = t3d3 - t3f1;
  tester.out() << "t3d4 = t3d3 - t3f1: " << t3d4 << std::endl;
  tester.check("t3d4", 
               (t3d4 == Tensor<3,double,Diagonal>(0,-4,-8)));


  // Test Tensor dot Tensor:

  // Full:
  double sum = 0.0;
  int i, j, k;
  t3f3 = dot(t3f1, t3f2);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      for (k = 0; k < 3; ++k) {
        t3f3(i,k) -= t3f1(i,j)*t3f2(j,k);
      }
    }
  }
  t3f3 = t3f3*t3f3;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      sum += t3f3(i,j);
    }
  }
  tester.check("dot(t3f1, t3f2)", (sum == 0));
  
  // Diagonal:
  sum = 0.0;
  t3f3 = dot(t3d1, t3d2);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      for (k = 0; k < 3; ++k) {
        t3f3(i,k) -= t3d1(i,j)*t3d2(j,k);
      }
    }
  }
  t3f3 = t3f3*t3f3;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      sum += t3f3(i,j);
    }
  }
  tester.check("dot(t3d1, t3d2)", (sum == 0));

  // Test Tensor dot Vector, and vice-versa:

  // Full:
  // Vector dot Tensor
  Vector<3> v31(1.0, 2.0, 3.0);
  tester.out() << "v31: " << v31 << std::endl;
  Vector<3> v32(9.0);
  v32 = dot(v31, t3f2);
  tester.out() << "v32 = dot(v31, t3f2): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(j) -= v31(i)*t3f2(i,j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(v31, t3f2)", (sum == 0));
  // Tensor dot Vector
  v32 = dot(t3f2, v31);
  tester.out() << "v32 = dot(t3f2, v31): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(i) -= t3f2(i,j)*v31(j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(t3f2, v31)", (sum == 0));
  
  // Diagonal:
  // Vector dot Tensor
  v32 = dot(v31, t3d2);
  tester.out() << "v32 = dot(v31, t3d2): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(j) -= v31(i)*t3d2(i,j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(v31, t3d2)", (sum == 0));
  // Tensor dot Vector
  v32 = dot(t3d2, v31);
  tester.out() << "v32 = dot(t3d2, v31): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(i) -= t3d2(i,j)*v31(j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(t3d2, v31)", (sum == 0));


  // --------------------------------------------------------------------------
  // 2D 
  // --------------------------------------------------------------------------

  Tensor<2,double,Full> t2f1(0.0, 2.0, 1.0, 3.0);
  tester.out() << "t2f1: " << t2f1 << std::endl;
  Tensor<2,double,Full> t2f2 = -t2f1;
  tester.out() << "t2f2: " << t2f2 << std::endl;

  Tensor<2,double,Diagonal> t2d1(1.0, 2.0);
  tester.out() << "t2d1: " << t2d1 << std::endl;
  Tensor<2,double,Diagonal> t2d2(-1.0, -2.0);
  tester.out() << "t2d2: " << t2d2 << std::endl;

  Tensor<2,double,Full> t2d1AsFull(1.0,0.0, 0.0,2.0);
  tester.out() << "t2d1AsFull: " << t2d1AsFull << std::endl;
  Tensor<2,double,Full> t2d2AsFull = -t2d1AsFull;
  tester.out() << "t2d2AsFull: " << t2d2AsFull << std::endl;

  Tensor<2,double,Diagonal> t2d3(9.0, 9.0), t2d4(9.0, 9.0);

  t2d3 = t2d1 + t2d2;
  tester.out() << "t2d3 = t2d1 + t2d2: " << t2d3 << std::endl;
  tester.check("t2d3", t2d3, Tensor<2,double,Diagonal>(0.0));
  tester.check("t2d3 against Full", 
               (t2d3 == Tensor<2,double,Diagonal>(0.0)));

  Tensor<2,double,Full> t2f3(99.9), t2f4(99.9), t2f5(99.9), t2f6(99.9), 
    t2f7(99.9);

  t2f3 = t2f1 + t2f2;
  tester.out() << "t2f3 = t2f1 + t2f2: " << t2f3 << std::endl;
  tester.check("t2f3", t2f3, Tensor<2,double,Full>(0.0));

  t2f4 = t2d1 + t2d2;
  tester.out() << "t2f4 = t2d1 + t2d2: " << t2f4 << std::endl;
  tester.check("t2f4", (t2f4 == t2d3));

  t2f5 = t2f1 + t2d2;
  tester.out() << "t2f5 = t2f1 + t2d2: " << t2f5 << std::endl;
  tester.check("t2f5", t2f5, t2f1 + t2d2AsFull);

  t2f6 = t2d2 + t2f1;
  tester.out() << "t2f6 = t2d2 + t2f1: " << t2f6 << std::endl;
  tester.check("t2f6", t2f6, t2f1 + t2d2AsFull);

  t2f6 -= t2f1;
  tester.out() << "t2f6 -= t2f1: " << t2f6 << std::endl;
  tester.check("t2f6", t2f6, t2d2AsFull);

  t2d4 = t2d3 - t2f1;
  tester.out() << "t2d4 = t2d3 - t2f1: " << t2d4 << std::endl;
  tester.check("t2d4", 
               (t2d4 == Tensor<2,double,Diagonal>(0, -3)));


  // Test Tensor dot Tensor:

  // Full:
  sum = 0.0;
  t2f3 = dot(t2f1, t2f2);
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      for (k = 0; k < 2; ++k) {
        t2f3(i,k) -= t2f1(i,j)*t2f2(j,k);
      }
    }
  }
  t2f3 = t2f3*t2f3;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      sum += t2f3(i,j);
    }
  }
  tester.check("dot(t2f1, t2f2)", (sum == 0));
  
  // Diagonal:
  sum = 0.0;
  t2f3 = dot(t2d1, t2d2);
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      for (k = 0; k < 2; ++k) {
        t2f3(i,k) -= t2d1(i,j)*t2d2(j,k);
      }
    }
  }
  t2f3 = t2f3*t2f3;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      sum += t2f3(i,j);
    }
  }
  tester.check("dot(t2d1, t2d2)", (sum == 0));

  // Test Tensor dot Vector, and vice-versa:

  // Full:
  // Vector dot Tensor
  Vector<2> v21(1.0, 2.0);
  tester.out() << "v21: " << v21 << std::endl;
  Vector<2> v22(9.0);
  v22 = dot(v21, t2f2);
  tester.out() << "v22 = dot(v21, t2f2): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(j) -= v21(i)*t2f2(i,j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(v21, t2f2)", (sum == 0));
  // Tensor dot Vector
  v22 = dot(t2f2, v21);
  tester.out() << "v22 = dot(t2f2, v21): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(i) -= t2f2(i,j)*v21(j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(t2f2, v21)", (sum == 0));
  
  // Diagonal:
  // Vector dot Tensor
  v22 = dot(v21, t2d2);
  tester.out() << "v22 = dot(v21, t2d2): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(j) -= v21(i)*t2d2(i,j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(v21, t2d2)", (sum == 0));
  // Tensor dot Vector
  v22 = dot(t2d2, v21);
  tester.out() << "v22 = dot(t2d2, v21): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(i) -= t2d2(i,j)*v21(j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(t2d2, v21)", (sum == 0));


  // --------------------------------------------------------------------------
  // 1D 
  // --------------------------------------------------------------------------

  Tensor<1,double,Full> t1f1(1.0);
  tester.out() << "t1f1: " << t1f1 << std::endl;
  Tensor<1,double,Full> t1f2 = -t1f1;
  tester.out() << "t1f2: " << t1f2 << std::endl;

  Tensor<1,double,Diagonal> t1d1(1.0);
  tester.out() << "t1d1: " << t1d1 << std::endl;
  Tensor<1,double,Diagonal> t1d2(-1.0);
  tester.out() << "t1d2: " << t1d2 << std::endl;

  Tensor<1,double,Full> t1d1AsFull(1.0);
  tester.out() << "t1d1AsFull: " << t1d1AsFull << std::endl;
  Tensor<1,double,Full> t1d2AsFull = -t1d1AsFull;
  tester.out() << "t1d2AsFull: " << t1d2AsFull << std::endl;

  Tensor<1,double,Diagonal> t1d3(9.0), t1d4(9.0);

  t1d3 = t1d1 + t1d2;
  tester.out() << "t1d3 = t1d1 + t1d2: " << t1d3 << std::endl;
  tester.check("t1d3", t1d3, Tensor<1,double,Diagonal>(0.0));
  tester.check("t1d3 against Full", 
               (t1d3 == Tensor<1,double,Diagonal>(0.0)));

  Tensor<1,double,Full> t1f3(99.9), t1f4(99.9), t1f5(99.9), t1f6(99.9), 
    t1f7(99.9);

  t1f3 = t1f1 + t1f2;
  tester.out() << "t1f3 = t1f1 + t1f2: " << t1f3 << std::endl;
  tester.check("t1f3", t1f3, Tensor<1,double,Full>(0.0));

  t1f4 = t1d1 + t1d2;
  tester.out() << "t1f4 = t1d1 + t1d2: " << t1f4 << std::endl;
  tester.check("t1f4", (t1f4 == t1d3));

  t1f5 = t1f1 + t1d2;
  tester.out() << "t1f5 = t1f1 + t1d2: " << t1f5 << std::endl;
  tester.check("t1f5", t1f5, t1f1 + t1d2AsFull);

  t1f6 = t1d2 + t1f1;
  tester.out() << "t1f6 = t1d2 + t1f1: " << t1f6 << std::endl;
  tester.check("t1f6", t1f6, t1f1 + t1d2AsFull);

  t1f6 -= t1f1;
  tester.out() << "t1f6 -= t1f1: " << t1f6 << std::endl;
  tester.check("t1f6", t1f6, t1d2AsFull);

  t1d4 = t1d3 - t1f1;
  tester.out() << "t1d4 = t1d3 - t1f1: " << t1d4 << std::endl;
  tester.check("t1d4", 
               (t1d4 == Tensor<1,double,Diagonal>(-1)));


  // Test Tensor dot Tensor:

  // Full:
  sum = 0.0;
  t1f3 = dot(t1f1, t1f2);
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      for (k = 0; k < 1; ++k) {
        t1f3(i,k) -= t1f1(i,j)*t1f2(j,k);
      }
    }
  }
  t1f3 = t1f3*t1f3;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      sum += t1f3(i,j);
    }
  }
  tester.check("dot(t1f1, t1f2)", (sum == 0));
  
  // Diagonal:
  sum = 0.0;
  t1f3 = dot(t1d1, t1d2);
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      for (k = 0; k < 1; ++k) {
        t1f3(i,k) -= t1d1(i,j)*t1d2(j,k);
      }
    }
  }
  t1f3 = t1f3*t1f3;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      sum += t1f3(i,j);
    }
  }
  tester.check("dot(t1d1, t1d2)", (sum == 0));

  // Test Tensor dot Vector, and vice-versa:

  // Full:
  // Vector dot Tensor
  Vector<1> v11(1.0);
  tester.out() << "v11: " << v11 << std::endl;
  Vector<1> v12(9.0);
  v12 = dot(v11, t1f2);
  tester.out() << "v12 = dot(v11, t1f2): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(j) -= v11(i)*t1f2(i,j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(v11, t1f2)", (sum == 0));
  // Tensor dot Vector
  v12 = dot(t1f2, v11);
  tester.out() << "v12 = dot(t1f2, v11): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(i) -= t1f2(i,j)*v11(j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(t1f2, v11)", (sum == 0));
  
  // Diagonal:
  // Vector dot Tensor
  v12 = dot(v11, t1d2);
  tester.out() << "v12 = dot(v11, t1d2): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(j) -= v11(i)*t1d2(i,j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(v11, t1d2)", (sum == 0));
  // Tensor dot Vector
  v12 = dot(t1d2, v11);
  tester.out() << "v12 = dot(t1d2, v11): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(i) -= t1d2(i,j)*v11(j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(t1d2, v11)", (sum == 0));


  int ret = tester.results("TestDiagonalTensors");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestDiagonalTensors.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
