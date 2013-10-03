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
// Various tests of Tensor<D,double,Symmetric>
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

  Tensor<3,double,Symmetric> t3s1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  tester.out() << "t3s1: " << t3s1 << std::endl;
  Tensor<3,double,Symmetric> t3s2(-1.0, -2.0, -3.0, -4.0, -5.0, -6.0);
  tester.out() << "t3s2: " << t3s2 << std::endl;

  Tensor<3,double,Full> t3s1AsFull(1.0,2.0,4.0, 2.0,3.0,5.0, 4.0,5.0,6.0);
  tester.out() << "t3s1AsFull: " << t3s1AsFull << std::endl;
  Tensor<3,double,Full> t3s2AsFull = -t3s1AsFull;
  tester.out() << "t3s2AsFull: " << t3s2AsFull << std::endl;

  Tensor<3,double,Symmetric> t3s3(9.0, 9.0, 9.0, 9.0, 9.0, 9.0), 
    t3s4(9.0, 9.0, 9.0, 9.0, 9.0, 9.0);

  t3s3 = t3s1 + t3s2;
  tester.out() << "t3s3 = t3s1 + t3s2: " << t3s3 << std::endl;
  tester.check("t3s3", t3s3, Tensor<3,double,Symmetric>(0.0));
  tester.check("t3s3 against Full", 
               (t3s3 == Tensor<3,double,Symmetric>(0.0)));

  Tensor<3,double,Full> t3f3(99.9), t3f4(99.9), t3f5(99.9), t3f6(99.9);

  t3f3 = t3f1 + t3f2; // No need to check results here; done in TestTensors

  t3f4 = t3s1 + t3s2;
  tester.out() << "t3f4 = t3s1 + t3s2: " << t3f4 << std::endl;
  tester.check("t3f4", (t3f4 == t3s3));

  t3f5 = t3f1 + t3s2;
  tester.out() << "t3f5 = t3f1 + t3s2: " << t3f5 << std::endl;
  tester.check("t3f5", t3f5, t3f1 + t3s2AsFull);

  t3f6 = t3s2 + t3f1;
  tester.out() << "t3f6 = t3s2 + t3f1: " << t3f6 << std::endl;
  tester.check("t3f6", t3f6, t3f1 + t3s2AsFull);

  t3f6 -= t3f1;
  tester.out() << "t3f6 -= t3f1: " << t3f6 << std::endl;
  tester.check("t3f6", t3f6, t3s2AsFull);

  t3s4 = t3s3 - t3f1;
  tester.out() << "t3s4 = t3s3 - t3f1: " << t3s4 << std::endl;
  tester.check("t3s4", 
               (t3s4 == Tensor<3,double,Symmetric>(0,-3,-4,-6,-7,-8)));


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
  
  // Symmetric:
  sum = 0.0;
  t3f3 = dot(t3s1, t3s2);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      for (k = 0; k < 3; ++k) {
        t3f3(i,k) -= t3s1(i,j)*t3s2(j,k);
      }
    }
  }
  t3f3 = t3f3*t3f3;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      sum += t3f3(i,j);
    }
  }
  tester.check("dot(t3s1, t3s2)", (sum == 0));

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
  
  // Symmetric:
  // Vector dot Tensor
  v32 = dot(v31, t3s2);
  tester.out() << "v32 = dot(v31, t3s2): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(j) -= v31(i)*t3s2(i,j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(v31, t3s2)", (sum == 0));
  // Tensor dot Vector
  v32 = dot(t3s2, v31);
  tester.out() << "v32 = dot(t3s2, v31): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(i) -= t3s2(i,j)*v31(j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(t3s2, v31)", (sum == 0));


  // --------------------------------------------------------------------------
  // 2D 
  // --------------------------------------------------------------------------

  Tensor<2,double,Full> t2f1(0.0, 2.0, 1.0, 3.0);
  tester.out() << "t2f1: " << t2f1 << std::endl;
  Tensor<2,double,Full> t2f2 = -t2f1;
  tester.out() << "t2f2: " << t2f2 << std::endl;

  Tensor<2,double,Symmetric> t2s1(1.0, 2.0, 3.0);
  tester.out() << "t2s1: " << t2s1 << std::endl;
  Tensor<2,double,Symmetric> t2s2(-1.0, -2.0, -3.0);
  tester.out() << "t2s2: " << t2s2 << std::endl;

  Tensor<2,double,Full> t2s1AsFull(1.0,2.0, 2.0,3.0);
  tester.out() << "t2s1AsFull: " << t2s1AsFull << std::endl;
  Tensor<2,double,Full> t2s2AsFull = -t2s1AsFull;
  tester.out() << "t2s2AsFull: " << t2s2AsFull << std::endl;

  Tensor<2,double,Symmetric> t2s3(9.0, 9.0, 9.0), t2s4(9.0, 9.0, 9.0);

  t2s3 = t2s1 + t2s2;
  tester.out() << "t2s3 = t2s1 + t2s2: " << t2s3 << std::endl;
  tester.check("t2s3", t2s3, Tensor<2,double,Symmetric>(0.0));
  tester.check("t2s3 against Full", 
               (t2s3 == Tensor<2,double,Symmetric>(0.0)));

  Tensor<2,double,Full> t2f3(99.9), t2f4(99.9), t2f5(99.9), t2f6(99.9), 
    t2f7(99.9);

  t2f3 = t2f1 + t2f2;
  tester.out() << "t2f3 = t2f1 + t2f2: " << t2f3 << std::endl;
  tester.check("t2f3", t2f3, Tensor<2,double,Full>(0.0));

  t2f4 = t2s1 + t2s2;
  tester.out() << "t2f4 = t2s1 + t2s2: " << t2f4 << std::endl;
  tester.check("t2f4", (t2f4 == t2s3));

  t2f5 = t2f1 + t2s2;
  tester.out() << "t2f5 = t2f1 + t2s2: " << t2f5 << std::endl;
  tester.check("t2f5", t2f5, t2f1 + t2s2AsFull);

  t2f6 = t2s2 + t2f1;
  tester.out() << "t2f6 = t2s2 + t2f1: " << t2f6 << std::endl;
  tester.check("t2f6", t2f6, t2f1 + t2s2AsFull);

  t2f6 -= t2f1;
  tester.out() << "t2f6 -= t2f1: " << t2f6 << std::endl;
  tester.check("t2f6", t2f6, t2s2AsFull);

  t2s4 = t2s3 - t2f1;
  tester.out() << "t2s4 = t2s3 - t2f1: " << t2s4 << std::endl;
  tester.check("t2s4", 
               (t2s4 == Tensor<2,double,Symmetric>(-0, -2, -3)));


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
  
  // Symmetric:
  sum = 0.0;
  t2f3 = dot(t2s1, t2s2);
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      for (k = 0; k < 2; ++k) {
        t2f3(i,k) -= t2s1(i,j)*t2s2(j,k);
      }
    }
  }
  t2f3 = t2f3*t2f3;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      sum += t2f3(i,j);
    }
  }
  tester.check("dot(t2s1, t2s2)", (sum == 0));

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
  
  // Symmetric:
  // Vector dot Tensor
  v22 = dot(v21, t2s2);
  tester.out() << "v22 = dot(v21, t2s2): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(j) -= v21(i)*t2s2(i,j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(v21, t2s2)", (sum == 0));
  // Tensor dot Vector
  v22 = dot(t2s2, v21);
  tester.out() << "v22 = dot(t2s2, v21): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(i) -= t2s2(i,j)*v21(j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(t2s2, v21)", (sum == 0));


  // --------------------------------------------------------------------------
  // 1D 
  // --------------------------------------------------------------------------

  Tensor<1,double,Full> t1f1(1.0);
  tester.out() << "t1f1: " << t1f1 << std::endl;
  Tensor<1,double,Full> t1f2 = -t1f1;
  tester.out() << "t1f2: " << t1f2 << std::endl;

  Tensor<1,double,Symmetric> t1s1(1.0);
  tester.out() << "t1s1: " << t1s1 << std::endl;
  Tensor<1,double,Symmetric> t1s2(-1.0);
  tester.out() << "t1s2: " << t1s2 << std::endl;

  Tensor<1,double,Full> t1s1AsFull(1.0);
  tester.out() << "t1s1AsFull: " << t1s1AsFull << std::endl;
  Tensor<1,double,Full> t1s2AsFull = -t1s1AsFull;
  tester.out() << "t1s2AsFull: " << t1s2AsFull << std::endl;

  Tensor<1,double,Symmetric> t1s3(9.0), t1s4(9.0);

  t1s3 = t1s1 + t1s2;
  tester.out() << "t1s3 = t1s1 + t1s2: " << t1s3 << std::endl;
  tester.check("t1s3", t1s3, Tensor<1,double,Symmetric>(0.0));
  tester.check("t1s3 against Full", 
               (t1s3 == Tensor<1,double,Symmetric>(0.0)));

  Tensor<1,double,Full> t1f3(99.9), t1f4(99.9), t1f5(99.9), t1f6(99.9), 
    t1f7(99.9);

  t1f3 = t1f1 + t1f2;
  tester.out() << "t1f3 = t1f1 + t1f2: " << t1f3 << std::endl;
  tester.check("t1f3", t1f3, Tensor<1,double,Full>(0.0));

  t1f4 = t1s1 + t1s2;
  tester.out() << "t1f4 = t1s1 + t1s2: " << t1f4 << std::endl;
  tester.check("t1f4", (t1f4 == t1s3));

  t1f5 = t1f1 + t1s2;
  tester.out() << "t1f5 = t1f1 + t1s2: " << t1f5 << std::endl;
  tester.check("t1f5", t1f5, t1f1 + t1s2AsFull);

  t1f6 = t1s2 + t1f1;
  tester.out() << "t1f6 = t1s2 + t1f1: " << t1f6 << std::endl;
  tester.check("t1f6", t1f6, t1f1 + t1s2AsFull);

  t1f6 -= t1f1;
  tester.out() << "t1f6 -= t1f1: " << t1f6 << std::endl;
  tester.check("t1f6", t1f6, t1s2AsFull);

  t1s4 = t1s3 - t1f1;
  tester.out() << "t1s4 = t1s3 - t1f1: " << t1s4 << std::endl;
  tester.check("t1s4", 
               (t1s4 == Tensor<1,double,Symmetric>(-1)));


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
  
  // Symmetric:
  sum = 0.0;
  t1f3 = dot(t1s1, t1s2);
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      for (k = 0; k < 1; ++k) {
        t1f3(i,k) -= t1s1(i,j)*t1s2(j,k);
      }
    }
  }
  t1f3 = t1f3*t1f3;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      sum += t1f3(i,j);
    }
  }
  tester.check("dot(t1s1, t1s2)", (sum == 0));

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
  
  // Symmetric:
  // Vector dot Tensor
  v12 = dot(v11, t1s2);
  tester.out() << "v12 = dot(v11, t1s2): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(j) -= v11(i)*t1s2(i,j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(v11, t1s2)", (sum == 0));
  // Tensor dot Vector
  v12 = dot(t1s2, v11);
  tester.out() << "v12 = dot(t1s2, v11): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(i) -= t1s2(i,j)*v11(j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(t1s2, v11)", (sum == 0));


  int ret = tester.results("TestSymmetricTensors");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestSymmetricTensors.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
