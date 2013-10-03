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
// Various tests of Tensor<D,double,[Full,Antisymmetric]>
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
  Tensor<3,double,Full> tf1(0.0, 3.0, 6.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0);
  tester.out() << "tf1: " << tf1 << std::endl;
  Tensor<3,double,Full> tf2 = -tf1;
  tester.out() << "tf2: " << tf2 << std::endl;

  Tensor<3,double,Antisymmetric> ta1(1.0, 2.0, 3.0);
  tester.out() << "ta1: " << ta1 << std::endl;
  Tensor<3,double,Antisymmetric> ta2(-1.0, -2.0, -3.0);
  tester.out() << "ta2: " << ta2 << std::endl;

  Tensor<3,double,Full> ta1AsFull(0.0,1.0,2.0, -1.0,0.0,3.0, -2.0,-3.0,0.0);
  tester.out() << "ta1AsFull: " << ta1AsFull << std::endl;
  Tensor<3,double,Full> ta2AsFull = -ta1AsFull;
  tester.out() << "ta2AsFull: " << ta2AsFull << std::endl;

  Tensor<3,double,Antisymmetric> ta3(9.0, 9.0, 9.0), ta4(9.0, 9.0, 9.0);

  ta3 = ta1 + ta2;
  tester.out() << "ta3 = ta1 + ta2: " << ta3 << std::endl;
  tester.check("ta3", ta3, Tensor<3,double,Antisymmetric>(0.0));
  tester.check("ta3 against Full", 
               (ta3 == Tensor<3,double,Antisymmetric>(0.0)));

  Tensor<3,double,Full> tf3(99.9), tf4(99.9), tf5(99.9), tf6(99.9), tf7(99.9);

  tf3 = tf1 + tf2;
  tester.out() << "tf3 = tf1 + tf2: " << tf3 << std::endl;
  tester.check("tf3", tf3, Tensor<3,double,Full>(0.0));

  tf4 = ta1 + ta2;
  tester.out() << "tf4 = ta1 + ta2: " << tf4 << std::endl;
  tester.check("tf4", (tf4 == ta3));

  tf5 = tf1 + ta2;
  tester.out() << "tf5 = tf1 + ta2: " << tf5 << std::endl;
  tester.check("tf5", tf5, tf1 + ta2AsFull);

  tf6 = ta2 + tf1;
  tester.out() << "tf6 = ta2 + tf1: " << tf6 << std::endl;
  tester.check("tf6", tf6, tf1 + ta2AsFull);

  tf6 -= tf1;
  tester.out() << "tf6 -= tf1: " << tf6 << std::endl;
  tester.check("tf6", tf6, ta2AsFull);

  ta4 = ta3 - tf1;
  tester.out() << "ta4 = ta3 - tf1: " << ta4 << std::endl;
  tester.check("ta4", 
               (ta4 == Tensor<3,double,Antisymmetric>(-3,-6,-7)));


  // Test Tensor dot Tensor:

  // Full:
  double sum = 0.0;
  int i, j, k;
  tf3 = dot(tf1, tf2);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      for (k = 0; k < 3; ++k) {
        tf3(i,k) -= tf1(i,j)*tf2(j,k);
      }
    }
  }
  tf3 = tf3*tf3;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      sum += tf3(i,j);
    }
  }
  tester.check("dot(tf1, tf2)", (sum == 0));
  
  // Antisymmetric:
  sum = 0.0;
  tf3 = dot(ta1, ta2);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      for (k = 0; k < 3; ++k) {
        tf3(i,k) -= ta1(i,j)*ta2(j,k);
      }
    }
  }
  tf3 = tf3*tf3;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      sum += tf3(i,j);
    }
  }
  tester.check("dot(ta1, ta2)", (sum == 0));

  // Test Tensor dot Vector, and vice-versa:

  // Full:
  // Vector dot Tensor
  Vector<3> v31(1.0, 2.0, 3.0);
  tester.out() << "v31: " << v31 << std::endl;
  Vector<3> v32(9.0);
  v32 = dot(v31, tf2);
  tester.out() << "v32 = dot(v31, tf2): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(j) -= v31(i)*tf2(i,j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(v31, tf2)", (sum == 0));
  // Tensor dot Vector
  v32 = dot(tf2, v31);
  tester.out() << "v32 = dot(tf2, v31): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(i) -= tf2(i,j)*v31(j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(tf2, v31)", (sum == 0));
  
  // Antisymmetric:
  // Vector dot Tensor
  v32 = dot(v31, ta2);
  tester.out() << "v32 = dot(v31, ta2): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(j) -= v31(i)*ta2(i,j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(v31, ta2)", (sum == 0));
  // Tensor dot Vector
  v32 = dot(ta2, v31);
  tester.out() << "v32 = dot(ta2, v31): " << v32 << std::endl;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      v32(i) -= ta2(i,j)*v31(j);
    }
  }
  v32 = v32*v32;
  sum = 0.0;
  for (i = 0; i < 3; ++i) {
    sum += v32(i);
  }
  tester.check("dot(ta2, v31)", (sum == 0));


  // --------------------------------------------------------------------------
  // 2D 
  // --------------------------------------------------------------------------

  Tensor<2,double,Full> t2f1(0.0, 2.0, 1.0, 3.0);
  tester.out() << "t2f1: " << t2f1 << std::endl;
  Tensor<2,double,Full> t2f2 = -t2f1;
  tester.out() << "t2f2: " << t2f2 << std::endl;

  Tensor<2,double,Antisymmetric> t2a1(1.0);
  tester.out() << "t2a1: " << t2a1 << std::endl;
  Tensor<2,double,Antisymmetric> t2a2(-1.0);
  tester.out() << "t2a2: " << t2a2 << std::endl;

  Tensor<2,double,Full> t2a1AsFull(0.0,1.0, -1.0,0.0);
  tester.out() << "t2a1AsFull: " << t2a1AsFull << std::endl;
  Tensor<2,double,Full> t2a2AsFull = -t2a1AsFull;
  tester.out() << "t2a2AsFull: " << t2a2AsFull << std::endl;

  Tensor<2,double,Antisymmetric> t2a3(9.0), t2a4(9.0);

  t2a3 = t2a1 + t2a2;
  tester.out() << "t2a3 = t2a1 + t2a2: " << t2a3 << std::endl;
  tester.check("t2a3", t2a3, Tensor<2,double,Antisymmetric>(0.0));
  tester.check("t2a3 against Full", 
               (t2a3 == Tensor<2,double,Antisymmetric>(0.0)));

  Tensor<2,double,Full> t2f3(99.9), t2f4(99.9), t2f5(99.9), t2f6(99.9), 
    t2f7(99.9);

  t2f3 = t2f1 + t2f2;
  tester.out() << "t2f3 = t2f1 + t2f2: " << t2f3 << std::endl;
  tester.check("t2f3", t2f3, Tensor<2,double,Full>(0.0));

  t2f4 = t2a1 + t2a2;
  tester.out() << "t2f4 = t2a1 + t2a2: " << t2f4 << std::endl;
  tester.check("t2f4", (t2f4 == t2a3));

  t2f5 = t2f1 + t2a2;
  tester.out() << "t2f5 = t2f1 + t2a2: " << t2f5 << std::endl;
  tester.check("t2f5", t2f5, t2f1 + t2a2AsFull);

  t2f6 = t2a2 + t2f1;
  tester.out() << "t2f6 = t2a2 + t2f1: " << t2f6 << std::endl;
  tester.check("t2f6", t2f6, t2f1 + t2a2AsFull);

  t2f6 -= t2f1;
  tester.out() << "t2f6 -= t2f1: " << t2f6 << std::endl;
  tester.check("t2f6", t2f6, t2a2AsFull);

  t2a4 = t2a3 - t2f1;
  tester.out() << "t2a4 = t2a3 - t2f1: " << t2a4 << std::endl;
  tester.check("t2a4", 
               (t2a4 == Tensor<2,double,Antisymmetric>(-2)));


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
  
  // Antisymmetric:
  sum = 0.0;
  t2f3 = dot(t2a1, t2a2);
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      for (k = 0; k < 2; ++k) {
        t2f3(i,k) -= t2a1(i,j)*t2a2(j,k);
      }
    }
  }
  t2f3 = t2f3*t2f3;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      sum += t2f3(i,j);
    }
  }
  tester.check("dot(t2a1, t2a2)", (sum == 0));

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
  
  // Antisymmetric:
  // Vector dot Tensor
  v22 = dot(v21, t2a2);
  tester.out() << "v22 = dot(v21, t2a2): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(j) -= v21(i)*t2a2(i,j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(v21, t2a2)", (sum == 0));
  // Tensor dot Vector
  v22 = dot(t2a2, v21);
  tester.out() << "v22 = dot(t2a2, v21): " << v22 << std::endl;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      v22(i) -= t2a2(i,j)*v21(j);
    }
  }
  v22 = v22*v22;
  sum = 0.0;
  for (i = 0; i < 2; ++i) {
    sum += v22(i);
  }
  tester.check("dot(t2a2, v21)", (sum == 0));


  // --------------------------------------------------------------------------
  // 1D 
  // --------------------------------------------------------------------------

  Tensor<1,double,Full> t1f1(1.0);
  tester.out() << "t1f1: " << t1f1 << std::endl;
  Tensor<1,double,Full> t1f2 = -t1f1;
  tester.out() << "t1f2: " << t1f2 << std::endl;

  Tensor<1,double,Antisymmetric> t1a1(1.0);
  tester.out() << "t1a1: " << t1a1 << std::endl;
  Tensor<1,double,Antisymmetric> t1a2(-1.0);
  tester.out() << "t1a2: " << t1a2 << std::endl;

  Tensor<1,double,Full> t1a1AsFull(0.0);
  tester.out() << "t1a1AsFull: " << t1a1AsFull << std::endl;
  Tensor<1,double,Full> t1a2AsFull = -t1a1AsFull;
  tester.out() << "t1a2AsFull: " << t1a2AsFull << std::endl;

  Tensor<1,double,Antisymmetric> t1a3(9.0), t1a4(9.0);

  t1a3 = t1a1 + t1a2;
  tester.out() << "t1a3 = t1a1 + t1a2: " << t1a3 << std::endl;
  tester.check("t1a3", t1a3, Tensor<1,double,Antisymmetric>(0.0));
  tester.check("t1a3 against Full", 
               (t1a3 == Tensor<1,double,Antisymmetric>(0.0)));

  Tensor<1,double,Full> t1f3(99.9), t1f4(99.9), t1f5(99.9), t1f6(99.9), 
    t1f7(99.9);

  t1f3 = t1f1 + t1f2;
  tester.out() << "t1f3 = t1f1 + t1f2: " << t1f3 << std::endl;
  tester.check("t1f3", t1f3, Tensor<1,double,Full>(0.0));

  t1f4 = t1a1 + t1a2;
  tester.out() << "t1f4 = t1a1 + t1a2: " << t1f4 << std::endl;
  tester.check("t1f4", (t1f4 == t1a3));

  t1f5 = t1f1 + t1a2;
  tester.out() << "t1f5 = t1f1 + t1a2: " << t1f5 << std::endl;
  tester.check("t1f5", t1f5, t1f1 + t1a2AsFull);

  t1f6 = t1a2 + t1f1;
  tester.out() << "t1f6 = t1a2 + t1f1: " << t1f6 << std::endl;
  tester.check("t1f6", t1f6, t1f1 + t1a2AsFull);

  t1f6 -= t1f1;
  tester.out() << "t1f6 -= t1f1: " << t1f6 << std::endl;
  tester.check("t1f6", t1f6, t1a2AsFull);

  t1a4 = t1a3 - t1f1;
  tester.out() << "t1a4 = t1a3 - t1f1: " << t1a4 << std::endl;
  tester.check("t1a4", 
               (t1a4 == Tensor<1,double,Antisymmetric>(-2)));


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
  
  // Antisymmetric:
  sum = 0.0;
  t1f3 = dot(t1a1, t1a2);
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      for (k = 0; k < 1; ++k) {
        t1f3(i,k) -= t1a1(i,j)*t1a2(j,k);
      }
    }
  }
  t1f3 = t1f3*t1f3;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      sum += t1f3(i,j);
    }
  }
  tester.check("dot(t1a1, t1a2)", (sum == 0));

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
  
  // Antisymmetric:
  // Vector dot Tensor
  v12 = dot(v11, t1a2);
  tester.out() << "v12 = dot(v11, t1a2): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(j) -= v11(i)*t1a2(i,j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(v11, t1a2)", (sum == 0));
  // Tensor dot Vector
  v12 = dot(t1a2, v11);
  tester.out() << "v12 = dot(t1a2, v11): " << v12 << std::endl;
  for (i = 0; i < 1; ++i) {
    for (j = 0; j < 1; ++j) {
      v12(i) -= t1a2(i,j)*v11(j);
    }
  }
  v12 = v12*v12;
  sum = 0.0;
  for (i = 0; i < 1; ++i) {
    sum += v12(i);
  }
  tester.check("dot(t1a2, v11)", (sum == 0));


  int ret = tester.results("TestTensors");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestTensors.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:17:12 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
