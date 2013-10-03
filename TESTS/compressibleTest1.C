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
// evaluation of compressible brick test #1
// (this test is really simple and only checks that the results are compressed
// when they should be)
//-----------------------------------------------------------------------------

// include files

#include "Pooma/Pooma.h"
#include "Array/Array.h"
#include "Engine/BrickEngine.h"
#include "Evaluator/CompressibleEval.h"
#include "Utilities/Tester.h"
#include <iostream>
#include <cmath>


const char* compressed(const Array<1,double,CompressibleBrick>& array)
{
  return (array.engine().compressed() ?
	  "compressed" : "uncompressed" );
}

int main(int argc, char* argv[])
{
  // initialize Pooma
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);

  int return_status = 0;

  int n = 20;
  double pi = 3.1415926535897932;
  Array<1> a(n),b(n),c(n),d(n);
  Array<1,double,CompressibleBrick> aa(n),bb(n),cc(n),dd(n);
  int i;
  for (i=0;i<n;++i)
  {
    a(i) = sin(0.1*pi*i);
    b(i) = cos(0.1*pi*i);
    c(i) = 1.0;
    d(i) = 2.0;
  }

  tester.out() << "Testing Compressile Bricks." << std::endl;

  aa = a;
  bb = b;
  cc = c;
  dd = d;

  Pooma::blockAndEvaluate();

  tester.out() << "  i    aa    bb   cc   dd " << std::endl;
  for (i=0; i<n; ++i)
  {
    tester.out() << i << " " << aa.read(i) << " " << bb.read(i) << " "
		 << cc.read(i) << " " << dd.read(i) << std::endl;
  }

  tester.out() << "aa: " << compressed(aa) << std::endl;
  tester.out() << "bb: " << compressed(bb) << std::endl;
  tester.out() << "cc: " << compressed(cc) << std::endl;
  tester.out() << "dd: " << compressed(dd) << std::endl;

  tester.check(!aa.engine().compressed());
  tester.check(!bb.engine().compressed());
  tester.check(cc.engine().compressed());
  tester.check(dd.engine().compressed());

  aa = aa*aa+bb*bb;
  bb = cc*dd+2.0*cc;
  bb += cc*dd+2.0*cc;
  cc = dd*c+d;

  Pooma::blockAndEvaluate();

  tester.out() << "  i    aa    bb   cc   dd " << std::endl;
  for (i=0; i<n; ++i)
  {
    tester.out() << i << " " << aa.read(i) << " " << bb.read(i) << " "
		 << cc.read(i) << " " << dd.read(i) << std::endl;
  }

  tester.out() << "aa: " << compressed(aa) << std::endl;
  tester.out() << "bb: " << compressed(bb) << std::endl;
  tester.out() << "cc: " << compressed(cc) << std::endl;
  tester.out() << "dd: " << compressed(dd) << std::endl;

  tester.check(!aa.engine().compressed());
  tester.check(bb.engine().compressed());
  tester.check(cc.engine().compressed());
  tester.check(dd.engine().compressed());

  a = b+dd*cc;

  Pooma::blockAndEvaluate();

  tester.out() << "a: ";
  for (i=0; i<n; ++i) {
    tester.out() << "(" << i << ")=" << a.read(i) << ",";
  }
  tester.out() << std::endl;

  tester.out() << "------------------------------------------------" << std::endl;
 
  int retval = tester.results("compressibleTest1");

  Pooma::finalize();  
  return return_status;

}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: compressibleTest1.cpp,v $   $Author: richi $
// $Revision: 1.17 $   $Date: 2004/11/10 22:13:03 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
