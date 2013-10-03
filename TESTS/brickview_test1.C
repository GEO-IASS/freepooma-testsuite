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
// BrickEngine test code
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"
#include "Domain/SliceInterval.h"
#include "Domain/SliceRange.h"
#include "Domain/AllDomain.h"
#include "Engine/BrickEngine.h"

typedef Engine<5,double,Brick> Brick5_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif
    tester.out() << "\n\nTesting 5D slice-range subsetting...\n\n";

    Interval<1> I(5);
    Interval<5> BD(I,I,I,I,I);

    Brick5_t b(BD);

    for (int i4 = 0; i4 < b.domain()[4].length(); ++i4)
      for (int i3 = 0; i3 < b.domain()[3].length(); ++i3)
	for (int i2 = 0; i2 < b.domain()[2].length(); ++i2)
	  for (int i1 = 0; i1 < b.domain()[1].length(); ++i1)
	    for (int i0 = 0; i0 < b.domain()[0].length(); ++i0)
	      b(i0,i1,i2,i3,i4) = i4+10*(i3+10*(i2+10*(i1+10*i0)));

    tester.out() << "b.domain()     = " << b.domain() << std::endl;
    tester.out() << std::endl;

    typedef NewDomain5<int, Range<1>, int, AllDomain<1>, Interval<1> > 
      NewDomain_t;
    typedef NewDomain_t::SliceType_t SliceType_t;
    SliceType_t VD;
    AllDomain<1> A;
    Interval<1> I1(1,3);
    Range<1> R(0,4,2);
    NewDomain_t::fillSlice(VD, b.domain(), 2, R, 1, A, I1);

    tester.out() << "VD = " << VD << std::endl;

    typedef NewEngine<Brick5_t, SliceType_t>::Type_t Engine_t;
    Engine_t v(b, VD);

    // v.domain() should be:     [0:2:1,0:4:1,0:2:1]

    tester.out() << "v.domain()     = " << v.domain() << std::endl;

    // v's values should be 2 | 0,2,4 | 1 | 0,1,2,3,4 | 1,2,3

    tester.out() << "v = \n";
    for (int i2 = 0; i2 < v.domain()[2].length(); ++i2)
      for (int i1 = 0; i1 < v.domain()[1].length(); ++i1)
	for (int i0 = 0; i0 < v.domain()[0].length(); ++i0)
	  tester.out() << v(i0,i1,i2) << std::endl;
            
    typedef NewDomain3<int, Range<1>, Interval<1> > NewDomain2_t;
    typedef NewDomain2_t::SliceType_t SliceType2_t;
    SliceType2_t VD2;
    Interval<1> I2(1,2);
    Range<1> R2(0,2,2);
    NewDomain2_t::fillSlice(VD2, v.domain(), 0, R2, I2);

    tester.out() << "VD2 = " << VD2 << std::endl;

    NewEngine<Engine_t, SliceType2_t>::Type_t v2(v, VD2);

    // v2.domain() should be:     [0:1:1,0:1:1]
      
    tester.out() << "v2.domain()     = " << v2.domain() << std::endl;

    // v2's values should be 2 | 0 | 1 | 0,2 | 2,3

    tester.out() << "v2 = \n";
    for (int i1 = 0; i1 < v2.domain()[1].length(); ++i1)
      for (int i0 = 0; i0 < v2.domain()[0].length(); ++i0)
	tester.out() << v2(i0,i1) << std::endl;
#if POOMA_EXCEPTIONS
  }
  catch(const char *err) 
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
  catch(const Pooma::Assertion &err)
    { 
      tester.exceptionHandler( err );
      tester.set( false );
    }
#endif    
  int ret = tester.results("brickview_test1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brickview_test1.cpp,v $   $Author: richard $
// $Revision: 1.16 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
