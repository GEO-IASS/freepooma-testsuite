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

#include <iomanip>

typedef Engine<1,double,Brick> Brick1_t;
typedef Engine<2,double,Brick> Brick2_t;

typedef Engine<1,double,BrickView> View1_t;
typedef Engine<2,double,BrickView> View2_t;

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif
    {
      // test 1

      tester.out() << "\n\nTesting 1D unit-stride subsetting...\n\n";

      Brick1_t b(Interval<1>(10));

      Interval<1> I(4,5);

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;
      tester.out() << "I = " << I << std::endl;

      for (int i = 0; i < 10; ++i)
	b(i) = i;

      tester.out() << "b = ";
      for (int i = 0; i < 10; ++i)
	tester.out() << b(i) << " ";
      tester.out() << std::endl;

      View1_t v(b,I);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain().length(); ++i)
	tester.out() << v(i) << " ";
      tester.out() << std::endl;

      v(0) = 0; v(1) = 0;

      tester.out() << "v = 0" << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain().length(); ++i)
	tester.out() << v(i) << " ";
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < 10; ++i)
	tester.out() << b(i) << " ";
      tester.out() << std::endl;

    }
    {
      // test 2

      tester.out() << "\n\nTesting 1D stride-2 subsetting...\n\n";

      Brick1_t b(Interval<1>(10));

      Range<1> I(3,7,2);

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;
      tester.out() << "I = " << I << std::endl;

      for (int i = 0; i < 10; ++i)
	b(i) = i;

      tester.out() << "b = ";
      for (int i = 0; i < 10; ++i)
	tester.out() << b(i) << " ";
      tester.out() << std::endl;

      View1_t v(b,I);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain().length(); ++i)
	tester.out() << v(i) << " ";
      tester.out() << std::endl;


      for (int i = 0; i < v.domain().length(); ++i)
	v(i) = 0;

      tester.out() << "After setting v(i) = 0" << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain().length(); ++i)
	tester.out() << v(i) << " ";
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < 10; ++i)
	tester.out() << b(i) << " ";
      tester.out() << std::endl;

    }
    {
      // test 3

      tester.out() << "\n\nTesting 2D unit-stride subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      Interval<2> VD;
      Interval<1> IV(1,3);
      VD[0] = IV; VD[1] = IV;

      tester.out() << "VD = " << I << std::endl;

      View2_t v(b,VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < v.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << v(i,j) << " ";
	}
      tester.out() << std::endl;

      for (int i = 0; i < v.domain()[0].length(); ++i)
	for (int j = 0; j < v.domain()[1].length(); ++j)
	  v(i,j) = 0;


      tester.out() << "After setting v(i,j) = 0" << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < v.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << v(i,j) << " ";
	}
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

    }
    {
      // test 4

      tester.out() << "\n\nTesting 2D stride-2 subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      Range<2> VD;
      Range<1> IV(1,3,2);
      VD[0] = IV; VD[1] = IV;

      tester.out() << "VD = " << I << std::endl;

      View2_t v(b,VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < v.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << v(i,j) << " ";
	}
      tester.out() << std::endl;

      for (int i = 0; i < v.domain()[0].length(); ++i)
	for (int j = 0; j < v.domain()[1].length(); ++j)
	  v(i,j) = 0;


      tester.out() << "After setting v(i,j) = 0" << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < v.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << v(i,j) << " ";
	}
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

    }
    {
      // test 5

      tester.out() << "\n\nTesting 2D wildcard subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      typedef NewDomain2<Range<1>, AllDomain<1> > NewDomain_t;
      typedef NewDomain_t::SliceType_t SliceType_t;
      SliceType_t VD;
      AllDomain<1> IV;
      Range<1> R(0,4,2);
      NewDomain_t::fillSlice(VD, b.domain(), R, IV);

      tester.out() << "VD = " << VD << std::endl;

      NewEngine<Brick2_t, SliceType_t>::Type_t v(b, VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < v.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << v(i,j) << " ";
	}
      tester.out() << std::endl;

      for (int i = 0; i < v.domain()[0].length(); ++i)
	for (int j = 0; j < v.domain()[1].length(); ++j)
	  v(i,j) = 0;

      tester.out() << "After setting v(i,j) = 0" << std::endl;

      tester.out() << "v = ";
      for (int i = 0; i < v.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < v.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << v(i,j) << " ";
	}
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

    }
    {
      // test 6

      tester.out() << "\n\nTesting 2D slice subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      typedef NewDomain2<Interval<1>, int> NewDomain_t;
      typedef NewDomain_t::SliceType_t SliceType_t;
      SliceType_t VD;
      Interval<1> IV(1,3);
      NewDomain_t::fillSlice(VD, b.domain(), IV, 3);

      tester.out() << "VD = " << VD << std::endl;

      NewEngine<Brick2_t, SliceType_t>::Type_t v(b, VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      for (int j = 0; j < v.domain()[0].length(); ++j)
	v(j) = 0;

      tester.out() << "After setting v(j) = 0" << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

    }
    {
      // test 7

      tester.out() << "\n\nTesting 2D slice subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      typedef NewDomain2<int, Interval<1> > NewDomain_t;
      typedef NewDomain_t::SliceType_t SliceType_t;
      SliceType_t VD;
      Interval<1> IV(1,3);
      NewDomain_t::fillSlice(VD, b.domain(), 2, IV);

      tester.out() << "VD = " << VD << std::endl;

      NewEngine<Brick2_t, SliceType_t>::Type_t v(b, VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      for (int j = 0; j < v.domain()[0].length(); ++j)
	v(j) = 0;

      tester.out() << "After setting v(j) = 0" << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;
    
    }
    {
      // test 8

      tester.out() << "\n\nTesting 2D slice subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      typedef NewDomain2<int, AllDomain<1> > NewDomain_t;
      typedef NewDomain_t::SliceType_t SliceType_t;
      SliceType_t VD;
      AllDomain<1> IV;
      NewDomain_t::fillSlice(VD, b.domain(), 2, IV);

      tester.out() << "VD = " << VD << std::endl;

      NewEngine<Brick2_t, SliceType_t>::Type_t v(b, VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      for (int j = 0; j < v.domain()[0].length(); ++j)
	v(j) = 0;

      tester.out() << "After setting v(j) = 0" << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

    }
    {
      // test 9

      tester.out() << "\n\nTesting 2D slice-range subsetting...\n\n";

      Interval<1> I(5);
      Interval<2> BD;
      BD[0] = I; BD[1] = I;

      Brick2_t b(BD);

      for (int j = 0; j < b.domain()[1].length(); ++j)
	for (int i = 0; i < b.domain()[0].length(); ++i)
	  b(i,j) = i + 10*j;

      tester.out() << "b.domain()     = " << b.domain() << std::endl;
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

      typedef NewDomain2<int, Range<1> > NewDomain_t;
      typedef NewDomain_t::SliceType_t SliceType_t;
      SliceType_t VD;
      Range<1> IV(0,4,2);
      NewDomain_t::fillSlice(VD, b.domain(), 2, IV);

      tester.out() << "VD = " << VD << std::endl;

      NewEngine<Brick2_t, SliceType_t>::Type_t v(b, VD);

      tester.out() << "v.domain()     = " << v.domain() << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      for (int j = 0; j < v.domain()[0].length(); ++j)
	v(j) = 0;

      tester.out() << "After setting v(j) = 0" << std::endl;

      tester.out() << "v = \n";
      for (int j = 0; j < v.domain()[0].length(); ++j)
	tester.out() << std::setw(2) << v(j) << " ";
      tester.out() << std::endl;

      tester.out() << "b = ";
      for (int i = 0; i < b.domain()[0].length(); ++i)
	{
	  tester.out() << "\n  ";
	  for (int j = 0; j < b.domain()[1].length(); ++j)
	    tester.out() << std::setw(2) << b(i,j) << " ";
	}
      tester.out() << std::endl;

    }
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
  int ret = tester.results("brickview_test2");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brickview_test2.cpp,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
