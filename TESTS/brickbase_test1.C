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
// BrickBase test code
//-----------------------------------------------------------------------------

// This had better stand on its own...

#include "Engine/BrickBase.h"

#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Loc.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif

    using Pooma::BrickBase;
    
    tester.out() << "\nTesting BrickBase." << std::endl;

// First we test BrickBase<1>'s constructors and resulting domains.

    Interval<1> I1(10);
    BrickBase<1> A1(I1);

    tester.check(A1.domain() == I1);
    tester.check(A1.strides()[0] == 1);
    tester.check(A1.offset() == 0);
    tester.out() << "A1's domain  = " << A1.domain() << std::endl;
    tester.out() << "A1's strides = " << A1.strides()[0] << std::endl;

// Next check the default constructor, copy constructor, and assignment operator:

    {
      BrickBase<1> D1(A1);
      tester.check(D1.domain() == I1);
      tester.check(D1.strides()[0] == 1);
      tester.check(D1.offset() == 0);
      tester.out() << "D1's domain  = " << D1.domain() << std::endl;
      tester.out() << "D1's strides = " << D1.strides()[0] << std::endl;
    }

    {
      BrickBase<1> D1;
      D1 = A1;
      tester.check(D1.domain() == I1);
      tester.check(D1.strides()[0] == 1);
      tester.check(D1.offset() == 0);
      tester.out() << "D1's domain  = " << D1.domain() << std::endl;
      tester.out() << "D1's strides = " << D1.strides()[0] << std::endl;
    
// Now check the offset(Domain) functions.

      for (int i = 0; i < 10; ++i)
      {
        Interval<1> Itest1(i,9);
        tester.check(D1.offset(Itest1)  == i);
        tester.check(D1.offset(i)       == i);
        tester.check(D1.offset0(Itest1) == i);
        tester.check(D1.offset0(i)      == i);
      }
    }

// Check for non-zero-based domains.

    Interval<1> J1(3,13);
    A1 = BrickBase<1>(J1);
    
    tester.check(A1.domain() == J1);
    tester.check(A1.strides()[0] == 1);
    tester.check(A1.offset() == -3);
    tester.out() << "A1's domain  = " << A1.domain() << std::endl;
    tester.out() << "A1's strides = " << A1.strides()[0] << std::endl;
    tester.out() << "A1's offset  = " << A1.offset() << std::endl;
    
    for (int i = 3; i < 10; ++i)
      {
        Range<1> Itest1(i,11,2);
        tester.check(A1.offset(Itest1)  == i-3);
        tester.check(A1.offset(i)       == i-3);
        tester.check(A1.offset0(Itest1) == i);
        tester.check(A1.offset0(i)      == i);
      }
      
    Interval<1> K1(-5,5);
    A1 = BrickBase<1>(K1);
    
    tester.check(A1.domain() == K1);
    tester.check(A1.strides()[0] == 1);
    tester.check(A1.offset() == 5);
    tester.out() << "A1's domain  = " << A1.domain() << std::endl;
    tester.out() << "A1's strides = " << A1.strides()[0] << std::endl;
    tester.out() << "A1's offset  = " << A1.offset() << std::endl;
    
    for (int i = -5; i < 6; ++i)
      {
        Loc<1> Itest1(i);
        tester.check(A1.offset(Itest1)  == i+5);
        tester.check(A1.offset(i)       == i+5);
        tester.check(A1.offset0(Itest1) == i);
        tester.check(A1.offset0(i)      == i);
      }

// Now repeat for 2D

    Interval<2> I2;
    I2[0] = I1; 
    I2[1] = I1;
    BrickBase<2> A2(I2);
    
    tester.check(A2.domain() == I2);
    tester.out() << "A2's domain = " << A2.domain() << std::endl;
    
    tester.check(A2.domain() == I2);
    tester.check(A2.strides()[0] == 1);
    tester.check(A2.strides()[1] == 10);
    tester.check(A2.offset() == 0);
    tester.out() << "A2's domain  = " << A2.domain() << std::endl;
    tester.out() << "A2's strides = " 
                 << A2.strides()[0] << " "
                 << A2.strides()[1] << std::endl;

// Next check the default constructor, copy constructor, and assignment operator:

    {
      BrickBase<2> D2(A2);
      tester.check(D2.domain() == I2);
      tester.check(D2.strides()[0] == 1);
      tester.check(D2.strides()[1] == 10);
      tester.check(D2.offset() == 0);
      tester.out() << "D2's domain  = " << D2.domain() << std::endl;
      tester.out() << "D2's strides = "
                   << D2.strides()[0] << " "
                   << D2.strides()[1] << std::endl;
    }

    {
      BrickBase<2> D2;
      D2 = A2;
      tester.check(D2.domain() == I2);
      tester.check(D2.strides()[0] == 1);
      tester.check(D2.strides()[1] == 10);
      tester.check(D2.offset() == 0);
      tester.out() << "D2's domain  = " << D2.domain() << std::endl;
      tester.out() << "D2's strides = "
                   << D2.strides()[0] << " "
                   << D2.strides()[1] << std::endl;
    
// Now check the offset(Domain) functions.

      for (int i = 0; i < 10; ++i)
      {
        Interval<1> Itest1(i,9);
        Interval<2> Itest2(Itest1,I1);
        tester.check(D2.offset(Itest2)        == i);
        tester.check(D2.offset(i,I1.first())  == i);
        tester.check(D2.offset0(Itest2)       == i);
        tester.check(D2.offset0(i,I1.first()) == i);
      }
    }

// Check for non-zero-based domains.

    Interval<1> JJ(3,13);
    Interval<2> J2(JJ,JJ);
    A2 = BrickBase<2>(J2);
    
    tester.check(A2.domain() == J2);
    tester.check(A2.strides()[0] == 1);
    tester.check(A2.offset() == -3 - 3 * A2.strides()[1]);
    tester.out() << "A2's domain  = " << A2.domain() << std::endl;
    tester.out() << "A2's strides = " 
                << A2.strides()[0] << " "
                << A2.strides()[1] << std::endl;
    tester.out() << "A2's offset  = " << A2.offset() << std::endl;
    
    for (int i = 3; i < 10; ++i)
      {
        Range<1> Itest1(i,11,2);
        Range<2> Itest2(Itest1,JJ);
        tester.check(A2.offset(Itest2)        == i - 3);
        tester.check(A2.offset(i,JJ.first())  == i - 3);
        tester.check(A2.offset0(Itest2)       == i - 3 - A2.offset());
        tester.check(A2.offset0(i,JJ.first()) == i - 3 - A2.offset());
      }
      
    Interval<1> KK(-5,5);
    Interval<2> K2(KK,KK);
    A2 = BrickBase<2>(K2);
    
    tester.check(A2.domain() == K2);
    tester.check(A2.strides()[0] == 1);
    tester.check(A2.strides()[1] == A2.domain()[0].length());
    tester.check(A2.offset() == 5 + 5 * A2.strides()[1]);
    tester.out() << "A2's domain  = " << A2.domain() << std::endl;
    tester.out() << "A2's strides = " 
                 << A2.strides()[0] << " "
                 << A2.strides()[1] << std::endl;
    tester.out() << "A2's offset  = " << A2.offset() << std::endl;
    
    int off = 0;
    for (int j = -5; j < 6; ++j)
      for (int i = -5; i < 6; ++i)
      {
        Loc<2> Itest2(i,j);
        tester.check(A2.offset(Itest2) == off);
        tester.check(A2.offset(i,j)    == off);
        ++off;
      }
    
    off = 0;
    for (int j = 0; j < A2.domain()[1].length(); ++j)
      for (int i = 0; i < A2.domain()[0].length(); ++i)
      {
        Loc<2> Itest2(i,j);
        tester.check(A2.offset0(Itest2) == off);
        tester.check(A2.offset0(i,j)    == off);
        ++off;
      }
    
    Interval<1> L1(-1,1);
    Interval<7> L7(L1,L1,L1,L1,L1,L1,L1);
    BrickBase<7> A7(L7);
    
    off = 0;
    for (int i7 = -1; i7 < 2; ++i7)
      for (int i6 = -1; i6 < 2; ++i6)
        for (int i5 = -1; i5 < 2; ++i5)
          for (int i4 = -1; i4 < 2; ++i4)
            for (int i3 = -1; i3 < 2; ++i3)
              for (int i2 = -1; i2 < 2; ++i2)
                for (int i1 = -1; i1 < 2; ++i1)
                    {
                      Loc<7> loc(i1,i2,i3,i4,i5,i6,i7);
                      tester.check(A7.offset(loc)                  == off);
                      tester.check(A7.offset(i1,i2,i3,i4,i5,i6,i7) == off);
                      tester.check(
                        A7.offset0(i1,i2,i3,i4,i5,i6,i7) + A7.offset() == off);
                      ++off;
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
  int ret = tester.results("brickbase_test1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brickbase_test1.cpp,v $   $Author: richi $
// $Revision: 1.4 $   $Date: 2004/11/29 16:21:34 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
