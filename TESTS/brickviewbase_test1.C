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

#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Domain/Loc.h"
#include "Engine/BrickBase.h"

int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
#if POOMA_EXCEPTIONS
  try {
#endif

    using Pooma::BrickBase;
    using Pooma::BrickViewBase;
    
    tester.out() << "\nTesting non-sliced BrickViewBase." << std::endl;

// First we test BrickBase<1>'s constructors and resulting domains.

    Interval<1> I1(10);
    BrickBase<1> A1(I1);
    
    Interval<1> IV1(5,7);
    BrickViewBase<1> AV1(A1,IV1);
    
    tester.check(AV1.domain() == Interval<1>(0,2));
    tester.check(AV1.strides()[0] == 1);
    tester.check(AV1.first(0) == 0);
    tester.out() << "AV1's domain      = " << AV1.domain() << std::endl;
    tester.out() << "AV1's strides     = " << AV1.strides()[0] << std::endl;

// Next check the copy constructor and assignment operator:

    {
      BrickViewBase<1> DV1(AV1);
      tester.check(DV1.domain() == Interval<1>(0,2));
      tester.check(DV1.strides()[0] == 1);
      tester.check(DV1.first(0) == 0);
      tester.out() << "DV1's domain      = " << DV1.domain() << std::endl;
      tester.out() << "DV1's strides     = " << DV1.strides()[0] << std::endl;
    }

    {
      BrickViewBase<1> DV1(A1,Interval<1>(5,6)); // dummy - no default
      DV1 = AV1;
      tester.check(DV1.domain() == Interval<1>(0,2));
      tester.check(DV1.strides()[0] == 1);
      tester.check(DV1.first(0) == 0);
      tester.out() << "DV1's domain      = " << DV1.domain() << std::endl;
      tester.out() << "DV1's strides     = " << DV1.strides()[0] << std::endl;
    
// Now check the offset(Domain) functions.

      for (int i = 0; i < 10; ++i)
      {
        Loc<1> Itest1(i);
        tester.check(DV1.offset(Itest1) == i);
        tester.check(DV1.offset(i)      == i);
      }
    }
    
// Check for non-zero-based domains.

    Interval<1> J1(3,13);
    A1 = BrickBase<1>(J1);
    AV1 = BrickViewBase<1>(A1,Interval<1>(4,12));
    
    tester.check(AV1.domain() == Interval<1>(9));
    tester.check(AV1.strides()[0] == 1);
    tester.check(AV1.first(0) == 0);
    tester.out() << "AV1's domain       = " << AV1.domain() << std::endl;
    tester.out() << "AV1's strides      = " << AV1.strides()[0] << std::endl;
    
    for (int i = 0; i < 9; ++i)
      {
        Range<1> Itest1(i,11,2);
        tester.check(AV1.offset(Itest1) == i);
        tester.check(AV1.offset(i)      == i);
      }
      
    Interval<1> K1(-5,5);
    A1 = BrickBase<1>(K1);
    AV1 = BrickViewBase<1>(A1,Interval<1>(-1,1));
    
    tester.check(AV1.domain() == Interval<1>(3));
    tester.check(AV1.strides()[0] == 1);
    tester.check(AV1.first(0) == 0);
    tester.out() << "AV1's domain       = " << AV1.domain() << std::endl;
    tester.out() << "AV1's strides      = " << AV1.strides()[0] << std::endl;
    
    for (int i = 0; i < 3; ++i)
      {
        Loc<1> Itest1(i);
        tester.check(AV1.offset(Itest1) == i);
        tester.check(AV1.offset(i)      == i);
      }
      
    A1 = BrickBase<1>(K1);
    BrickViewBase<1> AV1F(A1,Range<1>(-1,1,2));
    
    tester.check(AV1F.domain() == Interval<1>(2));
    tester.check(AV1F.strides()[0] == 2);
    tester.check(AV1F.first(0) == 0);
    tester.out() << "AV1F's domain       = " << AV1F.domain() << std::endl;
    tester.out() << "AV1F's strides      = " << AV1F.strides()[0] << std::endl;
    
    for (int i = 0; i < 2; ++i)
      {
        Loc<1> Itest1(i);
        tester.check(AV1F.offset(Itest1) == i*AV1F.strides()[0]);
        tester.check(AV1F.offset(i)      == i*AV1F.strides()[0]);
      }


// Now repeat for 2D

    Interval<2> I2;
    I2[0] = I1; 
    I2[1] = I1;
    BrickBase<2> A2(I2);
    
    Interval<2> IV2;
    IV2[0] = IV1;
    IV2[1] = IV1;    
    BrickViewBase<2> AV2(A2,IV2);
    
    Interval<2> IV2_0;
    IV2_0[0] = Interval<1>(IV1.length());
    IV2_0[1] = IV2_0[0];
    
    tester.check(AV2.domain() == IV2_0);
    tester.check(AV2.strides()[0] == 1);
    tester.check(AV2.strides()[1] == 10);
    tester.check(AV2.first(0) == 0);
    tester.check(AV2.first(1) == 0);
    tester.out() << "AV2's domain      = " << AV2.domain() << std::endl;
    tester.out() << "AV2's strides     = " 
                 << AV2.strides()[0] << " "
                 << AV2.strides()[1] << std::endl;

// Next check the default constructor, copy constructor, 
// and assignment operator:

    {
      BrickViewBase<2> DV2(AV2);
      tester.check(DV2.domain() == IV2_0);
      tester.check(DV2.strides()[0] == 1);
      tester.check(DV2.strides()[1] == 10);
      tester.check(DV2.first(0) == 0);
      tester.check(DV2.first(1) == 0);
      tester.out() << "DV2's domain      = " << DV2.domain() << std::endl;
      tester.out() << "DV2's strides     = "
                   << DV2.strides()[0] << " "
                   << DV2.strides()[1] << std::endl;
    }

    {
      BrickViewBase<2> DV2(A2,Interval<2>());
      DV2 = AV2;
      tester.check(DV2.domain() == IV2_0);
      tester.check(DV2.strides()[0] == 1);
      tester.check(DV2.strides()[1] == 10);
      tester.check(DV2.first(0) == 0);
      tester.check(DV2.first(1) == 0);
      tester.out() << "DV2's domain      = " << DV2.domain() << std::endl;
      tester.out() << "DV2's strides     = "
                   << DV2.strides()[0] << " "
                   << DV2.strides()[1] << std::endl;
    
// Now check the offset(Domain) functions.

      for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
        {
          Loc<2> loc(i,j);
          tester.check(DV2.offset(loc) == i + j * DV2.strides()[1]);
          tester.check(DV2.offset(i,j) == i + j * DV2.strides()[1]);
      }
    }

// Check for non-zero-based domains.

    Interval<1> JJ(3,13);
    Interval<2> J2(JJ,JJ);
    A2 = BrickBase<2>(J2);
 
    Interval<1> JJV(5,10);
    Interval<2> JV2(JJV,JJV);
    Interval<2> JV2_0(Interval<1>(JJV.length()), Interval<1>(JJV.length()));
    AV2 = BrickViewBase<2>(A2,JV2);
    
    tester.check(AV2.domain() == JV2_0);
    tester.check(AV2.strides()[0] == 1);
    tester.check(AV2.strides()[1] == 11);
    tester.check(AV2.first(0) == 0);
    tester.check(AV2.first(1) == 0);
    tester.out() << "AV2's domain      = " << AV2.domain() << std::endl;
    tester.out() << "AV2's strides     = " 
                << AV2.strides()[0] << " "
                << AV2.strides()[1] << std::endl;
    
    for (int i = 3; i < 10; ++i)
      for (int j = 3; j < 10; ++j)
      {
        Range<1> Itest1(i,11,2);
        Range<1> Jtest1(j,10,3);
        Range<2> Itest2(Itest1,Jtest1);
        tester.check(AV2.offset(Itest2) == i + j * AV2.strides()[1]);
        tester.check(AV2.offset(i,j)    == i + j * AV2.strides()[1]);
      }
      
    Interval<1> KK(-5,5);
    Interval<2> K2(KK,KK);
    A2 = BrickBase<2>(K2);
    
    Interval<1> KV1(-2,2);
    Interval<2> KV2(KV1,KV1);
    Interval<2> KV2_0(Interval<1>(KV1.length()),Interval<1>(KV1.length()));
    
    AV2 = BrickViewBase<2>(A2,KV2);
        
    tester.check(AV2.domain() == KV2_0);
    tester.check(AV2.strides()[0] == 1);
    tester.check(AV2.strides()[1] == A2.strides()[1]);
    tester.check(AV2.first(0) == 0);
    tester.check(AV2.first(1) == 0);
    tester.out() << "AV2's domain       = " << AV2.domain() << std::endl;
    tester.out() << "AV2's strides      = " 
                 << AV2.strides()[0] << " "
                 << AV2.strides()[1] << std::endl;
    
    for (int j = -2; j < 3; ++j)
      for (int i = -2; i < 3; ++i)
      {
        Loc<2> loc(i,j);
        tester.check(AV2.offset(loc) == i + j * AV2.strides()[1]);
        tester.check(AV2.offset(i,j) == i + j * AV2.strides()[1]);
      }
        
    A2 = BrickBase<2>(K2);
    
    Range<1> RV1(-1,1,2);
    Range<2> RV2(RV1,RV1);
    KV2_0 = Interval<2>(Interval<1>(RV1.length()),Interval<1>(RV1.length()));

    BrickViewBase<2> AV2F(A2,RV2);
    tester.check(AV2F.domain() == KV2_0);
    tester.check(AV2F.strides()[0] == 2);
    tester.check(AV2F.strides()[1] == 2 * A2.strides()[1]);
    tester.check(AV2F.first(0) == 0);
    tester.check(AV2F.first(1) == 0);
    tester.out() << "AV2F's domain       = " << AV2F.domain() << std::endl;
    tester.out() << "AV2F's strides      = " 
                 << AV2F.strides()[0] << " "
                 << AV2F.strides()[1] << std::endl;
    
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
      {
        Loc<2> loc(i,j);
        int off1 = AV2F.offset(loc);
        int off2 = AV2F.offset(i,j);
        int off3 = i*AV2F.strides()[0] + j*AV2F.strides()[1];
        tester.check(off1 == off3);
        tester.check(off2 == off3);
      }


    Interval<1> L1(-2,2);
    Interval<7> L7(L1,L1,L1,L1,L1,L1,L1);
    BrickBase<7> A7(L7);
    
    Interval<1> LV1(-1,1);
    Interval<7> LV7(LV1,LV1,LV1,LV1,LV1,LV1,LV1);
    Interval<1> LV0(0,2);
    Interval<7> LV7_0(LV0,LV0,LV0,LV0,LV0,LV0,LV0);

    BrickViewBase<7> AV7(A7,LV7);
    
    tester.check(AV7.domain() == LV7_0);
    
    for (int d = 0; d < 7; ++d)
    {
      tester.check(AV7.strides()[d] == A7.strides()[d]);
      tester.check(AV7.first(d) == 0);
    }
    tester.out() << "AV7's domain       = " << AV7.domain() << std::endl;
    tester.out() << "AV7's strides      = " ;
    for (int d = 0; d < 7; ++d)
    {
      tester.out() << AV7.strides()[d] << " ";
    }
    tester.out() << std::endl;
    
    for (int i7 = 0; i7 < 3; ++i7)
      for (int i6 = 0; i6 < 3; ++i6)
        for (int i5 = 0; i5 < 3; ++i5)
          for (int i4 = 0; i4 < 3; ++i4)
            for (int i3 = 0; i3 < 3; ++i3)
              for (int i2 = 0; i2 < 3; ++i2)
                for (int i1 = 0; i1 < 3; ++i1)
                    {
                      Loc<7> loc(i1,i2,i3,i4,i5,i6,i7);
                      int off1 =     i1 + 
                                     i2 * AV7.strides()[1] +
                                     i3 * AV7.strides()[2] +
                                     i4 * AV7.strides()[3] +
                                     i5 * AV7.strides()[4] +
                                     i6 * AV7.strides()[5] +
                                     i7 * AV7.strides()[6] ;
                      tester.check(AV7.offset(loc)                  == off1 );
                      tester.check(AV7.offset(i1,i2,i3,i4,i5,i6,i7) == off1);
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
  int ret = tester.results("brickviewbase_test1");
  Pooma::finalize();
  return ret;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: brickviewbase_test1.cpp,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
