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
// Array test 17: expression engine revisited.
//-----------------------------------------------------------------------------

// Turn on bounds checking.

#undef POOMA_BOUNDS_CHECK
#undef POOMA_BOUNDS_CHECK_DEFAULT
#define POOMA_BOUNDS_CHECK POOMA_YES
#define POOMA_BOUNDS_CHECK_DEFAULT POOMA_TRUE

// Include files

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/Range.h"
#include "Engine/BrickEngine.h"
#include "Engine/ExpressionEngine.h"
#include "Array/Array.h"

#include <iostream>

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  
  Interval<1> I(5), J(1,3);
  Array<1> a(I), b(I), c(I-10);
  
  for (int i = 0; i < 5; i++)
    {
      b(i) = i;
      c(i - 10) = -2.0 * i * i;
    }

#if POOMA_EXCEPTIONS    
  try
#endif
    {
      a = -4.0 * (b + c());
      tester.out() << a << std::endl;
      a = 0;

      a(J) = -4.0 * (b + c())(J);
      tester.out() << a << std::endl;
    }
#if POOMA_EXCEPTIONS    
  catch(Pooma::Assertion &err)
    {
      err.print(tester.out());
      tester.out() << std::endl;
    }
#endif

  Range<1> R(0,4,2);
  Interval<2> II(I, I);
  Array<2> aa(II), bb(II), cc(I-10,I);
  
  for (int j = 0; j < 5; j++)
    for (int k = 0; k < 5; k++)
      {
        bb(j, k) = j + k;
        cc(j - 10, k) = -j + k * k;
      }

#if POOMA_EXCEPTIONS    
  try
#endif
    {
      aa = -4.0 * (bb + cc());
      tester.out() << aa << std::endl;
      aa = 0;

      aa(I, J) = -4.0 * (bb + cc())(I, J);
      tester.out() << aa << std::endl;
      aa = 0;

      aa(2, J) = -4.0 * (bb + cc())(2, J);
      tester.out() << aa << std::endl;
    }
#if POOMA_EXCEPTIONS    
  catch(Pooma::Assertion &err)
    {
      err.print(tester.out());
      tester.out() << std::endl;
    }
#endif
  
  int retval = tester.results("array_test17");    
  Pooma::finalize();
  
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: array_test17.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:14 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
