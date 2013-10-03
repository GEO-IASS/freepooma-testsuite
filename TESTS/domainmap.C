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
#include "Pooma/Pooma.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/Touches.h"
#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/DomainMap.h"
#include <iostream>

#include "Utilities/Tester.h"

int main(int argc, char *argv[]) 
{

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  tester.out() << "Starting domain map test." << std::endl << std::endl;

  typedef DomainMap<Interval<2>,int> DMap_t;
  typedef DMap_t::Value_t DPair_t;

  Interval<1> x(100),y(100);
  Interval<2> xy(x,y);
  DMap_t domainMap(xy);

  int foo[]={30,31,32,40,41,42,33,43};
  int i,j;
  for (i=0;i<10;++i)
  {
    Interval<1> x1(i*10,(i+1)*10-1);
    for (j=0; j<10; ++j)
    {
      Interval<1> y1(j*10,(j+1)*10-1);
      domainMap.insert(DPair_t(Interval<2>(x1,y1),j+i*10));
    }
  }
  domainMap.update();

  Interval<1> x2(32,48),y2(2,38);
  Interval<2> xy2(x2,y2);

  typedef DMap_t::Touch_t DMTouch_t;
  DMTouch_t touch = domainMap.touch(xy2);

  tester.out() << "finding domains that touch domain " << xy2 << std::endl;

  typedef DMap_t::touch_iterator iterator;
  iterator a;
  i=0;
  for (a=touch.first;a!=touch.second;++a)
  {    
    tester.out()<<"touches "<<(*a);
    tester.check("  :", (*a) == foo[i] );
  ++i;
  }

  int retval = tester.results("Domain Map");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: domainmap.cpp,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:33 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
