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
// DomainRemoveOverlap test: Create Domains and use DomainRO function
//-----------------------------------------------------------------------------
#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Interval.h"
#include "Domain/DomainRemoveOverlap.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) 
{
  // Initialize POOMA and output stream, using Tester class
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv); 

  Interval<1> a(0,10),b(0,20),r(3,7);
  Interval<2> s(a,b),rr(r,r);
  typedef std::vector<Interval<2> > DomainList_t;
  DomainList_t res;

  tester.out() << " from " << s<< " remove "<<rr<<std::endl;

  res = DomainRemoveOverlap(s,rr);
 
  DomainList_t::iterator start = res.begin();
  DomainList_t::iterator end = res.end();
  for ( ; start!=end ; ++start)
    {
      tester.out() << *start << std::endl;
    }

  res.clear();
  Interval<2> k(Interval<1>(2,3),Interval<1>(-1,30));
  
  res = DomainRemoveOverlap(s,k);


  start = res.begin();
  end = res.end();
  for ( ; start!=end ; ++start)
    {
      tester.out() << *start << std::endl;
    }

  res.clear();

  Interval<2> k2(Interval<1>(2,3),Interval<1>(0,20));

  res = DomainRemoveOverlap(s,k2);
  tester.out() << " " <<std::endl;

  tester.out() <<"from " <<s<<" remove  "<<k2<<std::endl;

  tester.out() << " " <<std::endl;
  start = res.begin();
  end = res.end();
  for ( ; start!=end ; ++start)
    {
      tester.out() << *start << std::endl;
    }


  res.clear();

  Interval<2> k3(Interval<1>(-7,3),Interval<1>(-6,8));

  res = DomainRemoveOverlap(s,k3);
  tester.out() << " " <<std::endl;

  tester.out() <<"from "<< s<<"  remove "<<k3<<std::endl;

  tester.out() << " " <<std::endl;
  start = res.begin();
  end = res.end();
  for ( ; start!=end ; ++start)
    {
      tester.out() << *start << std::endl;
    }

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("DomainRO operations");
  Pooma::finalize();
  return retval;

}
		
		     
	      
	 
