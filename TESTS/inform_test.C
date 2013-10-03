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

#include "Pooma/Pooma.h"
#include "Domain/Interval.h"
#include "Utilities/Inform.h"
#include <iostream>
#include <iomanip>
#include <ios>
#include "Utilities/Tester.h"

int main(int argc, char *argv[])
{
  // initialize Pooma
  Pooma::initialize(argc,argv);
 
  Pooma::Tester tester(argc, argv);
  // create some Inform instances
  Inform A("POOMA-II:A");
  Inform B("POOMA-II:B");

  // add another connection to stderr for B
  Inform::ID_t connection = B.open(std::cerr);

  // write A's output also to a file
  Inform::ID_t connection2 = A.open("inform_test.dat", Inform::out);

  // simple test prints, which should have leading and trailing blank lines
  A << "------" << std::endl;
  A << std::endl
    << "This should have a leading and following blank line."
    << std::endl << std::endl;
  A << "------" << std::endl;
  B << "------" << std::endl;
  B << std::endl
    << "This should have a leading and following blank line."
    << std::endl << std::endl;
  B << "------" << std::endl;

  // print some domains to this, to test output to ostream:
  Interval<1> X(1,5);
  A << "Interval X = " << X << ", with no endl, just flush";
  A << std::flush;

  // use some manipulators
  int val = 2;
  int val2 = 1234;
  A << std::setw(4) << std::setfill('#') << val << ": should be ###2"
    << std::endl;
  B << val2 << " = " << std::hex << val2 << " (hex), ";
  B << std::oct << val2 << " (oct)" << std::endl;

  // close B's second connection
  B.close(connection);
  B << "This line should only appear once." << std::endl;

  // should be some blank lines
  A << std::endl << std::endl
    << "There should be two blank lines, then this." << std::endl;

  // inform about file to check
  A << std::endl
    << "The file 'inform_test.dat' should contain copies of all";
  A << " the lines written to the 'A' stream." << std::endl;

  int res = tester.results();

  // finalize Pooma
  Pooma::finalize();

  return res;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: inform_test.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:17:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
