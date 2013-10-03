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
// Paws test 1: Send and Receive an int and double set of scalars
//              in conjunction with test 2.
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"

#if POOMA_PAWS
#include "Pooma/Paws.h"
#endif // POOMA_PAWS

#include "Utilities/Tester.h"


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream, using Tester class

  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);
  tester.out() << argv[0] << ": Paws scalar send/receive test A" << std::endl;
  tester.out() << "--------------------------------------------" << std::endl;

#if POOMA_PAWS

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  std::string str1("Sender's Orig");
  std::string str2("Sender's Modified");
  std::string origstr1 = str1;

  // Create a Paws connection

  tester.out() << "Creating PawsConnection object ..." << std::endl;
  Connection<Paws> *paws = new Connection<Paws>("test1", argc, argv);
  tester.out() << "Finished creating PawsConnection object." << std::endl;

  // Establish connections for the two scalars

  tester.out() << "Connecting s1 = " << s1 << " for output ..." << std::endl;
  paws->connectScalar("s1", s1, ConnectionBase::out);
  tester.out() << "Connecting s2 = " << s2 << " for input ..." << std::endl;
  paws->connectScalar("s2", s2, ConnectionBase::in);

  // Establish connection for the string

  tester.out() << "Connecting str1 = '"<< str1 <<"' for output ..."<<std::endl;
  paws->connectScalar("str1", str1, ConnectionBase::out);

  // Wait for everything to be ready to proceed

  tester.out() << "Waiting for ready signal ..." << std::endl;
  paws->ready();
  tester.out() << "Ready complete, moving on." << std::endl;

  // Modify s2, str1, and update

  s2 *= 2;
  tester.out() << "Updating current s1 = " << s1 << "s2 = " << s2;
  tester.out() << " and str1 = '" << str1 << "' ..." << std::endl;
  paws->update();

  // Report the results

  tester.out() << "Received update.  New values:" << std::endl;
  tester.out() << "  s1 = " << s1 << " (should be " << origs1 << ")\n";
  tester.out() << "  s2 = " << s2 << " (should be " << origs2 << ")\n";
  tester.out() << "str1 = " << str1 << " (should be " << origstr1 << ")\n";
  tester.out() << std::endl;
  tester.check("s1 OK", s1 == origs1);
  tester.check("s2 OK", s2 == origs2);
  tester.check("str1 OK", str1 == origstr1);

  // Delete PAWS connection, disconnecting us from the other code.

  tester.out() << "Deleting Connection<Paws> object ..." << std::endl;
  delete paws;
 
#else // POOMA_PAWS

  tester.out() << "Please configure with --paws to use this test code!"
	       << std::endl;

#endif // POOMA_PAWS

  // Finish up and report results

  tester.out() << "-------------------------------------------" << std::endl;
  int retval = tester.results("Paws scalar send/receive test A");
  Pooma::finalize();
  return retval;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: paws_test1.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:23 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
