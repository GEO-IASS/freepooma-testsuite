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
// Lux demo: Display a series of byte-fields
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Pooma/Lux.h"
#include "Pooma/Arrays.h"
#include "Pooma/Domains.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[])
{
  // Initialize POOMA and output stream

  Pooma::initialize(argc, argv);
  Inform msg(argv[0]);

#if POOMA_LUX

  int i, sizeX = 0, sizeY = 0, sizeZ = 0;
  int numfiles = 0;

  // Get array size from command-line args

  for (i=1; i < argc; ++i)
    {
      std::string item(argv[i]);
      if (item == "-x" && i < (argc - 1))
	sizeX = atoi(argv[++i]);
      else if (item == "-y" && i < (argc - 1))
	sizeY = atoi(argv[++i]);
      else if (item == "-z" && i < (argc - 1))
	sizeZ = atoi(argv[++i]);
      else if (i < (argc - 1))
	numfiles++;
    }

  // Do some error-checking

  if (sizeX < 1 || sizeY < 1 || sizeZ < 1)
    {
      msg << "Bad size values, all must be > 0." << std::endl;
      exit(1);
    }

  if (numfiles < 1)
    {
      msg << "You must specify some files to display." << std::endl;
      exit(1);
    }

  // Create an array to display the data

  msg << "Initializing array ..." << std::endl;
  Interval<3> domain(sizeX, sizeY, sizeZ);
  Array<3, unsigned char, Brick> data(domain);
  data = 0;
  Pooma::blockAndEvaluate();

  // Create a Lux connection, and connect up the storage array

  msg << "Creating LuxConnection object ..." << std::endl;
  Connection<Lux> lux(argv[0]);

  msg << "Connecting data storage array ..." << std::endl;
  lux.connect("data", data);

  // Do, in a loop, updates of the arrays, and redisplay/interact.

  int count = 0;
  for (i=1; i < argc; ++i)
    {
      std::string item(argv[i]);
      if (item == "-x" || item == "-y" || item == "-z")
	{
	  ++i;
	}
      else
	{
	  msg << "Reading data from file '" << item << "' ..." << std::endl;

	  // Open the file

	  FILE *f = fopen(argv[i], "r");
	  if (f == 0)
	    {
	      msg << "Error opening file '" << item << "'." << std::endl;
	      break;
	    }

	  // Read in the data

	  if (fread(&(data(data.firsts())), domain.size(), 1, f) != 1)
	    {
	      msg << "Error reading from file '" << item << "'." << std::endl;
	      fclose(f);
	      break;
	    }

	  // Sanity check

	  int s = sum(data);
	  msg << "Sum of data read = " << s << std::endl;
	  msg << "Middle 1D slice of data:" << std::endl;
	  msg << data(Interval<1>(sizeX), sizeY/2, sizeZ/2) << std::endl;

	  // Update the display

	  msg << "Updating the display, for file " << ++count << " out of ";
	  msg << numfiles << " ..." << std::endl;
	  lux.ready();
	}
    }

  // Delete LUX connection, closing the window.

  msg << "Closing LUX connection ..." << std::endl;
  lux.close();

#else // POOMA_LUX

  msg << "Please configure with --lux to use this test code!"
      << std::endl;

#endif // POOMA_LUX

  // Finish up

  Pooma::finalize();
  return 0;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: bytefield.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:18 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
