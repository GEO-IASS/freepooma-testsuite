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
// Include files
//-----------------------------------------------------------------------------

#include "IO/DiskMeta.h"
#include "IO/MetaTokenIterator.h"
#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h> // atoi

//-----------------------------------------------------------------------------
// DiskMeta member function implementations
//-----------------------------------------------------------------------------

DiskMeta::DiskMeta(const char *basename) 
  : mycontext_m(Pooma::context()), iocontext_m(0),
    dim_m(-1), fieldsPerRecord_m(-1), numRecords_m(-1), numFileSets_m(-1)
{
  PInsist(basename != 0, "No filename supplied");
  filename_m = std::string(basename) + ".meta";
}

DiskMeta::~DiskMeta() 
{ }

// Macro used to implement error handling policy. 

#define IOInsist(c,str) \
if (!abortOnError && !(c)) { errorMsg_m = str; return false; } \
else PInsist(c,str)

bool DiskMeta::open(bool abortOnError)
{
  // Open the input file
  
  if (mycontext_m == iocontext_m)
    {
      fin_m.open(filename_m.c_str());
      if (!fin_m) IOInsist(0, "Couldn't open .meta file.");
    }
  return true;
}

bool DiskMeta::read(bool abortOnError)
{
  if (mycontext_m == iocontext_m)
    {
      using std::string;

      // Declare a string buffer for getting each line.
      // Give it lots of memory so reallocs are unlikely. 

      string line;
      line.reserve(128*1024);
      string keyword;
      keyword.reserve(32);

      int dl = 0; // which line of the Domain are we on

      while (std::getline(fin_m,line))
        {
          // Each line has the form
          //   keyword [=] values

          // First find the keyword

          Pooma::MetaTokenIterator wpos(line);
          Pooma::MetaTokenIterator wend;

          // If there are no words on this line, advance to the next line

          if (wpos == wend) continue;

          // Found a word - it must be a keyword.

          keyword = *wpos++;

          // Process the keyword

          if (keyword == "Type")
            {
              // Type = type

              IOInsist(wpos != wend, "Invalid line");
              type_m = *wpos++;
            }
          else if (keyword == "Dim")
            {
              // Dim = N

              IOInsist(wpos != wend, "Invalid line");
              dim_m = atoi(wpos++->c_str());
              IOInsist(dim_m > 0 && dim_m < 8, "Invalid dimension");
            }
          else if (keyword == "Domain")
            {
              // Domain = first last stride

              IOInsist(dl < dim_m, "Too many Domain entries");
              IOInsist(wpos != wend, "Invalid line");
              int first, last, stride;
              first = atoi(wpos++->c_str());
              IOInsist(wpos != wend, "Invalid line");
              last  = atoi(wpos++->c_str());
              IOInsist(wpos != wend, "Invalid line");
              stride = atoi(wpos++->c_str());
              domain_m[dl] = Interval<1>(first,last);
              ++dl;
            }
          else if (keyword == "Fields")
            {
              // Fields = N

              IOInsist(wpos != wend, "Invalid line");
              fieldsPerRecord_m = atoi(wpos++->c_str());
            }
          else if (keyword == "Records")
            {
              // Records = N

              IOInsist(wpos != wend, "Invalid line");
              numRecords_m = atoi(wpos++->c_str());
            }
          else if (keyword == "SMPs")
            {
              // SMPs = N

              IOInsist(wpos != wend, "Invalid line");
              numFileSets_m = atoi(wpos++->c_str());
            }
          else if (keyword == "VnodesInRecord")
            {
              // VnodesInRecord = N0 N1 ... NN

              IOInsist(wpos != wend, "Invalid line");
              while (wpos != wend)
                {
                  patchesPerRecord_m.push_back(atoi(wpos++->c_str()));
                }
            }
          else if (keyword == "VnodeTally")
            {
              // VnodeTally = N0 N1 ... NN

              IOInsist(wpos != wend, "Invalid line");
              while (wpos != wend)
                {
                  patchTally_m.push_back(atoi(wpos++->c_str()));
                }
            }
        } // while (getline...)

      // Some consistency checks

      IOInsist(dl == dim_m, "File did not specify Dim domains");
      IOInsist(patchesPerRecord_m.size() == numRecords_m, 
               "VnodesInRecord incomplete.");
      IOInsist(patchTally_m.size() == numRecords_m, 
               "VnodeTally incomplete.");
    }

  // If we get here, the file passed the sanity checks. Return
  // success.

  return true;
}

const Interval<1> &DiskMeta::domain(int d) const
{
  if (mycontext_m == iocontext_m)
    {  
      PAssert(d >= 0 && d < dim_m);
      return domain_m[d];
    }
  else
    {
      return domain_m[d]; // garbage, most likely
    }
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DiskMeta.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.4 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
