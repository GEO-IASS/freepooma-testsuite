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

// include files
#include "Array/PrintArray.h"
#include "Utilities/Inform.h"

// Forward declarations

// ----------------------------------------------------------------------------
// dataBrowserPrintArray
// 
// A static (global) PrintArray, used by DataBrowser<ArrayPrintDataBrowser<> >
// to store persistent values of PrintArray formatting parameter
// settings. Global functions allow the debugger-user to interactively set the
// parameters, and use them for a set of interactive print calls:
// ----------------------------------------------------------------------------

PrintArray dataBrowserPrintArray;

// ----------------------------------------------------------------------------
// Global functions for setting formatting parameters, stored in
// dataBrowserPrintArray:
// ----------------------------------------------------------------------------

//int dbDomainWidth() const { return dataBrowserPrintArraydomainWidth(); }
// All the non-set functions were const like this; can't do it here somewhy.
int dbDomainWidth() { return dataBrowserPrintArray.domainWidth(); }
void dbSetDomainWidth(int val) { dataBrowserPrintArray.setDomainWidth(val); }
int dbDataWidth() { return dataBrowserPrintArray.dataWidth(); }
void dbSetDataWidth(int val) { dataBrowserPrintArray.setDataWidth(val); }
int dbDataPrecision() { return dataBrowserPrintArray.dataPrecision(); }
void dbSetDataPrecision(int val) {dataBrowserPrintArray.setDataPrecision(val);}
int dbCarReturn() { return dataBrowserPrintArray.carReturn();}
void dbSetCarReturn(int val) { dataBrowserPrintArray.setCarReturn(val); }
bool dbScientific() { return dataBrowserPrintArray.scientific(); }
void dbSetScientific(bool val) { dataBrowserPrintArray.setScientific(val); }
int dbSpacing() { return dataBrowserPrintArray.spacing(); }
void dbSetSpacing(int val) { dataBrowserPrintArray.setSpacing(val); }

// ----------------------------------------------------------------------------
// dataBrowserInform
// 
// A global Inform*, used by the dbprint() functions; this allows setting up
// for output to a file, or whatever, when using
// DataBrowser<ArrayPrintDataBrowser<> > printing, interactively, or from the
// debugger prompt. The secondary pointer dataBrowserInformBackup is used as
// temporary pointer value storage in the no-argment dbSetInform() function
// below.
// ----------------------------------------------------------------------------

Inform dataBrowserDefaultInform(NULL, 0);
Inform *dataBrowserInform = &dataBrowserDefaultInform;
Inform *dataBrowserInformBackup = &dataBrowserDefaultInform;

// ----------------------------------------------------------------------------
// dbSetInform()
// 
// A global function for setting the desired Inform object to be used by
// dbprint() functions, either interactively from the debugger prompt, or from
// within a code. For example, allows using an Inform object that writes to a
// file. The no-argument prototype allows swapping back and forth from a
// non-default Inform object to the default one, using dataBrowserInformBackup
// as temporary pointer value storage.
// ----------------------------------------------------------------------------

void dbSetInform(Inform &inform) { dataBrowserInform = &inform; }

// ----------------------------------------------------------------------------
// dbSwapInform()
// 
// Allows swapping back and forth from a non-default Inform object to the
// default one, using dataBrowserInformBackup as temporary pointer value
// storage.
// ----------------------------------------------------------------------------

void dbSwapInform() 
{ 
  if (dataBrowserInform == &dataBrowserDefaultInform) {
    if (dataBrowserInformBackup == dataBrowserInform) {
      return;
    } else {
      dataBrowserInform = dataBrowserInformBackup;
      dataBrowserInformBackup = &dataBrowserDefaultInform;
      return;
    }
  } else {
    dataBrowserInformBackup = dataBrowserInform;
    dataBrowserInform = &dataBrowserDefaultInform;
  }
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DataBrowser.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:27 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
