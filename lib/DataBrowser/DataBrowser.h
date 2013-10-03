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
// Class:
//   DataBrowser
// Global Function Templates:
//   dbprint
//-----------------------------------------------------------------------------

#ifndef POOMA_DATABROWSER_DATABROWSER_H
#define POOMA_DATABROWSER_DATABROWSER_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
// 
// Classes:
//
// DataBrowser : ???
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 
// Global Function Templates:
// 
// dbprint() : ???
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------
#include "Array/PrintArray.h"
#include "DataBrowser/RangeMaker.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------
class Inform;

// ----------------------------------------------------------------------------
// Prototypes from DataBrowser.cmpl.cpp, needed for linking.

// Global PrintArray, for storing persistent formatting parameter set:
extern PrintArray dataBrowserPrintArray;

// Global functions for setting formatting parameters, stored in
// dataBrowserPrintArray:
int dbDomainWidth();
void dbSetDomainWidth(int val);
int dbDataWidth();
void dbSetDataWidth(int val);
int dbDataPrecision();
void dbSetDataPrecision(int val);
int dbCarReturn();
void dbSetCarReturn(int val);
bool dbScientific();
void dbSetScientific(bool val);
int dbSpacing();
void dbSetSpacing(int val);

// Global Inform*, for storing address of desired output Inform object:
extern Inform *dataBrowserInform;

// Global functions for setting dataBrowserInform:
void dbSetInform(Inform &inform);
void dbSwapInform();
// ----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
//
// Full Description:
//
// Classes:
//
// DataBrowser:
// 
// DataBrowser is a functor class ....
// 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 
// Global Function Templates:
// 
// dbprint() : ....
//-----------------------------------------------------------------------------

// General template: no implementation
template<class DataBrowserTag>
class DataBrowser
{
public:
  DataBrowser() {
    CTAssert(0);
  }
};

// For now, only support IJK-indexable array ASCII output to an Inform stream.
// In principle, DataBrowserTag classes could describe other browser types such
// as TK windows with spreadsheet-like or pixmap representations. Whether to
// use DataBrowser for that or use something else like DataConnect is not yet
// resolved.

// Tag Class to use for DataBrowserTag representing Array-type ASCII text
// output:

template<int Dim>
class ArrayPrintDataBrowser;


// Partial specializations of DataBrowser<ArrayPrintDataBrowser<Dim> > for
// specific values of Dim. These are useable for printing values from various
// POOMA container types: Arrays, DynamicArrays, and Fields centered on
// logically-rectilinear meshes using DiscreteGeometry.

template<int Dim>
class DataBrowser<ArrayPrintDataBrowser<Dim> >
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and constants

  // none?

  // --------------------------------------------------------------------------
  // Constructors.

  // Default ctor. Sets up typical Inform and PrintArray objects.
  DataBrowser() : inform_m(dataBrowserInform), pa_m() { }

  // Construct with a container:
  template<class Container>
  DataBrowser(const Container &c)
    : totalDomain_m(c.totalDomain()), view_m(c.totalDomain()), 
    inform_m(dataBrowserInform), pa_m()
  { }

  // Construct with a container and another domain to store as a view:
  template<class Container, class DomainType>
  DataBrowser(const Container &c, const DomainType &domain)
    : totalDomain_m(c.totalDomain()), view_m(domain), 
    inform_m(dataBrowserInform), pa_m()
  { }

  // --------------------------------------------------------------------------
  // Destructor.
  ~DataBrowser() { };

  // --------------------------------------------------------------------------
  // Methods.

  // Print the whole container:
  template<class Container>
  void print(const Container &c) { pa_m.print(*inform_m, c); }

  // Print a view of the whole container:
  template<class Container, class DomainType>
  void print(const Container &c, const DomainType &d) { 
    pa_m.print(*inform_m, c, d);
  }

  // Set all the formatting parameters from the example PrintArray's values:
  void setFormatParameters(PrintArray &pa) { pa_m.setFormatParameters(pa); }

  // Reset the Inform pointer:
  void setInform(Inform &inform) { inform_m = &inform; }

protected:

private:

  // POOMA PrintArray object used for the actual output:
  PrintArray pa_m;

  // (Current working) total domain:
  Range<Dim> totalDomain_m;
  
  // (Current working) view domain:
  Range<Dim> view_m;

  // POOMA Inform object for output:
  Inform *inform_m;
};


// ----------------------------------------------------------------------------
// 
// Global Function Templates:
//
// ----------------------------------------------------------------------------

// Print all elements in the container:
template<class Container>
void dbprint(const Container &c)
{
  DataBrowser<ArrayPrintDataBrowser<Container::dimensions> > 
    db(c, c.totalDomain());
  db.setFormatParameters(dataBrowserPrintArray);
  db.print(c);
}

// Print a specified view of  elements in the container:
template<class Container, class DomainType>
void dbprint(const Container &c, const DomainType &domain)
{
  DataBrowser<ArrayPrintDataBrowser<Container::dimensions> > 
    db(c, domain);
  db.setFormatParameters(dataBrowserPrintArray);
  db.print(c, domain);
}

// Print a specified Range<Dim> view of elements in the container, specified
// using a list of integers. This requires the RangeMaker monkey business. To
// support dimensionalitys 1-7, with sensible numbers of integer arguments for
// each, prototypes for 1-21 integer arguments are needed, excluding the
// numbers {11,13,16,17,19,20}.

template<class Container>
void dbprint(const Container &c, 
             const int  &i0)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 1>()(i0);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 2>()(i0,   i1);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 3>()(i0,   i1,  i2);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 4>()(i0,   i1,  i2,  i3);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 5>()(i0,   i1,  i2,  i3,  i4);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 6>()(i0,   i1,  i2,  i3,  i4,  i5);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 7>()(i0,   i1,  i2,  i3,  i4,  i5,  i6);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 8>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions, 9>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,10>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,12>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11,
             const int &i12, const int &i13)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,14>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11, i12, i13);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11,
             const int &i12, const int &i13, const int &i14)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,15>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11, i12, i13,
                                           i14);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11,
             const int &i12, const int &i13, const int &i14, const int &i15)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,16>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11, i12, i13,
                                           i14, i15);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11,
             const int &i12, const int &i13, const int &i14, const int &i15,
             const int &i16, const int &i17)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,18>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11, i12, i13,
                                           i14, i15, i16, i17);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11,
             const int &i12, const int &i13, const int &i14, const int &i15,
             const int &i16, const int &i17, const int &i18, const int &i19)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,20>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11, i12, i13,
                                           i14, i15, i16, i17, i18, i19);
  dbprint(c, domain);
}

template<class Container>
void dbprint(const Container &c, 
             const int  &i0, const int  &i1, const int  &i2, const int  &i3, 
             const int  &i4, const int  &i5, const int  &i6, const int  &i7, 
             const int  &i8, const int  &i9, const int &i10, const int &i11,
             const int &i12, const int &i13, const int &i14, const int &i15,
             const int &i16, const int &i17, const int &i18, const int &i19,
             const int &i20)
{
  Range<Container::dimensions> domain = 
    RangeMaker<Container::dimensions,21>()(i0,   i1,  i2,  i3,  i4,  i5,  i6,  
                                           i7,   i8,  i9, i10, i11, i12, i13,
                                           i14, i15, i16, i17, i18, i19, i20);
  dbprint(c, domain);
}

#endif     // POOMA_DATABROWSER_DATABROWSER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: DataBrowser.h,v $   $Author: richard $
// $Revision: 1.6 $   $Date: 2004/11/01 18:16:27 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
