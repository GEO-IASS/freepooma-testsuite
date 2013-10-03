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
// Classes: 
//   PrintArray - information and routines used to print out an
//                Array, or a view of one, to a stream.
//-----------------------------------------------------------------------------

#ifndef POOMA_ARRAY_PRINT_ARRAY_H
#define POOMA_ARRAY_PRINT_ARRAY_H

/** @file
 * @ingroup Array
 * @brief
 *    Extract the elements of an Array and print out the
 *    contents to a stream with some nice formatting.
 *
 *  The user can select
 *  how many elements to print per line, the precision, format, etc.
 *  This is done by constructing an PrintArray, and calling the
 *  'print' method with the stream to print to and the array to print.
 *  Also allows passing in a domain argument specifying a view (subset)
 *  of the Array elements; the output shows the global-index-space index
 *  values for the view.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Utilities/PAssert.h"
#include "Engine/IsValidLocation.h"
#include <iomanip>


//-----------------------------------------------------------------------------
// Forward References:
//-----------------------------------------------------------------------------

class PrintArray;


//-----------------------------------------------------------------------------
// Open POOMA namespace:
//-----------------------------------------------------------------------------

// namespace POOMA {

/**
 * PerformPrintArray struct: a simple wrapper around the templated 'print'
 * method of PrintArray.  This is here as a workaround to a Metrowerks
 * problem that prevents having templated member functions defined as
 * out-of-line.  This struct defines one static method 'print', which
 * prints an array to a stream (the two template parameters).  The Dim
 * parameter is used to specialize this to the case of a 1D domain.  This is
 * called by the 'print' templated member function of PrintArray.
 */

template<class S, class A, int Dim, class DomainType>
struct PerformPrintArray
{
  static void print(const PrintArray &, S &, const A &, const DomainType &);
};

template<class S, class A, class DomainType>
struct PerformPrintArray<S, A, 1, DomainType>
{
  static void print(const PrintArray &, S &, const A &, const DomainType &);
};


//-----------------------------------------------------------------------------
// PrintArray class declaration
//-----------------------------------------------------------------------------

/**
 * PrintArray is a utility program for Array objects, used
 * to print (nicely) the contents of an array to a provided stream.  All
 * the data in the array will be printed in ascii format; if you want to
 * display just a portion of an array, take a view of the array and give
 * that to the print() method. This will display the indexes as zero-based,
 * unit-stride. To avoid this, pass in both the Array and a domain object
 * (such as a Range<>) specifying the view, and the global-index-space
 * indices will be displayed.
 *
 * CONSTRUCTING A PrintArray:
 * --------------------------
 * When you construct a PrintArray, you can give it several format parameters
 * to control how to display the array.  These parameters are, in the
 * order they are given in the constructor (all have default values):
 *
 * domain width: the number of spaces that will be used to print out domain
 *               numbers.  If it is 3, say, then domains will be printed like
 *               [003:008]
 *    - query setting with "domainWidth()"
 *    - change setting with "setDomainWidth(int newval)"
 *
 * data width: the number of spaces, total, used to print out data values.
 *    - query setting with "dataWidth()"
 *    - change setting with "setDataWidth(int newval)"
 *
 * data precision: the number of digits past the decimal point displayed
 *                 when data values are printed
 *    - query setting with "dataPrecision()"
 *    - change setting with "setDataPrecision(int newval)"
 *
 * carriage return: if this value is < 0, then for each row of values from
 *                  the array, no carriage return is printed until the end
 *                  of the row.  If this number is > 0, it represents the
 *                  maximum number of values that will be printed before a
 *                  return.  Only the first set of numbers for a row of
 *                  the array will have a domain prefix included.
 *   - query setting with "carReturn()"
 *   - change setting with "setCarReturn(int newval)"
 *
 * scientific notation: a boolean flag, if true numbers are printed using
 *                      scientific notation, e.g. 10e-14
 *   - query setting with "scientific()"
 *   - change setting with "setScientific(bool newflag)"
 *
 * data spacing: the number of spaces to print between values
 *   - query setting with "spacing()"
 *   - change setting with "setSpacing(int newval)"
 *
 * PRINTING ARRAY OBJECTS WITH A PrintArray:
 * -----------------------------------------
 * PrintArray is not templated, so that you can reuse the same formatter
 * for different array's.  It has one templated member function 'print':
 *
 *   template<class S, class A>
 *   void print(S &s, const A &a) const
 *
 *   template<class S, class A, class DomainType>
 *   void print(S &s, const A &a, DomainType &d) const
 *
 * where 'S' must be a class with an ostream-like interface (such as cout, or
 * an Inform object), 'A' must be a class with an array interface, and
 * 'DomainType' must be a class with a Domain interface.  The first 'print'
 * prototype will take data from 'a' and print it to the stream using the
 * current format settings. The second prototype will print the view of 'a'
 * specified by 'd' and print it, showing the global index values from the
 * full domain of 'a'.
 *
 * 1-D arrays just have the one row printed, perhaps on multiple lines
 * if carReturn is > 0.  2-D arrays are printed as a table, with each
 * line prefixed by the domain it includes.  For example:
 *   [00:02][00] =          0          0          0
 *   [00:02][01] =          0          0          0
 * prints the values for [0:2][0:1] of a 2D array.  3-D and higher arrays
 * have a sequence of 2-D slices printed for them, each slice separated
 * by a line indicating which slice it is, and a separator.  Example:
 *   [0:2:1][0:3:1][2]:
 *   ----------------------------------------------------
 *   [00:02][00][02] =          0          0          0
 *   [00:02][01][02] =          0          0          0
 *   [00:02][02][02] =          0          0          0
 *   [00:02][03][02] =          0          0          0
 */ 

class PrintArray
{
public:
  //============================================================
  // Constructors.
  //============================================================

  //------------------------------------------------------------
  // Construct an PrintArray object with an array, and the settings for
  // how the array should be printed.

  PrintArray(int domainWidth = 3, int dataWidth = 10,
	     int dataPrecision = 4, int carReturn = -1,
	     bool scientific = false, int spacing = 1)
    : domainwidth_m(domainWidth), datawidth_m(dataWidth), 
      dataprecision_m(dataPrecision), carreturn_m(carReturn),
      spacing_m(spacing), scientific_m(scientific)
    {
      PAssert(domainwidth_m > 0);
      PAssert(datawidth_m > 0);
      PAssert(dataprecision_m > 0);
      PAssert(spacing_m >= 0);
    }

  //------------------------------------------------------------
  // Copy constructor.

  PrintArray(const PrintArray &a)
    : domainwidth_m(a.domainwidth_m), datawidth_m(a.datawidth_m),
      dataprecision_m(a.dataprecision_m), carreturn_m(a.carreturn_m),
      scientific_m(a.scientific_m)
    {
    }


  //============================================================
  // Destructor.
  //============================================================

  ~PrintArray()
    {
    }


  //============================================================
  // PrintArray print methods.
  //============================================================

  //------------------------------------------------------------
  // print this array to the given stream.  This just invokes a
  // static method in an external class in order to get around a
  // Metrowerks bug.

  template<class S, class A, class DomainType>
  void print(S &s, const A &a, const DomainType &d) const
  {
    Pooma::blockAndEvaluate();
    PerformPrintArray<S,A,A::dimensions,DomainType>::print(*this, s, a, d);
  }

  //------------------------------------------------------------
  // Print the specified view of the array:
  template<class S, class A>
  void print(S &s, const A &a) const
  {
    Pooma::blockAndEvaluate();
    PerformPrintArray<S, A, A::dimensions, typename A::Domain_t>::
      print(*this, s, a, a.totalDomain());
  }


  //============================================================
  // PrintArray format settings accessors.
  //============================================================

  //------------------------------------------------------------
  // get/set the number of places used to print out domain numbers

  int domainWidth() const
    {
      return domainwidth_m;
    }

  void setDomainWidth(int val)
    {
      domainwidth_m = val;
      PAssert(domainwidth_m > 0);
    }

  //------------------------------------------------------------
  // get/set the number of places used to print out array data values

  int dataWidth() const
    {
      return datawidth_m;
    }

  void setDataWidth(int val)
    {
      datawidth_m = val;
      PAssert(datawidth_m > 0);
    }

  //------------------------------------------------------------
  // get/set the precision of the array data values

  int dataPrecision() const
    {
      return dataprecision_m;
    }

  void setDataPrecision(int val)
    {
      dataprecision_m = val;
      PAssert(dataprecision_m > 0);
    }

  //------------------------------------------------------------
  // get/set the number of carriage returns used

  int carReturn() const
    {
      return carreturn_m;
    }

  void setCarReturn(int val)
    {
      carreturn_m = val;
    }

  //------------------------------------------------------------
  // get/set the flag indicating whether to use scientific notation

  bool scientific() const
    {
      return scientific_m;
    }

  void setScientific(bool val)
    {
      scientific_m = val;
    }

  //------------------------------------------------------------
  // get/set the number of spaces between numbers

  int spacing() const
    {
      return spacing_m;
    }

  void setSpacing(int val)
    {
      spacing_m = val;
      PAssert(spacing_m >= 0);
    }

  // -----------------------------------------------------------
  // Set all the formatting parameters from the example PrintArray's values:
  void setFormatParameters(PrintArray &pa) {
    domainwidth_m = pa.domainWidth();
    datawidth_m = pa.dataWidth();
    dataprecision_m = pa.dataPrecision();
    carreturn_m = pa.carReturn();
    spacing_m = pa.spacing();
    scientific_m = pa.scientific();
  }

private:
  //------------------------------------------------------------
  // The width for domain numbers.

  int domainwidth_m;

  //------------------------------------------------------------
  // The width for array element values.

  int datawidth_m;

  //------------------------------------------------------------
  // The precision of the array element values.

  int dataprecision_m;

  //------------------------------------------------------------
  // How long before a carriage return is printed.

  int carreturn_m;

  //------------------------------------------------------------
  // The number of spaces between values.

  int spacing_m;

  //------------------------------------------------------------
  // Should scientific notation be used?

  bool scientific_m;
};



/// PerformPrintArray print method definition.
/// print takes data from an array, and prints it nicely to a stream.
/// S is a template parameter for an ostream-like object.
/// A is a template parameter for a array-like object.
/// This is the 1-D specialized case.

template<class S, class A, class DomainType>
void
PerformPrintArray<S,A,1,DomainType>::print(const PrintArray &p, S &s, 
                                           const A &a, const DomainType &d)
{
  // make sure this is the right function

  CTAssert(A::dimensions == 1);

  // determine the domain and domain-iterator type in the given array

  typedef          DomainType               Domain_t;
  typedef typename Domain_t::const_iterator Iterator_t;

  Domain_t domain(d);

  // create an iterator over the domain of the array

  Iterator_t griditer = domain.begin();
  Iterator_t  enditer = domain.end();

  // For a domain of total extent 1 (single Array element), treat specially:
  if (domain.size() == 1) {

    // Print the index value:
    s <<  "(";
    if (domain[0].first() < 0)
      s.fill(' ');
    else
      s.fill('0');
    s.width(p.domainWidth());
    s << domain[0].first();
    s << ")";
    s.fill(' ');
    s << " = ";

    // print the number
    if (p.scientific())
      s.setf(std::ios::scientific);
    s.precision(p.dataPrecision());
    s.width(p.dataWidth());
    s << a.read(*griditer);
    s << std::endl;

    return;
  }

  // For non-single-element-domain, proceed as always:

  // print out the prefix

  s << "(";
  if (domain[0].first() < 0)
    s.fill(' ');
  else
    s.fill('0');
  s.width(p.domainWidth());
  s << domain[0].first() << ":";
  if (domain[0].last() < 0)
    s.fill(' ');
  else
    s.fill('0');
  s.width(p.domainWidth());
  s << domain[0].last() << ":";
  if (domain[0].stride() < 0)
    s.fill(' ');
  else
    s.fill('0');
  s.width(p.domainWidth());
  s << domain[0].stride() << ") = ";
  s.fill(' ');

  // loop over the elements, printing out values as necessary

  int i, printed = 0;
  while (griditer != enditer)
    {
      // determine the number of spaces to print first
      int spacing = 0;
      if (printed > 0)
	{
	  spacing = p.spacing();
	  if (p.carReturn() >= 0 && printed >= p.carReturn())
	    {
	      s << std::endl;
	      spacing = 3*p.domainWidth() + 7;
	      printed = 0;
	    }
	}

      // print out spaces
      for (i=0; i < spacing; ++i)
	s << " ";

      // print the number
      if (p.scientific())
	s.setf(std::ios::scientific);
      s.precision(p.dataPrecision());
      s.width(p.dataWidth());
      s << a.read(*griditer);

      // increment iterator and counter
      ++griditer;
      ++printed;
    }

  // print final newline when done

  s << std::endl;
}


/// PerformPrintArray print method definition.
/// print takes data from an array, and prints it nicely to a stream.
/// S is a template parameter for an ostream-like object.
/// A is a template parameter for a array-like object.
/// This is the N-D general case, for N > 1.  It prints out 2D 'slices' for
/// the first two dimensions, and loops over the other dimensions.

template<class S, class A, int Dim, class DomainType>
void
PerformPrintArray<S,A,Dim,DomainType>::print(const PrintArray &p, S &s, 
                                             const A &a, const DomainType &d)
{
  int i, j, k;

  // make sure this is the right function

  CTAssert(A::dimensions == Dim && Dim > 1);

  // determine the domain and domain-iterator type in the given array

  typedef          DomainType               Domain_t;
  typedef typename Domain_t::Element_t      Element_t;
  typedef typename Domain_t::const_iterator Iterator_t;

  Domain_t domain(d);

  // create an iterator over the domain of the array

  Iterator_t griditer = domain.begin();
  Iterator_t  enditer = domain.end();

  // For a domain of total extent 1 (single Array element), treat specially:
  if (domain.size() == 1) {

    // Print the index values:
    s <<  "(";
    if (domain[0].first() < 0)
      s.fill(' ');
    else
      s.fill('0');
    s.width(p.domainWidth());
    s << domain[0].first();
    for (int d = 1; d < Dim; d++) {
      s << ",";
      if (domain[d].first() < 0)
        s.fill(' ');
      else
        s.fill('0');
      s.width(p.domainWidth());
      s << domain[d].first();
    }
    s << ")";
    s.fill(' ');
    s << " = ";

    // print the number
    if (p.scientific())
      s.setf(std::ios::scientific);
    s.precision(p.dataPrecision());
    s.width(p.dataWidth());
    s << a.read(*griditer);
    s << std::endl;

    return;
  }

  // For non-single-element-domain, proceed as always:

  // get 1-D domains info for the first two dimensions, and use these sizes
  // to determine how to do the inner dimensional loops

  Element_t x0 = domain[0].first();
  Element_t x1 = domain[0].last();
  Element_t xs = domain[0].stride();
  Element_t y0 = domain[1].first();
  Element_t y1 = domain[1].last();
  Element_t ys = domain[1].stride();

  // Start looping over all the elements.  We can stop printing when
  // we hit the end of the griditer.  We print out 2D slices, and for
  // higher dimensions, we print out a line saying which slice is coming
  // up next, for example:
  //   [1:5:1][2:8:1][2][4][0]:
  //   -------------------------------------------

  // For 3D and higher, also print out the entire domain specification:
  if (Dim > 2) {
    s << std::endl << "~~~~~~~~~~~~~~ ";
    s << "(" << domain[0].first() << ":" << domain[0].last()
      << ":" << domain[0].stride();
    for (int d = 1; d < Dim; d++) {
      s << "," << domain[d].first() << ":" << domain[d].last() << ":" 
        << domain[d].stride();
    }
    s << ")" << " ~~~~~~~~~~~~~~" << std::endl;
  }

  while (griditer != enditer)
    {
      // print out the higher-dim size statement, if necessary
      if (Dim > 2) {
        s << std::endl << "(" << domain[0].first() << ":" << domain[0].last()
          << ":" << domain[0].stride() << "," << domain[1].first() << ":"
          << domain[1].last() << ":" << domain[1].stride();
        for (i=2; i < Dim; ++i)
          s << "," << (*griditer)[i].first();
        s << "):" << std::endl;
        s << "----------------------------------------------------" 
          << std::endl;
      }

      // loop over all the elements of the next 2D slice now
      for (j=y0; j <= y1; j += ys)
	{
	  // print out the prefix for the next 2D slice line
	  s <<  "(";
          if (x0 < 0)
            s.fill(' ');
          else
            s.fill('0');
          s.width(p.domainWidth());
          s << x0 << ":";
          if (x1 < 0)
            s.fill(' ');
          else
            s.fill('0');
          s.width(p.domainWidth());
	  s << x1 << ":";
          if (xs < 0)
            s.fill(' ');
          else
            s.fill('0');
          s.width(p.domainWidth());
	  s << xs;
	  for (i=1; i < Dim; ++i)
	    {
              s.fill(' ');
              s << ",";
              if ((*griditer)[i].first() < 0)
                s.fill(' ');
              else
                s.fill('0');
              s.width(p.domainWidth());
	      s << (*griditer)[i].first();
	    }
          s << ")";
	  s.fill(' ');
          s << " = ";

	  // print all the values along this 1-D strip
	  int printed = 0;
	  for (i=x0; i <= x1; i += xs)
	    {
	      // determine the number of spaces to print first
	      int spacing = 0;
	      if (printed > 0)
		{
		  spacing = p.spacing();
		  if (p.carReturn() >= 0 && printed >= p.carReturn())
		    {
		      s << std::endl;
// 		      spacing = (Dim + 1)*(p.domainWidth() + 2) + 4;
		      spacing = (Dim + 2)*p.domainWidth() + 2 + Dim + 4;
		      printed = 0;
		    }
		}

	      // print out spaces
	      for (k=0; k < spacing; ++k)
		s << ' ';

	      // print the number
	      if (p.scientific())
		s.setf(std::ios::scientific);
	      s.precision(p.dataPrecision());
	      s.width(p.dataWidth());

	      typedef typename A::Engine_t::Tag_t Etag_t;

	      Etag_t tag;
 	      if(isValidLocation(a,*griditer,tag)) 
		s << a.read(*griditer);
 	      else
 		s<< ".";

	      // increment iterator and counter
	      ++griditer;
	      ++printed;
	    }

	  // print final newline when done with this line
	  s << std::endl;
	}
    }
}

// } // namespace Pooma

#endif // POOMA_ARRAY_PRINT_ARRAY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PrintArray.h,v $   $Author: richard $
// $Revision: 1.23 $   $Date: 2004/11/01 18:16:13 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
