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
//   PrintField - information and routines used to print out a
//                Field to a stream.
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_PRINTFIELD_H
#define POOMA_FIELD_PRINTFIELD_H

/** @file
 * @ingroup Field
 * @brief
 * Extract the elements of a Field and print out the
 * contents to a stream with some nice formatting.
 *
 * The user can select
 * how many elements to print per line, the precision, format, etc.
 * This is done by constructing an PrintField, and calling the
 * 'print' method with the stream to print to and the field to print.
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"         // Pooma::blockAndEvaluate().
#include "Utilities/PAssert.h"
#include <iomanip>


//-----------------------------------------------------------------------------
// Forward References:
//-----------------------------------------------------------------------------

class PrintField;



//-----------------------------------------------------------------------------
// Open POOMA namespace:
//-----------------------------------------------------------------------------

// namespace POOMA {

/**
 * PerformPrintField struct: a simple wrapper around the templated 'print'
 * method of PrintField.
 * This struct defines one static method 'print', which
 * prints an field to a stream (the two template parameters).  The Dim
 * parameter is used to specialize this to the case of a 1D domain.  This is
 * called by the 'print' templated member function of PrintField.
 */

template<int Dim>
struct PerformPrintField
{
  template <class S, class A>
  static void print(const PrintField &, S &, const A &);
};

template<>
struct PerformPrintField<1>
{
  template<class S, class A>
  static void print(const PrintField &, S &, const A &);
};


//-----------------------------------------------------------------------------
// PrintField class declaration
//-----------------------------------------------------------------------------

/**
 * PrintField is a utility program for ConstField and Field object, used
 * to print (nicely) the contents of an field to a provided stream.  All
 * the data in the field will be printed in ascii format; if you want to
 * display just a portion of an field, take a view of the field and give
 * that to the print() method.
 *
 * CONSTRUCTING A PrintField:
 * --------------------------
 * When you construct a PrintField, you can give it several format parameters
 * to control how to display the field.  These parameters are, in the
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
 *                  the field, no carriage return is printed until the end
 *                  of the row.  If this number is > 0, it represents the
 *                  maximum number of values that will be printed before a
 *                  return.  Only the first set of numbers for a row of
 *                  the field will have a domain prefix included.
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
 * PRINTING FIELD OBJECTS WITH A PrintField:
 * -----------------------------------------
 * PrintField is not templated, so that you can reuse the same formatter
 * for different fields.  It has one templated member function 'print':
 *
 *   template<class S, class A>
 *   void print(S &s, const A &a) const
 *
 * where 'S' must be an object with an ostream-like interface (such as
 * cout, or an Inform object), and 'A' must be an object with a Field
 * interface.  'print' will take data from A and print it to the stream
 * using the current format settings.
 *
 * 1-D fields just have the one row printed, perhaps on multiple lines
 * if carReturn is > 0.  2-D fields are printed as a table, with each
 * line prefixed by the domain it includes.  For example:
 * <PRE>
 *   [00:02][00] =          0          0          0
 *   [00:02][01] =          0          0          0
 * </PRE>
 * prints the values for [0:2][0:1] of a 2D field.  3-D and higher arrays
 * have a sequence of 2-D slices printed for them, each slice separated
 * by a line indicating which slice it is, and a separator.  Example:
 * <PRE>
 *   [0:2:1][0:3:1][2]:
 *   ----------------------------------------------------
 *   [00:02][00][02] =          0          0          0
 *   [00:02][01][02] =          0          0          0
 *   [00:02][02][02] =          0          0          0
 *   [00:02][03][02] =          0          0          0
 * </PRE>
 */ 

class PrintField
{
public:
  //============================================================
  // Constructors.
  //============================================================

  //------------------------------------------------------------
  // Construct an PrintField object with an field, and the settings for
  // how the field should be printed.

  PrintField(int domainWidth = 3, int dataWidth = 10,
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

  PrintField(const PrintField &a)
    : domainwidth_m(a.domainwidth_m), datawidth_m(a.datawidth_m),
      dataprecision_m(a.dataprecision_m), carreturn_m(a.carreturn_m),
      scientific_m(a.scientific_m)
    {
    }


  //============================================================
  // Destructor.
  //============================================================

  ~PrintField()
    {
    }


  //============================================================
  // PrintField print methods.
  //============================================================

  //------------------------------------------------------------
  // print this field to the given stream.  This just invokes a
  // static method in an external class in order to get around a
  // Metrowerks bug.

  template<class S, class A>
  void print(S &s, const A &a) const
    {
      forEach(a, PerformUpdateTag(), NullCombine());
      Pooma::blockAndEvaluate();

      for (int m = 0; m < a.numMaterials(); m++)
        for (int c = 0; c < a.centeringSize(); c++)
          {
            s << "Material #" << m << ", Centering #" << c << " " << a.centering(c) 
              << "\n"<< "-------------\n";
            PerformPrintField<A::dimensions>::print(*this, s, a.subField(m, c));
          }
    }    


  //============================================================
  // PrintField format settings accessors.
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
  // get/set the number of places used to print out field data values

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
  // get/set the precision of the field data values

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

private:
  //------------------------------------------------------------
  // The width for domain numbers.

  int domainwidth_m;

  //------------------------------------------------------------
  // The width for field element values.

  int datawidth_m;

  //------------------------------------------------------------
  // The precision of the field element values.

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


//-----------------------------------------------------------------------------
// PerformPrintField print method definition.
// print takes data from an field, and prints it nicely to a stream.
// S is a template parameter for an ostream-like object.
// A is a template parameter for a Field-like object.
// This is the 1-D specialized case.
//-----------------------------------------------------------------------------

template<class S, class A>
void
PerformPrintField<1>::print(const PrintField &p, S &s, const A &a)
{
  // make sure this is the right function

  CTAssert(A::dimensions == 1);

  // determine the domain and domain-iterator type in the given field

  typedef typename A::Domain_t              Domain_t;
  typedef typename Domain_t::const_iterator Iterator_t;

  // create an iterator over the domain of the field

  Iterator_t griditer = a.domain().begin();
  Iterator_t  enditer = a.domain().end();

  // print out the prefix

  s << "[";
  if (a.domain()[0].first() < 0)
    s.fill(' ');
  else
    s.fill('0');
  s.width(p.domainWidth());
  s << a.domain()[0].first() << ":";
  if (a.domain()[0].last() < 0)
    s.fill(' ');
  else
    s.fill('0');
  s.width(p.domainWidth());
  s << a.domain()[0].last() << "] = ";
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
              s << "\n";
	          spacing = 2*p.domainWidth() + 6;
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

  s << "\n";
}


//-----------------------------------------------------------------------------
// PerformPrintField print method definition.
// print takes data from an field, and prints it nicely to a stream.
// S is a template parameter for an ostream-like object.
// A is a template parameter for a Field-like object.
// This is the N-D general case, for N > 1.  It prints out 2D 'slices' for
// the first two dimensions, and loops over the other dimensions.
//-----------------------------------------------------------------------------

template<int Dim>
template<class S, class A>
void
PerformPrintField<Dim>::print(const PrintField &p, S &s, const A &a)
{
  int i, j, k;

  // make sure this is the right function

  CTAssert(A::dimensions == Dim && Dim > 1);

  // determine the domain and domain-iterator type in the given field

  typedef typename A::Domain_t              Domain_t;
  typedef typename Domain_t::Element_t      Element_t;
  typedef typename Domain_t::const_iterator Iterator_t;

  // create an iterator over the domain of the field

  Iterator_t griditer = a.domain().begin();
  Iterator_t  enditer = a.domain().end();

  // get 1-D domains info for the first two dimensions, and use these sizes
  // to determine how to do the inner dimensional loops

  Element_t x0 = a.domain()[0].first();
  Element_t x1 = a.domain()[0].last();
  Element_t xs = a.domain()[0].stride();
  Element_t y0 = a.domain()[1].first();
  Element_t y1 = a.domain()[1].last();
  Element_t ys = a.domain()[1].stride();

  // Start looping over all the elements.  We can stop printing when
  // we hit the end of the griditer.  We print out 2D slices, and for
  // higher dimensions, we print out a line saying which slice is coming
  // up next, for example:
  //   [1:5:1][2:8:1][2][4][0]:
  //   -------------------------------------------

  while (griditer != enditer)
    {
      // print out the higher-dim size statement, if necessary
      if (Dim > 2)
	    {
	      s << '\n' << a.domain()[0] << a.domain()[1];
	      for (i=2; i < Dim; ++i)
	        s << "[" << (*griditer)[i].first() << "]";
	      s << ":" << '\n';
	      s << "----------------------------------------------------\n";
        }

      // loop over all the elements of the next 2D slice now
      for (j=y0; j <= y1; j += ys)
	    {
	      // print out the prefix for the next 2D slice line
	      s <<  "[";
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
	      s << x1 << "]";
	      for (i=1; i < Dim; ++i)
	        {
	          s << "[";
              if ((*griditer)[i].first() < 0)
                s.fill(' ');
              else
                s.fill('0');
              s.width(p.domainWidth());
	          s << (*griditer)[i].first() << "]";
	        }
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
		              s << '\n';
		              spacing = (Dim + 1)*(p.domainWidth() + 2) + 4;
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
	          s << a.read(*griditer);

	          // increment iterator and counter
	          ++griditer;
	          ++printed;
	        }

	      // print final newline when done with this line
	      s << '\n';
        }
    }
}

// } // namespace Pooma

#endif // POOMA_FIELD_PRINTFIELD_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PrintField.h,v $   $Author: richard $
// $Revision: 1.7 $   $Date: 2004/11/01 18:16:43 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
