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

#ifndef POOMA_CONNECT_PAWS_RESIZE_H
#define POOMA_CONNECT_PAWS_RESIZE_H

//-----------------------------------------------------------------------------
// Classes:
// Resize<T>
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Overview:
//
// Resize<T> is a simple class with one static method 'resize' that is
// specialized on the type T to do an array resize operation for various
// engine types.  Some cannot be resized and will throw a runtime error,
// others will get their sizes changed (which may or may not preserve the
// existing contents).
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Pooma/Arrays.h"
#include "Utilities/PAssert.h"


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

//-----------------------------------------------------------------------------
//
// Resize<T> : general version just throws a runtime assertion failure
//
//-----------------------------------------------------------------------------

template<class T>
struct Resize
{
  template<class Dom>
  static void resize(const T &, const Dom &)
  {
    PInsist(false, "Resize<T>::resize(): Cannot resize the given type.");
  }
};


//-----------------------------------------------------------------------------
//
// Resize<Array<T>> : For Array, by default try an initialize
//
//-----------------------------------------------------------------------------

template<int Dim, class T, class E>
struct Resize< Array<Dim, T, E> >
{
  template<class Dom>
  static void resize(Array<Dim, T, E> &array, const Dom &domain)
  {
    PInsist(false,
	    "Resize<Array<D,T,E>::resize(): Cannot resize the given type.");
  }
};


//-----------------------------------------------------------------------------
//
// Resize< Array<1,T,SharedBrick> : For special 1D engines,
// resize each patch separately to get an equal size in all patches.
//
//-----------------------------------------------------------------------------

template<class T>
struct Resize< Array<1, T, SharedBrick> >
{
  template<class Dom>
  static void resize(Array<1, T, SharedBrick> &array,
		     const Dom &domain)
  {
    // This will only work with 1D domains

    CTAssert(DomainTraits<Dom>::dimensions == 1);

    // Get the number of elements from the domain, and do resize

    int patchsize = DomainTraits<Dom>::getSize(domain);
    int currsize = array.domain().size();
    if (currsize > patchsize)
      {
	array.engine().destroy(Interval<1>(patchsize, currsize - 1), 0,
			       ShiftUp());
      }
    else if (currsize < patchsize)
      {
	array.engine().create(patchsize - currsize, 0);
      }

    // At the end, do a sync to get the layout updated.

    array.engine().sync();
  }
};


//-----------------------------------------------------------------------------
//
// Resize< Array<1,T,MultiBrick<Grid,Brick> > : For special 1D engines,
// resize each patch separately to get an equal size in all patches.
//
//-----------------------------------------------------------------------------

template<class T>
struct Resize< Array<1, T, MultiPatch<GridTag, Brick> > >
{
  template<class Dom>
  static void resize(Array<1, T, MultiPatch<GridTag, Brick> > &array,
		     const Dom &domain)
  {
    // This will only work with 1D domains

    CTAssert(DomainTraits<Dom>::dimensions == 1);

    // We must have some patches in the destination array

    int patches = array.numPatchesLocal();
    PAssert(patches > 0);

    // Get the number of elements from the domain

    int newsize = DomainTraits<Dom>::getSize(domain);

    // Compute the number of elements to put in each patch

    int basicpatchsize = newsize / patches;
    int extra = newsize % patches;

    // Resize each patch ... this could be made parallel since each resize
    // is independent.

    for (int p=0; p < patches; ++p)
      {
	int patchsize = basicpatchsize + (p < extra ? 1 : 0);
	int currsize = array.patch(p).domain().size();
	if (currsize > patchsize)
	  {
	    array.engine().destroy(Interval<1>(patchsize, currsize - 1), p,
				   ShiftUp());
	  }
	else if (currsize < patchsize)
	  {
	    array.engine().create(patchsize - currsize, p);
	  }
      }

    // At the end, do a sync to get the layout updated.

    array.engine().sync();
  }
};


// }; // namespace Pooma
///////////////////////////////////////////////////////////////////////////////


#endif // POOMA_CONNECT_PAWS_RESIZE_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Resize.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:19 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
