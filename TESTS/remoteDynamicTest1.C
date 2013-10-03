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
// Remote Dynamic Engine test code
//-----------------------------------------------------------------------------

#include "Pooma/Pooma.h"
#include "Utilities/Tester.h"
#include "Domain/Loc.h"
#include "Domain/Interval.h"
#include "Domain/IndirectionList.h"
#include "Engine/DynamicEngine.h"
#include "Engine/RemoteDynamicEngine.h"
#include "Array/Array.h"

#include <iostream>
#include <vector>

#if POOMA_MESSAGING

struct PackObject
{
  typedef Array<1, double, Remote<Dynamic> > Array_t;
  typedef IndirectionList<int> List_t;
  PackObject(Array_t &array, List_t &list)
    : array_m(&array), list_m(&list)
  {
  }

  PackObject() { }

  Array_t *array_m;
  List_t  *list_m;
  char    *buffer_m;
};

namespace Cheetah {

template<>
class Serialize<CHEETAH, PackObject>
{
public:
  static inline int
  size(const PackObject &pack)
  {
    int dummy;
    return Serialize<CHEETAH, int>::size(dummy) +
      pack.array_m->engine().packSize(*(pack.list_m));
  }

  static inline int
  pack(const PackObject &pack, char *buffer)
  {
    int length;
    char *data = buffer + Serialize<CHEETAH, int>::size(length);
    length = pack.array_m->engine().pack(*(pack.list_m), data, false);
    return (length + Serialize<CHEETAH, int>::pack(length, buffer));
  }

  static inline int
  unpack(PackObject* &pack, char *buffer)
  {
    pack = new PackObject();
    int *length;
    int sizeInt = Serialize<CHEETAH, int>::unpack(length, buffer);
    pack->buffer_m = buffer + sizeInt;

    return sizeInt + *length;
  }

  static inline void
  cleanup(PackObject *pack)
  {
    delete pack;
  }
};

} // namespace Cheetah

static bool ready_g = false;

void unpackFunction(Array<1, double, Remote<Dynamic> > *b, PackObject &pack)
{
  b->engine().unpack(Interval<1>(7, 11), pack.buffer_m, false);
  ready_g = true;
}

#endif


int main(int argc, char *argv[])
{
  Pooma::initialize(argc,argv);
  Pooma::Tester tester(argc,argv);
  tester.out().setOutputContext(-1);

  int myContext = Pooma::context();
  int numContexts = Pooma::contexts();

  // Create the total domain.

  Interval<1> x(12);
  
  Array<1, double, Remote<Dynamic> > a(x);

  // Store some stuff.
  
  int i;
  for (i = 0; i < 12; ++i)
  {
     a(i) = double(i);
  }

  tester.out() << " Array a = " << a << std::endl;

#if POOMA_CHEETAH

  IndirectionList<int> list(5);

  list(0) = 3;
  list(1) = 5;
  list(2) = 7;
  list(3) = 8;
  list(4) = 9;

  if (myContext == 0)
  {
    int size = a.engine().packSize(list);

    char *buffer = new char [size];

    a.engine().pack(list, buffer, false);

    a.engine().unpack(Interval<1>(7, 11), buffer, false);

    delete [] buffer;
  }

  tester.out() << " Array a = " << a << std::endl;

  if (numContexts > 1)
  {
    // Create an engine on another context...
    Engine<1, double, Remote<Dynamic> > eng(1, x);
    Array<1, double, Remote<Dynamic> > b(eng);

    b = 0;
    Pooma::blockAndEvaluate();

    if (myContext == 0)
    {
      PackObject pack(a, list);
      int toContext = 1;
      int tag = Pooma::sendTag(toContext);
      tester.out() << "Sending data to context " << toContext
                   << " with tag " << tag << std::endl;
      Pooma::particleSwapHandler()->send(toContext, tag, pack);
    }
    if (myContext == 1)
    {
      int fromContext = 0;
      int tag = Pooma::receiveTag(fromContext);
      tester.out() << "Receiving data from context " << fromContext
                   << " with tag " << tag << std::endl;
      Pooma::particleSwapHandler()->request(fromContext, tag,
					    unpackFunction, &b);

      while (!ready_g)
      {
	Pooma::poll();
      }
    }

    tester.out() << b << std::endl;
  }

#endif

  int ret = tester.results("remoteDynamicTest1");
  Pooma::finalize();
  return 0;
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: remoteDynamicTest1.cpp,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:38 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
