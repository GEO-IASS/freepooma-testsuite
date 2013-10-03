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
// MultiArg1<A1>
// MultiArg2<A1,A2>
// etc.
// Functions:
// applyMultiArg
//-----------------------------------------------------------------------------

#ifndef POOMA_FUNCTIONS_MULTIARG_H
#define POOMA_FUNCTIONS_MULTIARG_H

//////////////////////////////////////////////////////////////////////

/** @file
 * @ingroup Utilities
 * @brief
 * The MultiArg1...N classes are intended to be used to wrap multiple arrays,
 * fields, or particles where a common set of operations need to be performed
 * on the set.
 *
 *  Typical operations would be: taking views, getting locks,
 * performing intersections.
 *
 * It would be nicer in some sense to have an inhomogenous container with
 * common interfaces through a base class.  That approach would avoid the
 * following huge number of MultiArg classes.
 * The difficulty arises in extracting
 * the arrays at the end and handing them off to a function.  To extract the
 * arrays requires knowing their types, requiring you to know the complete
 * signature of the function you are calling.  At that point, you've already
 * accumulated all the type information in something at least as complicated
 * as the MultiArg<> classes here.
 */

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Pooma/View.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

template<class A1> struct MultiArg1;
template<class A1, class A2> struct MultiArg2;
template<class A1, class A2, class A3> struct MultiArg3;
template<class A1, class A2, class A3, class A4> struct MultiArg4;
template<class A1, class A2, class A3, class A4, class A5> struct MultiArg5;
template<class A1, class A2, class A3, class A4, class A5, class A6> struct MultiArg6;
template<class A1, class A2, class A3, class A4, class A5, class A6, class A7> struct MultiArg7;

// These MultiArgView structs are workarounds for Compaq's
// c++ compiler.  It had some odd problems with the definition
// of View1<MA<A>> in terms of View1<A>.

template<class A1, class Dom>
struct MultiArgView1
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef MultiArg1<A1_t> Type_t;
};

template<class A1, class A2, class Dom>
struct MultiArgView2
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef typename View1<A2, Dom>::Type_t A2_t;
  typedef MultiArg2<A1_t, A2_t> Type_t;
};

template<class A1, class A2, class A3, class Dom>
struct MultiArgView3
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef typename View1<A2, Dom>::Type_t A2_t;
  typedef typename View1<A3, Dom>::Type_t A3_t;
  typedef MultiArg3<A1_t, A2_t, A3_t> Type_t;
};

template<class A1, class A2, class A3, class A4, class Dom>
struct MultiArgView4
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef typename View1<A2, Dom>::Type_t A2_t;
  typedef typename View1<A3, Dom>::Type_t A3_t;
  typedef typename View1<A4, Dom>::Type_t A4_t;
  typedef MultiArg4<A1_t, A2_t, A3_t, A4_t> Type_t;
};

template<class A1, class A2, class A3, class A4, class A5, class Dom>
struct MultiArgView5
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef typename View1<A2, Dom>::Type_t A2_t;
  typedef typename View1<A3, Dom>::Type_t A3_t;
  typedef typename View1<A4, Dom>::Type_t A4_t;
  typedef typename View1<A5, Dom>::Type_t A5_t;
  typedef MultiArg5<A1_t, A2_t, A3_t, A4_t, A5_t> Type_t;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class Dom>
struct MultiArgView6
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef typename View1<A2, Dom>::Type_t A2_t;
  typedef typename View1<A3, Dom>::Type_t A3_t;
  typedef typename View1<A4, Dom>::Type_t A4_t;
  typedef typename View1<A5, Dom>::Type_t A5_t;
  typedef typename View1<A6, Dom>::Type_t A6_t;
  typedef MultiArg6<A1_t, A2_t, A3_t, A4_t, A5_t, A6_t> Type_t;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class Dom>
struct MultiArgView7
{
  typedef typename View1<A1, Dom>::Type_t A1_t;
  typedef typename View1<A2, Dom>::Type_t A2_t;
  typedef typename View1<A3, Dom>::Type_t A3_t;
  typedef typename View1<A4, Dom>::Type_t A4_t;
  typedef typename View1<A5, Dom>::Type_t A5_t;
  typedef typename View1<A6, Dom>::Type_t A6_t;
  typedef typename View1<A7, Dom>::Type_t A7_t;
  typedef MultiArg7<A1_t, A2_t, A3_t, A4_t, A5_t, A6_t, A7_t> Type_t;
};

template<class A1, class Dom>
struct View1<MultiArg1<A1>, Dom>
{
  typedef typename MultiArgView1<A1, Dom>::Type_t Type_t;
};

template<class A1, class A2, class Dom>
struct View1<MultiArg2<A1, A2>, Dom>
{
  typedef typename MultiArgView2<A1, A2, Dom>::Type_t Type_t;
};

template<class A1, class A2, class A3, class Dom>
struct View1<MultiArg3<A1, A2, A3>, Dom>
{
  typedef typename MultiArgView3<A1, A2, A3, Dom>::Type_t Type_t;
};

template<class A1, class A2, class A3, class A4, class Dom>
struct View1<MultiArg4<A1, A2, A3, A4>, Dom>
{
  typedef typename MultiArgView4<A1, A2, A3, A4, Dom>::Type_t Type_t;
};

template<class A1, class A2, class A3, class A4, class A5, class Dom>
struct View1<MultiArg5<A1, A2, A3, A4, A5>, Dom>
{
  typedef typename MultiArgView5<A1, A2, A3, A4, A5, Dom>::Type_t Type_t;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class Dom>
struct View1<MultiArg6<A1, A2, A3, A4, A5, A6>, Dom>
{
  typedef typename MultiArgView6<A1, A2, A3, A4, A5, A6, Dom>::Type_t Type_t;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class Dom>
struct View1<MultiArg7<A1, A2, A3, A4, A5, A6, A7>, Dom>
{
  typedef typename MultiArgView7<A1, A2, A3, A4, A5, A6, A7, Dom>::Type_t Type_t;
};

template<class A1>
struct MultiArg1
{
  enum { size = 1 };

  MultiArg1(const A1 &a1)
    : a1_m(a1)
  {
  }

  template<class Dom>
  typename View1<MultiArg1<A1>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg1<A1>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom));
  }

  A1 a1_m;
};

template<class A1, class A2>
struct MultiArg2
{
  enum { size = 2 };

  MultiArg2(const A1 &a1, const A2 &a2)
    : a1_m(a1), a2_m(a2)
  {
  }

  template<class Dom>
  typename View1<MultiArg2<A1, A2>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg2<A1, A2>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom), a2_m(dom));
  }

  A1 a1_m;
  A2 a2_m;
};

template<class A1, class A2, class A3>
struct MultiArg3
{
  enum { size = 3 };

  MultiArg3(const A1 &a1, const A2 &a2, const A3 &a3)
    : a1_m(a1), a2_m(a2), a3_m(a3)
  {
  }

  template<class Dom>
  typename View1<MultiArg3<A1, A2, A3>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg3<A1, A2, A3>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom), a2_m(dom), a3_m(dom));
  }

  A1 a1_m;
  A2 a2_m;
  A3 a3_m;
};

template<class A1, class A2, class A3, class A4>
struct MultiArg4
{
  enum { size = 4 };

  MultiArg4(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4)
    : a1_m(a1), a2_m(a2), a3_m(a3), a4_m(a4)
  {
  }

  template<class Dom>
  typename View1<MultiArg4<A1, A2, A3, A4>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg4<A1, A2, A3, A4>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom), a2_m(dom), a3_m(dom), a4_m(dom));
  }

  A1 a1_m;
  A2 a2_m;
  A3 a3_m;
  A4 a4_m;
};

template<class A1, class A2, class A3, class A4, class A5>
struct MultiArg5
{
  enum { size = 5 };

  MultiArg5(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5)
    : a1_m(a1), a2_m(a2), a3_m(a3), a4_m(a4), a5_m(a5)
  {
  }

  template<class Dom>
  typename View1<MultiArg5<A1, A2, A3, A4, A5>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg5<A1, A2, A3, A4, A5>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom), a2_m(dom), a3_m(dom), a4_m(dom), a5_m(dom));
  }

  A1 a1_m;
  A2 a2_m;
  A3 a3_m;
  A4 a4_m;
  A5 a5_m;
};

template<class A1, class A2, class A3, class A4, class A5, class A6>
struct MultiArg6
{
  enum { size = 6 };

  MultiArg6(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6)
    : a1_m(a1), a2_m(a2), a3_m(a3), a4_m(a4), a5_m(a5), a6_m(a6)
  {
  }

  template<class Dom>
  typename View1<MultiArg6<A1, A2, A3, A4, A5, A6>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg6<A1, A2, A3, A4, A5, A6>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom), a2_m(dom), a3_m(dom), a4_m(dom), a5_m(dom), a6_m(dom));
  }

  A1 a1_m;
  A2 a2_m;
  A3 a3_m;
  A4 a4_m;
  A5 a5_m;
  A6 a6_m;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7>
struct MultiArg7
{
  enum { size = 7 };

  MultiArg7(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7)
    : a1_m(a1), a2_m(a2), a3_m(a3), a4_m(a4), a5_m(a5), a6_m(a6), a7_m(a7)
  {
  }

  template<class Dom>
  typename View1<MultiArg7<A1, A2, A3, A4, A5, A6, A7>, Dom>::Type_t
  operator()(const Dom &dom) const
  {
    typedef typename View1<MultiArg7<A1, A2, A3, A4, A5, A6, A7>, Dom>::Type_t Ret_t;
    return Ret_t(a1_m(dom), a2_m(dom), a3_m(dom), a4_m(dom), a5_m(dom), a6_m(dom), a7_m(dom));
  }

  A1 a1_m;
  A2 a2_m;
  A3 a3_m;
  A4 a4_m;
  A5 a5_m;
  A6 a6_m;
  A7 a7_m;
};

template<class A1, class Function>
void applyMultiArg(const MultiArg1<A1> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
}

template<class A1, class A2, class Function>
void applyMultiArg(const MultiArg2<A1, A2> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
  f(node.a2_m, condition[1]);
}

template<class A1, class A2, class A3, class Function>
void applyMultiArg(const MultiArg3<A1, A2, A3> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
  f(node.a2_m, condition[1]);
  f(node.a3_m, condition[2]);
}

template<class A1, class A2, class A3, class A4, class Function>
void applyMultiArg(const MultiArg4<A1, A2, A3, A4> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
  f(node.a2_m, condition[1]);
  f(node.a3_m, condition[2]);
  f(node.a4_m, condition[3]);
}

template<class A1, class A2, class A3, class A4, class A5, class Function>
void applyMultiArg(const MultiArg5<A1, A2, A3, A4, A5> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
  f(node.a2_m, condition[1]);
  f(node.a3_m, condition[2]);
  f(node.a4_m, condition[3]);
  f(node.a5_m, condition[4]);
}

template<class A1, class A2, class A3, class A4, class A5, class A6, class Function>
void applyMultiArg(const MultiArg6<A1, A2, A3, A4, A5, A6> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
  f(node.a2_m, condition[1]);
  f(node.a3_m, condition[2]);
  f(node.a4_m, condition[3]);
  f(node.a5_m, condition[4]);
  f(node.a6_m, condition[5]);
}

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class Function>
void applyMultiArg(const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &node,
		   const Function &f,
		   const std::vector<bool> &condition)
{
  f(node.a1_m, condition[0]);
  f(node.a2_m, condition[1]);
  f(node.a3_m, condition[2]);
  f(node.a4_m, condition[3]);
  f(node.a5_m, condition[4]);
  f(node.a6_m, condition[5]);
  f(node.a7_m, condition[6]);
}

template<class A1, class Function>
void applyMultiArg(const MultiArg1<A1> &node,
		   const Function &f)
{
  f(node.a1_m);
}

template<class A1, class A2, class Function>
void applyMultiArg(const MultiArg2<A1, A2> &node,
		   const Function &f)
{
  f(node.a1_m);
  f(node.a2_m);
}

template<class A1, class A2, class A3, class Function>
void applyMultiArg(const MultiArg3<A1, A2, A3> &node,
		   const Function &f)
{
  f(node.a1_m);
  f(node.a2_m);
  f(node.a3_m);
}

template<class A1, class A2, class A3, class A4, class Function>
void applyMultiArg(const MultiArg4<A1, A2, A3, A4> &node,
		   const Function &f)
{
  f(node.a1_m);
  f(node.a2_m);
  f(node.a3_m);
  f(node.a4_m);
}

template<class A1, class A2, class A3, class A4, class A5, class Function>
void applyMultiArg(const MultiArg5<A1, A2, A3, A4, A5> &node,
		   const Function &f)
{
  f(node.a1_m);
  f(node.a2_m);
  f(node.a3_m);
  f(node.a4_m);
  f(node.a5_m);
}

template<class A1, class A2, class A3, class A4, class A5, class A6, class Function>
void applyMultiArg(const MultiArg6<A1, A2, A3, A4, A5, A6> &node,
		   const Function &f)
{
  f(node.a1_m);
  f(node.a2_m);
  f(node.a3_m);
  f(node.a4_m);
  f(node.a5_m);
  f(node.a6_m);
}

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class Function>
void applyMultiArg(const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &node,
		   const Function &f)
{
  f(node.a1_m);
  f(node.a2_m);
  f(node.a3_m);
  f(node.a4_m);
  f(node.a5_m);
  f(node.a6_m);
  f(node.a7_m);
}

template<class A1, class Function>
void applyMultiArgIf(const MultiArg1<A1> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);
}

template<class A1, class A2, class Function>
void applyMultiArgIf(const MultiArg2<A1, A2> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);

  if (condition[1])
    f(node.a2_m);
}

template<class A1, class A2, class A3, class Function>
void applyMultiArgIf(const MultiArg3<A1, A2, A3> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);

  if (condition[1])
    f(node.a2_m);

  if (condition[2])
    f(node.a3_m);
}

template<class A1, class A2, class A3, class A4, class Function>
void applyMultiArgIf(const MultiArg4<A1, A2, A3, A4> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);

  if (condition[1])
    f(node.a2_m);

  if (condition[2])
    f(node.a3_m);

  if (condition[3])
    f(node.a4_m);
}

template<class A1, class A2, class A3, class A4, class A5, class Function>
void applyMultiArgIf(const MultiArg5<A1, A2, A3, A4, A5> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);

  if (condition[1])
    f(node.a2_m);

  if (condition[2])
    f(node.a3_m);

  if (condition[3])
    f(node.a4_m);

  if (condition[4])
    f(node.a5_m);
}

template<class A1, class A2, class A3, class A4, class A5, class A6, class Function>
void applyMultiArgIf(const MultiArg6<A1, A2, A3, A4, A5, A6> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);

  if (condition[1])
    f(node.a2_m);

  if (condition[2])
    f(node.a3_m);

  if (condition[3])
    f(node.a4_m);

  if (condition[4])
    f(node.a5_m);

  if (condition[5])
    f(node.a6_m);
}


template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class Function>
void applyMultiArgIf(const MultiArg7<A1, A2, A3, A4, A5, A6, A7> &node,
		 const Function &f,
		 const std::vector<bool> &condition)
{
  if (condition[0])
    f(node.a1_m);

  if (condition[1])
    f(node.a2_m);

  if (condition[2])
    f(node.a3_m);

  if (condition[3])
    f(node.a4_m);

  if (condition[4])
    f(node.a5_m);

  if (condition[5])
    f(node.a6_m);

  if (condition[6])
    f(node.a7_m);
}

#endif     // POOMA_FUNCTIONS_MULTIARG_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: MultiArg.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:49 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
