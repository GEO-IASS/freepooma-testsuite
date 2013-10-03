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
// BisectionMapper
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Partition
 * @brief
 * BisectionMapper is a ContextMapper implementation.
 */

#ifndef POOMA_BISECTIONMAPPER_H
#define POOMA_BISECTIONMAPPER_H

#include <list>



/** 
 * BisectionMapper is a ContextMapper implementation.
 *
 * It assigns contexts to nodes by recursively bisecting the partition.
 */

template <int Dim>
class BisectionMapper
  : public ContextMapper<Dim>
{ 
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================
  typedef Interval<Dim>                       Domain_t;
  typedef Node<Domain_t>                      Value_t;
  typedef std::vector<Value_t *>              List_t;

  template<class Partitioner>
  BisectionMapper(const Partitioner & gp,
		   const Loc<Dim> &nblocks) 
    : blocks_m(gp.blocks())
  {
  }

  template<class Partitioner>  
  BisectionMapper(const Partitioner & gp)
     : blocks_m(gp.blocks())
  {
  }

  BisectionMapper(const Loc<Dim>& blocks)
     : blocks_m(blocks)
  {
  }

  // 

  void map(const List_t & templist) const;

  // Memeber Data
  Loc<Dim> blocks_m;

};

template <int Dim>
void BisectionMapper<Dim>::map(const List_t & templist) const 
{
  int ncontexts = Pooma::contexts();
  int npatch = 1; 
  for (int i =0;i<Dim; ++i) 
    npatch*=blocks_m[i].first();

  std::list<Domain_t> bvec;

  Domain_t allb;
  for (int i = 0; i<Dim; ++i)
    allb[i]=Interval<1>(0,blocks_m[i].first()-1);
  bvec.push_back(allb);

  while ( bvec.size() < ncontexts )
    {
      int s = 0;
      typename std::list<Domain_t>::iterator bstart = bvec.begin();
      typename std::list<Domain_t>::iterator bend = bvec.end();
      typename std::list<Domain_t>::iterator bpatch;
      // find the largest patch.
      for ( ; bstart != bend ; ++bstart)
	{
	  if (s < (*bstart).size() )
	    {
	      bpatch = bstart; 
	      s = (*bstart).size();
	    }
	}
      // now find the largest dimension on the largest patch
      int d = 0;
      int sd = 0;
      for (int i = 0; i<Dim; ++i)
	{
	  if ( sd < (*bpatch)[i].size() )
	    {
	      d = i;
	      sd = (*bpatch)[i].size();
	    }
	}
      Domain_t hi(*bpatch),lo(*bpatch);
      int lopoint = hi[d].first();
      int hipoint = hi[d].last();
      int mid     = lopoint + ( (hipoint - lopoint)/2);

      if (lopoint<=mid)
	lo[d] = Interval<1>(lopoint,mid);
      else
	lo[d] = Interval<1>(lopoint,lopoint);
      if ( hipoint>=mid+1)
	hi[d] = Interval<1>(mid+1,hipoint);
      else
	hi[d] = Interval<1>(hipoint,hipoint);
      bvec.erase(bpatch++);	
      bvec.insert(bpatch,lo);
      bvec.insert(bpatch,hi);
    }
  // now step through the intervals, using their elements as indexes into
  // all_m;
  int strides[Dim];
  strides[0] = 1;
  for ( int i=1; i<Dim; ++i)
    strides[i] = strides[i-1]*blocks_m[i-1].first();

  typename std::list<Domain_t>::iterator start = bvec.begin();
  typename std::list<Domain_t>::iterator end = bvec.end();
  int pcontext = 0;
  for ( ; start != end ; ++start)
    {
      int idx[Dim],mi[Dim],mx[Dim];
      for ( int  i = 0 ; i < Dim ; ++i)
	{
	  idx[i] = mi[i] = (*start)[i].first();
	  mx[i]  = (*start)[i].last();
	}

      while ( idx[Dim-1] <= mx[Dim-1] )
	{
	  int allIdx = 0;
	  for ( int i = 0 ; i < Dim ; ++i)
	    allIdx += idx[i]*strides[i];
	  (*templist[allIdx]).context() = pcontext;
	  ++idx[0];
	  for ( int i = 0 ; i < Dim ; ++i)
	    {
	      if ( idx[i] > mx[i] )
		{
		  if ( i!=(Dim-1) ) 
		    {
		      ++idx[i+1];
		      idx[i]=mi[i];
		    }
		  else
		    ++idx[i];
		}
	      else
		break;
	    }
	}
      ++pcontext;
    }
  // set the affinity and local ID values
  this->setAffinity(templist);
}


#endif   // POOMA_BISECTIONMAPPER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: BisectionMapper.h,v $   $Author: richi $
// $Revision: 1.10 $   $Date: 2004/11/10 22:17:08 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
