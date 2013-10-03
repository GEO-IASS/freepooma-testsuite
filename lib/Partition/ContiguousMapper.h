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

#ifndef POOMA_CONTIGUOUSMAPPER_H
#define POOMA_CONTIGUOUSMAPPER_H

/** @file
 * @ingroup Partition
 * @brief
 * ContiguousMapper is a ContextMapper implementation.
 */



/**
 * ContiguousMapper is a ContextMapper implementation.
 *
 * It assigns contexts to nodes in a contiguous fashion.
 */

template<int Dim>
class ContiguousMapper
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
  ContiguousMapper(const Partitioner & gp,
		   const Loc<Dim> &nblocks) 
    : blocks_m(gp.blocks())
  {
  }

  template<class Partitioner>
  ContiguousMapper(const Partitioner & gp)
     : blocks_m(gp.blocks())
  {
  }
  
  ContiguousMapper(const Loc<Dim> & blocks)
     : blocks_m(blocks)
  {
  }

  void map(const List_t & templist) const;

  // Member Data
  Loc<Dim> blocks_m;
};

template<int Dim>
void ContiguousMapper<Dim>::map(const List_t & templist) const 
{  
  int idx[Dim];
  for (int i = 0;i<Dim;++i) 
    idx[i]=0;
  int strides[Dim];
  strides[Dim-1] = 1;
  for ( int i=Dim-2; i>=0; --i)
    strides[i] = strides[i+1]*blocks_m[i+1].last();

  int npatch = 1;
  for (int i=0; i<Dim; ++i)
    npatch *= blocks_m[i].first();

  int ncontexts = Pooma::contexts();
  int npc = npatch/ncontexts;

  int remainder = npatch - (npc*ncontexts);

  int pcontext = 0;
  int c = 0; 
  int patchdone = 0;
  int patchleft = npatch;

  int incriment[Dim];
  for (int i =0 ; i<Dim; ++i) incriment[i] = 1;
  while ( true  )
    {
      int allIdx = 0;
      for ( int i = 0 ; i < Dim ; ++i) 
	allIdx += idx[i]*strides[i];
      (*templist[allIdx]).context() = pcontext;

      ++c;
      ++patchdone;
      --patchleft;

      if(c >= npc )
	{
	  // if we are at the end of a context, and there are still
	  // remainder patches left over, extend the number of
	  // patches per context by one, if we are not at an edge
	  // in the first dimension. However, if we have extra
	  // patches left over, and are in danger of running out
	  // of contexts to add patches to, add them now. 

	  if (c == npc && remainder >0 && 
	      ((idx[0]-1 >= 0 && idx[0]+1<=(blocks_m[0].first()-1)) || 
	       (patchleft - ((npc+1)*(ncontexts-(pcontext+1))) >=0 )) )
	    --remainder;
	  else
	    {
	      c = 0;
	      ++pcontext;
	    }
	}

      bool t = true;
      for ( int i = 0 ; i < Dim ; ++i)
	{
	  t = t && (
		    idx[i] == (blocks_m[i]-1) && incriment[i] == 1
		    || 
		    idx[i] == 0 && incriment[i] == -1);
	}
      if (t)
	break;

      idx[0] += incriment[0];
      for ( int i = 0 ; i < Dim ; ++i)
	{
	  if ( idx[i] > blocks_m[i].last()-1)
	    {
	      idx[i+1]+=incriment[i+1];
	      idx[i]=blocks_m[i].last()-1;
	      incriment[i] *= -1;
	    }
	  else if (idx[i]<0)
	    {
	      idx[i+1]+=incriment[i+1];
	      idx[i]=0;
	      incriment[i] *= -1;
	    }
	  else
	    break;
	}
    }

  ContextMapper<Dim>::setAffinity(templist);
}


#endif
