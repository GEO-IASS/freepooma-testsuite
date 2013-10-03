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

/** @file
 * @ingroup Engine
 * @brief
 * Managing intersections of engines:
 *  - Intersector           
 *    class for managing intersections of engines.
 *  - IntersectEngineTag
 *    a functor for intersecting UMP engines inside an expression
 */

#ifndef POOMA_ENGINE_INTERSECTOR_H
#define POOMA_ENGINE_INTERSECTOR_H

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include "Domain/Range.h"
#include "Utilities/RefCounted.h"
#include "Utilities/RefCountedPtr.h"
#include "Utilities/Unique.h"
#include "Layout/INode.h"
#include "Layout/GlobalIDDataBase.h"
#include "Layout/GuardLayers.h"
#include "Layout/TouchesConstruct.h"

#include <vector>


//-----------------------------------------------------------------------------
// Forward declaratations
//-----------------------------------------------------------------------------

template<int Dim>
class IntersectorData
  : public RefCounted
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef IntersectorData<Dim>                          This_t;
  typedef std::vector<int>                              IDContainer_t;
  typedef Range<7>                                      BaseDomain_t;
  typedef std::vector<BaseDomain_t>                     BaseDomainContainer_t;
  typedef INode<Dim>                                    INode_t;
  typedef std::vector<INode_t>                          INodeContainer_t;
  typedef typename INodeContainer_t::const_iterator     const_iterator;
  typedef Unique::Value_t                               LayoutID_t;
  
  enum { dimensions = Dim };
  
  
  //===========================================================================
  // Constructors
  //===========================================================================

  // Default constructor is trival.
  
  inline IntersectorData() { }

  //===========================================================================
  // Destructor
  //===========================================================================

  // Trival since members manage their own data.
  
  inline ~IntersectorData() { }
  
  template<class Engine>
  void intersect(const Engine &engine) 
  {
    // First, we need to check through our list of layout IDs and see if we've
    // either seen this layout or another layout with the same baseID before.
    
    typedef typename Engine::Layout_t Layout_t;
    const Layout_t &layout(engine.layout());

    int n = ids_m.size();
    for (int i = 0; i < n; ++i)
    {
      // If we've seen this ID before, we're done.
 
      if (ids_m[i] == layout.ID())
	return;
          
      // If we've seen the base ID before and the base domain is the same
      // we're done.
        
      if (baseIDs_m[i] == layout.baseID()
	  && sameBaseDomain(i, layout.baseDomain()))
      {
	shared(layout.ID(),ids_m[i]);
	return;
      }
    }
          
    touches(layout);
  }

  template<class Engine, int Dim2>
  bool intersect(const Engine &engine, const GuardLayers<Dim2> &guard,
		 GuardLayers<Dim2> &usedGuards) 
  {
    CTAssert(Engine::dimensions == Dim);

    // First, we need to check through our list of layout IDs and see if we've
    // either seen this layout or another layout with the same baseID before.
    
    typedef typename Engine::Layout_t Layout_t;
    const Layout_t &layout(engine.layout());

    int n = ids_m.size();
    for (int i = 0; i < n; ++i)
    {
      // If we've seen this ID before, we're done.
        
      if (ids_m[i] == layout.ID())
	return false;
   
      // If we've seen the base ID before and the base domain is the same
      // we're done.
        
      if (baseIDs_m[i] == layout.baseID()
	  && sameBaseDomain(i, layout.baseDomain(), guard))
      {
	shared(layout.ID(),ids_m[i]);

	// was: return (!sameBaseDomain(i,layout.baseDomain()));

        // We should be able to find out the actual shape of the
	// used internal guards here, rather than just returning bool.
	// Something like:

	// But what do, if Dim2 > baseDims_m[i]!?
	if (baseDims_m[i] < Dim2)
	  return true;

	bool used = false;
	for (int j = 0; j < Dim2; j++)
	{
	  usedGuards.lower(j) = std::max(0, baseDomains_m[i][j].first() - layout.baseDomain()[j].first());
	  if (usedGuards.lower(j) != 0)
	    used = true;
	  usedGuards.upper(j) = std::max(0, layout.baseDomain()[j].last() - baseDomains_m[i][j].last());
	  if (usedGuards.upper(j) != 0)
	    used = true;
	}
	return used;
      }
    }
          
    // current touches operation works on the owned region, so we don't
    // use the guard cells.  If we start using touchesAlloc, then you
    // need to return true here, and the bypass calculation above
    // becomes somewhat more complicated.

    touches(layout);
    return false;
  }

  // Check to see whether or not we've seen the specified domain in the
  // base domain list.
  
  template<int Dim2>
  bool sameBaseDomain(int i, const Range<Dim2> &domain,
		      const GuardLayers<Dim2> &guard)
  {
    // If they're not the same dimensionality, there is no match.
        
    if (baseDims_m[i] != Dim2)
      return false;

    // Check each 1D domain and report failure if any one is not equal.
        
    for (int j = 0; j < Dim2; j++)
    {
      if (baseDomains_m[i][j].stride() != domain[j].stride()) return false;
      if (baseDomains_m[i][j].first()  > domain[j].first() + guard.lower(j))
	return false;
      if (baseDomains_m[i][j].last()  < domain[j].last() - guard.upper(j))
	return false;
    }
      
    return true;
  }

  template<int Dim2>
  bool sameBaseDomain(int i, const Range<Dim2> &domain)
  {
    // If they're not the same dimensionality, there is no match.
        
    if (baseDims_m[i] != Dim2)
      return false;

    // Check each 1D domain and report failure if any one is not equal.
        
    for (int j = 0; j < Dim2; j++)
      if (baseDomains_m[i][j] != domain[j]) return false;
      
    return true;
  }

  template<int Dim2>
  bool sameBaseDomain(int i, const Interval<Dim2> &domain,
		      const GuardLayers<Dim2> & guard)
  {
    // If they're not the same dimensionality, there is no match.
        
    if (baseDims_m[i] != Dim2)
      return false;

    // Check each 1D domain and report failure if any one is not equal.
        
    for (int j = 0; j < Dim2; j++)
    {
      if (baseDomains_m[i][j].stride() != 1) return false;
      if (baseDomains_m[i][j].first()  > domain[j].first() + guard.lower(j))
	return false;
      if (baseDomains_m[i][j].last()  < domain[j].last() - guard.upper(j))
	return false;
    }
      
    return true;
  }

  template<int Dim2>
  void pushBaseDomain(const Range<Dim2> &domain)
  {
    // Grab the index we will be working with.
        
    int i = baseDomains_m.size();

    // Push the number of active dimensions and an empty BaseDomain_t
    // onto our lists.
    
    baseDims_m.push_back(Dim2);    
    baseDomains_m.push_back(BaseDomain_t());

    // Reset the domains.
    
    for (int j = 0; j < Dim2; j++)
      baseDomains_m[i][j] = 
        Range<1>(domain[j].first(), domain[j].last(), domain[j].stride());
  }
  
  // The following are exact copies of the above with Range replaced
  // by Interval. These could be done with Domain, but CodeWarrior
  // has a problem with this. 

  template<int Dim2>
  bool sameBaseDomain(int i, const Interval<Dim2> &domain)
  {
    // If they're not the same dimensionality, there is no match.
        
    if (baseDims_m[i] != Dim2)
      return false;

    // Check each 1D domain and report failure if any one is not equal.
        
    for (int j = 0; j < Dim2; j++)
      if (baseDomains_m[i][j] != domain[j]) return false;
      
    return true;
  }

  template<int Dim2>
  void pushBaseDomain(const Interval<Dim2> &domain)
  {
    // Grab the index we will be working with.
        
    int i = baseDomains_m.size();

    // Push the number of active dimensions and an empty BaseDomain_t
    // onto our lists.
    
    baseDims_m.push_back(Dim2);    
    baseDomains_m.push_back(BaseDomain_t());

    // Reset the domains.
    
    for (int j = 0; j < Dim2; j++)
      baseDomains_m[i][j] = 
        Range<1>(domain[j].first(), domain[j].last(), domain[j].stride());
  }
  
  template<class Layout>
  void touches(const Layout &l)
  {
    int n = ids_m.size();

    // This is a new layout that will contribute unique intersections, save
    // the data.
        
    ids_m.push_back(l.ID());
    baseIDs_m.push_back(l.baseID());
    pushBaseDomain(l.baseDomain());

    // If we previously had no layouts stored, we simply need to fill our list
    // with INodes constructed using the nodes from the iterators.
    // If this is the second (or later) layout we're being asked to intersect,
    // we need to loop over the existing INodes and have the 
    // layout's touches() function insert new versions at the end of the list.
    
    if (n == 0)
    {
      typename Layout::const_iterator p = l.beginGlobal();
      while (p != l.endGlobal())
      {
	if (! (*p).domain().empty())
	  inodes_m.push_back(INode_t(*p,l.ID(),
				     &(gidStore_m)));
	++p;
      }
    }
    else
    {
      int ni = inodes_m.size();
      for (int i = 0; i < ni; i++)
	l.touches(inodes_m[i].domain(),
		  std::back_inserter(inodes_m),
		  inodes_m[i].touchesConstructINode(l.ID())
		  );
        
      // The INodes we just pushed supercede the previous set,
      // which we now erase.
        
      inodes_m.erase(inodes_m.begin(),
			     inodes_m.begin() + ni);
    }  
  }

  inline
  void shared(LayoutID_t id1, LayoutID_t id2)
  {
    gidStore_m.shared(id1,id2);
  }

  //private:  

  // Don't ever want to copy one of these.
  IntersectorData(const This_t &);
  This_t &operator=(const This_t &);
  
  //===========================================================================
  // Data members
  //===========================================================================
  
  IDContainer_t ids_m, baseIDs_m, baseDims_m;
  BaseDomainContainer_t baseDomains_m;
  INodeContainer_t inodes_m;
  GlobalIDDataBase gidStore_m;
  
};


template<int Dim>
class Intersector
{
public:

  //===========================================================================
  // Exported typedefs and constants
  //===========================================================================

  typedef IntersectorData<Dim>                          IntersectorData_t;
  typedef Intersector<Dim>                              This_t;
  typedef typename IntersectorData_t::IDContainer_t     IDContainer_t;
  typedef typename IntersectorData_t::BaseDomain_t      BaseDomain_t;
  typedef typename IntersectorData_t::BaseDomainContainer_t
                                                        BaseDomainContainer_t;
  typedef typename IntersectorData_t::INode_t           INode_t;
  typedef typename IntersectorData_t::INodeContainer_t  INodeContainer_t;
  typedef typename IntersectorData_t::const_iterator    const_iterator;
  typedef RefCountedPtr<IntersectorData_t>              DataPtr_t;
  
  enum { dimensions = Dim };
  
  Intersector()
    : pdata_m(new IntersectorData_t())
  { }

  Intersector(const This_t &model)
    : pdata_m(model.pdata_m)
  { }

  This_t &operator=(const This_t &model)
  {
    if (this != &model)
      pdata_m = model.pdata_m;
    return *this;
  }

  ~Intersector() { }

  inline DataPtr_t &data() { return pdata_m; }
  inline const DataPtr_t &data() const { return pdata_m; }

  //===========================================================================
  // Accessors
  //===========================================================================

  // STL iterator support.
  
  inline const_iterator begin() const { return data()->inodes_m.begin(); }
  
  inline const_iterator end() const { return data()->inodes_m.end(); }

  inline int size() const { return data()->inodes_m.size(); }

  //===========================================================================
  // Intersect routines
  //===========================================================================

  // All domains.
  
  template<class Engine>
  inline
  void intersect(const Engine &l) 
  {
    data()->intersect(l);
  }

  template<class Engine, int Dim2>
  inline
  bool intersect(const Engine &l, const GuardLayers<Dim2> &guard, GuardLayers<Dim2> &usedGuards) 
  {
    return (data()->intersect(l,guard,usedGuards));
  }

private:
  DataPtr_t pdata_m;
};

#endif // POOMA_ENGINE_INTERSECTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Intersector.h,v $   $Author: richard $
// $Revision: 1.17 $   $Date: 2004/11/01 18:16:37 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
