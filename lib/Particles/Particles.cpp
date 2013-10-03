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

// Include files

#include "Domain/Contains.h"
#include "Domain/Intersect.h"
#include "Domain/IteratorPairDomain.h"

//-----------------------------------------------------------------------------
// Classes: 
//   Particles<Traits> template definitions.
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//
// Default constructor: this should be used if
//   o The particle layout has a useful default constructor, and
//   o The attribute layout has a useful default constructor, and
//   o The initialize() method will be used later.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
Particles<ParticleTraits>::Particles()
  : destroyMethod_m(BackFill::code), destroyList_m(0)
{
}


//-----------------------------------------------------------------------------
//
// Construct a Particles object, using the given particle layout object.
// Let the particle layout object initialize the attribute data layout.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
Particles<ParticleTraits>::Particles(const ParticleLayout_t& pl)
  : particleLayout_m(pl),
    destroyMethod_m(BackFill::code)
{
  // The particle layout object may need to have the attribute layout
  // organized in a specific way.  Thus, we defer attribute layout
  // initialization to the particle layout object.

  particleLayout_m.initializeAttributeLayout(attributeLayout_m);

  // Also initialize the destroy list.

  int npatches = attributeLayout_m.sizeLocal();
  typedef DynamicArray<int,Dynamic> dlist_t;
  destroyList_m = new dlist_t[npatches];
  for (int i=0; i < npatches; ++i)
    destroyList_m[i].initialize(Interval<1>());
}


//-----------------------------------------------------------------------------
//
// Copy constructor: only copy the layout data and destroy flag, do
// NOT copy the attributes.  The new particle object starts out with
// just layouts and no attributes (except for the internal destroy list).
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
Particles<ParticleTraits>::Particles(const This_t &p)
  : particleLayout_m(p.particleLayout_m),
    attributeLayout_m(p.attributeLayout_m),
    destroyMethod_m(p.destroyMethod_m)
{
  // Also initialize the destroy list.

  int npatches = attributeLayout_m.sizeLocal();
  typedef DynamicArray<int,Dynamic> dlist_t;
  destroyList_m = new dlist_t[npatches];
  for (int i=0; i < npatches; ++i)
    destroyList_m[i].initialize(Interval<1>());
}


//-----------------------------------------------------------------------------
//
// Initialize this object with the given particle layout.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
void Particles<ParticleTraits>::initialize(const ParticleLayout_t& pl)
{
  // Initialize the particle layout, by copying info from the provided
  // particle layout.

  particleLayout_m.initialize(pl);

  // The particle layout object may need to have the attribute layout
  // organized in a specific way.  Thus, we defer attribute layout
  // initialization to the particle layout object.

  particleLayout_m.initializeAttributeLayout(attributeLayout_m);

  // Also initialize the destroy list.

  int npatches = attributeLayout_m.sizeLocal();
  typedef DynamicArray<int,Dynamic> dlist_t;
  destroyList_m = new dlist_t[npatches];
  for (int i=0; i < npatches; ++i)
    destroyList_m[i].initialize(Interval<1>());
}


//-----------------------------------------------------------------------------
//
// Destructor
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
Particles<ParticleTraits>::~Particles()
{
  delete [] destroyList_m;
}


//-----------------------------------------------------------------------------
//
// Return the number of particles we have queued to destroy in the
// next performDestroy call, for all the patches (if pid < 0), or just
// for the specified patch.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
typename Particles<ParticleTraits>::Size_t
Particles<ParticleTraits>::deferredDestroyAmount(PatchID_t pid) const
{
  Size_t retval = 0;

  // Add in all relevant sizes

  if (pid < 0)
    {
      PatchID_t i, npatch = attributeLayout_m.sizeLocal();
      for (i = 0; i < npatch; ++i)
	retval += destroyList(i).domain().size();
    }
  else
    {
      retval = destroyList(pid).domain().size();
    }

  return retval;
}


//-----------------------------------------------------------------------------
//
// Destroy set of particles specified by DomainType, which may be a
// 1D Range of particle index numbers or a list of index numbers.
// This version does the destroy immediately, instead of deferring
// it until later.  If the patchID is < 0, the domain should be
// a global domain.  If the patchID is >= 0, the domain should be
// a "local" one, where index = 0 refers to the first particle in
// that patch.  If the final argument is true (the default), the
// domain is renumbered after the destroy is complete.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
template <class DomainType>
void Particles<ParticleTraits>::destroy(const DomainType& domain,
                                        PatchID_t pid, bool renum)
{
  // Ask the attribute layout to ask all registered attributes
  // to destroy these elements.

  if (destroyMethod_m == DynamicEvents::backfill)
    {
      if (pid < 0)
        attributeLayout_m.destroy(domain, BackFill());
      else
        attributeLayout_m.destroy(domain, pid, BackFill());
    }
  else if (destroyMethod_m == DynamicEvents::shiftup)
    {
      if (pid < 0)
        attributeLayout_m.destroy(domain, ShiftUp());
      else
        attributeLayout_m.destroy(domain, pid, ShiftUp());
    }
  else
    {
      PInsist(false, "Unknown destroy method in Particles::destroy!");
    }

  // If requested, renumber the domain now.
  
  if (renum)
    renumber();
}

//-----------------------------------------------------------------------------
//
// Same as above, but use begin/end iterator pair to specify domain.
//
//-----------------------------------------------------------------------------

template <class ParticleTraits>
template <class Iter>
void Particles<ParticleTraits>::destroy(Iter begin, Iter end,
                                        PatchID_t pid, bool renum)
{
  // Ask the attribute layout to ask all registered attributes
  // to destroy these elements.

  if (destroyMethod_m == DynamicEvents::backfill)
    {
      if (pid < 0) {
	Pooma::IteratorPairDomain<Iter> domain(begin,end);
        attributeLayout_m.destroy(domain, BackFill());
      }
      else {
	Pooma::IteratorPairDomain<const int*> domain(begin,end);
        attributeLayout_m.destroy(domain, pid, BackFill());
      }
    }
  else if (destroyMethod_m == DynamicEvents::shiftup)
    {
      if (pid < 0) {
	Pooma::IteratorPairDomain<Iter> domain(begin,end);
        attributeLayout_m.destroy(domain, ShiftUp());
      }
      else {
	Pooma::IteratorPairDomain<const int*> domain(begin,end);
        attributeLayout_m.destroy(domain, pid, ShiftUp());
      }
    }
  else
    {
      PInsist(false, "Unknown destroy method in Particles::destroy!");
    }

  // If requested, renumber the domain now.
  
  if (renum)
    renumber();
}

//-----------------------------------------------------------------------------
//
// Destroy the particles in the given global domain by storing them
// a list for use in a later performDestroy call.  This lets you
// perform several different destroy operations without going through
// all the shifting-around work until destroys are done.  This does
// NOT renumber the domain.  If the patchID is < 0, the domain should
// be a global domain. If the patchID is >= 0, the domain should
// be a "local" one, where index = 0 refers to the first particle in
// that patch.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
template <class DomainType>
void Particles<ParticleTraits>::deferredDestroy(const DomainType& domain,
                                                PatchID_t pid)
{
  if (pid >= 0)
    {
      PAssert(pid < attributeLayout_m.sizeLocal());

      // This is only for one patch request, and the domain values are
      // already relative, so we can just stuff them at the end of the
      // patch destroy list for the specified patch.

      int i = 0;
      int destroys = domain[0].size();
      int next = destroyList(pid).domain().size();
      destroyList(pid).create(destroys);
      while (i < destroys)
        destroyList(pid)(next++) = domain[0](i++);
    }
  else
    {
      // Check that the destroy list is contained by the global domain.

      PInsist(contains(attributeLayout_m.domain(),domain),
	      "Destroy request outside of global domain!");

      // We assume the domain is sorted in ascending order.
      // Convert this list into a set of patch-specific destroy requests
      // by intersecting the list with the domain of each of the 
      // separate patches.  This code was borrowed from the destroy()
      // method in DynamicLayout.

      // Find pieces of this total destroy domain in each subdomain,
      // and destroy them.

      int is = 0, ie = 0;
      PatchID_t patch = 0, numPatches = attributeLayout_m.sizeLocal();

      // Skip to the first non-empty local patch

      bool found = false;
      while (!found && patch < numPatches) {
        found = !attributeLayout_m.ownedDomain(patch).empty();
        if (!found) patch++;
      }
      if (!found) return;

      // Some portion of the destroy domain may precede all of the 
      // domain controlled by this context, so skip that part.

      while (domain(is) < attributeLayout_m.ownedDomain(patch).first() && 
             is < domain.size()) is++;
      ie = is;

      while (patch < numPatches && ie < domain.size())
        {
	  // find last item on this patch
	  for ( ; patch < numPatches && ie < domain.size(); ++ie)
	    if (domain(ie)>attributeLayout_m.ownedDomain(patch).last()) break;

	  // if we didn't find any intersection, go to next patch
	  if (ie == is)
	    {
	      ++patch;
	      continue;
	    }

	  // put the intersection points into this patch's destroy list
	  int patchOffset = attributeLayout_m.ownedDomain(patch).first();
	  int currSize = destroyList(patch).domain().size();
	  int newDestroys = ie - is;
	  destroyList(patch).create(newDestroys);

          // insert local offset index values into destroy list patch
          for (int ii = 0; ii < newDestroys; ++ii)
            destroyList(patch)(currSize+ii) = domain(is+ii) - patchOffset;

	  // Move on to next local patch
	  ++patch;
	  is = ie;
        }

    }
}

//-----------------------------------------------------------------------------
//
// Same as above, but use begin/end iterator pair to specify domain.
//
//-----------------------------------------------------------------------------

template <class ParticleTraits>
template <class Iter>
void Particles<ParticleTraits>::deferredDestroy(Iter begin, Iter end,
                                                PatchID_t pid)
{
  if (pid >= 0)
    {
      Pooma::IteratorPairDomain<const int*> domain(begin,end);

      PAssert(pid < attributeLayout_m.sizeLocal());

      // This is only for one patch request, and the domain values are
      // already relative, so we can just stuff them at the end of the
      // patch destroy list for the specified patch.

      int i = 0;
      int destroys = domain[0].size();
      int next = destroyList(pid).domain().size();
      destroyList(pid).create(destroys);
      while (i < destroys)
        destroyList(pid)(next++) = domain[0](i++);
    }
  else
    {
      Pooma::IteratorPairDomain<Iter> domain(begin,end);

      // Check that the destroy list is contained by the global domain.

      PInsist(contains(attributeLayout_m.domain(),domain),
	      "Destroy request outside of global domain!");

      // We assume the domain is sorted in ascending order.
      // Convert this list into a set of patch-specific destroy requests
      // by intersecting the list with the domain of each of the 
      // separate patches.  This code was borrowed from the destroy()
      // method in DynamicLayout.

      // Find pieces of this total destroy domain in each subdomain,
      // and destroy them.

      int is = 0, ie = 0;
      PatchID_t patch = 0, numPatches = attributeLayout_m.sizeLocal();

      // Skip to the first non-empty local patch

      bool found = false;
      while (!found && patch < numPatches) {
        found = !attributeLayout_m.ownedDomain(patch).empty();
        if (!found) patch++;
      }
      if (!found) return;

      // Some portion of the destroy domain may precede all of the 
      // domain controlled by this context, so skip that part.

      while (domain(is) < attributeLayout_m.ownedDomain(patch).first() && 
             is < domain.size()) is++;
      ie = is;

      while (patch < numPatches && ie < domain.size())
        {
	  // find last item on this patch
	  for ( ; patch < numPatches && ie < domain.size(); ++ie)
	    if (domain(ie)>attributeLayout_m.ownedDomain(patch).last()) break;

	  // if we didn't find any intersection, go to next patch
	  if (ie == is)
	    {
	      ++patch;
	      continue;
	    }

	  // put the intersection points into this patch's destroy list
	  int patchOffset = attributeLayout_m.ownedDomain(patch).first();
	  int currSize = destroyList(patch).domain().size();
	  int newDestroys = ie - is;
	  destroyList(patch).create(newDestroys);

          // insert local offset index values into destroy list patch
          for (int ii = 0; ii < newDestroys; ++ii)
            destroyList(patch)(currSize+ii) = domain(is+ii) - patchOffset;

	  // Move on to next local patch
	  ++patch;
	  is = ie;
        }

    }
}

//-----------------------------------------------------------------------------
//
// Carry out all the previously requested destroy commands, and
// renumber the domain if requested.  If the requested patch ID is
// < 0, do all deferred destroys, otherwise do it for just the requested
// patch.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
void Particles<ParticleTraits>::performDestroy(PatchID_t pid, bool renum)
{
  PatchID_t i, npatch = attributeLayout_m.sizeLocal();

  if (pid < 0)
    {
      for (i = 0; i < npatch; ++i)
	{
	  // skip this patch if there are no destroy requests for it.

	  if (destroyList(i).domain().empty())
	    continue;

	  // do the destroys on this patch

	  if (destroyMethod_m == DynamicEvents::backfill)
	    {
	      attributeLayout_m.destroy(IndirectionList<int>(destroyList(i)),
					i, BackFill());
	    }
	  else if (destroyMethod_m == DynamicEvents::shiftup)
	    {
	      attributeLayout_m.destroy(IndirectionList<int>(destroyList(i)),
					i, ShiftUp());
	    }
	  else
	    {
	      PInsist(false, "Unknown destroy method in Particles::destroy!");
	    }

	  // clear the destroy list for this patch

	  destroyList(i).destroy(destroyList(i).domain());
	}
    }
  else
    {
      // Just destroy for the given patch

      PAssert(pid < npatch);
      i = pid;

      // Do the destroy if we have something to do

      if (!destroyList(i).domain().empty())
	{
	  if (destroyMethod_m == DynamicEvents::backfill)
	    {
	      attributeLayout_m.destroy(IndirectionList<int>(destroyList(i)),
					i, BackFill());
	    }
	  else if (destroyMethod_m == DynamicEvents::shiftup)
	    {
	      attributeLayout_m.destroy(IndirectionList<int>(destroyList(i)),
					i, ShiftUp());
	    }
	  else
	    {
	      PInsist(false, "Unknown destroy method in Particles::destroy!");
	    }

	  // clear the destroy list for this patch
	  
	  destroyList(i).destroy(destroyList(i).domain());
	}
    }

  // if requested, recompute the global domain of the Attributes

  if (renum)
    renumber();
}


//-----------------------------------------------------------------------------
//
// Add a new boundary condition of given type, acting on info from the
// Subject and affecting values in the Object.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
template <class Subject, class Object, class BCType>
void Particles<ParticleTraits>::addBoundaryCondition(const Subject &s,
						     const Object &o,
						     const BCType &bc)
{
  // create BC with this subject and object and add to our list
  
  bcList_m.addBC(s, o, bc);
}


//-----------------------------------------------------------------------------
//
// Add a new boundary condition of given type, acting on info from the
// Subject and affecting values in the same Subject.
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
template <class Subject, class BCType>
void Particles<ParticleTraits>::addBoundaryCondition(const Subject &s,
						     const BCType &bc)
{
  // create BC with this subject and add to our list

  bcList_m.addBC(s, bc);
}


//-----------------------------------------------------------------------------
//
// Remove the Nth boundary condition in this Particles object.
//
//-----------------------------------------------------------------------------

template <class ParticleTraits>
void Particles<ParticleTraits>::
removeBoundaryCondition(ParticleBCList::Size_t i)
{
  // delete ith BC and remove from our list
  
  bcList_m.removeBC(i);
}


//-----------------------------------------------------------------------------
//
// Remove all the boundary conditions
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
void Particles<ParticleTraits>::removeBoundaryConditions()
{
  // Just keep removing the first one until none are left.

  while (bcList_m.size() > 0)
    bcList_m.removeBC(0);
}


//-----------------------------------------------------------------------------
//
// Apply the boundary conditions to the existing particles for all
// patches, or just one (if all, the value for the patch ID is < 0).
//
//-----------------------------------------------------------------------------

template<class ParticleTraits>
void Particles<ParticleTraits>::applyBoundaryConditions(PatchID_t pid)
{
  // loop over the list of particle BC's and apply each one to all
  // the different patches.

  ParticleBCList::Size_t i, n = bcList_m.size();
  for (i=0; i<n; ++i)
    bcList_m(i)->applyBoundaryCondition(pid);

  // wait for them to finish ... this could be replaced with some
  // waits on a counting semaphore in the future.  This is only
  // needed if we're doing this for all patches.

  if (pid < 0)
    Pooma::blockAndEvaluate();
}


//-----------------------------------------------------------------------------
//
// Create the specified number of particles, total, spread across
// all patches.  An equal (or close to) number of particles are added
// to each patch on each context.  This must be called in SPMD-style.
// The domain is renumbered at the end if the second argument is true.
//
//-----------------------------------------------------------------------------

template <class ParticleTraits>
void Particles<ParticleTraits>::globalCreate(Size_t np, bool renum)
{
  if (np > 0)
    {
      // Divide np by the number of patches

      Size_t perpatch = np / attributeLayout_m.sizeGlobal();
      Size_t extra    = np % attributeLayout_m.sizeGlobal();

      // Find out how many extra particles go in our local context

      Size_t myextra  = extra / Pooma::contexts();

      if ( Pooma::context() < (extra % Pooma::contexts()) )
        myextra++;

      // Create elements in each local patch that has particles
      // Do not renumber now; we'll do this at the end if requested

      int npatches = attributeLayout_m.sizeLocal();
      for (int i=0; i < npatches; ++i)
        create(perpatch + (i < myextra ? 1 : 0), i, false);

      // Renumber the domain now

      if (renum)
        renumber();
    }
}


//-----------------------------------------------------------------------------
//
// Print the contents of a Particles object to the given stream.
//
//-----------------------------------------------------------------------------

template <class ParticleTraits>
template <class Out>
void Particles<ParticleTraits>::print(Out& o) const
{
  o << "Particles:\n";
  o << "  Particle layout     = " << particleLayout_m << "\n";
  o << "  Attribute layout    = " << attributeLayout_m << "\n";
  if (destroyMethod_m == DynamicEvents::backfill)
    o << "  Destroy Method      = BackFill\n";
  else if (destroyMethod_m == DynamicEvents::shiftup)
    o << "  Destroy Method      = ShiftUp\n";
  else
    o << "  Destroy Method      = Unknown\n";
  o << "  Cached Destroy Cmds =\n";
  int i, npatches = attributeLayout_m.sizeLocal();
  for (i = 0; i < npatches; ++i)
    o << "    Local Patch " << i << ": " << destroyList(i);
  o << "\n";
  o << "  Boundary conditions = " << bcList_m << "\n";
  o << "  Attribute values:\n";
  o << attributes_m << "\n";
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Particles.cpp,v $   $Author: richard $
// $Revision: 1.31 $   $Date: 2004/11/01 18:16:59 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
