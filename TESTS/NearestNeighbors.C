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
// Test computing the set of nearest neighbors.
//-----------------------------------------------------------------------------

#include "Pooma/Fields.h"
#include "Utilities/Tester.h"
#include <vector>
#include <algorithm>


// Instead of checking all the nearest neighbor's FieldOffset values,
// we check "random" values.

// Check if a FieldOffset is in the FieldOffsetList.
template <int Dim>
inline bool
checkForFieldOffset(const FieldOffsetList<Dim> &lst,
		    const FieldOffset<Dim> &offset)
{
  for (typename FieldOffsetList<Dim>::size_type index = 0;
       index < lst.size();
       ++index)
    if (lst[index] == offset)
      return true;
  return false;
}


// Check for a particular FieldOffset within a vector of FieldOffsetLists.
// The arguments should be:
//   a tester object,
//   a C-string describing the test(s),
//   a vector of FieldOffsetLists,
//   the vector's correct size,
//   the index of the particular FieldOffsetList to check,
//   the FieldOffsetList correct size,
//   a FieldOffset that should be (or not be) present in the list,
//   whether the FieldOffset should be present.

template <int Dim>
inline bool
checkFieldOffset(Pooma::Tester &tester,
		 const char *testExplanation,
		 const std::vector<FieldOffsetList<Dim> > &nn,
		 const typename std::vector<FieldOffsetList<Dim> >::size_type nnSize,
		 const typename std::vector<FieldOffsetList<Dim> >::size_type listNum,
		 const typename FieldOffsetList<Dim>::size_type listSize,
		 const FieldOffset<Dim> &offset,
		 const bool offsetPresent = true)
{
  PInsist(listNum < nnSize,
	  "Incorrect std::vector<FieldOffsetList> index.");

  return 
    tester.check(testExplanation, nn.size() == nnSize) &&
    tester.check(testExplanation, nn[listNum].size() == listSize) &&
    tester.check(testExplanation,
		 checkForFieldOffset(nn[listNum], offset) == offsetPresent);
}


// Check that the distances are the same for all input values.

// Compute the Manhattan distance of a difference between positions.

template <int Dim, class T, class E>
inline double
manhattanDistance(const Vector<Dim, T, E> &difference)
{
  double answer = 0.0;
  for (int coordinate = Dim-1; coordinate >= 0; --coordinate)
    answer += std::abs(difference(coordinate));
  return answer;
}

// Compute the Manhattan distance between an input centering's value
// shifted by an FieldOffset and an output centering's value.

template <int Dim>
inline double
manhattanDistance(const Centering<Dim> &inputCentering,
		  const FieldOffset<Dim> &offset,
		  const Centering<Dim> &outputCentering,
		  const int outputIndex)
{
  // Compute the actual input position.
  return manhattanDistance(outputCentering.position(outputIndex) -
			   inputPosition(inputCentering, offset));
}

// Check that the distance between the input and output values are the
// same for all the input values.

template <int Dim>
inline bool
sameDistances(const std::vector<FieldOffsetList<Dim> > &nn,
	      const Centering<Dim> &inputCentering,
	      const Centering<Dim> &outputCentering)
{
  typedef typename std::vector<FieldOffsetList<Dim> >::size_type nn_size_type;
  typedef typename FieldOffsetList<Dim>::size_type fol_size_type;
  PInsist(nn.size() == outputCentering.size(),
	  "Nearest neighbors and output centering must have the same length.");

  for (nn_size_type outputIndex = 0; outputIndex < nn.size(); ++outputIndex)
    {
      const double distance =
	manhattanDistance(inputCentering, nn[outputIndex][0],
			  outputCentering, outputIndex);
      for (fol_size_type inputIndex = 1;
	   inputIndex < nn[outputIndex].size();
	   ++inputIndex)
	if (std::abs(distance -
		     manhattanDistance(inputCentering,
				       nn[outputIndex][inputIndex],
				       outputCentering, outputIndex))
	    > 1.0e-08)
	  return false;
    }

  return true;
}


int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  Centering<2> inputCenteringTwo, outputCenteringTwo;
  Centering<3> inputCenteringThree, outputCenteringThree;
  FieldOffsetList<2> fieldOffsetListTwo;

  // Test 2D Continuous Cell -> Continuous Cell.

  inputCenteringTwo = canonicalCentering<2>(CellType, Continuous, AllDim);
  outputCenteringTwo = canonicalCentering<2>(CellType, Continuous, AllDim);
  checkFieldOffset(tester, "cell->cell intracell", 
		   nearestNeighbors(inputCenteringTwo, outputCenteringTwo,
				    true),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 1,
		   FieldOffset<2>(Loc<2>(0)));

  checkFieldOffset(tester, "cell->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 1,
		   FieldOffset<2>(Loc<2>(0)));

  // Test 2D Continuous Vertex -> Continuous Cell.
  
  inputCenteringTwo = canonicalCentering<2>(VertexType, Continuous, AllDim);
  outputCenteringTwo = canonicalCentering<2>(CellType, Continuous, AllDim);
  checkFieldOffset(tester, "vertex->cell intracell", 
		   nearestNeighbors(inputCenteringTwo, outputCenteringTwo,
				    true),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 1,
		   FieldOffset<2>(Loc<2>(0)));
  tester.check("vertex->cell intracell distances",
	       sameDistances(nearestNeighbors(inputCenteringTwo,
					      outputCenteringTwo,
					      true),
			     inputCenteringTwo,
			     outputCenteringTwo));

  checkFieldOffset(tester, "vertex->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 4,
		   FieldOffset<2>(Loc<2>(0)));
  checkFieldOffset(tester, "vertex->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 4,
		   FieldOffset<2>(Loc<2>(1,1)));
  tester.check("vertex->cell intercell distances",
	       sameDistances(nearestNeighbors(inputCenteringTwo,
					      outputCenteringTwo),
			     inputCenteringTwo,
			     outputCenteringTwo));

  fieldOffsetListTwo =
    nearestNeighbors(inputCenteringTwo,
		     FieldOffset<2>(Loc<2>(0,0)), outputCenteringTwo);
  tester.check("vertex->cell intercell",
	       fieldOffsetListTwo.size() == 4 && 
	       checkForFieldOffset(fieldOffsetListTwo,
				   FieldOffset<2>(Loc<2>(0,0))));

  fieldOffsetListTwo =
    nearestNeighbors(inputCenteringTwo,
		     FieldOffset<2>(Loc<2>(0,0)), outputCenteringTwo, true);
  tester.check("vertex->cell intracell",
	       fieldOffsetListTwo.size() == 1 && 
	       checkForFieldOffset(fieldOffsetListTwo,
				   FieldOffset<2>(Loc<2>(0,0))));

  // Test 2D Discontinuous Vertex -> Continuous Cell.
  
  inputCenteringTwo = canonicalCentering<2>(VertexType, Discontinuous, AllDim);
  outputCenteringTwo = canonicalCentering<2>(CellType, Continuous, AllDim);
  checkFieldOffset(tester, "discontinuous vertex->cell intracell", 
		   nearestNeighbors(inputCenteringTwo, outputCenteringTwo,
				    true),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 4,
		   FieldOffset<2>(Loc<2>(0), 0));
  checkFieldOffset(tester, "discontinuous vertex->cell intracell", 
		   nearestNeighbors(inputCenteringTwo, outputCenteringTwo,
				    true),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 4,
		   FieldOffset<2>(Loc<2>(0), 3));
  tester.check("discontinuous vertex->cell intracell distances",
	       sameDistances(nearestNeighbors(inputCenteringTwo,
					      outputCenteringTwo,
					      true),
			     inputCenteringTwo,
			     outputCenteringTwo));

  checkFieldOffset(tester, "discontinuous vertex->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 16,
		   FieldOffset<2>(Loc<2>(0), 0));
  checkFieldOffset(tester, "discontinuous vertex->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 16,
		   FieldOffset<2>(Loc<2>(0), 3));
  checkFieldOffset(tester, "discontinuous vertex->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 16,
		   FieldOffset<2>(Loc<2>(0), 3));
  checkFieldOffset(tester, "discontinuous vertex->cell intercell", 
		   nearestNeighbors(inputCenteringTwo,
				    outputCenteringTwo),
		   /* vector: */ 1, 0, 
		   /* FieldOffsetList: */ 16,
		   FieldOffset<2>(Loc<2>(-1,0), 3), false);
  tester.check("discontinuous vertex->cell intercell distances",
	       sameDistances(nearestNeighbors(inputCenteringTwo,
					      outputCenteringTwo),
			     inputCenteringTwo,
			     outputCenteringTwo));


  // Test 3D Continuous Face -> Continuous Edge.
  
  inputCenteringThree = canonicalCentering<3>(FaceType, Continuous, AllDim);
  outputCenteringThree = canonicalCentering<3>(EdgeType, Continuous, AllDim);
  checkFieldOffset(tester, "face->edge intracell", 
		   nearestNeighbors(inputCenteringThree,
				    outputCenteringThree,
				    true),
		   /* vector: */ 3, 1,
		   /* FieldOffsetList: */ 2,
		   FieldOffset<3>(Loc<3>(0), 2));
  tester.check("face->edge intracell distances",
	       sameDistances(nearestNeighbors(inputCenteringThree,
					      outputCenteringThree,
					      true),
			     inputCenteringThree,
			     outputCenteringThree));

  checkFieldOffset(tester, "face->edge intercell", 
		   nearestNeighbors(inputCenteringThree,
				    outputCenteringThree),
		   /* vector: */ 3, 1,
		   /* FieldOffsetList: */ 4,
		   FieldOffset<3>(Loc<3>(-1,0,0), 2));
  checkFieldOffset(tester, "face->edge intercell", 
		   nearestNeighbors(inputCenteringThree,
				    outputCenteringThree),
		   /* vector: */ 3, 2,
		   /* FieldOffsetList: */ 4,
		   FieldOffset<3>(Loc<3>(-1,0,0), 1));
  checkFieldOffset(tester, "face->edge intercell", 
		   nearestNeighbors(inputCenteringThree,
				    outputCenteringThree),
		   /* vector: */ 3, 2,
		   /* FieldOffsetList: */ 4,
		   FieldOffset<3>(Loc<3>(-1,-1,-1), 1), false);
  tester.check("face->edge intercell distances",
	       sameDistances(nearestNeighbors(inputCenteringThree,
					      outputCenteringThree),
			     inputCenteringThree,
			     outputCenteringThree));

  int ret = tester.results("NearestNeighbors");
  Pooma::finalize();
  return ret; 
}

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NearestNeighbors.cpp,v $   $Author: richi $
// $Revision: 1.5 $   $Date: 2004/11/24 10:15:01 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
