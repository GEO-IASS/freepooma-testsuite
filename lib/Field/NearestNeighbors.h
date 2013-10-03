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
// Function:
//   nearestNeighbors
// Class:
//   NearestNeighborClass
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_NEAREST_NEIGHBORS_H
#define POOMA_FIELD_NEAREST_NEIGHBORS_H

/** @file
 * @ingroup Field
 * @brief
 * yields FieldOffsets corresponding to the given parameters
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Domain/Loc.h"
#include "Field/FieldCentering.h"
#include "Field/FieldOffset.h"
#include <vector>
#include <numeric>
#include <utility>
#include <functional>
#include <algorithm>


//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------


//----------------------------------------------------------------------
// NearestNeighborClass.
//----------------------------------------------------------------------

/**
 * This class implements the work of nearestNeighbors().
 */

template <int Dim, bool IntraCellOnly = false>
class NearestNeighborClass {
public:

  typedef FieldOffset<Dim> FieldOffset_t;
  typedef FieldOffsetList<Dim> FieldOffsetList_t;
  typedef std::vector<FieldOffset_t> FieldOffset_vt;
  typedef std::vector<FieldOffsetList_t> Answer_t;
  typedef Centering<Dim> Center;
  typedef typename Center::Positions Positions;
  typedef typename Center::Position Position;

  // To compute the set of input values, we maintain a set of input
  // values and their differences from the output value.  In fact, we
  // maintain the input value index and its difference.

  typedef std::pair<int, Position> MinimumPair;
  typedef std::vector<MinimumPair> MinimumSet;


  /// The constructor performs no work.  The function operators do all
  /// the work.

  NearestNeighborClass() {}


  /// Return nearest neighbors for all value in an output centering.

  inline Answer_t
  operator()(const Center &inputCentering, const Center &outputCentering)
  {
    PInsist(inputCentering.size() > 0,
	    "The input centering must be non-empty.");

    Answer_t answer;
    answer.resize(outputCentering.size());
    const Positions inputPositions = inputCentering.positions();
    const Positions outputPositions = outputCentering.positions();

    // Determine nearest neighbors for each output value.

    for (typename Answer_t::size_type outputIndex = 0;
	 outputIndex < outputCentering.size();
	 ++outputIndex)
      answer[outputIndex] = nearestNeighbors(inputPositions,
					     outputPositions[outputIndex]);
    return answer;
  }


  /// Return the nearest neighbors for one output position, specified
  /// by a FieldOffset.

  inline FieldOffsetList_t
  operator()(const Center &inputCentering,
	     const FieldOffset_t &fieldOffset,
	     const Center &outputCentering)
  {
    PInsist(inputCentering.size() > 0,
	    "The input centering must be non-empty.");
    PInsist(fieldOffset.subFieldNumber() < outputCentering.size(),
	    "The field offset must correspond to the output centering.");

    return nearestNeighbors(inputCentering.positions(),
			    outputCentering.position(fieldOffset.subFieldNumber()));
  }


  /// Return the nearest neighbors for multiple output positions, specified
  /// by a FieldOffsetList.

  inline std::vector<FieldOffsetList_t>
  operator()(const Center &inputCentering,
	     const FieldOffsetList_t &fieldOffsetList,
	     const Center &outputCentering)
  {
    PInsist(inputCentering.size() > 0,
	    "The input centering must be non-empty.");

    Answer_t answer;
    answer.resize(fieldOffsetList.size());
    const Positions inputPositions = inputCentering.positions();

    // Determine nearest neighbors for each field offset.

    for (typename FieldOffsetList_t::size_type folIndex = 0;
	 folIndex < outputCentering.size();
	 ++folIndex) {
      PInsist(fieldOffsetList[folIndex].subFieldNumber() < outputCentering.size(),
	    "The field offset must correspond to the output centering.");
      answer[folIndex] =
	nearestNeighbors(inputPositions,
			 outputCentering.position(fieldOffsetList[folIndex].subFieldNumber()));
    }

    return answer;
  }


private:

  /// Given an input centering and a position in logical
  /// coordinate space, return a FieldOffsetList of the nearest
  /// neighbors according to Manhattan distance.

  inline FieldOffsetList_t
  nearestNeighbors(const Positions &inputPositions,
		   const Position outputValue)
  {

    // Compute all input values in the first shell.

    MinimumSet minimumSet;	// all input values in first shell
    typename Positions::size_type inputIndex = 0;

    // Use the first input value to start computing the minimum.

    Position positionDifference = inputPositions[inputIndex] - outputValue;
    double minimumDistance = 
      (IntraCellOnly ?
       manhattanDistance<Manhattan>(positionDifference) :
       manhattanDistance<ManhattanGrid>(positionDifference));
    minimumSet.push_back(std::make_pair(inputIndex, positionDifference));

    // Compute the minimum over the rest of the input values.

    for (++inputIndex;
	 inputIndex < inputPositions.size();
	 ++inputIndex) {
      positionDifference = inputPositions[inputIndex] - outputValue;
      const double distance = 
	(IntraCellOnly ?
	 manhattanDistance<Manhattan>(positionDifference) :
	 manhattanDistance<ManhattanGrid>(positionDifference));
      if (distance < minimumDistance + epsilon) {
	if (distance < minimumDistance) {
	  minimumSet.clear();
	  minimumDistance = distance;
	}
	minimumSet.push_back(std::make_pair(inputIndex,
					    positionDifference));
      }
    }


    // Convert the minimum set to a set of FieldOffsets.
    // minimumSet has all the minimum distance locations.

    FieldOffset_vt answerHolder;
    if (IntraCellOnly) {
      for (typename MinimumSet::size_type minIndex = 0;
	   minIndex < minimumSet.size();
	   ++minIndex)
	answerHolder.push_back(FieldOffset_t(Loc<Dim>(0),
					     minimumSet[minIndex].first));
    }
    else {
      FieldOffset_vt partialAnswer;
      for (typename MinimumSet::size_type minIndex = 0;
	   minIndex < minimumSet.size();
	   ++minIndex)
	{
	  // Compute the cell offsets, appending to the set of answers.

	  partialAnswer = computeCellOffsets(minimumSet[minIndex].first,
					     minimumSet[minIndex].second);
	  answerHolder.insert(answerHolder.end(),
			      partialAnswer.begin(), partialAnswer.end());
	}

      // Remove all duplicates from the answer set.

      std::sort(answerHolder.begin(), answerHolder.end(),
		CompareFieldOffset());
      answerHolder.erase(std::unique(answerHolder.begin(),
				     answerHolder.end(),
				     EqualFieldOffset()),
			 answerHolder.end());
    }

    return answerHolder;
  }

  /// Given a difference between two positions in logical coordinate
  /// space, return the Manhattan norm distance taking into account
  /// that input values are repeated in every grid cell.

  struct ManhattanGrid : public std::binary_function<double, double, double>
  {
    double operator()(const double totalSoFar, double coordinate) const {
      const double absCoordinate = std::abs(coordinate);
      return totalSoFar + std::min(absCoordinate, 1-absCoordinate);
    }
  };

  /// Given a difference between two positions in logical coordinate
  /// space, return the Manhattan norm distance not taking into account
  /// that input values are repeated in every grid cell.

  struct Manhattan : public std::binary_function<double, double, double>
  {
    double operator()(const double totalSoFar, double coordinate) const {
      return totalSoFar + std::abs(coordinate);
    }
  };

  template <class Distance>
  inline static
  double manhattanDistance(const Position &difference)
  {
    double answer = 0.0;;
    for (int coordinate = Dim-1; coordinate >= 0; --coordinate)
      answer = Distance()(answer, difference(coordinate));
    return answer;
  }

  /// Given an input value in the first shell and its position
  /// difference from the given output value, return a vector of
  /// FieldOffsets of input values in the first shell, taking into
  /// account the repetition of input values throughout the grid.  This
  /// is non-trivial because
  /// -# input values are replicated and multiple values may
  ///     be the same distance away
  /// -# the closest input location may be in a different cell

  inline static const FieldOffset_vt
  computeCellOffsets(const int inputValueIndex, const Position &difference)
  {
    // Start with one empty tuple.

    FieldOffset_vt answer(1);
    int numTuples = 1;

    // Store the cell offsets for input values.

    int cellOffsetCoordinates[2];	// For our problem, there can
					// be at most two cell offsets for
					// each dimension.
    int numOffsets;			// number of cell offsets in
					// cellOffsetCoordinates

    for (int dimension = 0; dimension < Dim; ++dimension) {
      numOffsets =
	convertDifferenceToCellOffsets(difference(dimension),
				       cellOffsetCoordinates);
      PInsist(numOffsets >= 1 && numOffsets <= 2,
	      "Incorrect number of cell offsets");

      if (numOffsets == 2)
	// Duplicate the tuples.
	answer.insert(answer.end(), answer.begin(), answer.end());

      for (int coc = 0; coc < numOffsets; ++coc)
	for (int tuple = 0; tuple < numTuples; ++tuple)
	  answer[numTuples * coc + tuple].modifyCellOffset()[dimension] =
	    cellOffsetCoordinates[coc];

      numTuples *= numOffsets;
    }

    // Set the subField numbers.

    for (int i = numTuples-1; i >= 0; --i)
      answer[i].setSubFieldNumber(inputValueIndex);

    return answer;
  }

  /// Given one coordinate of a difference between two coordinates,
  /// return the corresponding cell offset(s), either one or two.

  inline static
  int convertDifferenceToCellOffsets(const double difference,
				     int cellOffsetCoordinate[])
  {
    if (difference < -0.5 - epsilon) {
      cellOffsetCoordinate[0] = +1;
      return 1;
    }
    else if (std::abs(difference + 0.5) < epsilon) {
      cellOffsetCoordinate[0] = +1;
      cellOffsetCoordinate[1] = 0;
      return 2;
    }
    else if (difference <= 0.5 - epsilon) {
      cellOffsetCoordinate[0] = 0;
      return 1;
    }
    else if (std::abs(difference - 0.5) < epsilon) {
      cellOffsetCoordinate[0] = 0;
      cellOffsetCoordinate[1] = -1;
      return 2;
    }
    else if (difference < 1.0 + epsilon) {
      cellOffsetCoordinate[0] = -1;
      return 1;
    }
    else {
      PInsist(0, "Out of range difference");
      return 0;	// Keep the compiler quiet.
    }
  }

  /// Specify a partial order for FieldOffsets to use when removing
  /// duplicates.

  struct CompareFieldOffset :
    public std::binary_function<FieldOffset_t, FieldOffset_t, bool> {
    bool operator()(const FieldOffset_t &op1, const FieldOffset_t &op2) {
      return (op1.cellOffset() < op2.cellOffset()) ||
	(op1.cellOffset() == op2.cellOffset() &&
	 op1.subFieldNumber() < op2.subFieldNumber());
    }
  };
  
  /// Specify an equality operator for FieldOffsets to use when removing
  /// duplicates.

  struct EqualFieldOffset :
    public std::binary_function<FieldOffset_t, FieldOffset_t, bool> {
    bool operator()(const FieldOffset_t &op1, const FieldOffset_t &op2) {
      return op1.cellOffset() == op2.cellOffset() &&
	op1.subFieldNumber() == op2.subFieldNumber();
    }
  };
  
  /// Use epsilon when comparing floating-point numbers, which cannot
  /// be represented precisely.

  static const double epsilon;

};


template <int Dim, bool IntraCellOnly>
const double
NearestNeighborClass<Dim, IntraCellOnly>::epsilon = 1.0e-08;


//-----------------------------------------------------------------------------
// nearestNeighbors.
//-----------------------------------------------------------------------------

///@name nearestNeighbors()
/// Given input and output centerings, this function computes the
/// "first shell" of nearest neighbors for each output value.  That is,
/// for each output value, it computes the FieldOffsetList containing
/// the input values that are closest, with respect to the Manhattan
/// norm (l_1 norm).  The FieldOffsetLists are stored in a std::vector
/// in the same order as the output values occur in the output
/// centering.
///
/// If only values within the corresponding input cell are desired,
/// specify a third parameter.

//@{

template <int Dim>
inline
std::vector<FieldOffsetList<Dim> >
nearestNeighbors(const Centering<Dim> &inputCentering,
		 const Centering<Dim> &outputCentering)
{
  return NearestNeighborClass<Dim>()(inputCentering, outputCentering);
}

template <int Dim>
inline
std::vector<FieldOffsetList<Dim> >
nearestNeighbors(const Centering<Dim> &inputCentering,
		 const Centering<Dim> &outputCentering,
		 const bool)
{
  return NearestNeighborClass<Dim, true>()(inputCentering, outputCentering);
}

template <int Dim>
inline
std::vector<FieldOffsetList<Dim> >
nearestNeighbors(const Centering<Dim> &inputCentering,
		 const FieldOffsetList<Dim> &fOL,
		 const Centering<Dim> &outputCentering)
{
  return NearestNeighborClass<Dim>()(inputCentering, fOL, outputCentering);
}

template <int Dim>
inline
std::vector<FieldOffsetList<Dim> >
nearestNeighbors(const Centering<Dim> &inputCentering,
		 const FieldOffsetList<Dim> &fOL,
		 const Centering<Dim> &outputCentering,
		 const bool)
{
  return NearestNeighborClass<Dim, true>()
    (inputCentering, fOL, outputCentering);
}

template <int Dim>
inline
FieldOffsetList<Dim>
nearestNeighbors(const Centering<Dim> &inputCentering,
		 const FieldOffset<Dim> &fieldOffset,
		 const Centering<Dim> &outputCentering)
{
  return NearestNeighborClass<Dim>()
    (inputCentering, fieldOffset, outputCentering);
}

template <int Dim>
inline
FieldOffsetList<Dim>
nearestNeighbors(const Centering<Dim> &inputCentering,
		 const FieldOffset<Dim> &fieldOffset,
		 const Centering<Dim> &outputCentering,
		 const bool)
{
  return NearestNeighborClass<Dim, true>()
    (inputCentering, fieldOffset, outputCentering);
}

//@}

/// Given an input centering and a field offset obtained from a nearest
/// neighbor calculation, this function computes the position of the
/// corresponding point in cell logical coordinates.

template<int Dim>
inline Vector<Dim, double>
inputPosition(const Centering<Dim> &inputCentering,
              const FieldOffset<Dim> &fieldOffset)
{
  Vector<Dim, double> ret =
    inputCentering.position(fieldOffset.subFieldNumber());
  for (int i = 0; i < Dim; ++i)
    ret(i) += fieldOffset.cellOffset()[i].first();
  return ret;
}

#endif // POOMA_FIELD_NEAREST_NEIGHBORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: NearestNeighbors.h,v $   $Author: richard $
// $Revision: 1.5 $   $Date: 2004/11/01 18:16:43 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
