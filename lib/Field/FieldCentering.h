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
//   Centering
//-----------------------------------------------------------------------------

#ifndef POOMA_FIELD_CENTERING_H
#define POOMA_FIELD_CENTERING_H


/** @file
 * @ingroup Field
 * @brief
 * specifies value locations within a field's cell
 *
 * Centering specifies value locations within a field's cell.
 *
 * CanonicalCentering yields some canonical centerings.
 *
 * canonicalCentering<Dim>(type, discontinuous, dimension)
 * yields the specified canonical centering
 */

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "Tiny/Vector.h"
#include "Utilities/PAssert.h"
#include <vector>
#include <ostream>

//-----------------------------------------------------------------------------
// Forward declarations:
//-----------------------------------------------------------------------------

template <int Dim>
class Centering;

template <int Dim>
class CanonicalCentering;



//---------------------------------------------------------------------------
// Global enumerations.

/// Indicate a cell's centering type.
enum CenteringType {
	VertexType, 	///< Vertex type
	EdgeType, 	///< Edge type
	FaceType, 	///< Face type
	CellType 	///< Cell type
};

/// Indicate whether a cell's values are shared xor not shared.
enum ContinuityType {
	Continuous = 0, ///< Shared xor
	Discontinuous 	///< Not shared
};

/// Permit specifying various dimensions.
enum {
	XDim = 1, 			///< X components
	YDim = XDim << 1, 		///< Y components
	ZDim = YDim << 1, 		///< Z components
	AllDim = XDim | YDim | ZDim 	///< All components
}; 


/**
 * Centering indicates the positions of values within a field's cell.
 * It is specified using the following fields:
 * - centering type: an enumeration for vertex, edge, face, or cell
 *   centering type
 * - discontinuous: a boolean indicating whether, for a value located
 *   on a cell boundary, each of the neighboring cells has its own value
 *   (discontinuous) xor one value for all neighboring cells is provided
 *   (continuous)
 * - list of values: Each list element is a pair of an orientation and
 *   a position.  The position, a Vector<dim>, specifies the value's
 *   position with respect to the cell's logical coordinate system,
 *   which is either [0.0,1.0)^{dim} or [0.0,1.0]^{dim} depending on whether
 *   values are continuous or discontinuous, respectively.  The
 *   orientation in Z^{dim}_2, represented using a Loc<dim>, indicates
 *   which zeroes (or ones if discontinuous) in the position must be
 *   zero (or one) because of the centering type.  For example, a
 *   continuous face centering for an x-face must have a 0 in the
 *   x-component.  Other coordinates can be zero but need not be.
 *
 *   In practice, we use two lists of values, one for orientations and
 *   one for positions, with elements at the same position related to
 *   each other.
 *
 * For a dim-D cell, a face is a (dim-1)-D object.  An edge
 * centering is always a 1-D object.  For two dimensions, edge and
 * face centerings are the same although it is sometimes useful to
 * distinguish them when writing programs that work for various
 * dimensions.
 * 
 * Adjacent cells can share values.  For example, a vertex-centered value
 * in a 3-D field is shared by eight cells.  A program might require that
 * each cell maintain its own value at the point.  To do so, specify
 * discontinuous values.  To be a valid discontinuous centering, values
 * on cell boundaries must be arranged so that every adjacent cell also
 * has a value at that same position in space.  For example, any
 * discontinuous edge value in three-dimensional space must have related
 * values in the three adjacent cells.  (Why do we have this restriction?)
 * 
 * Each value in the list should be specified exactly once.  For example,
 * the canonical continuous vertex-centered value is at position (0.0,
 * ..., 0.0).  This cell specifies a value at its origin.  Neighboring
 * cells specify values at their origins.  Collectively, all the field's
 * vertex values are specified.  For continuous values, positions should
 * be in the range [0.0,1.0)^{dim}.  For discontinuous values, positions
 * should be in the range [0.0,1.0]^{dim}.  To implement the field, each
 * specified position corresponds to a subfield.
 * 
 * Orientations are used when creating storage for field values.  When
 * creating storage, the number of values does not necessarily match
 * the number of cells.  For example, to create a field with n^2 cells
 * requires (n+1)^2 vertices.  Although orientations might appear to
 * be redundant since nonzero values indicate an orientation, an
 * x-face centered value may have position (0,0,0.5), which does not
 * indicate whether it is an x- or y-face.
 * 
 * For three-dimensions, example orientations include
 * - vertex type: 000 since vertex positions are completely determined
 * - x-edge: 100 since the x-coordinate varies along the edge
 * - x-face: 011 since the x-coordinate is fixed at zero (or one) but the other
 *         values can vary
 * - cell type: 111 since values can be placed at any location.
 *
 * Example centerings include:
 * 
 * - continuous vertex centering: Values occur at every vertex in the
 *   field, but each cell is responsible only for the vertex located at its
 *   origin.  Thus, the centering need only declare one
 *   orientation-position pair:
 *   <PRE>VertexType,
 *   false / * continuous * /,
 *   {((0, 0, 0), (0.0, 0.0, 0.0))}</PRE>
 * - discontinuous vertex centering: Each cell maintains its own set of
 *   2^dim vertex values:
 *   <PRE>VertexType,
 *   true / * discontinuous * /,
 *   {((0, 0, 0), (0.0, 0.0, 0.0)),
 *    ((0, 0, 0), (1.0, 0.0, 0.0)),
 *    ((0, 0, 0), (0.0, 1.0, 0.0)),
 *    ((0, 0, 0), (0.0, 0.0, 1.0)),
 *    ((0, 0, 0), (1.0, 1.0, 0.0)),
 *    ((0, 0, 0), (1.0, 0.0, 1.0)),
 *    ((0, 0, 0), (0.0, 1.0, 1.0)),
 *    ((0, 0, 0), (1.0, 1.0, 1.0))}</PRE>
 * - continuous face centering: Since the faces are shared, each cell is
 *   responsible for only the faces intersecting the origin.  Here we
 *   specify one value on each face.
 *   <PRE>FaceType,
 *   false / * continuous * /,
 *   {((0, 1, 1), (0.0, 0.5, 0.5)),
 *    ((1, 0, 1), (0.5, 0.0, 0.5)),
 *    ((1, 1, 0), (0.5, 0.5, 0.0))}</PRE>
 */

template <int Dim>
class Centering 
{
public:

  //---------------------------------------------------------------------------
  // Exported typedefs and enumerations.
    
  /// An orientation.
  typedef Loc<Dim> Orientation;

  /// A position.
  typedef Vector<Dim> Position;

  /// A list of value orientations.
  typedef std::vector<Orientation> Orientations;

  /// A list of value positions. 
  typedef std::vector<Position> Positions;
  
  //---------------------------------------------------------------------------
  //@name User-callable constructors. These ctors are meant to be called by users.
  //@{

  /// Create a centering without any values.

  Centering(CenteringType cent = CellType, ContinuityType cont = Continuous)
    : centering_type_m(cent), discontinuous_m(cont),
      orientations_m(0), positions_m(0)
    { }

  /// Create a centering with values specified in two vectors.  The two
  /// vectors should have the same length and corresponding entries
  /// specifying values.

  Centering(CenteringType cent, ContinuityType cont,
	    const Orientations &orientations, const Positions &positions)
    : centering_type_m(cent), discontinuous_m(cont),
      orientations_m(orientations), positions_m(positions)
  {
#if POOMA_BOUNDS_CHECK
    PInsist(orientations_m.size() == positions_m.size(),
	    "Centering differing vector length error.");
#endif
      return;
  }

  /// Sub-centering constructor.
  Centering(const Centering<Dim>& model, int c)
    : centering_type_m(model.centering_type_m),
      discontinuous_m(model.discontinuous_m),
      orientations_m(1, model.orientations_m[c]),
      positions_m(1, model.positions_m[c])
    { }

  //@}

  //---------------------------------------------------------------------------
  /// Sub-centering creation function.  This is not really meant to be
  /// called by users.
  /// Return a centering with one specified value.
  Centering<Dim> operator[](int iSubField) const
  {
#if POOMA_BOUNDS_CHECK
    PInsist(iSubField >= 0 && iSubField < size(),
	    "Illegal attempt to extract a non-existent centering.");
#endif
    return Centering<Dim>(*this, iSubField);
  }


  //---------------------------------------------------------------------------
  /// Empty destructor is fine for us.
  
  ~Centering() { }

  //---------------------------------------------------------------------------
  //@name Accessors.
  //@{

  inline CenteringType centeringType() const
    {
      return centering_type_m;
    }

  inline ContinuityType continuityType() const
    {
      return discontinuous_m;
    }

  inline bool discontinuous() const
  {
    return discontinuous_m == Discontinuous;
  }

  inline bool continuous() const
  {
    return discontinuous_m == Continuous;
  }

  inline const Orientations &orientations() const
  {
    return orientations_m;
  }

  inline const Positions &positions() const
  {
    return positions_m;
  }

  inline const Orientation &orientation(int i) const
  {
    return orientations_m[i];
  }

  inline const Position &position(int i) const
  {
    return positions_m[i];
  }

  inline int size() const
  {
#if POOMA_BOUNDS_CHECK
    PInsist(orientations_m.size() == positions_m.size(),
	    "In a centering, the number of orientations must match the number of positions.");
#endif
    return orientations_m.size();
  }

  //@}

  //---------------------------------------------------------------------------
  /// Add a value to a centering.  There is no check that the value is
  /// not already present.

  inline void
  addValue(const Orientation &orientation,
	   const Position &position)
    {
      orientations_m.push_back(orientation);
      positions_m.push_back(position);
      return;
    }

private:

  /// The type of cell centering.
  CenteringType centering_type_m;

  /// Discontinuous (or true or 1) indicates that values on cell
  /// boundaries are not shared with neighboring cells.  Continuous (or
  /// false or 0) indicates that values on cell boundaries are shared
  /// with neighboring cells.

  ContinuityType discontinuous_m;

  /// The list of value orientations.

  Orientations orientations_m;

  /// The list of value positions.
  /// This must have the same length as orientations_m.

  Positions positions_m;
};


/**
 * This object makes available some canonical centerings.  By calling with
 * -# a centering type, e.g., CellType or VertexType,
 * -# whether the centering should be discontinuous or not,
 * -# a dimension \in [0,Dim),
 *
 * the corresponding centering is returned.  Some parameters do not
 * make sense for some centerings.
 *
 * The canonical centerings include:
 * - (CellType, / * ignored * /, / * ignored * /) :
 *    a cell centering with one value at the cell's center
 * - (VertexType, true/false, / * ignored * /) :
 *    a vertex centering with values at all cell vertices
 * - (EdgeType, true/false, dimension) :
 *    an edge centering with values on the specified edges
 * - (FaceType, true/false, dimension) :
 *    a face centering with values on the specified faces

 * The dimension field should be the bitwise-or of X, Y, Z, and All,
 * where the All value equals X | Y | Z.  For example, using `Y | Z'
 * yields the edges along the y- and z-axes or yields the y- and
 * z-faces.
 */

template <int Dim>
class CanonicalCentering {
public:
  
  //---------------------------------------------------------------------------
  /// Internal POOMA constructors. These ctors are used internally by POOMA.
  /// They are not really meant to be called by users.

  CanonicalCentering();

  //---------------------------------------------------------------------------
  /// Deallocate centering_table_m.

  ~CanonicalCentering () {
    if (--class_count_m == 0) {
      for (int i = 0; i <= CellType; ++i) {
	for (int j = 0; j < 2; ++j)
	  delete [] centering_table_m[i][j];
	delete [] centering_table_m[i];
      }
      delete [] centering_table_m;
    }
  }

  //---------------------------------------------------------------------------
  /// Return the desired centering.

  inline Centering<Dim> operator()
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     int dimension = 0) const
  {
    if (dimension == 0)
      dimension = AllDim;
    // Permit the user to specify higher dimension values but return
    // the ones appropriate for this dimension.
    dimension %= (1<<Dim);
    return centering_table_m[type][discontinuous][dimension];
  }

private:

  /// Add the given orientation and position to the corresponding
  /// vectors.
  inline static void addValue(typename Centering<Dim>::Orientations &os,
			      typename Centering<Dim>::Positions &pos,
			      const typename Centering<Dim>::Orientation &o,
			      const typename Centering<Dim>::Position &p)
  {
    os.push_back(o);
    pos.push_back(p);
    return;
  }

  /// Return a container Orientations containing all the values in the operands.
  template <class T>
  inline static
  T combine(const T &op1, const T &op2)
  {
    T answer(op1);
    answer.insert(answer.end(), op2.begin(), op2.end());
    return answer;
  }

  /// Number of extent objects.

  static int class_count_m;

  /// Table containing the centerings, which are the cross product of
  /// the centering type, discontinuous xor continuous, all the
  /// possible combinations of dimensions.  We ignore illegal entries
  /// such as centering_table_m[x][y][0].

  static Centering<Dim>*** centering_table_m;
};

//-----------------------------------------------------------------------------
// Overload the << operator to print a Centering to a stream.
//-----------------------------------------------------------------------------

template <int Dim>
std::ostream &operator<<(std::ostream &o, 
                         const Centering<Dim> &centering)
{
  switch (centering.centeringType())
  {
  case VertexType:
    o << "Vertex";
    break;
  case EdgeType:
    o << "Edge";
    break;
  case FaceType:
    o << "Face";
    break;
  case CellType:
    o << "Cell";
    break;
  }

  o << "," << (centering.continuous() ? "Continuous" : "Discontinuous")
    << ",{";
  for (int i = 0; i < centering.size();)
  {
    o << "[" << centering.orientation(i)
      << "," << centering.position(i) << "]";
    ++i;
    if (i < centering.size())
      o << ",";
  }
  o << "}";

  return o;
}


//-----------------------------------------------------------------------------
// Provide == and != comparison operators for centerings.
//-----------------------------------------------------------------------------

template <int Dim>
bool operator==(const Centering<Dim> &centering1, const Centering<Dim> &centering2) 
{
  return
    centering1.centeringType() == centering2.centeringType() &&
    centering1.discontinuous() == centering2.discontinuous() &&
    centering1.orientations() == centering2.orientations() &&
    centering1.positions() == centering2.positions();
}

template <int Dim>
bool operator!=(const Centering<Dim> &centering1, const Centering<Dim> &centering2) 
{
  return !(centering1 == centering2);
}


//-----------------------------------------------------------------------------
// Define CanonicalCentering's static members.
//-----------------------------------------------------------------------------

/// Number of extent objects.

template <int Dim>
int CanonicalCentering<Dim>::class_count_m = 0;

/// Table containing the centerings, which are the cross product of
/// the centering type, discontinuous xor continuous, all the
/// possible combinations of dimensions.  We ignore illegal entries
/// such as centering_table_m[x][y][0].

template <int Dim>
Centering<Dim>*** CanonicalCentering<Dim>::centering_table_m = 0;


//-----------------------------------------------------------------------------
// canonicalCentering.
//-----------------------------------------------------------------------------

///@name canonicalCentering
///
/// canonicalCentering<Dim>(type, discontinuous, dimension) is a
/// functional wrapper around a CanonicalCentering object.
//@{

// The canonical centering objects
extern const CanonicalCentering<1> canonicalCenteringOne_g;
extern const CanonicalCentering<2> canonicalCenteringTwo_g;
extern const CanonicalCentering<3> canonicalCenteringThree_g;

#if POOMA_NO_TEMPLATEFUNC_DEFAULTARGS
template <int Dim>
const Centering<Dim> canonicalCentering
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension);
#else
template <int Dim>
const Centering<Dim> canonicalCentering
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension = 0);
#endif

template <>
const Centering<1> canonicalCentering<1>
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension);

template <>
const Centering<2> canonicalCentering<2>
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension);

template <>
const Centering<3> canonicalCentering<3>
    (const enum CenteringType type,
     const enum ContinuityType discontinuous,
     const int dimension);

#if POOMA_NO_TEMPLATEFUNC_DEFAULTARGS
template <int Dim>
inline const Centering<Dim> canonicalCentering
    (const enum CenteringType type,
     const enum ContinuityType discontinuous)
{
	return canonicalCentering<Dim>(type, discontinuous, 0);
}
#endif
//@}


///@name Functions for translating domains based on centerings.
//@{

/// cellDomainToCenteringDomain(cellDom, centering, i)
///   computes the domain of the i'th subfield for a field that has the given
///   cell domain.

template<int Dim>
inline Interval<Dim>
cellDomainToCenteringDomain(const Interval<Dim> &cellDom,
                            const Centering<Dim> &centering, int i)
{
  if (centering.discontinuous())
  {
    return cellDom;
  }
  else
  {
    return shrinkRight(growRight(cellDom, 1), centering.orientation(i));
  }
}
                                          
/// centeringDomainToCellDomain(cDom, centering, i)
/// the inverse function.

template<int Dim>
inline Interval<Dim>
centeringDomainToCellDomain(const Interval<Dim> &cDom,
                            const Centering<Dim> &centering, int i)
{
  if (centering.discontinuous())
  {
    return cDom;
  }
  else
  {
    return shrinkRight(growRight(cDom, centering.orientation(i)), 1);
  }
}

//@}

#endif // POOMA_FIELD_CENTERING_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FieldCentering.h,v $   $Author: richard $
// $Revision: 1.10 $   $Date: 2004/11/01 18:16:42 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
