/*! \file PolyIntersect.hpp
  \brief Finds the intersection points between two polygons
  \author Elad Steinberg
 */

#include "geotests.hpp"
#include "Edge.hpp"

//! \brief Flags for if the segment is inner or outer
enum InFlags {UnKnown,Pi,Qi};
//! \brief Flags for there is an intersction or not or if the segments are parallel
enum IntersectFlags {True,False,Par};

/*!
\brief Calculates the intersection between two convex polygons
\param poly0 The first polygon
\param poly1 The second polygon
\return The vertices of the intersecting polygon, can be empty
*/
vector<Vector2D> ConvexIntersect(vector<Vector2D> const& poly0,vector<Vector2D>
	const& poly1);
/*!
\brief Checks if two segments intersect
\param p0 The first vertice of the first segment
\param p1 The second vertice of the first segment
\param q0 The first vertice of the second segment
\param q1 The second vertice of the second segment
\param Intersection The location of the intersection (if there is) given as output
\return A flag stating if thereis an intersection no intersection or the two segements parallel
*/
IntersectFlags SegmentIntersection(Vector2D const& p0,Vector2D const& p1,
	Vector2D const& q0,Vector2D const& q1,Vector2D &Intersection);

/*
vector<Vector2D> GetParEdge(Vector2D const& p0,Vector2D const& p1,
	Vector2D const& q0,Vector2D const& q1);
	*/
