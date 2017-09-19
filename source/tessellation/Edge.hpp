/*! \file Edge.hpp
  \brief Edge between cells
  \author Elad Steinberg
*/

#ifndef EDGE_HPP
#define EDGE_HPP

#include <vector>
#include "geometry.hpp"

//! \brief Interface between two cells
class Edge
{
public:

  //! \brief Points at the ends of the edge
  std::pair<Vector2D, Vector2D> vertices;

  //! \brief Neighboring cells
  std::pair<int, int> neighbors;

  /*! \brief Class constructor
    \param p1 Position of first vertex
    \param p2 Position of second vertex
    \param neighbor1 Index of first neighbor cell
    \param neighbor2 Index of second neighbor cell
  */
  Edge(Vector2D const& p1, Vector2D const& p2,
       int neighbor1, int neighbor2);

  Edge& operator=(const Edge& other);

  Edge(void);

  ~Edge(void);

  /*! \brief Copy constructor
    \param other Source edge
  */
  Edge(Edge const& other);

  /*! \brief Returns the length of the edge
    \return Length
  */
  double GetLength(void) const;
};

/*! \brief Calculates a unit vector parallel to an edge
  \param edge Edge
  \return Unit vector parallel to the edge
*/
Vector2D Parallel(Edge const& edge);

/*!
  \brief Calculates the distance of a point to an edge
  \param point The point
  \param edge The edge
  \returns The minimum distance between the point and the edge
*/
double DistanceToEdge(Vector2D const& point,Edge const& edge);

/*!
  \brief Calculates the intersection of two edges
  \param edge1 The first edge
  \param edge2 The second edge
  \param Intersection Gets the value of the intersection point
  \param eps The relative tolerance in units of the smaller edge
  \returns true if the intersection is on the edges and false if it is outside them
*/
bool SegmentIntersection(Edge const& edge1,Edge const& edge2,
			 Vector2D &Intersection,double eps=1e-8);

/*! \brief Calculates the centroid of an edge
  \param edge An edge
  \return Centroid
 */
Vector2D calc_centroid(const Edge& edge);

#endif	// EDGE_HPP
