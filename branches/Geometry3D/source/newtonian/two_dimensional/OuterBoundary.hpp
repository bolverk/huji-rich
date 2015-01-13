/*! \file OuterBoundary.hpp
  \brief Outer Boundary Conditions
  \author Elad Steinberg
 */

#ifndef OUTERBOUNDARY_HPP
#define OUTERBOUNDARY_HPP 1

#include "../../tessellation/tessellation.hpp"
#include <cmath>

//! \brief Directions of boundaries of the computational domain
enum Directions {Left, Right, Up, Down};

//! \brief Type of boundary
enum BoundaryType{Rectengular, Periodic,HalfPeriodic};

//! \brief Abstract class for geometric boundary conditions for the tessellation
class OuterBoundary
{
public:
	/*!
	\brief Returns the boundary type
	\return The boundary type
	*/
  virtual BoundaryType GetBoundaryType(void) const = 0;
  /*!
	\brief Returns the boundary coordinate
	\param dir The direction of the boundary
	\return The boundary coordinate
	*/
  virtual double GetGridBoundary(Directions dir) const = 0;

	/*!
	\brief Return wheter an edge is reflective or not
	\param edge The edge to check
	\returns Is the edge reflective
	*/
  virtual bool AreWeReflective(Edge const& edge)const=0;

  /*!
  \brief Checks if the point is a reflected point outside the domain
  \param point The point to check
  \return Whether the point is reflected or not
  */
  virtual bool PointIsReflective(Vector2D const& point)const=0;

  virtual ~OuterBoundary(void);
  /*!
  \brief Returns the outer box as a set of edges in the order: Right, Up, Left and Down. All neighbors are set to zero.
  \return The edges of the boundary box.
  */
  vector<Edge> GetBoxEdges(void) const;
};

#endif // OUTERBOUNDARY_HPP
