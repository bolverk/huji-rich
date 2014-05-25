/*! \brief Outer Boundary Conditions
  \author Elad Steinberg
 */

#ifndef OUTERBOUNDARY_HPP
#define OUTERBOUNDARY_HPP 1

#include "../../tessellation/tessellation.hpp"
#include <cmath>

enum Directions {Left, Right, Up, Down};
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

  virtual ~OuterBoundary(void);
};
  

#endif // OUTERBOUNDARY_HPP
