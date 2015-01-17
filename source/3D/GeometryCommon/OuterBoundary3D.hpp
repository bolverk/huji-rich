/*! \file OuterBoundary3D.hpp
  \brief Outer Boundary Conditions
  \author Elad Steinberg
 */

#ifndef OUTERBOUNDARY3D_HPP
#define OUTERBOUNDARY3D_HPP 1

#include <cmath>

//! \brief Type of boundary
enum BoundaryType{Rectengular, Periodic,HalfPeriodic};

//! \brief Type of boundary point
enum BoundaryPoint{BackLowerLeft,FrontUpperRight};


//! \brief Abstract class for geometric boundary conditions for the tessellation
class OuterBoundary3D
{
public:
	/*!
	\brief Returns the boundary type
	\return The boundary type
	*/
  virtual BoundaryType GetBoundaryType(void) const = 0;
  /*!
	\brief Returns the boundary point
	\param point The point of the boundary
	\return The boundary point
	*/
  virtual Vector3D const& GetGridBoundary(BoundaryPoint point) const = 0;

	/*!
	\brief Return whether an face is reflective or not
	\param face The Face to check
	\returns Is the Face reflective
	*/
  virtual bool AreWeReflective(Face const& face)const=0;
  //! \brief Virtual destructor
   virtual ~OuterBoundary3D(void);
};

#endif // OUTERBOUNDARY3D_HPP
