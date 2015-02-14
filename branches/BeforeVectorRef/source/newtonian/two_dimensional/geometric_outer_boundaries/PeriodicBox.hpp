/*! \file PeriodicBox.hpp
\brief Square box outer boundary conditions
\author Elad Steinberg
*/

#ifndef PERIODICBOX_HPP
#define PERIODICBOX_HPP 1

#include "../OuterBoundary.hpp"

//! \brief Square box outer boundary conditions
class PeriodicBox: public OuterBoundary
{
public:
	bool AreWeReflective(Edge const& edge)const;

	BoundaryType GetBoundaryType(void) const;

	bool PointIsReflective(Vector2D const& point)const;

	double GetGridBoundary(Directions dir) const;
	/*!
	\brief Class constructor
	\param left The left coordinate
	\param right The right coordinate
	\param up The yp coordinate
	\param down The down coordinate
	*/
	PeriodicBox(double left, double right,double up, double down);

  /*! \brief Class constructor
    \param lower_left Lower left corner
    \param upper_right Upper right corner
   */
  PeriodicBox(const Vector2D& lower_left,
	      const Vector2D& upper_right);

  /*! \brief Returns the lower left and upper right corners
    \return std::pair. first is lower left corner and second is upper right corner 
   */
  std::pair<Vector2D, Vector2D> getBoundaries(void) const;

private:
  const double _left,_right,_up,_down;
};

#endif // PERIODICBOX_HPP
