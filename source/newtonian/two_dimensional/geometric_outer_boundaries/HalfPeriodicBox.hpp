/*! \file HalfPeriodicBox.hpp
\brief Square box outer boundary conditions with two sides reflective and two periodic. The x direction is taken to be periodic.
\author Elad Steinberg
*/

#ifndef HALFPERIODICBOX_HPP
#define HALFPERIODICBOX_HPP 1

#include "../OuterBoundary.hpp"

//! \brief Square box outer boundary conditions with two sides reflective and two periodic. The x direction is taken to be periodic.
class HalfPeriodicBox: public OuterBoundary
{
public:
	bool AreWeReflective(Edge const& edge)const override;

	BoundaryType GetBoundaryType(void) const override;

	bool PointIsReflective(Vector2D const& point)const override;

	double GetGridBoundary(Directions dir) const override;
	/*!
	\brief Class constructor
	\param left The left coordinate
	\param right The right coordinate
	\param up The up coordinate
	\param down The down coordinate
	*/
	HalfPeriodicBox(double left, double right,double up, double down);

	/*! \brief Returns the lower left and upper right corners
	\return std::pair. first is lower left corner and second is upper right corner
	*/
	std::pair<Vector2D, Vector2D> getBoundaries(void) const;

private:
	double _left,_right,_up,_down;
};

#endif // HALFPERIODICBOX_HPP
