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

	double GetGridBoundary(Directions dir) const;
	/*!
	\brief Class constructor
	\param left The left coordinate
	\param right The right coordinate
	\param up The yp coordinate
	\param down The down coordinate
	*/
	PeriodicBox(double left, double right,double up, double down);
private:
	double _left,_right,_up,_down;
};

#endif // PERIODICBOX_HPP
