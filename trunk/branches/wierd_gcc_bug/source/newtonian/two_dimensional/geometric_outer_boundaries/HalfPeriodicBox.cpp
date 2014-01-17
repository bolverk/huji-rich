#include "HalfPeriodicBox.hpp"
#include "../../../misc/universal_error.hpp"

HalfPeriodicBox::HalfPeriodicBox(double left, double right,double up, double down)
	:_left(left),_right(right),_up(up),_down(down)
{
	if(left>=right||down>=up)
	  throw UniversalError("Invalid values for grid boundaries");
}

BoundaryType HalfPeriodicBox::GetBoundaryType(void) const
{
	return HalfPeriodic;
}

double HalfPeriodicBox::GetGridBoundary(Directions dir) const
{
	if(dir==Left)
		return _left;
	else if(dir==Right)
		return _right;
	else if(dir==Up)
		return _up;
	else if(dir==Down)
		return _down;
	else
	  throw UniversalError("Unknown direction");
}

bool HalfPeriodicBox::AreWeReflective(Edge const& edge) const
{
	double length=edge.GetLength();
	// Are we periodic or reflective?
	if(abs((edge.get_y(0)-edge.get_y(1)))<1e-5*length)
		if(abs(edge.get_y(0)-_up)<1e-5*length
			||abs(edge.get_y(0)-_down)<1e-5*length)
			return true;
		else
			return false;
	else
		return false;
}

