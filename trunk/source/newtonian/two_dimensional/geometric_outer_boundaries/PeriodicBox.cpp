#include "PeriodicBox.hpp"
#include "../../../misc/universal_error.hpp"

bool PeriodicBox::PointIsReflective(Vector2D const& point)const
{
	return false;
}

bool PeriodicBox::AreWeReflective(Edge const& /*edge*/)const
{
	return false;
}

PeriodicBox::PeriodicBox(double left, double right,double up, double down)
	:_left(left),_right(right),_up(up),_down(down)
{
	if(left>=right||down>=up)
	  throw UniversalError("Invalid values for grid boundaries");
}

BoundaryType PeriodicBox::GetBoundaryType(void) const
{
	return Periodic;
}

double PeriodicBox::GetGridBoundary(Directions dir) const
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