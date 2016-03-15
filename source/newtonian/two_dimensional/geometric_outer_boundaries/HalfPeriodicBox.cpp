#include "HalfPeriodicBox.hpp"
#include "../../../misc/universal_error.hpp"

bool HalfPeriodicBox::PointIsReflective(Vector2D const& point)const
{
	if(point.y<_down||point.y>_up)
		return true;
	return false;
}

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
	if(fabs((edge.vertices.first.y-edge.vertices.second.y))<1e-5*length)
		if(fabs(edge.vertices.first.y-_up)<1e-5*length
			||fabs(edge.vertices.first.y-_down)<1e-5*length)
			return true;
		else
			return false;
	else
		return false;
}
