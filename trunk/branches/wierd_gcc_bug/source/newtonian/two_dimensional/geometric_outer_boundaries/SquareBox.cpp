#include "SquareBox.hpp"
#include "../../../misc/universal_error.hpp"

bool SquareBox::AreWeReflective(Edge const& /*edge*/)const
{
	return true;
}

SquareBox::SquareBox(double left, double right,double up, double down)
	:_left(left),_right(right),_up(up),_down(down)
{
	if(left>=right||down>=up)
	  throw UniversalError("Invalid values for grid boundaries");
}

BoundaryType SquareBox::GetBoundaryType(void) const
{
	return Rectengular;
}

double SquareBox::GetGridBoundary(Directions dir) const
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
