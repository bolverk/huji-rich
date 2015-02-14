#include "PeriodicBox.hpp"
#include "../../../misc/universal_error.hpp"

bool PeriodicBox::PointIsReflective(Vector2D const& /*point*/)const
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

PeriodicBox::PeriodicBox(const Vector2D& lower_left,
			 const Vector2D& upper_right):
  _left(lower_left.x),
  _right(upper_right.x),
  _up(upper_right.y),
  _down(lower_left.y)
{
  assert((_up>_down) && (_right>_left));
}

std::pair<Vector2D, Vector2D> PeriodicBox::getBoundaries(void) const
{
  return std::pair<Vector2D,Vector2D>(Vector2D(_left,_down),
				      Vector2D(_right,_up));
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
