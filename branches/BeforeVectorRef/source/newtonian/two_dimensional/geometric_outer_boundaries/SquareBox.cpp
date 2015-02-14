#include "SquareBox.hpp"
#include "../../../misc/universal_error.hpp"

bool SquareBox::PointIsReflective(Vector2D const& point)const
{
	if(point.x<_left||point.x>_right)
		return true;
	if(point.y<_down||point.y>_up)
		return true;
	return false;
}

bool SquareBox::AreWeReflective(Edge const& /*edge*/)const
{
  return true;
}

namespace {
  void assert_directions(double left, double right,
			 double up, double down)
  {
    if(left>=right||down>=up)
      throw UniversalError("Invalid values for grid boundaries");
  }
}

SquareBox::SquareBox(double left, double right,double up, double down):
  _left(left),_right(right),_up(up),_down(down)
{
  assert_directions(left,right,up,down);
}

SquareBox::SquareBox(Vector2D const& bottom_left,
		     Vector2D const& top_right):
  _left(bottom_left.x),
  _right(top_right.x),
  _up(top_right.y),
  _down(bottom_left.y)
{
  assert_directions(_left,_right,_up,_down);
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

pair<Vector2D, Vector2D> SquareBox::getBoundary(void) const
{
  return pair<Vector2D,Vector2D>(Vector2D(_left,_down),
				 Vector2D(_right,_up));
}

SquareBox::~SquareBox(void) {}
