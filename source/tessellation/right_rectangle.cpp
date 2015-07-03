#include "right_rectangle.hpp"

RightRectangle::RightRectangle(const Vector2D& lower_left,
			       const Vector2D& upper_right):
  lower_left_(lower_left), upper_right_(upper_right) 
{
  assert(upper_right.x>lower_left.x &&
	 upper_right.y>lower_left.y &&
	 "vertices are reversed");
}

RightRectangle::RightRectangle(const pair<Vector2D,Vector2D>& ll_ur):
  lower_left_(ll_ur.first), upper_right_(ll_ur.second)
{
  assert(upper_right_.x>lower_left_.x &&
	 upper_right_.y>lower_left_.y &&
	 "vertices are reversed");
}

namespace
{
  bool is_between(double x_cand,
		  double x_low,
		  double x_high)
  {
    return (x_high>x_cand)&&(x_cand>x_low);
  }
}

bool RightRectangle::operator()(const Vector2D& point) const
{
  return is_between(point.x,lower_left_.x,upper_right_.x)&&
    is_between(point.y,lower_left_.y,upper_right_.y);
}
