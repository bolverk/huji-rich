#include "polar2cart.hpp"

Vector2D polar2cart(double radius, double angle)
{
  return Vector2D(radius*cos(angle),
		  radius*sin(angle));
}
