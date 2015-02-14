#include "shape_2d.hpp"
#include <cmath>

Shape2D::~Shape2D(void) {}

Circle::Circle(Vector2D const& center,
	       double radius):
  center_(center),
  radius_(radius) {}

const Vector2D& Circle::getCenter(void) const
{
  return center_;
}

double Circle::getRadius(void) const
{
  return radius_;
}

bool Circle::operator()(Vector2D const& r) const
{
  return ScalarProd(r-center_,r-center_)<pow(radius_,2);
}

Outside::Outside(Shape2D const& shape):
  shape_(shape) {}

bool Outside::operator()(Vector2D const& r) const
{
  return !shape_(r);
}
