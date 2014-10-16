#include "piecewise.hpp"

Piecewise::Piecewise(const Shape2D& shape,
		     const SpatialDistribution& inside,
		     const SpatialDistribution& outside):
  shape_(shape), inside_(inside), outside_(outside) {}

double Piecewise::operator()(const Vector2D& point) const
{
  return shape_(point) ? inside_(point) : outside_(point);
}
