#include "uniform2d.hpp"

Uniform2D::Uniform2D(double val): _val(val) {}

double Uniform2D::operator()(Vector2D const& /*r*/) const
{
  return _val;
}
