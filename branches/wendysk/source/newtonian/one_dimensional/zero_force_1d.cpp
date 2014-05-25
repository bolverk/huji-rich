#include "zero_force_1d.hpp"

Conserved ZeroForce1D::calc
(vector<double> const& /*vertices*/,
 vector<Primitive> const& /*cells*/,
 int /*point*/,
 double /*t*/,
 double /*dt*/) const
{
  return Conserved();
}
		  
