#include "zero_force_1d.hpp"

Conserved ZeroForce1D::operator()
(vector<double> const& /*vertices*/,
 vector<Primitive> const& /*cells*/,
 size_t /*point*/,
 double /*t*/,
 double /*dt*/) const
{
  return Conserved();
}
		  
