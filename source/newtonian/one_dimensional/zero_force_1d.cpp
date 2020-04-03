#include "zero_force_1d.hpp"

Conserved ZeroForce1D::operator()
  (const SimulationState1D& /*state*/,
   size_t /*point*/,
   double /*t*/,
   double /*dt*/) const
{
  return Conserved();
}
		  
