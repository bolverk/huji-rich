#include "zero_force_1d.hpp"

Extensive ZeroForce1D::operator()
  (const SimulationState1D& state,
   size_t /*point*/,
   const vector<Extensive>& /*fluxes*/,
   const PhysicalGeometry1D& /*pg*/,
   double /*t*/,
   double /*dt*/) const
{
  Extensive res;
  res.mass = 0;
  res.momentum = Vector2D(0,0);
  res.energy = 0;
  res.tracers = vector<double>(state.getCells().at(0).tracers.size(),0);
  return res;
}
		  
