#include "zero_force.hpp"

vector<Extensive> ZeroForce::operator()
(const Tessellation& tess,
 const PhysicalGeometry& /*pg*/,
 const vector<ComputationalCell>& /*cells*/,
 const vector<Extensive>& /*fluxes*/,
 const vector<Vector2D>& /*point_velocities*/,
 const double /*t*/) const
{
  return vector<Extensive>(static_cast<size_t>(tess.GetPointNo()));
}
