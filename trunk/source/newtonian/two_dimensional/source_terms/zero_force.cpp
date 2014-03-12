#include "zero_force.hpp"

Conserved ZeroForce::Calculate
(Tessellation const& /*tess*/,
 vector<Primitive> const& /*cells*/,
 int /*point*/,vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 vector<vector<double> > const &/*tracer_extensive*/,
 vector<double> &/*dtracer*/,
 double /*t*/,
 double /*dt*/)
{
  return Conserved();
}

Vector2D ZeroAcceleration::Calculate
(Tessellation const& /*tess*/,
 vector<Primitive> const& /*cells*/,
 int /*point*/,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 double /*time*/,
 double /*dt*/)
{
  return Vector2D(0,0);
}
