#include "coriolis.hpp"

Coriolis::Coriolis(double angular_velocity):
  omega_(angular_velocity) {}

Conserved Coriolis::Calculate
(Tessellation const* tess,
 vector<Primitive> const& cells,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const* /*hbc*/,
 vector<vector<double> > const& /*tracers*/,
 vector<double>& /*dtracer*/,
 double /*time*/,
 double /*dt*/)
{
  const Vector2D velocity = cells[point].Velocity;
  const Vector2D acceleration = 2*omega_*zcross(velocity);
  const double mass = tess->GetVolume(point)*cells[point].Density;
  const Vector2D force = mass*acceleration;
  return Conserved(0,force,0);
}
