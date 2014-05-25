#include "centripetal.hpp"

Centripetal::Centripetal(double angular_velocity,
			 Vector2D const& center):
  omega_(angular_velocity),
  center_(center) {}

Vector2D Centripetal::Calculate
(Tessellation const* tess,
 vector<Primitive> const& /*cells*/,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const* /*hbc*/,
 double /*time*/,
 double /*dt*/)
{
  const Vector2D rvec(tess->GetCellCM(point) - center_);
  return pow(omega_,2)*rvec;
}
