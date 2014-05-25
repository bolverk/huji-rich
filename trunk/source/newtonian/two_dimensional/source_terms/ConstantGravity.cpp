#include "ConstantGravity.hpp"

ConstantGravity::ConstantGravity(Vector2D const& force):
force_(force){}

Vector2D ConstantGravity::Calculate
(Tessellation const& /*tess*/,
 vector<Primitive> const& /*cells*/,
 int /*point*/,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 double /*time*/,
 double /*dt*/)
{
	return force_;
}
