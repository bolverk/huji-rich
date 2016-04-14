#include "ConstantGravity.hpp"

ConstantGravity::ConstantGravity(Vector2D const& force):
force_(force){}

Vector2D ConstantGravity::operator()
(const Tessellation& /*tess*/,
 const vector<ComputationalCell>& /*cells*/,
 const vector<Extensive>& /*fluxes*/,
 const double /*time*/,
 const int /*point*/,
 TracerStickerNames const& /*tracerstickernames*/) const
{
  return force_;
}
