#include "CenterGravity.hpp"

CenterGravity::CenterGravity
(double M, double Rmin, const Vector2D& center):
  M_(M),Rmin_(Rmin),center_(center){}

Vector2D CenterGravity::operator()
(const Tessellation& tess,
 const vector<ComputationalCell>& /*cells*/,
 const vector<Extensive>& /*fluxes*/,
 const double /*time*/,
 const int point,
 TracerStickerNames const& /*tracerstickernames*/) const
{
  const Vector2D pos(tess.GetCellCM(point)-center_);
  const double r = abs(pos);
  return (-1)*pos*M_/(r*r*r+Rmin_*Rmin_*Rmin_);
}
