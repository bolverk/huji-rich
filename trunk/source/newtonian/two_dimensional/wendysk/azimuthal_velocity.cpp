#include "azimuthal_velocity.hpp"

AzimuthalVelocity::AzimuthalVelocity(Func1Var const& radial,
				     Vector2D const& center,
				     char comp):
  radial_(radial),
  center_(center),
  comp_(comp) {}

double AzimuthalVelocity::EvalAt(Vector2D const& p) const
{
  const Vector2D rvec = p - center_;
  const double radius = abs(rvec);
  const double vq = radial_.eval(radius);
  if(comp_=='x')
    return -vq*rvec.y/radius;
  else if(comp_=='y')
    return vq*rvec.x/radius;
  else
    throw "Unknown component name "+comp_;
}
