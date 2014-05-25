#include "eulerian.hpp"

Vector2D Eulerian::CalcVelocity
(int /*index*/, Tessellation const* /*tessellation*/,
 vector<Primitive> const& /*primitives*/,double /*time*/)
{
  return Vector2D(0,0);
}
