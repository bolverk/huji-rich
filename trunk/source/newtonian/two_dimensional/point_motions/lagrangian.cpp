#include "lagrangian.hpp"

Vector2D Lagrangian::CalcVelocity
(int index, Tessellation const& /*tessellation*/,
 vector<Primitive> const& primitives,double /*time*/)
{
  return primitives[index].Velocity;
}
