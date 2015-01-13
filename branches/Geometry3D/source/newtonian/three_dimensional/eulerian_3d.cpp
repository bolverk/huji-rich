#include "eulerian_3d.hpp"

Eulerian3D::Eulerian3D(void) {}

Vector3D Eulerian3D::operator()(const Vector3D& /*pos*/) const
{
  return Vector3D();
}
