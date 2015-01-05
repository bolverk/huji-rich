#ifndef POINT_MOTION_HPP
#define POINT_MOTION_HPP 1

#include "../../3D/Tessellation/Vector3D.hpp"

class PointMotion3D
{
public:

  virtual Vector3D operator()(const Vector3D& pos) const = 0;

  virtual ~PointMotion3D(void);
};

#endif // POINT_MOTION_HPP
