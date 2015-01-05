#ifndef POINT_MOTION_HPP
#define POINT_MOTION_HPP 1

class PointMotion3D
{
public:

  virtual Vector3D operator()(const Vector3D& pos) const = 0;
};

#endif // POINT_MOTION_HPP
